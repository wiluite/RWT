//
// Created by root on 8/25/22.
//
#pragma once

#include <transwarp.h>
#include <iostream>
#include <boost/math/constants/constants.hpp>
#ifndef __clang__
#include <parallel/numeric>
#include <parallel/algorithm>
#else
#include <omp.h>
#include <numeric>
#endif

namespace detail
{
    namespace rwt
    {
        int ricker_filter_size(double f_, double rate_, double width = 5.0) noexcept
        {
            constexpr double ricker_coefficient = 2.2508; // About "world-wide" constant do not touch it.
            return static_cast<int>(lround (width / ricker_coefficient / f_ * rate_ + 0.5));
        }

        // TODO: Improve make_odd(), align_by_ten(), ricker_filter_size() functions to do a simmetrical filter without sly tricks
        int align_by_ten(int sz) noexcept
        {
            return sz / 10 * 10 + 10;
        }

        int make_odd(int sz) noexcept
        {
            return (sz % 2) ? sz : (sz + 1);
        }

        constexpr size_t ricker_filter_max_size = 400'000;
        using ricker_filter_data_type = std::array<double, ricker_filter_max_size + 1>;

        struct ricker_filter {
            ricker_filter_data_type data;
            int actual_sz;
        };

        ricker_filter create_filter(int dots) noexcept(false)
        {
            if (ricker_filter_max_size <= dots)
                throw std::runtime_error ("Insufficient filter buffer size (ricker_filter_max_size). Enlarge it, please.");

            auto const tricky_dots = (dots - 1) / 10;
            double const c = 2 / ((sqrt (3)) * pow (boost::math::constants::pi<double>(), 0.25));
            double const one_per_mh_dots = 1.0 / tricky_dots;
            double t = -5.0;

            ricker_filter rf {};
            for (int i = 0; i < dots; ++i)
            {
                const double pow_t_2 = pow (t, 2);
                rf.data[i] = c * exp(-pow_t_2 / 2) * (1 - pow_t_2);
                t += one_per_mh_dots;
            }
            rf.actual_sz = dots;
            return rf;
        }

        using sample_type = int32_t;
        using signal_sequence_type = std::vector<sample_type>;
        using transform_result_type = double;

        using sample_quantity_type = int;
        using prepare_type = std::tuple<sample_quantity_type, const signal_sequence_type &>;
        using output_type = std::shared_ptr<std::vector<transform_result_type>>;

        struct use_local_filter_strategy final
        {
            explicit use_local_filter_strategy (ricker_filter const & rfd) noexcept : filter(rfd.actual_sz)
            {
                std::copy (std::begin(rfd.data), std::begin(rfd.data) + rfd.actual_sz, std::begin(filter));
            };
            auto operator()() const noexcept
            {
                return filter.begin();
            };
        private:
            std::vector <ricker_filter_data_type::value_type> filter;
        };

        std::once_flag filt_flag;
        static ricker_filter_data_type::value_type filter_[ricker_filter_max_size];
        struct use_static_filter_strategy
        {
            explicit use_static_filter_strategy (ricker_filter const & rfd) noexcept
            {
                std::call_once (filt_flag, [&]() {std::copy(std::begin(rfd.data), std::begin(rfd.data)+rfd.actual_sz, filter_);});
            };
            auto operator()() const noexcept
            {
                return &filter_[0];
            };
        };

        template <class FILTER_LOCATION_STRATEGY = use_static_filter_strategy>
        output_type partial_task(prepare_type input, uint8_t part_index, uint8_t parts, ricker_filter const & rfd) noexcept
        {
            auto const sample_quantity = std::get<0>(input);
            auto const & signal_sequence = std::get<1>(input);

            auto const default_iteration_number = sample_quantity / parts;
            auto const begin_iter = std::begin(signal_sequence) + ((part_index - 1) * default_iteration_number);
            auto const end_iter = (part_index != parts) ? (std::begin(signal_sequence) + (part_index * default_iteration_number + rfd.actual_sz)) : std::end(signal_sequence);

            alignas (64) sample_type signal_subsequence [std::distance(begin_iter, end_iter)];
            std::copy (begin_iter, end_iter, signal_subsequence);

            auto const total_iterations = (part_index != parts) ? (sample_quantity / parts) : (sample_quantity / parts + (sample_quantity % parts));
            auto data = std::make_shared<std::vector<transform_result_type>>(total_iterations);

            auto dest = std::begin(*data);

            auto signal_iter = signal_subsequence;

            FILTER_LOCATION_STRATEGY filt_strat (rfd);

            std::transform(signal_iter, signal_iter + total_iterations, dest, [&](auto & elem)
            {
                return std::inner_product(filt_strat(), filt_strat() + rfd.actual_sz, std::addressof(elem), 0.0) / rfd.actual_sz;
            });

            return data;
        }

        namespace tw = transwarp;
        constexpr int hardware_threads = 8; /*std::thread::hardware_concurrency();*/
        tw::parallel & get_executor()
        {
            static tw::parallel executor(hardware_threads);
            return executor;
        }

        static std::shared_ptr<tw::task<output_type>> make_task_by_index (uint8_t index, signal_sequence_type const & ss, size_t from, size_t to, ricker_filter const & fd)
        {
            return tw::make_task(tw::consume, [&] (prepare_type input, uint8_t index)
                                 {
                                     return partial_task(input, index, hardware_threads, fd);
                                 },
                                 tw::make_value_task (std::make_tuple(to - from, ss)), tw::make_value_task(index)
            );
        }

        static std::shared_ptr<tw::task<void>> gather_task(signal_sequence_type const & ss, size_t from, size_t to,
                                                           transform_result_type * const result, ricker_filter const & filter_data) {
            std::vector<std::shared_ptr<tw::task<output_type>>> vec (hardware_threads);
            uint8_t task_idx = 1;
            for (auto & elem : vec)
            {
                elem = make_task_by_index(static_cast<uint8_t>(task_idx++), ss, from, to, filter_data);
            }

            // capture by value
            return tw::make_task(tw::consume, [=](std::vector<output_type> const & parents) {
                size_t copy_pos = 0;
                for (auto & elem : parents)
                {
                    std::copy (elem->begin(), elem->begin() + static_cast<int>(elem->size()), &result[copy_pos]);
                    copy_pos += elem->size();
                }
            }, vec);

        }

        void transform_tw(signal_sequence_type const & sig_seq, transform_result_type * const result, ricker_filter const & filter_data)
        {
            auto const sample_rate = sig_seq.size()-filter_data.actual_sz;
            auto final_task = gather_task(sig_seq, 0, static_cast<size_t>(sample_rate), result, filter_data);
            final_task->schedule_all(rwt::get_executor());
            final_task->get();
        }

    } //rwt

    auto inner_product_lambda = [](auto&&... args){return std::inner_product(decltype(args)(args)...);};
#if !defined(__clang__)
    auto parallel_inner_product_lambda = [](auto&&... args){return __gnu_parallel::inner_product(decltype(args)(args)...);};
#else
    auto parallel_inner_product_lambda = [](auto&&... args){return 0.0;};
#endif

    template <typename L>
    inline void seq_transform(rwt::signal_sequence_type const & sig_seq, rwt::transform_result_type * const result, rwt::ricker_filter const & filter_data, L const & lambda)
    {
        auto const sample_rate = static_cast<int>(sig_seq.size())-filter_data.actual_sz;
        std::transform(std::begin(sig_seq), std::begin(sig_seq) + sample_rate, result, [&](auto & elem){
            return lambda(std::begin(filter_data.data), std::begin(filter_data.data) + filter_data.actual_sz, std::addressof(elem), 0.0) / filter_data.actual_sz;
        });
    }

#if 0
    template <typename L>
    inline void seq_transform_like_for_each(rwt::signal_sequence_type const & sig_seq, rwt::transform_result_type * const result, rwt::ricker_filter const & filter_data, L const & lambda)
    {
        auto const sample_rate = sig_seq.size()-filter_data.actual_sz;
        auto const iterations = static_cast<size_t>(sample_rate);
        for (int curr_iteration = 0; curr_iteration != iterations; ++curr_iteration)
        {
            result[curr_iteration] = static_cast<rwt::transform_result_type> (
                    lambda(std::begin(filter_data.data), std::begin(filter_data.data) + filter_data.actual_sz, std::begin(sig_seq) + curr_iteration, 0.0) / filter_data.actual_sz
            );
        }
    }
#endif

    void transform_uni(rwt::signal_sequence_type const & sig_seq, rwt::transform_result_type * const result, rwt::ricker_filter const & filter_data)
    {
        seq_transform (sig_seq, result, filter_data, inner_product_lambda);
    }

    void transform_omp (rwt::signal_sequence_type const & sig_seq, rwt::transform_result_type * const result, rwt::ricker_filter const & filter_data)
    {
        seq_transform (sig_seq, result, filter_data, parallel_inner_product_lambda);
    }

#if 0
    template <typename L>
    inline void par_transform(rwt::signal_sequence_type const & sig_seq, rwt::transform_result_type * const result, rwt::ricker_filter const & filter_data, L const & lambda)
    {
        auto const sample_rate = static_cast<int>(sig_seq.size()) - filter_data.actual_sz;
        __gnu_parallel::transform(std::begin(sig_seq), std::begin(sig_seq) + sample_rate, result, [&](auto & elem){
            return lambda(std::begin(filter_data.data), std::begin(filter_data.data) + filter_data.actual_sz, std::addressof(elem), 0.0) / filter_data.actual_sz;
        });
    }
#endif
    void transform_honest_omp (rwt::signal_sequence_type const & sig_seq, rwt::transform_result_type * const result, rwt::ricker_filter const & filter_data)
    {
#if 0
        par_transform (sig_seq, result, filter_data, parallel_inner_product_lambda);
#else   // to ensure nothing for compiler rest non-optimized
        auto const sample_rate = static_cast<int>(sig_seq.size())-filter_data.actual_sz;
        __gnu_parallel::transform(std::begin(sig_seq), std::begin(sig_seq) + sample_rate, result, [&](auto & elem){
            return __gnu_parallel::inner_product(std::begin(filter_data.data), std::begin(filter_data.data) + filter_data.actual_sz, std::addressof(elem), 0.0) / filter_data.actual_sz;
        });
#endif
    }

    template <typename F, typename ...T>
    auto measure_it(char const* f_name, F const & fun, T &&... arguments)
    {
        auto const time_point1 = std::chrono::high_resolution_clock::now();
        fun (std::forward<T>(arguments)...);
        auto const time_point2 = std::chrono::high_resolution_clock::now();
        std::cout << f_name << std::chrono::duration_cast<std::chrono::milliseconds>(time_point2-time_point1).count() << '\n';
        return std::chrono::duration_cast<std::chrono::milliseconds>(time_point2-time_point1).count();
        static_assert (std::is_same_v<decltype(std::chrono::high_resolution_clock::now()) const, decltype(time_point1)>);
    }

}
