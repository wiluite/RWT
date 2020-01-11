#include <transwarp.h>
#include <iostream>
#include <boost/math/constants/constants.hpp>
#include <parallel/numeric>
#include <boost/shared_array.hpp>

namespace daq
{
    constexpr double f = 15.0;       // Hz
    constexpr double rate = 95000.0; // Hz

    int ricker_filter_size(double f_, double rate_, double width = 5.0) noexcept
    {
        constexpr double ricker_coefficient = 2.2508;
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

    constexpr size_t ricker_filter_max_size = 20000;
    using ricker_filter_data_type = std::array<double, ricker_filter_max_size + 1>;

    struct ricker_filter_data {
        ricker_filter_data_type& data;
        size_t& sz;
    };

    void create_filter(int dots, ricker_filter_data & rfd) noexcept
    {
        assert (rfd.data.size() > dots);

        auto const tricky_dots = (dots - 1) / 10;
        double const c = 2 / ((sqrt (3)) * pow (boost::math::constants::pi<double>(), 0.25));
        double const one_per_mh_dots = 1.0 / tricky_dots;
        double t = -5.0;

        for (int i = 0; i < dots; ++i)
        {
            const double pow_t_2 = pow (t, 2);
            rfd.data[i] = c * exp(-pow_t_2 / 2) * (1 - pow_t_2);
            t += one_per_mh_dots;
        }

        rfd.sz = dots;
    }

    using sample_type = int32_t;
    using signal_sequence_type = std::vector<sample_type>;
    using transform_result_type = double;

    using sample_quantity_type = int;
    using prepare_type = std::tuple<sample_quantity_type, const signal_sequence_type &>;
    using output_type = std::shared_ptr<std::vector<transform_result_type>>;

    output_type partial_task(prepare_type input, uint8_t part_index, uint8_t part_number, ricker_filter_data const & rfd)
    {
        auto const sample_quantity = std::get<0>(input);
        auto const & signal_sequence = std::get<1>(input);

        auto const default_iteration_number = sample_quantity / part_number;
        auto const begin_iter = std::begin(signal_sequence) + ((part_index - 1) * default_iteration_number);
        auto const end_iter = (part_index != part_number) ? (std::begin(signal_sequence) + ((part_index) * default_iteration_number + rfd.sz)) : std::end(signal_sequence);

        sample_type signal_subsequence [std::distance(begin_iter, end_iter)];
        ricker_filter_data_type::value_type filter [rfd.sz];

        std::copy (begin_iter, end_iter, signal_subsequence);
        std::copy (std::begin(rfd.data), std::begin(rfd.data) + rfd.sz, filter);

        auto const total_iterations = (part_index != part_number) ? (sample_quantity / part_number) : (sample_quantity / part_number + (sample_quantity % part_number));
        auto const data = std::make_shared<std::vector<transform_result_type>>(total_iterations);

        auto curr = std::begin(*data);
        auto const end = curr + total_iterations;

        auto pure_s_iter = signal_subsequence;

        while (curr != end)
        {
            *curr++ = static_cast<transform_result_type>(std::inner_product(filter, filter + rfd.sz, pure_s_iter++, 0.0) / rfd.sz);
        }

        return data;
    }

    namespace tw = transwarp;
    constexpr int hardware_threads = 8; /*std::thread::hardware_concurrency();*/
    tw::parallel & get_executor()
    {
        static tw::parallel executor(hardware_threads);
        return executor;
    }

    static std::shared_ptr<tw::task<output_type>> make_task_by_index (uint8_t index, signal_sequence_type const & ss, size_t from, size_t to, ricker_filter_data const & fd)
    {
        return tw::make_task(tw::consume, [&] (prepare_type input, uint8_t index)
                             {
                                 return partial_task(input, index, hardware_threads, fd);
                             },
                             tw::make_value_task (std::make_tuple(to - from, ss)), tw::make_value_task(index)
        );
    }

    static std::shared_ptr<tw::task<void>> gather_task(signal_sequence_type const & ss, size_t from, size_t to,
                                                       transform_result_type * const result, ricker_filter_data const & filter_data) {

        std::vector<std::shared_ptr<tw::task<output_type>>> vec (hardware_threads);
        uint8_t curr_idx = 0;
        std::for_each(std::begin(vec), std::end(vec), [&](std::shared_ptr<tw::task<output_type>> & elem)
                      {
                          elem = make_task_by_index(static_cast<uint8_t>(++curr_idx), ss, from, to, filter_data);
                      }
        );

        // capture by value
        return tw::make_task(tw::consume, [=](std::vector<output_type> const & parents) {
            size_t copy_pos = 0;
            for (int i = 0; i < hardware_threads; ++i)
            {
                std::copy (parents[i]->begin(), parents[i]->begin()+parents[i]->size(), &result[copy_pos]);
                copy_pos += parents[i]->size();
            }}, vec);
    }

    auto inner_product_lambda = [](auto&&... args){return std::inner_product(decltype(args)(args)...);};
    auto parallel_inner_product_lambda = [](auto&&... args){return __gnu_parallel::inner_product(decltype(args)(args)...);};

    template <typename L>
    inline void seq_transform(signal_sequence_type const & ss, transform_result_type * const result, ricker_filter_data const & filter_data, L const & lambda)
    {
        auto const iterations = static_cast<size_t>(rate);
        int curr_iteration = 0;
        while (curr_iteration < iterations)
        {
            result[curr_iteration] = static_cast<transform_result_type> (
                    lambda(std::begin(filter_data.data), std::begin(filter_data.data) + filter_data.sz, std::begin(ss) + curr_iteration, 0.0) / filter_data.sz
            );
            ++curr_iteration;
        }
    }

    void transform_uni(signal_sequence_type const & ss, transform_result_type * const result, ricker_filter_data const & filter_data)
    {
        seq_transform (ss, result, filter_data, inner_product_lambda);
    }

    void transform_omp (signal_sequence_type const & ss, transform_result_type * const result, ricker_filter_data const & filter_data)
    {
        seq_transform (ss, result, filter_data, parallel_inner_product_lambda);
    }

    void transform_tw(signal_sequence_type const & ss, transform_result_type * const result, ricker_filter_data const & filter_data)
    {
        auto final_task = gather_task(ss, 0, static_cast<size_t>(rate), result, filter_data);
        final_task->schedule_all(get_executor());
        final_task->get();
    }
}

template <typename F, typename ...T>
void measure_it(char const* f_name, F const & f, T &&... arguments)
{
    auto const time_point1 = std::chrono::high_resolution_clock::now();
    f (std::forward<T>(arguments)...);
    auto const time_point2 = std::chrono::high_resolution_clock::now();
    std::cout << f_name << std::chrono::duration_cast<std::chrono::milliseconds>(time_point2-time_point1).count() << std::endl;
}

int main()
{
    using namespace daq;

    ricker_filter_data_type filter_data {};
    size_t filter_size {};
    ricker_filter_data filter {filter_data, filter_size};

    create_filter (make_odd(align_by_ten(ricker_filter_size(f, rate))), filter);

    signal_sequence_type ss (static_cast<size_t>(rate) + filter.sz);
    std::iota(std::begin(ss), std::end(ss), -rate/2);

    boost::shared_array<transform_result_type> wavelet_transform {new transform_result_type[static_cast<size_t>(rate)]{0}};

    measure_it("transform_uni: ", transform_uni, ss, wavelet_transform.get(), filter);
    measure_it("transform_omp: ", transform_omp, ss, wavelet_transform.get(), filter);
    measure_it("transform_tw:  ", transform_tw, ss, wavelet_transform.get(), filter);

    return 0;
}

