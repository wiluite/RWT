//
// Created by wiluite on 18.01.2020-2022
//

#define BOOST_TEST_MODULE boost_test_module_
#include <boost/test/unit_test.hpp> // UTF ??
#include "detail/rwt_omp.h"
#include <random>

using namespace detail;
using namespace detail::rwt;

constexpr double f = 15.0;      // Hz
constexpr double rate = 95000;  // Hz

std::array<detail::rwt::transform_result_type, static_cast<size_t>(rate)> transform_uni_result, transform_omp_result, transform_tw_result, transform_honest_omp_result;

template<class T, size_t N>
struct epsilon_compare_test
{
    explicit epsilon_compare_test(std::array<T, N> & arr)
    {
        constexpr double epsilon_ = 0.00000015; //0.000000001;
        std::for_each(std::begin(transform_uni_result), std::end(transform_uni_result), [&](auto && elem)
        {
            auto const curr_index = std::addressof(elem)-std::addressof(arr[0]);
            if (elem != arr[curr_index])
                BOOST_REQUIRE (abs (elem - arr[curr_index]) < epsilon_);
        });
    }
};

BOOST_AUTO_TEST_CASE(test1) {
    auto const filter = create_filter (make_odd(align_by_ten(ricker_filter_size(f, rate))));

    signal_sequence_type ss (static_cast<size_t>(rate) + filter.actual_sz);
    std::mt19937 rng;
    std::generate_n(std::begin(ss), std::distance(std::begin(ss), std::end(ss)), rng);

    transform_uni (ss, std::addressof(transform_uni_result[0]), filter);
    transform_omp (ss, std::addressof(transform_omp_result[0]), filter);
    transform_tw (ss, std::addressof(transform_tw_result[0]), filter);
    transform_honest_omp (ss, std::addressof(transform_honest_omp_result[0]), filter);
    BOOST_REQUIRE (transform_uni_result == transform_tw_result);

    // transform_omp results differ from the others for a very small epsilon (don't ask why)
    epsilon_compare_test test1 (transform_omp_result);
    epsilon_compare_test test2 (transform_honest_omp_result);
}
