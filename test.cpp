//
// Created by wiluite on 18.01.2020.
//

#define BOOST_TEST_MODULE boost_test_module_
#include <boost/test/unit_test.hpp> // UTF ??
#define MAIN_IS_ABSENT
#include "main.cpp"

BOOST_AUTO_TEST_CASE(test1) {

    using namespace daq;

    ricker_filter_data_type filter_data {};
    size_t filter_size {};
    ricker_filter filter {filter_data, filter_size};

    create_filter (make_odd(align_by_ten(ricker_filter_size(f, rate))), filter);

    signal_sequence_type ss (static_cast<size_t>(rate) + filter.sz);
    std::iota(std::begin(ss), std::end(ss), -rate/2);

    transform_uni (ss, std::addressof(transform_uni_result[0]), filter);
    transform_omp (ss, std::addressof(transform_omp_result[0]), filter);
    transform_tw (ss, std::addressof(transform_tw_result[0]), filter);

    // transform_uni always equals transform_tw
    BOOST_REQUIRE (transform_uni_result == transform_tw_result);

    // transform_omp_result differ from the others for a very small epsilon (dont ask why)
    auto cnt = 0;
    constexpr double epsilon_ = 0.0000000001;
    for (auto elem : transform_uni_result)
    {
        if (elem != transform_omp_result[cnt])
        {
            BOOST_REQUIRE ((abs (elem-transform_omp_result[cnt]) < epsilon_));
        }
        ++cnt;
    }
}
