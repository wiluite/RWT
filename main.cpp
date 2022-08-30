#include "detail/rwt_omp.h"
#include <random>

constexpr double f = 15.0;
constexpr double rate = 95000;

std::array<detail::rwt::transform_result_type, static_cast<size_t>(rate)> transform_uni_result, transform_omp_result, transform_tw_result, transform_honest_omp_result;

using namespace detail;
using namespace detail::rwt;

 struct compare_repeat_strategy{
    compare_repeat_strategy(size_t iter_cnt, size_t iterations, ricker_filter const & filter, signal_sequence_type const & ss) :
    iter_cnt_(iter_cnt), iterations_(iterations), filter_(filter), ss_(ss) {}

    void process()
    {
        while (iter_cnt_--)
        {
            if (iter_cnt_ % 2) {
                times_sums[0] += measure_it("transform-tw:  ", transform_tw, ss_,
                                            std::addressof(transform_tw_result[0]), filter_);
                times_sums[1] += measure_it("transform-honest-omp: ", transform_honest_omp, ss_,
                                            std::addressof(transform_honest_omp_result[0]), filter_);
            } else {
                times_sums[1] += measure_it("transform-honest-omp: ", transform_honest_omp, ss_,
                                            std::addressof(transform_honest_omp_result[0]), filter_);
                times_sums[0] += measure_it("transform-tw:  ", transform_tw, ss_,
                                            std::addressof(transform_tw_result[0]), filter_);
            }
        }
        std::cout << "Average:\n";
        std::cout << "tw:  " << times_sums[0] / iterations_ << '\n';
        std::cout << "honest_omp: " << times_sums[1] / iterations_ << '\n';
    }
 private:
    size_t iter_cnt_;
    size_t iterations_;
    signal_sequence_type const & ss_;
    ricker_filter const & filter_;
    std::array<decltype(std::chrono::duration_cast<std::chrono::milliseconds>
             (std::chrono::high_resolution_clock::now() - std::chrono::high_resolution_clock::now()).count()), 2> times_sums {};

};

struct vacation_tw_repeat_strategy {
    vacation_tw_repeat_strategy (size_t iter_cnt, ricker_filter const & filter, signal_sequence_type const & ss) :
    iter_cnt_(iter_cnt), filter_(filter), ss_(ss) {}
    void process()
    {
        while (iter_cnt_--)
        {
            if (iter_cnt_ % 2) {
                auto const newtime = measure_it("transform-tw:  ", transform_tw,  ss_, std::addressof(transform_tw_result[0]),  filter_);
                if (newtime < 1000)
                    std::this_thread::sleep_for(std::chrono::milliseconds(1000-newtime));
            } else {
                auto const newtime = measure_it("transform-tw:  ", transform_tw,  ss_, std::addressof(transform_tw_result[0]),  filter_);
                if (newtime < 1000)
                    std::this_thread::sleep_for(std::chrono::milliseconds(1000-newtime));
            }
        }
    }
private:
    size_t iter_cnt_;
    signal_sequence_type const & ss_;
    ricker_filter const & filter_;
};

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << argv[0] << " <number>\n";
        exit(1);
    }

    size_t iter_cnt, iterations;
    try {
        iter_cnt = iterations = std::stoi(std::string(argv[1]));
    } catch (std::invalid_argument const & e)
    {
        std::cout << e.what() << '\n';
        exit(1);
    }

    ricker_filter filter{};
    try {
        filter = create_filter (make_odd(align_by_ten(ricker_filter_size(f, rate))));
    } catch (std::exception const & e)
    {
        std::cout << e.what() << '\n';
        return 1;
    }

    signal_sequence_type ss (static_cast<size_t>(rate) + filter.actual_sz);
    std::mt19937 rng;
    std::generate_n(std::begin(ss), std::distance(std::begin(ss), std::end(ss)), rng);

#if 1
    measure_it("transform-uni:         ", transform_uni, ss, std::addressof(transform_uni_result[0]), filter);
#endif
#if 0
    measure_it("transform-omp:         ", transform_omp, ss, std::addressof(transform_omp_result[0]), filter);
    measure_it("transform-tw:          ", transform_tw, ss, std::addressof(transform_tw_result[0]), filter);
    measure_it("transform-honest-omp:  ", transform_honest_omp, ss, std::addressof(transform_honest_omp_result[0]), filter);
#else
    compare_repeat_strategy rs (iter_cnt, iterations, filter, ss);
    //vacation_tw_repeat_strategy rs (iter_cnt, filter, ss);
    rs.process();
#endif
    return 0;
}

