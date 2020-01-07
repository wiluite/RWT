#include <transwarp.h>
#include <iostream>
#include <boost/math/constants/constants.hpp>

namespace daq
{
    using sample_type_ = int32_t;

    constexpr size_t ricker_filter_max_size = 20000;
    using ricker_filter_data_type = std::array<double, ricker_filter_max_size + 1>;

    struct ricker_filter_data {
        ricker_filter_data_type& data;
        size_t& sz;
    };

    /**
        * Вычисление длины вейвлет-фильтра по частоте искомого сигнала и частоте дискретизации АЦП на канал
        * @param f Искомая частота сигнала
        * @param rate Частота дискретизации АЦП
        * @param width Константа 5
        * @return Предварительный размер фильтра
        */
    int ricker_filter_size(double f, double rate, double width = 5.0)
    {
        constexpr double ricker_coefficient = 2.2508;
        return static_cast<int>(lround (width / ricker_coefficient / f * rate + 0.5));
    }

    int align_by_ten(int sz)
    {
        return sz / 10 * 10 + 10;
    }

    int make_odd(int sz)
    {
        return (sz % 2) ? sz : (sz + 1);
    }

    void create_filter(int dots, ricker_filter_data & mh_filter)
    {
        auto const tricky_dots = (dots - 1) / 10;
        double const c = 2 / ((sqrt (3)) * pow (boost::math::constants::pi<double>(), 0.25));
        double const one_per_mh_dots = 1.0 / tricky_dots;
        double t = -5.0;

        for (int i = 0; i < dots; ++i)
        {
            const double pow_t_2 = pow (t, 2);
            mh_filter.data[i] = c * exp(-pow_t_2 / 2) * (1 - pow_t_2);
            t += one_per_mh_dots;
            std::cout << (double)mh_filter.data[i] << std::endl;
        }

        mh_filter.sz = dots;
    }
}



int main()
{
    using namespace daq;

    ricker_filter_data_type filter_data {};
    size_t filter_size {};
    ricker_filter_data filter {filter_data, filter_size};

    create_filter (make_odd(align_by_ten(ricker_filter_size(11.0, 19531.0))), filter);

    std::cout << "test" << std::endl;
    std::cout << filter.sz << std::endl;
    return 0;
}

