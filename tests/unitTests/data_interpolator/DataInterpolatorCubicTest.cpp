#include <cmath>
#include <cstdlib>
#include <vector>

#include <core/Environment.h>
#include <core/debug/TestSubsystem.h>
#include <core/DataTypes.h>

#include "data_interpolator/DataInterpolator.h"
#include "wind_turbine_core/math/Vector3.h"

namespace data_interpolator_test {

    using namespace walberla;

    using CubicScalar_T = turbine_core::ScalarDataInterpolator<true>;
    using CubicVector_T = turbine_core::Vector3DataInterpolator<true>;

    real_t polynomial( real_t x ) {
        return real_t(0.5) * x * x * x - real_t(2) * x * x + real_t(3) * x - real_t(1);
    }

    bool close( real_t a, real_t b ) {
        return std::fabs( a - b ) <= real_t(1e-3) * ( real_t(1) + std::fabs( b ) );
    }

    void cubicScalarTest() {

        const uint_t n = 9;
        std::vector<real_t> x(n), y(n);
        for( uint_t i = 0; i < n; ++i ) {
            x[i] = real_t(0.3) * real_t(i) + real_t(0.02) * real_t(i * i);
            y[i] = polynomial( x[i] );
        }

        CubicScalar_T interpolator;
        interpolator.setPoints( n, x.data(), y.data() );

        for( uint_t i = 0; i < n; ++i ) {
            WALBERLA_CHECK( std::isfinite( interpolator.weight()[i] ) );
            WALBERLA_CHECK( interpolator.weight()[i] != real_t(0) );
        }

        for( uint_t i = 0; i < n; ++i ) {
            WALBERLA_CHECK( close( interpolator( x[i] ), y[i] ) );
        }

        for( uint_t i = 0; i + 1 < n; ++i ) {
            const real_t mid = real_t(0.5) * ( x[i] + x[i+1] );
            WALBERLA_CHECK( close( interpolator( mid ), polynomial( mid ) ) );
        }
    }

    void cubicVectorTest() {

        const uint_t n = 7;
        std::vector<real_t> x(n);
        std::vector< turbine_core::Vector3<real_t> > y(n);
        for( uint_t i = 0; i < n; ++i ) {
            x[i] = real_t(i);
            y[i] = turbine_core::Vector3<real_t>( polynomial( x[i] ), real_t(2) * x[i], real_t(1) );
        }

        CubicVector_T interpolator;
        interpolator.setPoints( n, x.data(), y.data() );

        const real_t q = real_t(3.5);
        const auto v = interpolator( q );
        WALBERLA_CHECK( close( v[0], polynomial( q ) ) );
        WALBERLA_CHECK( close( v[1], real_t(2) * q ) );
        WALBERLA_CHECK( close( v[2], real_t(1) ) );
    }

    int main( int argc, char** argv ) {
        Environment walberlaEnv( argc, argv );

        debug::enterTestMode();

        cubicScalarTest();
        cubicVectorTest();

        return EXIT_SUCCESS;
    }

}

int main( int argc, char** argv ) {
    return data_interpolator_test::main(argc, argv);
}
