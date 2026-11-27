#include <math.h>

#include "aerodynamics/AirfoilPolar.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using turbine_core::real_t;
using turbine_core::uint_t;
using turbine_core::Vector2;
using turbine_core::aerodynamics::AirfoilPolar;

namespace {

AirfoilPolar makePolar(const std::vector<real_t> & aoa, real_t cl0, real_t clAmp, real_t cdBase) {
    std::vector<Vector2<real_t>> polar;
    for (real_t a : aoa)
        polar.emplace_back(cl0 + clAmp * std::sin(real_t(3) * a), cdBase + real_t(0.1) * std::cos(real_t(5) * a) * std::cos(real_t(5) * a));
    return AirfoilPolar(aoa.size(), aoa.data(), polar.data());
}

std::vector<real_t> uniformGrid(uint_t n) {
    const real_t pi = real_t(3.14159265358979323846);
    std::vector<real_t> aoa(n);
    for (uint_t i = 0; i < n; ++i)
        aoa[i] = -pi + real_t(2) * pi * real_t(i) / real_t(n - 1);
    return aoa;
}

real_t maxBlendError(real_t weight, const AirfoilPolar & lower, const AirfoilPolar & upper, const AirfoilPolar & blend) {
    const real_t pi = real_t(3.14159265358979323846);
    real_t err = 0;
    for (int i = 0; i <= 4000; ++i) {
        const real_t a = -pi + real_t(2) * pi * real_t(i) / real_t(4000);
        const auto lo = lower(a);
        const auto up = upper(a);
        const auto b = blend(a);
        for (int k = 0; k < 2; ++k) {
            const real_t exact = (real_t(1) - weight) * lo[k] + weight * up[k];
            err = std::max(err, std::abs(b[k] - exact));
        }
    }
    return err;
}

int failures = 0;

void check(bool ok, const char * what) {
    std::printf("%s: %s\n", ok ? "PASS" : "FAIL", what);
    if (!ok) ++failures;
}

}

int main() {
    const real_t tol = real_t(2e-5);

    {
        const AirfoilPolar coarse = makePolar(uniformGrid(3), real_t(0), real_t(0), real_t(0.5));
        const AirfoilPolar fine = makePolar(uniformGrid(73), real_t(0.3), real_t(0.8), real_t(0.02));
        const AirfoilPolar blend(real_t(0.5), coarse, fine);
        const real_t err = maxBlendError(real_t(0.5), coarse, fine, blend);
        std::printf("coarse 3 points, fine 73 points: blend has %u points, max error %g\n", unsigned(blend.length()), double(err));
        check(err < tol, "blend of 3-point and 73-point polars equals the exact blend");
        check(blend.length() >= fine.length(), "blend keeps every point of the fine polar");
    }

    {
        const AirfoilPolar coarse = makePolar(uniformGrid(3), real_t(0), real_t(0), real_t(0.5));
        const AirfoilPolar fine = makePolar(uniformGrid(73), real_t(0.3), real_t(0.8), real_t(0.02));
        const AirfoilPolar blend(real_t(0.5), fine, coarse);
        const real_t err = maxBlendError(real_t(0.5), fine, coarse, blend);
        std::printf("fine 73 points, coarse 3 points: blend has %u points, max error %g\n", unsigned(blend.length()), double(err));
        check(err < tol, "blend with the fine polar as lower station equals the exact blend");
    }

    {
        const AirfoilPolar a = makePolar(uniformGrid(37), real_t(0.1), real_t(0.5), real_t(0.03));
        const AirfoilPolar b = makePolar(uniformGrid(37), real_t(0.4), real_t(0.2), real_t(0.06));
        const AirfoilPolar blend(real_t(0.25), a, b);
        check(blend.length() == a.length(), "identical grids give a table of the same length");
        check(maxBlendError(real_t(0.25), a, b, blend) < tol, "identical grids give the exact blend");
    }

    {
        const AirfoilPolar coarse = makePolar(uniformGrid(3), real_t(0), real_t(0), real_t(0.5));
        const AirfoilPolar fine = makePolar(uniformGrid(73), real_t(0.3), real_t(0.8), real_t(0.02));
        check(maxBlendError(real_t(0), coarse, fine, AirfoilPolar(real_t(0), coarse, fine)) < tol, "weight 0 reproduces the lower polar");
        check(maxBlendError(real_t(1), coarse, fine, AirfoilPolar(real_t(1), coarse, fine)) < tol, "weight 1 reproduces the upper polar");
    }

#ifdef AIRFOIL_DATA_DIR
    {
        const std::string dir = AIRFOIL_DATA_DIR;
        const AirfoilPolar cylinder(dir + "/nrel_5mw_Cylinder2.foil");
        const AirfoilPolar du40(dir + "/DU40_A17.foil");
        const AirfoilPolar blend(real_t(0.5), cylinder, du40);
        const real_t err = maxBlendError(real_t(0.5), cylinder, du40, blend);
        std::printf("Cylinder2 (%u points) and DU40 (%u points): blend has %u points, max error %g\n",
                    unsigned(cylinder.length()), unsigned(du40.length()), unsigned(blend.length()), double(err));
        check(err < tol, "NREL 5 MW Cylinder2 to DU40 blend equals the exact blend");
        check(blend.length() >= du40.length(), "NREL 5 MW Cylinder2 to DU40 blend keeps every DU40 point");
    }
#endif

    return failures == 0 ? 0 : 1;
}
