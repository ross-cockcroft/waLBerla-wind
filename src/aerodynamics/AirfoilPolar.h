
#pragma once

#ifndef TURBINECORE_AIRFOILPOLAR_H
#define TURBINECORE_AIRFOILPOLAR_H

#include "wind_turbine_core/math/Vector2.h"

#include "wind_turbine_core/ProjectDefines.h"
#include "data_interpolator/DataInterpolator.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <core/math/Constants.h>
#include <core/Abort.h>

namespace turbine_core {

    namespace aerodynamics {

        struct AirfoilPolar {

            HOST_DEVICE_PREFIX AirfoilPolar() {}

            HOST_PREFIX explicit AirfoilPolar(const std::string & filename) {
                std::ifstream file(filename);

                std::string line;

                std::vector<real_t> aoa;
                std::vector<Vector2<real_t>> polar;

                if (!file.is_open()) {
                    WALBERLA_ABORT("Error opening file " << filename << ", please check file existence.");
                } else {
                    std::getline(file, line); // skip header
                    std::getline(file, line); // skip separation line

                    std::getline(file, line);

                    std::getline(file, line); // skip separation line

                    while (std::getline(file, line)) {

                        bool is_empty = true;
                        for (char ch : line) {
                            is_empty = is_empty && isspace(ch);
                        }
                        if(is_empty)
                            continue;

                        std::stringstream in(line);
                        real_t a, cl, cd;
                        in >> a >> cl >> cd;

                        aoa.push_back(a / real_t(180) * walberla::math::pi);
                        polar.emplace_back(cl, cd);
                    }

                    file.close();
                }

                if (aoa.empty()) {
                    WALBERLA_ABORT("Airfoil polar file " << filename << " holds no data.");
                }

                interpolator.setPoints(aoa.size(), aoa.data(), polar.data());

                nPoints = aoa.size();
                aMin = aoa.front();
                aMax = aoa.back();

            }

            HOST_DEVICE_PREFIX AirfoilPolar(const uint_t & _nPoints, const real_t * aoa,
                                            const Vector2<real_t> * polar)
            {
                nPoints = _nPoints;
                aMin = aoa[0];
                aMax = aoa[_nPoints-1];
                interpolator.setPoints(_nPoints, aoa, polar);
            }

            HOST_DEVICE_PREFIX AirfoilPolar( const real_t & weight, const AirfoilPolar & lower, const AirfoilPolar & upper ) {

                const uint_t nLo = lower.interpolator.length();
                const uint_t nUp = upper.interpolator.length();
                const real_t * xLo = lower.interpolator.x();
                const real_t * xUp = upper.interpolator.x();

                auto * aoa = (real_t*) malloc((nLo + nUp) * sizeof(real_t));
                auto * polar = (Vector2<real_t>*) malloc((nLo + nUp) * sizeof(Vector2<real_t>));

                const real_t tol = real_t(1e-7);
                uint_t iLo = 0, iUp = 0, n = 0;
                while (iLo < nLo || iUp < nUp) {
                    real_t a;
                    if (iUp >= nUp || (iLo < nLo && xLo[iLo] <= xUp[iUp]))
                        a = xLo[iLo++];
                    else
                        a = xUp[iUp++];

                    if (n > 0 && a - aoa[n-1] <= tol)
                        continue;

                    auto lowerPolar = lower.interpolator(a);
                    auto upperPolar = upper.interpolator(a);

                    aoa[n] = a;
                    polar[n][0] = (static_cast<real_t>(1) - weight) * lowerPolar[0] + weight * upperPolar[0];
                    polar[n][1] = (static_cast<real_t>(1) - weight) * lowerPolar[1] + weight * upperPolar[1];
                    ++n;
                }

                interpolator.setPoints(n, aoa, polar);

                nPoints = n;
                aMin = aoa[0];
                aMax = aoa[n-1];

                free(aoa);
                free(polar);

            }

            HOST_DEVICE_PREFIX AirfoilPolar( const AirfoilPolar & ap )
                    : interpolator(ap.interpolator), nPoints(ap.nPoints), aMin(ap.aMin), aMax(ap.aMax)
            {}

            HOST_DEVICE_PREFIX AirfoilPolar & operator=(const AirfoilPolar & ap) {
                if(this != &ap) {
                    interpolator = ap.interpolator;
                    nPoints = ap.nPoints;
                    aMin = ap.aMin;
                    aMax = ap.aMax;
                }

                return *this;
            }

            HOST_DEVICE_PREFIX uint_t length() const {
                return interpolator.length();
            }

            HOST_DEVICE_PREFIX Vector2<real_t> interpolate( const real_t & x ) const {
                return interpolator(x);
            }

            HOST_DEVICE_PREFIX Vector2<real_t> operator()( const real_t & x ) const {
                return interpolator(x);
            }

            using Interpolator_T = data_interpolation::DataInterpolator<real_t, Vector2<real_t>, false>;
            Interpolator_T interpolator;

            uint_t nPoints{0};
            real_t aMin{0.};
            real_t aMax{0.};
        };

    } // namespace aerodynamics

} // namespace turbine_core

#endif //TURBINECORE_AIRFOILPOLAR_H
