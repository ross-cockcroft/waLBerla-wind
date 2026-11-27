#include <blockforest/all.h>

#include <core/debug/TestSubsystem.h>
#include <core/mpi/Environment.h>
#include <core/math/Vector3.h>

#include <field/AddToStorage.h>
#include <field/GhostLayerField.h>

#include <cmath>

#include "walberla_helper/blockforest/BlockInfo.h"
#include "walberla_helper/field/Field.h"
#include "walberla_helper/field/interpolators/KernelFieldInterpolator.h"
#include "walberla_helper/field/interpolators/TrilinearFieldInterpolator.h"
#include "walberla_helper/field/projector_kernels/SmoothedDiracDeltaKernel.h"

namespace turbine_core {

    namespace actuator_sampling_tests {

        using namespace walberla;

        typedef GhostLayerField<real_t, 1> ScalarField_T;

        real_t linearField(real_t x, real_t y, real_t z) {
            return real_t(2) + real_t(0.3) * x - real_t(0.2) * y + real_t(0.5) * z;
        }

        template<typename Function_T>
        real_t sampleAtPoints(const shared_ptr<StructuredBlockForest> & blocks, const BlockDataID & fieldID, Function_T exactField, const char * name) {

            const walberla::Vector3<real_t> points[] = {
                walberla::Vector3<real_t>(real_t(3.5), real_t(3.5), real_t(3.5)),
                walberla::Vector3<real_t>(real_t(3.9), real_t(4.2), real_t(2.8)),
                walberla::Vector3<real_t>(real_t(4.25), real_t(3.1), real_t(5.7)),
                walberla::Vector3<real_t>(real_t(1.6), real_t(6.4), real_t(4.9))
            };

            real_t maxTrilinearError{0};
            real_t maxKernelError{0};

            for (auto &p : points) {
                auto block = blocks->getBlock(p);
                WALBERLA_CHECK_NOT_NULLPTR(block);
                blockforest::BlockInfo blockInfo(block, blocks);
                field::Field<real_t> field(block->getData<ScalarField_T>(fieldID));

                real_t exact = exactField(p[0], p[1], p[2]);

                real_t trilinear{0};
                projectors::TrilinearFieldInterpolator<field::Field<real_t>> trilinearInterpolator(blockInfo, &field);
                trilinearInterpolator.get(p, &trilinear);

                real_t kernel{0};
                projectors::KernelFieldInterpolator<field::Field<real_t>, projectors::SmoothedDiracDeltaKernel> kernelInterpolator(blockInfo, &field);
                kernelInterpolator.get(p, &kernel);

                WALBERLA_LOG_RESULT(name << " point (" << p[0] << ", " << p[1] << ", " << p[2] << ") exact " << exact
                                    << " trilinear " << trilinear << " kernel " << kernel
                                    << " trilinear error " << trilinear - exact << " kernel error " << kernel - exact);

                maxTrilinearError = std::max(maxTrilinearError, std::abs(trilinear - exact));
                maxKernelError = std::max(maxKernelError, std::abs(kernel - exact));
            }

            WALBERLA_LOG_RESULT(name << " max trilinear error " << maxTrilinearError << " max kernel error " << maxKernelError);
            return maxTrilinearError;
        }

        void fillField(const shared_ptr<StructuredBlockForest> & blocks, const BlockDataID & fieldID, real_t (*function)(real_t, real_t, real_t)) {
            for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
                auto f = blockIt->getData<ScalarField_T>(fieldID);
                WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(f,
                    walberla::Vector3<real_t> c = blocks->getBlockLocalCellCenter(*blockIt, Cell(x, y, z));
                    f->get(x, y, z) = function(c[0], c[1], c[2]);
                );
            }
        }

        int main(int argc, char **argv) {

            mpi::Environment mpiEnv(argc, argv);
            debug::enterTestMode();

            auto blocks = walberla::blockforest::createUniformBlockGrid(1, 1, 1, 8, 8, 8, real_t(1), 0, false, false,
                                                                        false, false, false, false);

            BlockDataID fieldID = walberla::field::addToStorage<ScalarField_T>(blocks, "scalar field", real_t(0), walberla::field::fzyx, uint_t(1));

            fillField(blocks, fieldID, linearField);
            real_t trilinearErrorLinear = sampleAtPoints(blocks, fieldID, linearField, "linear");

            WALBERLA_CHECK_LESS(trilinearErrorLinear, real_t(1e-5));

            return 0;
        }

    }
}

int main(int argc, char **argv) {
    return turbine_core::actuator_sampling_tests::main(argc, argv);
}
