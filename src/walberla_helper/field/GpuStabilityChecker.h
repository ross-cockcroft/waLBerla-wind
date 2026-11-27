#pragma once

#ifndef TURBINECORE_GPU_STABILITY_CHECKER_H
#define TURBINECORE_GPU_STABILITY_CHECKER_H

#include <functional>
#include <string>

#include <core/config/Config.h>
#include <field/StabilityChecker.h>
#include <gpu/FieldCopy.h>
#include <gpu/GPUWrapper.h>

#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace gpu {

        template<typename ScalarField_T, typename GPUField_T, typename FlagField_T, typename BlockStorage_T>
        std::function<void()> makeDensityStabilityChecker(const walberla::shared_ptr<walberla::Config> &config,
                                                          const walberla::shared_ptr<BlockStorage_T> &blocks,
                                                          const walberla::BlockDataID &densityCpuID,
                                                          const walberla::BlockDataID &densityGpuID,
                                                          const walberla::ConstBlockDataID &flagFieldID,
                                                          const walberla::Set<walberla::FlagUID> &fluidFlags,
                                                          const std::string &domainName) {

            const auto block = config->getGlobalBlock().getBlock("StabilityChecker");
            const walberla::uint_t frequency = block ? block.template getParameter<walberla::uint_t>("checkFrequency", 0) : walberla::uint_t(0);
            const bool streamOutput = block ? block.template getParameter<bool>("streamOutput", true) : true;
            const bool vtkOutput = block ? block.template getParameter<bool>("vtkOutput", true) : true;

            auto checker = walberla::field::makeStabilityChecker<ScalarField_T, FlagField_T>(blocks, densityCpuID, flagFieldID, fluidFlags,
                                                                                            walberla::uint_t(1), streamOutput, vtkOutput);
            if (block) {
                if (block.isDefined("vtkBaseFolder")) { checker->setVTKBaseFolder(block.template getParameter<std::string>("vtkBaseFolder")); }
                if (block.isDefined("vtkExecutionFolder")) { checker->setVTKExecutionFolder(block.template getParameter<std::string>("vtkExecutionFolder")); }
                if (block.isDefined("vtkBinary")) { checker->setVTKBinary(block.template getParameter<bool>("vtkBinary")); }
                if (block.isDefined("vtkLittleEndian")) { checker->setVTKLittleEndian(block.template getParameter<bool>("vtkLittleEndian")); }
                if (block.isDefined("vtkMPIIO")) { checker->setVTKMPIIO(block.template getParameter<bool>("vtkMPIIO")); }
                if (block.isDefined("vtkForcePVTU")) { checker->setVTKForcePVTU(block.template getParameter<bool>("vtkForcePVTU")); }
            }
            const std::string identifier = block ? block.template getParameter<std::string>("vtkIdentifier", "error_field") : std::string("error_field");
            checker->setVTKIdentifier(domainName.empty() ? identifier : domainName + "_" + identifier);

            return [=, counter = walberla::uint_t(0)]() mutable {
                if (frequency == 0 || counter++ % frequency != 0) {
                    return;
                }
                gpuDeviceSynchronize();
                walberla::gpu::fieldCpy<ScalarField_T, GPUField_T>(blocks, densityCpuID, densityGpuID);
                (*checker)();
            };
        }

    } // namespace gpu

} // namespace turbine_core

#endif //TURBINECORE_GPU_STABILITY_CHECKER_H
