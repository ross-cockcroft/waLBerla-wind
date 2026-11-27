//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file SnapshotHandler.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_SNAPSHOTHANDLER_H
#define TURBINECORE_SNAPSHOTHANDLER_H

#pragma once

#include <string>
#include <vector>

#include <blockforest/SetupBlockForest.h>
#include <blockforest/StructuredBlockForest.h>
#include <core/config/Config.h>
#include <core/mpi/Broadcast.h>

#include <vtk/VTKOutput.h>

#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace domain {

        struct SnapshotHandler {

            explicit SnapshotHandler(std::shared_ptr<walberla::Config> & globalConfig,
                                     const std::vector<std::string> & outputBlockNames = {"Output"}) {

                const auto config = globalConfig->getOneBlock("DomainSetup");

                // load parameters
                loadSnapshot           = config.getParameter< bool >("loadSnapshot", false);
                storeSnapshot          = config.getParameter< bool >("storeSnapshot", false);
                snapshotFrequency      = config.getParameter< uint_t >("snapshotFrequency", uint_t(1));
                if(storeSnapshot) {
                    WALBERLA_CHECK(snapshotFrequency > uint_t(0), "Snapshot saving frequency \"snapshotFrequency\" must be greater than 0, "
                                                                  "please check your input file.")
                }

                baseFolder             = config.getParameter< std::string >("snapshotBaseFolder", "snapshot");
                configFile             = config.getParameter< std::string >("snapshotConfigFile", "");
                welfordFile            = config.getParameter< std::string >("welfordFile", "");
                blockforestFile        = config.getParameter< std::string >("snapshotBlockforestFile", "");
                pdfFile                = config.getParameter< std::string >("snapshotPdfFile", "");
                forceFile              = config.getParameter< std::string >("snapshotForceFile", "");
                meanVelFile            = config.getParameter< std::string >("snapshotMeanVelFile", "");
                sosFile                = config.getParameter< std::string >("snapshotSOSFile", "");
                socFile                = config.getParameter< std::string >("snapshotSOCFile", "");

                // early out
                if(!loadSnapshot && !storeSnapshot)
                    return;

                if(welfordFile.empty() && !configFile.empty()) {
                    welfordFile = configFile + ".welford";
                    WALBERLA_LOG_INFO_ON_ROOT("No \"welfordFile\" given, using " << welfordFile << " in the snapshot folder for the Welford counter.")
                }

                if(storeSnapshot) {
                    if(pdfFile.empty()) WALBERLA_ABORT("No pdf filename given for snapshot creation.")
                    if(forceFile.empty()) WALBERLA_ABORT("No force filename given for snapshot creation.")
                    if(configFile.empty()) WALBERLA_ABORT("No config filename given for snapshot creation.")

                    // create directories if they do not exist already
                    WALBERLA_ROOT_SECTION() {
                        std::filesystem::path snapshotPath(baseFolder);
                        if (!std::filesystem::exists(snapshotPath)) {
                            std::filesystem::create_directory(snapshotPath);
                        }
                    }
                }

                if(loadSnapshot) {
                    if(configFile.empty()) {
                        WALBERLA_ABORT("No config filename given, \"snapshotConfigFile\" is needed to restore the time step.")
                    }
                    // read snapshot config file
                    WALBERLA_ROOT_SECTION()
                    {
                        auto filePath = std::filesystem::path(baseFolder) / std::filesystem::path(configFile);
                        {
                            std::ifstream file(filePath);
                            if (!file) {
                                WALBERLA_ABORT("Error: " << configFile << " could not be opened!")
                            } else {
                                file >> startingTimestep;
                                if (file.fail()) {
                                    WALBERLA_ABORT("Error: no time step could be read from " << configFile << ".")
                                }
                            }
                        }

                    }
                    walberla::mpi::broadcastObject(startingTimestep);

                    WALBERLA_LOG_INFO_ON_ROOT("Successfully read config parameters from snapshot config file:")
                    WALBERLA_LOG_INFO_ON_ROOT(" - startingTimestep = " << startingTimestep)

                    for(const auto & outputBlockName : outputBlockNames) {
                        setInitialExecutionCount(globalConfig, outputBlockName);
                    }

                }

            }

            template< typename Timeloop_T >
            void setTimeloopStep(Timeloop_T & timeloop) const {
                if(loadSnapshot)
                    timeloop.setCurrentTimeStep(startingTimestep + 1);
            }

            uint_t readWelfordCounter() {
                if(welfordFile.empty()) {
                    WALBERLA_ABORT("No welford filename given, \"welfordFile\" is needed to restore the Welford counter.")
                }
                WALBERLA_ROOT_SECTION()
                {
                    std::ifstream file(std::filesystem::path(baseFolder) / std::filesystem::path(welfordFile));
                    if (!file) {
                        WALBERLA_ABORT("Error: " << welfordFile << " could not be opened!")
                    }
                    file >> welfordCounter;
                    if (file.fail()) {
                        WALBERLA_ABORT("Error: no Welford counter could be read from " << welfordFile << ".")
                    }
                }
                walberla::mpi::broadcastObject(welfordCounter);
                WALBERLA_LOG_INFO_ON_ROOT("Read welford counter from snapshot welford file:")
                WALBERLA_LOG_INFO_ON_ROOT(" - welfordCounter = " << welfordCounter)
                return welfordCounter;
            }

            void createSnapshot(const uint_t t, const std::shared_ptr<walberla::StructuredBlockForest> & blocks,
                                const walberla::BlockDataID & pdfFieldID, const walberla::BlockDataID & forceFieldID) {
                if(storeSnapshot) {
                    if (t % snapshotFrequency == 0 && t > 0) {
                        const auto pdfFileName = std::string(std::filesystem::path(baseFolder) / std::filesystem::path(pdfFile));
                        const auto forceFileName = std::string(std::filesystem::path(baseFolder) / std::filesystem::path(forceFile));

                        blocks->saveBlockData(pdfFileName, pdfFieldID);
                        blocks->saveBlockData(forceFileName, forceFieldID);

                        WALBERLA_ROOT_SECTION()
                        {
                            const auto checkpointConfigFile = std::string(std::filesystem::path(baseFolder) / std::filesystem::path(configFile));
                            std::ofstream file;
                            file.open(checkpointConfigFile.c_str());

                            if(!file.is_open()) {
                                WALBERLA_ABORT("Error: " << configFile << " could not be opened for writing!")
                            }
                            file << std::setprecision(16);
                            file << t << "\n";
                            file.close();
                        }
                    }

                }
            }

            void createSnapshot(const uint_t t, const std::shared_ptr<walberla::StructuredBlockForest> & blocks,
                                const walberla::BlockDataID & pdfFieldID, const walberla::BlockDataID & forceFieldID,
                                const uint_t _welfordCounter,
                                const walberla::BlockDataID & meanVelocityFieldID ) {
                createSnapshot(t, blocks, pdfFieldID, forceFieldID);

                if(storeSnapshot) {
                    if(meanVelFile.empty()) WALBERLA_ABORT("No mean velocity filename given for snapshot creation.")

                    if (t % snapshotFrequency == 0 && t > 0) {
                        const auto meanVelocityFileName = std::string(std::filesystem::path(baseFolder) / std::filesystem::path(meanVelFile));

                        blocks->saveBlockData(meanVelocityFileName, meanVelocityFieldID);

                        WALBERLA_ROOT_SECTION()
                        {
                            std::ofstream file;

                            const auto welfordConfigFile = std::string(std::filesystem::path(baseFolder) / std::filesystem::path(welfordFile));
                            file.open(welfordConfigFile.c_str());

                            if(!file.is_open()) {
                                WALBERLA_ABORT("Error: " << welfordFile << " could not be opened for writing!")
                            }
                            file << std::setprecision(16);
                            file << _welfordCounter << "\n";
                            file.close();
                        }
                    }

                }
            }

            void createSnapshot(const uint_t t, const std::shared_ptr<walberla::StructuredBlockForest> & blocks,
                                const walberla::BlockDataID & pdfFieldID, const walberla::BlockDataID & forceFieldID,
                                const uint_t _welfordCounter,
                                const walberla::BlockDataID & meanVelocityFieldID,
                                const walberla::BlockDataID & sosFieldID) {
                createSnapshot(t, blocks, pdfFieldID, forceFieldID, _welfordCounter, meanVelocityFieldID);

                if(storeSnapshot) {
                    if(sosFile.empty()) WALBERLA_ABORT("No SoS filename given for snapshot creation.")

                    if (t % snapshotFrequency == 0 && t > 0) {
                        const auto sosFileName = std::string(std::filesystem::path(baseFolder) / std::filesystem::path(sosFile));

                        blocks->saveBlockData(sosFileName, sosFieldID);
                    }
                }
            }

            void createSnapshot(const uint_t t, const std::shared_ptr<walberla::StructuredBlockForest> & blocks,
                                const walberla::BlockDataID & pdfFieldID, const walberla::BlockDataID & forceFieldID,
                                const uint_t _welfordCounter,
                                const walberla::BlockDataID & meanVelocityFieldID,
                                const walberla::BlockDataID & sosFieldID,
                                const walberla::BlockDataID & socFieldID) {
                createSnapshot(t, blocks, pdfFieldID, forceFieldID, _welfordCounter, meanVelocityFieldID, sosFieldID);

                if(storeSnapshot) {
                    if(socFile.empty()) WALBERLA_ABORT("No SoC filename given for snapshot creation.")

                    if (t % snapshotFrequency == 0 && t > 0) {
                        const auto socFileName = std::string(std::filesystem::path(baseFolder) / std::filesystem::path(socFile));

                        blocks->saveBlockData(socFileName, socFieldID);
                    }
                }
            }

            bool loadSnapshot;
            bool storeSnapshot;

            uint_t snapshotFrequency;
            uint_t startingTimestep{0};
            uint_t welfordCounter{0};

            std::string baseFolder;
            std::string configFile;
            std::string blockforestFile;
            std::string pdfFile;
            std::string forceFile;
            std::string meanVelFile;
            std::string sosFile;
            std::string socFile;
            std::string welfordFile;

        private:

            void setInitialExecutionCount(std::shared_ptr<walberla::Config> & globalConfig, const std::string & outputBlockName) const {
                std::vector< walberla::config::Config::Block* > configOutputBlock;
                globalConfig->getWritableGlobalBlock().getWritableBlocks(outputBlockName, configOutputBlock, 0, 1);
                for(auto * outputBlock : configOutputBlock) {
                    std::vector< walberla::config::Config::Block* > configVTKBlock;
                    outputBlock->getWritableBlocks("VTK", configVTKBlock);
                    for(auto * vtkBlock : configVTKBlock) {
                        const uint_t skipTimestep = vtkBlock->getParameter< uint_t >("skipTimestep", uint_t(0));
                        const uint_t executionCount = (startingTimestep + 1 > skipTimestep) ? startingTimestep + 1 - skipTimestep : uint_t(0);
                        vtkBlock->setOrAddParameter("initialExecutionCount", std::to_string(executionCount));
                    }
                }
            }
        };

    } // domain

} // turbine_core 

#endif // TURBINECORE_SNAPSHOTHANDLER_H
