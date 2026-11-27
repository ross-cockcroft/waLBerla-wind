waLBerla_add_test_executable( AirfoilPolarBlendTest AirfoilPolarBlendTest.cpp )
target_link_libraries( AirfoilPolarBlendTest walberla::core )
target_compile_definitions( AirfoilPolarBlendTest PRIVATE AIRFOIL_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../../../apps/GeneratedApp_cpu_uniform/airfoil_data" )
waLBerla_execute_test( NAME AirfoilPolarBlendTest )
