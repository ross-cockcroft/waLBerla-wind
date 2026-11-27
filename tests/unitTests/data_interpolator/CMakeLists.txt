########################
### DataInterpolator ###
########################

waLBerla_add_test_executable( DataInterpolatorCubicTest DataInterpolatorCubicTest.cpp )
target_link_libraries( DataInterpolatorCubicTest walberla::core )
waLBerla_execute_test( NAME DataInterpolatorCubicTest )
