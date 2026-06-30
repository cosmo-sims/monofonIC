# CMake configuration for monofonIC tests

include(cmake/TestSetup.cmake)
include(cmake/RegressionTests.cmake)
include(cmake/PlaneWaveTests.cmake)
include(cmake/MPIConsistencyTests.cmake)
include(cmake/TransferFunctionTests.cmake)

if(ENABLE_CLASS)
    include(cmake/BackscalingTests.cmake)
endif()
