file(GLOB TEST_TRANSFER_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)
list(FILTER TEST_TRANSFER_SOURCES EXCLUDE REGEX ".*main\\.cc$")
file(GLOB TEST_TRANSFER_PLUGINS
    ${PROJECT_SOURCE_DIR}/src/plugins/transfer_*.cc
    ${PROJECT_SOURCE_DIR}/src/plugins/random_*.cc
)

add_executable(test_transfer_functions
    ${CMAKE_CURRENT_SOURCE_DIR}/test_transfer_functions.cc
    ${TEST_TRANSFER_SOURCES}
    ${TEST_TRANSFER_PLUGINS}
)

if(ENABLE_PANPHASIA AND ENABLE_MPI)
    target_sources(test_transfer_functions PRIVATE
        ${PROJECT_SOURCE_DIR}/external/panphasia/panphasia_routines.f
        ${PROJECT_SOURCE_DIR}/external/panphasia/generic_lecuyer.f90
        ${PROJECT_SOURCE_DIR}/external/panphasia_ho/high_order_panphasia_routines.c
        ${PROJECT_SOURCE_DIR}/external/panphasia_ho/pan_mpi_routines.c
        ${PROJECT_SOURCE_DIR}/external/panphasia_ho/uniform_rand_threefry4x64.c
    )
endif()

target_include_directories(test_transfer_functions PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/external
)
target_link_libraries(test_transfer_functions PRIVATE GSL::gsl)

if(MPI_CXX_FOUND)
    if(CODE_PRECISION STREQUAL "FLOAT")
        if(FFTW3_SINGLE_MPI_FOUND)
            target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_SINGLE_MPI)
            target_compile_definitions(test_transfer_functions PRIVATE "USE_FFTW_MPI")
        endif()
    elseif(CODE_PRECISION STREQUAL "DOUBLE")
        if(FFTW3_DOUBLE_MPI_FOUND)
            target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_DOUBLE_MPI)
            target_compile_definitions(test_transfer_functions PRIVATE "USE_FFTW_MPI")
        endif()
    elseif(CODE_PRECISION STREQUAL "LONGDOUBLE")
        if(FFTW3_LONGDOUBLE_MPI_FOUND)
            target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_LONGDOUBLE_MPI)
            target_compile_definitions(test_transfer_functions PRIVATE "USE_FFTW_MPI")
        endif()
    endif()
    target_link_libraries(test_transfer_functions PRIVATE MPI::MPI_CXX)
    target_compile_definitions(test_transfer_functions PRIVATE "USE_MPI")
endif()

if(CODE_PRECISION STREQUAL "FLOAT")
    if(FFTW3_SINGLE_THREADS_FOUND)
        target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_SINGLE_THREADS)
        target_compile_definitions(test_transfer_functions PRIVATE "USE_FFTW_THREADS")
    endif()
    target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_SINGLE_SERIAL)
elseif(CODE_PRECISION STREQUAL "DOUBLE")
    if(FFTW3_DOUBLE_THREADS_FOUND)
        target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_DOUBLE_THREADS)
        target_compile_definitions(test_transfer_functions PRIVATE "USE_FFTW_THREADS")
    endif()
    target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_DOUBLE_SERIAL)
elseif(CODE_PRECISION STREQUAL "LONGDOUBLE")
    if(FFTW3_LONGDOUBLE_THREADS_FOUND)
        target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_LONGDOUBLE_THREADS)
        target_compile_definitions(test_transfer_functions PRIVATE "USE_FFTW_THREADS")
    endif()
    target_link_libraries(test_transfer_functions PRIVATE FFTW3::FFTW3_LONGDOUBLE_SERIAL)
endif()

if(HDF5_FOUND)
    target_link_libraries(test_transfer_functions PRIVATE ${HDF5_LIBRARIES})
    target_include_directories(test_transfer_functions PRIVATE ${HDF5_INCLUDE_DIRS})
    target_compile_definitions(test_transfer_functions PRIVATE "USE_HDF5")
endif()
if(ENABLE_PANPHASIA)
    target_compile_definitions(test_transfer_functions PRIVATE "USE_PANPHASIA")
endif()
if(ENABLE_CLASS)
    target_link_libraries(test_transfer_functions PRIVATE class::libclass)
    target_compile_definitions(test_transfer_functions PRIVATE "USE_CLASS")
endif()

function(add_transfer_test PLUGIN_NAME CONFIG_NAME)
    set(TEST_NAME "test_transfer_${PLUGIN_NAME}_${CONFIG_NAME}")
    set(CONFIG_FILE "${TEST_CONFIG_DIR}/transfer/${CONFIG_NAME}.conf")
    set(REFERENCE_FILE "${TEST_REFERENCE_DIR}/transfer/${PLUGIN_NAME}/${CONFIG_NAME}.txt")
    add_test(
        NAME ${TEST_NAME}
        COMMAND $<TARGET_FILE:test_transfer_functions>
            --test ${PLUGIN_NAME} ${CONFIG_FILE} ${REFERENCE_FILE} 1e-6
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(${TEST_NAME} PROPERTIES
        TIMEOUT 120
        LABELS "transfer;regression"
    )
    if(NOT EXISTS ${REFERENCE_FILE})
        set_tests_properties(${TEST_NAME} PROPERTIES DISABLED TRUE)
        message(STATUS "Test ${TEST_NAME}: Reference file not found - test will be disabled")
    else()
        message(STATUS "Test ${TEST_NAME}: Registered")
    endif()
endfunction()

set(TF_CONFIGS fiducial dark_energy massive_nu low_omega_m high_z)
if(ENABLE_CLASS)
    foreach(CONFIG ${TF_CONFIGS})
        add_transfer_test("CLASS" ${CONFIG})
    endforeach()
endif()
foreach(CONFIG ${TF_CONFIGS})
    add_transfer_test("eisenstein" ${CONFIG})
endforeach()

message(STATUS "")
message(STATUS "Transfer function regression tests configured")
message(STATUS "  Generate references with: cd build && bash ../tests/scripts/generate_transfer_references.sh")
message(STATUS "  Run tests with: ctest -R test_transfer")
