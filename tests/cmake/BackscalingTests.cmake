add_executable(test_class_reference ${class_SOURCE_DIR}/main/class.c)
set_target_properties(test_class_reference PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(test_class_reference PRIVATE class::libclass)

function(add_backscaling_class_test TEST_NAME BACKSCALING_RADIATION BACKSCALING_NEUTRINOS M_NU1 M_NU2 M_NU3)
    set(TEST_CONFIG_DIR_LOCAL "${CMAKE_CURRENT_BINARY_DIR}/backscaling_configs")
    set(TEST_WORK_DIR "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}")
    file(MAKE_DIRECTORY "${TEST_CONFIG_DIR_LOCAL}")
    file(MAKE_DIRECTORY "${TEST_WORK_DIR}")
    set(CONFIG_FILE "${TEST_CONFIG_DIR_LOCAL}/${TEST_NAME}.conf")
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/configs/backscaling/base.inc"
        "${CONFIG_FILE}"
        @ONLY
    )
    add_test(
        NAME ${TEST_NAME}
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_backscaling_to_class.py"
            --monofonic $<TARGET_FILE:monofonIC>
            --class $<TARGET_FILE:test_class_reference>
            --config "${CONFIG_FILE}"
            --workdir "${TEST_WORK_DIR}"
        WORKING_DIRECTORY "${TEST_WORK_DIR}"
    )
    set_tests_properties(${TEST_NAME} PROPERTIES
        TIMEOUT 120
        LABELS "backscaling;class;regression"
    )
endfunction()

add_backscaling_class_test(test_backscaling_r1_n1_one true true 0.06 0 0)
add_backscaling_class_test(test_backscaling_r1_n0_one true false 0.06 0 0)
add_backscaling_class_test(test_backscaling_r0_n1_one false true 0.06 0 0)
add_backscaling_class_test(test_backscaling_r0_n0_one false false 0.06 0 0)
add_backscaling_class_test(test_backscaling_zero_mass true true 0 0 0)
add_backscaling_class_test(test_backscaling_three_mass true true 0.02 0.02 0.02)

set(TARGET_CLASS_INPUT_WORK_DIR "${CMAKE_CURRENT_BINARY_DIR}/test_target_class_input_invariance")
file(MAKE_DIRECTORY "${TARGET_CLASS_INPUT_WORK_DIR}")
add_test(
    NAME test_target_class_input_invariance
    COMMAND ${Python3_EXECUTABLE}
        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_target_class_inputs.py"
        --monofonic $<TARGET_FILE:monofonIC>
        --workdir "${TARGET_CLASS_INPUT_WORK_DIR}"
        "${CMAKE_CURRENT_BINARY_DIR}/backscaling_configs/test_backscaling_r1_n1_one.conf"
        "${CMAKE_CURRENT_BINARY_DIR}/backscaling_configs/test_backscaling_r1_n0_one.conf"
        "${CMAKE_CURRENT_BINARY_DIR}/backscaling_configs/test_backscaling_r0_n1_one.conf"
        "${CMAKE_CURRENT_BINARY_DIR}/backscaling_configs/test_backscaling_r0_n0_one.conf"
)
set_tests_properties(test_target_class_input_invariance PROPERTIES
    TIMEOUT 300
    LABELS "backscaling;class;regression"
    RUN_SERIAL TRUE
)
