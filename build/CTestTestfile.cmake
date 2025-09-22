# CMake generated Testfile for 
# Source directory: /home/kevin/meshcut
# Build directory: /home/kevin/meshcut/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(basic_test "/home/kevin/meshcut/build/meshcut_test")
set_tests_properties(basic_test PROPERTIES  _BACKTRACE_TRIPLES "/home/kevin/meshcut/CMakeLists.txt;59;add_test;/home/kevin/meshcut/CMakeLists.txt;0;")
add_test(advanced_test "/home/kevin/meshcut/build/meshcut_advanced_test")
set_tests_properties(advanced_test PROPERTIES  _BACKTRACE_TRIPLES "/home/kevin/meshcut/CMakeLists.txt;60;add_test;/home/kevin/meshcut/CMakeLists.txt;0;")
add_test(visualization_test "/home/kevin/meshcut/build/meshcut_visualization_test")
set_tests_properties(visualization_test PROPERTIES  _BACKTRACE_TRIPLES "/home/kevin/meshcut/CMakeLists.txt;61;add_test;/home/kevin/meshcut/CMakeLists.txt;0;")
add_test(performance_test "/home/kevin/meshcut/build/meshcut_performance_test")
set_tests_properties(performance_test PROPERTIES  _BACKTRACE_TRIPLES "/home/kevin/meshcut/CMakeLists.txt;62;add_test;/home/kevin/meshcut/CMakeLists.txt;0;")
