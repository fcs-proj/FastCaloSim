# Function to add a google test to the project
function(create_test_suite suite_name)
  # Glob the source file 
  # Create the executable from the provided source file
  add_executable(${suite_name}_tests source/${suite_name}_tests.cxx)
  # Link the necessary libraries
  target_link_libraries(${suite_name}_tests PRIVATE FastCaloSim::FastCaloSim GTest::gtest_main)
  # Enable GTest discovery
  gtest_discover_tests(${suite_name}_tests)
  # Deactivate compiler warnings and static analyzers
  # Note: this is necessary until we have fixed all FastCaloSim related warnings
  deactivate_checks(${suite_name}_tests)
endfunction()

# Function to create multiple tests, e.g. create_tests(test1 test2 test3)
function(create_test_suites)
  foreach(suite IN LISTS ARGN)
    create_test_suite(${suite})
  endforeach()
endfunction()