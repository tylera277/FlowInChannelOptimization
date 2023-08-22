add_executable(test_simulation )#test_simulation.cpp)

target_sources(test_simulation PRIVATE 
test_simulation.cpp
../../src/simulation/simulation.cpp
../../src/simulation/simulation.hpp 
../catch.hpp
)

find_package (Eigen3 3.3 REQUIRED)

target_link_libraries(test_simulation PUBLIC Eigen3::Eigen)
