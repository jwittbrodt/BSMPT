add_library(BounceCalculation ModelWrapper.cpp)
target_link_libraries(BounceCalculation PUBLIC Models
                                               Vevacious::libVevaciousPlusPlus)

add_executable(bounceTest test.cpp)
target_link_libraries(bounceTest PRIVATE BounceCalculation)
