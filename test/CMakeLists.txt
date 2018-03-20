add_executable(test_reg test_registration.cpp ${ANTs_SOURCE_DIR}/Examples/antsRegistration.cxx )
target_link_libraries(test_reg 
      ${ITK_LIBRARIES}  
      ${Boost_LIBRARIES}
      ${ANTS_LIBS}
      glog 
      nlohmann_json
    )

install(TARGETS test_reg DESTINATION bin)

add_test(NAME test_reg 
  COMMAND test_reg )

file(GLOB SRCS 
  test_main.cpp 
  dicom_tests.cpp
  resolute_tests.cpp
)

add_executable(testRESOLUTE ${SRCS})

target_link_libraries(testRESOLUTE
  ${ITK_LIBRARIES}  
  ${Boost_LIBRARIES}
  ${ANTS_LIBS}
  gtest
  gmock
  glog 
  nlohmann_json
)

install(TARGETS testRESOLUTE DESTINATION bin)

add_test(NAME testRESOLUTE 
  COMMAND testRESOLUTE )

