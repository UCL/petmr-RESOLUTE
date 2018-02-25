add_executable(resolute Resolute.cpp 
  ${ANTs_SOURCE_DIR}/Examples/antsRegistration.cxx
  ${ANTs_SOURCE_DIR}/Examples/antsApplyTransforms.cxx
  )
target_link_libraries(resolute 
      ${ITK_LIBRARIES}  
      ${Boost_LIBRARIES}
      ${ANTS_LIBS}
      glog 
      nlohmann_json
    )

install(TARGETS resolute DESTINATION bin)

add_executable(test_mask TestMaskExtract.cpp  )
target_link_libraries(test_mask 
      ${ITK_LIBRARIES}  
      ${Boost_LIBRARIES}
      glog 
      nlohmann_json
    )

install(TARGETS test_mask DESTINATION bin)

add_executable(test_reg TestRegistration.cpp ${ANTs_SOURCE_DIR}/Examples/antsRegistration.cxx )
target_link_libraries(test_reg 
      ${ITK_LIBRARIES}  
      ${Boost_LIBRARIES}
      ${ANTS_LIBS}
      glog 
      nlohmann_json
    )

install(TARGETS test_reg DESTINATION bin)