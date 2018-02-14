add_executable(resolute Resolute.cpp  )
target_link_libraries(resolute 
      ${ITK_LIBRARIES}  
      ${Boost_LIBRARIES}
      glog 
      nlohmann_json
    )

add_executable(test_mask TestMaskExtract.cpp  )
target_link_libraries(test_mask 
      ${ITK_LIBRARIES}  
      ${Boost_LIBRARIES}
      glog 
      nlohmann_json
    )

install(TARGETS test_mask DESTINATION bin)