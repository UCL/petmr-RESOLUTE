add_executable(resolute Resolute.cpp  )
target_link_libraries(resolute 
      ${ITK_LIBRARIES}  
      ${Boost_LIBRARIES}
      glog 
      nlohmann_json
    )

install(TARGETS resolute DESTINATION bin)