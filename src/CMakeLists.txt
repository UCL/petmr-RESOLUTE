add_executable(resolute Resolute.cpp 
  ${ANTs_SOURCE_DIR}/Examples/antsRegistration.cxx
  ${ANTs_SOURCE_DIR}/Examples/antsApplyTransforms.cxx
  )
target_link_libraries(resolute 
      ${ANTS_LIBS}
      ${ITK_LIBRARIES}  
      ${Boost_LIBRARIES}
      glog 
      nlohmann_json
    )

install(TARGETS resolute DESTINATION bin)
