file(GLOB SRCS *.cpp)

add_executable(testRESOLUTE ${SRCS})

target_link_libraries(testRESOLUTE
  ${ITK_LIBRARIES}  
  ${Boost_LIBRARIES}
  ${ANTS_LIBS}
  gtest
  gmock
)

install(TARGETS testRESOLUTE DESTINATION bin)

add_test(NAME testRESOLUTE 
  COMMAND testRESOLUTE )
