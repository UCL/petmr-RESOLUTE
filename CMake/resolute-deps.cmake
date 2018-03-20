find_package(Boost 1.54.0 MODULE
  COMPONENTS
    chrono
    date_time
    filesystem
    program_options
    regex
    system
    thread
  REQUIRED
)

include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(glog REQUIRED)

if(BUILD_TESTING)
  find_package(GTest REQUIRED)
endif()

find_package(nlohmann_json)

#ANTS
# find ANTS includes

set(ANTs_SOURCE_DIR "" CACHE PATH "ANTs source location")
message("ANTs_SOURCE_DIR = ${ANTs_SOURCE_DIR}")

set(ANTs_LIBRARY_DIR "" CACHE PATH "ANTs library location")
message("ANTs_LIBRARY_DIR = ${ANTs_LIBRARY_DIR}")

include_directories(${ANTs_SOURCE_DIR}/Temporary)
include_directories(${ANTs_SOURCE_DIR}/Tensor)
include_directories(${ANTs_SOURCE_DIR}/Utilities)
include_directories(${ANTs_SOURCE_DIR}/Examples)
include_directories(${ANTs_SOURCE_DIR}/Examples/include)
include_directories(${ANTs_SOURCE_DIR}/ImageRegistration)

link_directories(${ANTs_LIBRARY_DIR})
set(ANTS_LIBS antsUtilities l_ANTS)
