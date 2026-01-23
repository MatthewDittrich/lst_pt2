cmake_minimum_required(VERSION 3.10)
project(LSTAnalysis)

# 1. Find ROOT Package
find_package(ROOT REQUIRED COMPONENTS RIO Tree Hist Matrix MathCore Minuit2 Gpad Graf)

# 2. Setup ROOT flags and include directories
include(${ROOT_USE_FILE})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${ROOT_INCLUDE_DIRS})

# 3. Define Source Files
set(SOURCES 
    main.cpp 
    LSTReader.cpp 
    PhysicsAlgo.cpp 
    HistManager.cpp 
    Plotter.cpp
    Gator.cpp
)

# 4. Create the Executable
add_executable(LSTAnalysis ${SOURCES})

# 5. Link Libraries (ROOT + Minuit2)
target_link_libraries(LSTAnalysis ${ROOT_LIBRARIES} Minuit2)
