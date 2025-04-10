# Set the minimum required CMake version
cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON) # Ensure C++17 is enforced

include(FetchContent)

FetchContent_Declare(
  tensorstore
  URL "https://github.com/google/tensorstore/archive/refs/tags/v0.1.71.tar.gz"
)

FetchContent_MakeAvailable(tensorstore)

# Define the project name
project(nTracer_cdn)

# Add an executable target
add_executable(nTracer_cdn src/app.cpp)

target_link_libraries(nTracer_cdn pthread)
target_link_libraries(nTracer_cdn sqlite3)

target_link_libraries(nTracer_cdn avcodec)
target_link_libraries(nTracer_cdn avformat)
target_link_libraries(nTracer_cdn avutil)
target_link_libraries(nTracer_cdn swscale)

find_library(ZSTD_LIB zstd PATH ./zstd/lib/)
target_link_libraries(nTracer_cdn "${ZSTD_LIB}")

find_library(X264_LIB x264 PATH ./x264/)
target_link_libraries(nTracer_cdn "${X264_LIB}")

target_compile_options(nTracer_cdn PUBLIC -O3 -march=native -Wl,--no-as-needed -ldl)

# Define a target that depends on TensorStore...
target_link_libraries(
  nTracer_cdn
#  PRIVATE
    tensorstore::tensorstore tensorstore::all_drivers
)