cmake_minimum_required(VERSION 3.10)
project(inferAdapter)
set(EXECUTABLE_OUTPUT_PATH ./bin)
include_directories(include 
                    include/tools)
aux_source_directory(./ SOURCE)
add_executable(inferAdapter ${SOURCE})
target_link_libraries(inferAdapter hts)
