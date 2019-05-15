cmake_minimum_required(VERSION 2.6)
project(squarefull)

add_executable(squarefull main.cpp sqf.cpp sqfmarmet.cpp primetester.cpp sqfsorter.cpp factor.cpp checkresult.cpp)

install(TARGETS squarefull RUNTIME DESTINATION bin)

target_link_libraries(squarefull libntl.so)
