add_executable(Plotter)

target_include_directories(Plotter PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../
    SYSTEM PUBLIC ${ROOT_INCLUDE_DIRS}
)

target_sources(Plotter PRIVATE Plotter.h Plotter.cpp main.cpp)

target_link_libraries(Plotter PRIVATE SimDict ${ROOT_LIBRARIES})

set_target_properties(Plotter PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${ASIM_BINARY_DIR})