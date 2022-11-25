add_library(SimDict SHARED)
target_include_directories(SimDict
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    SYSTEM PUBLIC ${ROOT_INCLUDE_DIRS}
)

ROOT_GENERATE_DICTIONARY(sim_dict Nucleus.h LINKDEF LinkDef_SimDict.h MODULE SimDict)

target_sources(SimDict PRIVATE Nucleus.h Nucleus.cpp)

target_link_libraries(SimDict PUBLIC ${ROOT_LIBRARIES})

set_target_properties(SimDict PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${ASIM_LIBRARY_DIR})
add_custom_command(TARGET SimDict POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_BINARY_DIR}/libSimDict_rdict.pcm
    ${ASIM_LIBRARY_DIR}/libSimDict_rdict.pcm
)

