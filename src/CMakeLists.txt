add_subdirectory(Dict)
add_executable(AnasenSim)

target_include_directories(AnasenSim
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    SYSTEM PRIVATE ${ROOT_INCLUDE_DIRS}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../vendor/catima/
)

target_sources(AnasenSim PRIVATE
    Sim/SimBase.h
    Sim/Application.h
    Sim/Application.cpp
    Sim/MassLookup.h
    Sim/MassLookup.cpp
    Sim/Target.h
    Sim/Target.cpp
    Sim/RxnType.h
    Sim/Reaction.h
    Sim/Reaction.cpp
    Sim/ReactionSystem.h
    Sim/ReactionSystem.cpp
    Sim/DecaySystem.h
    Sim/DecaySystem.cpp
    Sim/OneStepSystem.h
    Sim/OneStepSystem.cpp
    Sim/TwoStepSystem.h
    Sim/TwoStepSystem.cpp
)

target_link_libraries(AnasenSim PRIVATE catima SimDict ${ROOT_LIBS})
set_target_properties(AnasenSim PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${ASIM_BINARY_DIR})