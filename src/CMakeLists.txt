set(CMAKE_POSITION_INDEPENDENT_CODE ON)

SET(lib_name rheocell_common)

IF(CMAKE_BUILD_TYPE MATCHES Debug)
	SET(exe_name rheoCell_debug)
ELSEIF(MPI)
	SET(exe_name rheoCell_mpi)
ELSE()
	SET(exe_name rheoCell)
ENDIF()

SET(observables_SOURCES
	Observables/BaseObservable.cpp
	Observables/ObservableFactory.cpp
	Observables/ObservableOutput.cpp
	Observables/Configuration.cpp
	Observables/VelocityField.cpp
	Observables/ForceField.cpp
	Observables/NematicField.cpp
	Observables/StressField.cpp
	Observables/ShapeField.cpp
)

SET(backends_SOURCES
	Backends/BackendFactory.cpp
	Backends/SimBackend.cpp
	Backends/FDBackend.cpp
	Backends/FD_CPUBackend.cpp 
	Backends/FD_CHBackend.cpp 
)

SET(interactions_SOURCES
	Interactions/InteractionFactory.cpp
	Interactions/BaseInteraction.cpp
	Interactions/SimpleMultiField.cpp
	Interactions/ActiveMultiField.cpp
	Interactions/ActiveNematic.cpp
	Interactions/LEBcActiveNematic.cpp
	Interactions/WetModel.cpp
	Interactions/LEBcWetModel.cpp
	Interactions/GeneralizedWetModel.cpp
	Interactions/DifferentialAdhesion.cpp
)

SET(forces_SOURCES
	Forces/ForceFactory.cpp
	Forces/BaseForce.cpp
	Forces/ChannelWalls.cpp
	Forces/ShearFlowChannel.cpp
)

SET(box_SOURCES
	Boxes/BaseBox.cpp
	Boxes/SquareBox.cpp
	Boxes/Channel.cpp
	Boxes/MovingWalls.cpp
	Boxes/SquareWalls.cpp
	Boxes/OrthogonalBox.cpp
	Boxes/LeesEdwardsBox.cpp
	Boxes/BoxFactory.cpp
)

SET(common_SOURCES
	Fields/BaseField.cpp
	Fields/MultiPhaseField.cpp
	Fields/LEBcMultiPhaseField.cpp
	Managers/SimManager.cpp
	Utilities/LR_vector.cpp
	Utilities/Utils.cpp
	Utilities/RCexception.cpp
	Utilities/Logger.cpp
	Utilities/parse_input/parse_input.cpp
	Utilities/time_scales/time_scales.cpp
	Utilities/SignalManager.cpp
	Utilities/ConfigInfo.cpp
	Utilities/FlattenedConfigInfo.cpp
	${forces_SOURCES}
	${observables_SOURCES}
	${lists_SOURCES}
	${backends_SOURCES}
	${interactions_SOURCES}
	${move_SOURCES}
	${box_SOURCES}
)


SET(confGenerator_SOURCES
	confGenerator.cpp
	Managers/GeneratorManager.cpp
)

SET(analyser_SOURCES
	analyser.cpp
)

SET(rheoCell_SOURCES
	rheoCell.cpp 
)
	
ADD_DEFINITIONS(-DNOCUDA)
SET(common_SOURCES
	${common_SOURCES}
	Utilities/Timings.cpp
)

ADD_EXECUTABLE(${exe_name} ${rheoCell_SOURCES})
ADD_LIBRARY(${lib_name} ${common_SOURCES})
ADD_EXECUTABLE(confGenerator ${confGenerator_SOURCES})
ADD_EXECUTABLE(analyser ${analyser_SOURCES})

TARGET_LINK_LIBRARIES(${lib_name} ${CMAKE_DL_LIBS})

TARGET_LINK_LIBRARIES(${exe_name} ${lib_name})
TARGET_LINK_LIBRARIES(confGenerator ${lib_name})
TARGET_LINK_LIBRARIES(analyser ${lib_name})
