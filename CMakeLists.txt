CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

if(POLICY CMP0054)
	cmake_policy(SET CMP0054 NEW)
endif()

if(POLICY CMP0037)
        cmake_policy(SET CMP0037 OLD)
endif()

# this is required when working with newer gcc versions (> 5.4)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FORCE_INLINES")

OPTION(CUDA "Set to ON to compile with CUDA support" OFF)
OPTION(MPI "Set to ON to compile with MPI support" OFF)
OPTION(Debug "Set to ON to compile with debug symbols" OFF)
OPTION(G "Set to ON to compile with optimisations and debug symbols" OFF)
OPTION(INTEL "Use the Intel compiler" OFF)
OPTION(MOSIX "Make oxDNA compatible with MOSIX" OFF)
OPTION(SIGNAL "Enable SignalManager - set to OFF for OSX compatibility" OFF)
OPTION(DOUBLE "Set the numerical precision to double" ON)
OPTION(CUDA_DOUBLE "Set the numerical precision for the CUDA backend to double" OFF)
OPTION(CUDA_COMMON_ARCH "Set to OFF to autodetect the GPU and compile for its architecture, set to ON (default) to compile for the most common architectures" ON)
OPTION(NATIVE_COMPILATION "Set to OFF to compile without the -march=native flag. This may be required when compiling binaries to be used elsewhere" ON)

OPTION(Python "Enable Python bindings" OFF)
OPTION(OxpySystemInstall "Install oxpy as a system-wide package" OFF)

# these operations have to be performed before PROJECT(oxDNA) or we will have problems at linking time
if(INTEL)
	FIND_PROGRAM(CMAKE_C_COMPILER NAMES icc)
	FIND_PROGRAM(CMAKE_CXX_COMPILER NAMES icpc)
	FIND_PROGRAM(CMAKE_AR NAMES xiar)
	FIND_PROGRAM(CMAKE_LINKER NAMES xild)

	if(CMAKE_C_COMPILER MATCHES CMAKE_C_COMPILER-NOTFOUND OR
	CMAKE_CXX_COMPILER MATCHES CMAKE_CXX_COMPILER-NOTFOUND OR
	CMAKE_AR MATCHES CMAKE_AR-NOTFOUND OR
	CMAKE_LINKER MATCHES CMAKE_LINKER-NOTFOUND)
	    MESSAGE(FATAL_ERROR "Cannot find Intel compiler")
	ENDIF()

	if(CUDA)
		MESSAGE(FATAL_ERROR "oxDNA with CUDA support cannot be compiled with the Intel compiler")
	ENDIF(CUDA)
ENDIF(INTEL) 

PROJECT(rheoCell)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmakeHelp)
find_package(Eigen3 REQUIRED)

# The default behaviour is to compile with optimisations
SET(CMAKE_BUILD_TYPE Release)

# c++-14 support is mandatory
IF(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.1)
	SET(CMAKE_CXX_STANDARD 14)
	SET(CMAKE_CXX_STANDARD_REQUIRED ON)
ELSE()
	ADD_COMPILE_OPTIONS(-std=c++14)
ENDIF()

IF(Debug)
	SET(CMAKE_BUILD_TYPE "Debug")
	MESSAGE(STATUS "Compiling in Debug mode")
ENDIF(Debug)

IF(G)
	SET(CMAKE_BUILD_TYPE RelWithDebInfo)
	MESSAGE(STATUS "Compiling in Release+Debug symbols mode")
ENDIF(G)

IF(MOSIX)
	ADD_DEFINITIONS(-DMOSIX)
	MESSAGE(STATUS "Disabling timings")
ENDIF(MOSIX)

IF(SIGNAL)
	ADD_DEFINITIONS(-DSIGNAL)
ENDIF(SIGNAL)

IF(NOT DOUBLE)
	ADD_DEFINITIONS(-DFLOAT_PRECISION)
ENDIF()

include_directories(src/extern 
	${EIGEN3_INCLUDE_DIR})

ADD_DEFINITIONS(-Wall -Wshadow -Wextra -Wno-unused-parameter)

IF(NOT Debug)
	IF(INTEL)
		ADD_DEFINITIONS(-fast)
	ELSE()
		ADD_DEFINITIONS(-ffast-math)
		if(NATIVE_COMPILATION)
			include(CheckCXXCompilerFlag)
            unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
            CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
            if(COMPILER_SUPPORTS_MARCH_NATIVE)
                    ADD_DEFINITIONS(-march=native)
			else()
				message(STATUS "The '-march=native' flag is not supported, compiling without it")
            endif()
		else()
			message(STATUS "Compiling without the '-march=native' flag")
		endif()
	ENDIF(INTEL)
ENDIF()

IF(MPI)
	ADD_DEFINITIONS(-DHAVE_MPI)
ENDIF(MPI)

# return the current date as day/month/year
MACRO(CURRENT_DATE result)
    IF(WIN32)
        EXECUTE_PROCESS(COMMAND "cmd" " /C date /T" OUTPUT_VARIABLE ${result})
        string(REGEX REPLACE "(..)/(..)/..(..).*" "\\1/\\2/\\3" ${result} ${${result}})
    ELSEIF(UNIX)
        EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE ${result})
        string(REGEX REPLACE "(..)/(..)/..(..).*" "\\1/\\2/\\3" ${result} ${${result}})
    ELSE(WIN32)
        MESSAGE(SEND_ERROR "date not implemented")
        SET(${result} 000000)
    ENDIF(WIN32)
ENDMACRO(CURRENT_DATE)

CURRENT_DATE(today)
ADD_DEFINITIONS(-D'RELEASE="${tag}"' -D'GIT_COMMIT="${revision}"' -D'BUILD_TIME="${today}"')

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

ADD_SUBDIRECTORY(src)
