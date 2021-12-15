cmake_minimum_required(VERSION 3.9.0)
project(zdPlasmaChem
	VERSION 1.0
	DESCRIPTION "Compiling the zero-dimensional plasma"
	LANGUAGES Fortran)

enable_language(Fortran)


add_executable(${PROJECT_NAME} main.f90)
add_library(config m_config.f90)
add_dependencies()

set(SOURCES
	main.f90
	m_config.f90
	m_spline_interp.f90
	m_lookup_table.f90
	m_types.f90
	m_units_constants.f90
	m_gas.f90
	m_table_data.f90
	m_transport_data.f90
	m_chemistry.f90
	)


set(CMAKE_Fortran_FLAGS "-O0 -Wall  -fcheck=array-temps,bounds,do,mem,pointer -pg -ffpe-trap=invalid,zero,overflow -pedantic -finit-real=snan")

set(HDF5_INCLUDE_DIRS "/home/hemaditya/hdf5_config_fortran_build/fortran/src")
set(HDF5_LIB_DIRS "/home/hemaditya/hdf5_config_fortran_build/fortran/src")
set(HDF5_LIBRARIES "-lhdf5 -lhdf5_fortran")
MESSAGE(STATUS "HDF5_INCLUDE_DIRS set to:${HDF5_INCLUDE_DIRS}")
target_include_directories(hdf5_test PUBLIC ${HDF5_INCLUDE_DIRS})
MESSAGE(STATUS "HDF5_LIBRARIES set to:${HDF5_LIBRARIES}")
target_link_directories(hdf5_test PUBLIC ${HDF5_LIB_DIRS})
target_link_libraries(hdf5_test PUBLIC ${HDF5_LIBRARIES})

