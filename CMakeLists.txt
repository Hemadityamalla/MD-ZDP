cmake_minimum_required(VERSION 3.9.0)
project(zdPlasmaChem
	VERSION 1.0
	DESCRIPTION "Compiling the zero-dimensional plasma"
	LANGUAGES Fortran)

enable_language(Fortran)

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

add_executable(zdPlasmaChem ${SOURCES})

set(CMAKE_Fortran_FLAGS "-O0 -Wall  -fcheck=array-temps,bounds,do,mem,pointer -pg -ffpe-trap=invalid,zero,overflow -pedantic -finit-real=snan")


