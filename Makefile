FC := gfortran
FFLAGS := -O3 -Wall  -fcheck=array-temps,bounds,do,mem,pointer -pg -ffpe-trap=invalid,zero,overflow -pedantic -finit-real=snan -fbounds-check
LDFLAGS := 
TARGET := zdPlasmaChem
OBJ := m_chemistry.o m_config.o m_gas.o m_lookup_table.o m_spline_interp.o m_table_data.o m_transport_data.o m_types.o m_units_constants.o main.o

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) -D NDIM=0


m_config.o: m_config.f90
	$(FC) -c $^
m_spline_interp.o: m_spline_interp.f90
	$(FC) -c $^
m_lookup_table.o: m_lookup_table.f90
	$(FC) -c $^
m_types.o: m_types.f90
	$(FC) -c $^
main.o: m_types.o m_config.o m_chemistry.o main.f90 m_transport_data.o 
	$(FC) -c $^
m_units_constants.o: m_units_constants.f90
	$(FC) -c $^
m_gas.o: m_gas.f90 m_types.o m_units_constants.o m_config.o m_table_data.o
	$(FC) -c $^
m_table_data.o: m_table_data.f90 m_types.o m_spline_interp.o m_config.o m_lookup_table.o
	$(FC) -c $^

m_transport_data.o: m_transport_data.f90 m_lookup_table.o m_spline_interp.o m_types.o m_config.o m_table_data.o m_units_constants.o m_gas.o
	$(FC) -c $^

m_chemistry.o: m_chemistry.f90 m_types.o m_lookup_table.o m_table_data.o m_units_constants.o m_transport_data.o m_config.o
	$(FC) -c $^



.PHONY: clean

clean:
	rm -f $(OBJ) *.mod $(TARGET) *.o gmon.out
