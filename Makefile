FC := /usr/bin/h5fc
FFLAGS := -O3 -Wall -fcheck=all -cpp
LDFLAGS :=
TARGET := main
# Required due to bug: https://bugzilla.redhat.com/show_bug.cgi?id=1971826
INCDIRS := /usr/lib64/gfortran/modules

OBJ := m_chemistry.o m_config.o m_gas.o m_lookup_table.o m_spline_interp.o m_output.o	\
m_table_data.o m_transport_data.o m_types.o m_units_constants.o main.o

all: $(TARGET)

# How to get .o object files from .f90 source files
%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get .mod files from .f90 source files (remake only if they have been
# removed, otherwise assume they are up to date)
%.mod: %.f90 %.o
	@test -f $@ || $(FC) -c -o $(@:.mod=.o) $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get executables from .o object files
%: %.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: clean

clean:
	rm -f $(OBJ) *.mod $(TARGET) *.o gmon.out

# Dependencies
$(TARGET): $(OBJ)

# Generated with list_dependencies.sh

example_1.o: m_config.mod
example_2.o: m_config.mod
main.o: m_chemistry.mod
main.o: m_config.mod
main.o: m_gas.mod
main.o: m_lookup_table.mod
main.o: m_output.mod
main.o: m_table_data.mod
main.o: m_transport_data.mod
main.o: m_types.mod
main.o: m_units_constants.mod
m_chemistry.o: m_config.mod
m_chemistry.o: m_gas.mod
m_chemistry.o: m_lookup_table.mod
m_chemistry.o: m_table_data.mod
m_chemistry.o: m_transport_data.mod
m_chemistry.o: m_types.mod
m_chemistry.o: m_units_constants.mod
m_gas.o: m_config.mod
m_gas.o: m_table_data.mod
m_gas.o: m_types.mod
m_gas.o: m_units_constants.mod
m_output.o: m_chemistry.mod
m_output.o: m_config.mod
m_output.o: m_types.mod
m_table_data.o: m_config.mod
m_table_data.o: m_lookup_table.mod
m_table_data.o: m_spline_interp.mod
m_table_data.o: m_types.mod
m_transport_data.o: m_config.mod
m_transport_data.o: m_gas.mod
m_transport_data.o: m_lookup_table.mod
m_transport_data.o: m_spline_interp.mod
m_transport_data.o: m_table_data.mod
m_transport_data.o: m_types.mod
m_transport_data.o: m_units_constants.mod
