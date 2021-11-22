!> Module that stores parameters related to the gas
module m_output
  use m_types

  implicit none
  character(len=string_len), public :: output_name


  public :: output_initialize
  public :: output_HDF5_solution

contains

  !> Initialize this module
  subroutine output_initialize(cfg)
    use m_config

    type(CFG_t), intent(inout) :: cfg

    call CFG_add_get(cfg, "output%name", output_name, &
      "Name of the output file")

  end subroutine output_initialize

  subroutine output_HDF5_solution(odes_s, filename, t, idx)
     use m_chemistry
     use m_types
     use HDF5
     implicit none
     type(ode_sys), intent(in) :: odes_s
     character(len=*), intent(in) :: filename
     character(len=50), save :: fmt, fmt_header
     character(len=50) :: test_fmt
     integer :: my_unit, n, i, n_vars, error, space_rank
     integer, intent(in) :: idx
     real(dp), intent(in) :: t
     logical, save :: first_time = .true.
     
     integer(HSIZE_T) :: data_dims(1)
     integer(HID_T) :: file_id, dspace_id
     integer(HID_T), allocatable :: dset_id(:)
     integer :: max_rows = 5

     

     n_vars = odes_s%n_vars
     allocate(dset_id(n_vars))
     if (first_time) then
        first_time = .false.

        call h5open_f(error)

        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error)

        space_rank = 1
        data_dims(1) = max_rows

        call h5screate_simple_f(space_rank, data_dims, dspace_id, error)
        do i=1,n_vars
               call h5dcreate_f(file_id, trim(odes_s%var_names(i)), H5T_NATIVE_DOUBLE, dspace_id, dset_id(i), error)
               call h5dclose_f(dset_id(i), error)

        end do

        




        call h5fclose_f(file_id, error)

        call h5close_f(error)
        

     end if

     !write(fmt, "A, I0, A"), 
     !fmt = "(E16.8, E16.8)"
     !open(newunit=my_unit, file=trim(filename), action="write", &
     !   position="append")

     !write(my_unit, fmt) t, odes_s%vars(idx)
     !close(my_unit)
   end subroutine output_HDF5_solution

end module m_output
