!> Module that stores parameters related to the gas
module m_output
  use m_types

  implicit none
  character(len=string_len), public, protected :: output_name


  public :: output_initialize
  public :: output_HDF5_solution
  private :: write_var_names
  private :: write_time_steps

contains

  !> Initialize this module
  subroutine output_initialize(cfg)
    use m_config
    use HDF5

    type(CFG_t), intent(inout) :: cfg
    integer :: error
    integer(HID_T) :: file_id, group_id1, group_id2
    character(len=string_len) :: outname

    call CFG_add_get(cfg, "output%name", outname, &
      "Name of the output file")

    output_name = trim(outname)

    call h5open_f(error)
    call h5fcreate_f(trim(outname) // ".h5", H5F_ACC_TRUNC_F, file_id, error)

    call h5gcreate_f(file_id, "times", group_id1, error)
    call h5gclose_f(group_id1, error)

    !call h5gcreate_f(file_id, "tsteps", group_id2, error)
    !call h5gclose_f(group_id2, error)

    call h5fclose_f(file_id, error)
    call h5close_f(error)


  end subroutine output_initialize

  subroutine output_HDF5_solution(odes_s, filename, t, iter)
     use m_chemistry
     use m_types
     use HDF5
     implicit none
     type(ode_sys), intent(in) :: odes_s
     character(len=*), intent(in) :: filename
     real(dp), intent(in) :: t
     integer, intent(in) :: iter 
     !character(len=50), save :: fmt, fmt_header
     !character(len=50) :: test_fmt
     logical, save :: first_time = .true.
     integer :: my_unit, n, i, n_vars 
     integer :: error
     integer :: odeVar_rank, odeVar_dspace_i, odeVar_dset_i, dvar_i
     
     integer(HSIZE_T) :: data_dims(1)
     integer(HID_T) :: file_id, dspace_id, dset_id
     integer(HID_T) :: grp_time_id
     character(len=name_len) :: dset_name
     !integer(HID_T), allocatable :: dset_id(:)

     

     n_vars = sum(odes_s%var_matrix)
     !allocate(dset_id(n_vars))


     !write(dset_name, "(A,E15.6)") "time_", t
     write(dset_name, "(A,I0.6)") "time_", iter 

     if (first_time) then
        first_time = .false.
        call write_var_names(odes_s, filename)
     end if
     call write_time_steps(t, filename)

     call h5open_f(error)
     call h5fopen_f(trim(filename) // ".h5", H5F_ACC_RDWR_F, file_id, error)
     call h5gopen_f(file_id, "/times", grp_time_id, error)
     data_dims = n_vars !We need to do this for HDF5 for integer type conversion
     call h5screate_simple_f(1, data_dims, dspace_id, error)
     call h5dcreate_f(grp_time_id, trim(dset_name), H5T_NATIVE_DOUBLE, &
        dspace_id, dset_id, error)
     
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
      pack(odes_s%vars, odes_s%var_matrix==1), data_dims, error)



     call h5sclose_f(dspace_id, error)
     call h5dclose_f(dset_id, error)
     call h5gclose_f(grp_time_id, error)
     call h5fclose_f(file_id, error)
     call h5close_f(error)

     !write(fmt, "A, I0, A"), 
     !fmt = "(E16.8, E16.8)"
     !open(newunit=my_unit, file=trim(filename), action="write", &
     !   position="append")

     !write(my_unit, fmt) t, odes_s%vars(idx)
     !close(my_unit)
   end subroutine output_HDF5_solution

   subroutine write_var_names(o_s, filename)
     use m_chemistry
     use m_types
     implicit none
     type(ode_sys), intent(in) :: o_s
     character(len=*), intent(in) :: filename
     character(len=name_len), allocatable :: var_names(:)
     integer :: opVarNum
     integer :: i, iter, f_unit
     character(len=string_len) :: varFname

     varFname = trim(filename) // "_var_names.txt"
     opVarNum = sum(o_s%var_matrix)
     allocate(var_names(1:opVarNum))

     open(newunit=f_unit, file=trim(varFname), action="write")
     iter = 1
     do i=1, o_s%n_vars
      if (o_s%var_matrix(i) == 0) then
        cycle
      end if
      var_names(iter) = o_s%var_names(i)
      !print *, "write_var_names", trim(var_names(iter))
      iter = iter + 1

      write(f_unit, "(A)") trim(o_s%var_names(i))
     end do
     write(f_unit, "(A)")"-----------------------"
     close(f_unit)
     
   
   end subroutine write_var_names
   subroutine write_time_steps(time, filename)
     use m_types
     implicit none
     real(dp), intent(in) :: time
     character(len=*), intent(in) :: filename
     integer :: opVarNum
     integer :: i, iter, f_unit
     character(len=string_len) :: varFname

     varFname = trim(filename) // "_var_names.txt"

     open(newunit=f_unit, file=trim(varFname), action="write", &
        position="append")
      write(f_unit, "(E16.8)") time
     close(f_unit)
     
   
   end subroutine write_time_steps


end module m_output
