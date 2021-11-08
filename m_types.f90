!> Module with basic types
module m_types

  implicit none
  public

  character(len=*), parameter :: undefined_str = "UNDEFINED"

  !> Default length of strings
  integer, parameter :: string_len = 200

  !> Default length of names
  integer, parameter :: name_len = 20

  !> Default length of component names
  integer, parameter :: comp_len = 20

  !> Maximum size of the ODE system
  integer, parameter :: max_var_lim = 100

  !> Double precision kind
  integer, parameter :: dp = kind(0.0d0)

  type ode_sys
    integer :: n_vars
    character(len=name_len) :: var_names(max_var_lim)
    integer :: var_num_copies(max_var_lim) = 1
    real(dp), allocatable :: vars(:)
  end type ode_sys

  contains

  subroutine add_ode_var(ode_system, name, n_copies, ix)
    type(ode_sys), intent(inout) :: ode_system
    character(len=*), intent(in) :: name
    integer, intent(in), optional :: n_copies
    integer, intent(out), optional :: ix

    integer :: n, ncpy
    ncpy = 1; if (present(n_copies)) ncpy = n_copies
    if (ncpy < 1) error stop "variable copies < 1"

    do n=1, ncpy
        ode_system%n_vars = ode_system%n_vars + 1
        if (n==1) then
            if (present(ix)) ix = ode_system%n_vars
            ode_system%var_names(ode_system%n_vars) = name
            ode_system%var_num_copies = ncpy
        else
            write(ode_system%var_names(ode_system%n_vars), "(A,I0)") trim(name) // "_", n
            ode_system%var_num_copies(ode_system%n_vars) = 0
        end if
    end do
  end subroutine add_ode_var
 
  integer function find_ode_var(ode_s, name)
    type(ode_sys), intent(in) :: ode_s
    character(len=*), intent(in) :: name
    integer :: n

    do n=1,ode_s%n_vars
      if (ode_s%var_names(n) == name) exit
    end do

    if (n == ode_s%n_vars+1) then
      print *, "Variable name: ", trim(name)
      error stop "find_ode_var: variable not found"
    end if
  end function find_ode_var

end module m_types
