#include "cpp_macros.h"
!> Module that stores parameters related to the gas
module m_gas
  use m_types
  !use m_af_types
  use m_units_constants

  implicit none
  private
  ! Check for constant gas density
  logical, public, protected  :: gas_constant_density = .true.
  ! Pressure of the gas in bar
  real(dp), public, protected :: gas_pressure = 1.0_dp

  ! Gas temperature in Kelvin
  real(dp), public, protected :: gas_temperature = 300.0_dp

  ! List of gas fractions (the last one is 1.0 for the total density)
  real(dp), allocatable, public, protected :: gas_fractions(:)

  ! List of gas densities (the last one is the total density)
  real(dp), allocatable, public, protected :: gas_densities(:)

  ! List of gas components (the last one is M, the total density)
  character(len=comp_len), allocatable, public, protected :: gas_components(:)

  ! Gas number density (1/m3)
  real(dp), public, protected :: gas_number_density

  ! Inverse gas number density (1/m3)
  real(dp), public, protected :: gas_inverse_number_density

  ! Convert V/m to Townsend
  real(dp), parameter, public :: SI_to_Townsend = 1e21_dp

  ! Convert Townsend to V/m
  real(dp), parameter, public :: Townsend_to_SI = 1e-21_dp

  ! Index of the gas number density
  integer, public, protected :: i_gas_dens = -1
  integer, public, protected :: i_gas_pres = -1
  integer, public, protected :: i_gas_temp = -1

  ! Gas mean molecular weight (kg)
  real(dp), public, protected :: gas_molecular_weight = 28.8_dp * UC_atomic_mass

  ! Ratio of heat capacities (polytropic index)
  real(dp), public, protected :: gas_euler_gamma = 1.4_dp


  public :: gas_initialize
  public :: gas_index

contains

  !> Initialize this module
  subroutine gas_initialize(odes, cfg)
    use m_config
    use m_units_constants
    !use m_user_methods
    !use m_dt
    use m_table_data

    type(CFG_t), intent(inout) :: cfg
    type(ode_sys), intent(inout) :: odes
    integer :: n

    call add_ode_var(odes, "gas_density", ix=i_gas_dens)
    call add_ode_var(odes, "gas_pressure", ix=i_gas_pres)
    call add_ode_var(odes, "gas_temperature", ix=i_gas_temp)


    call CFG_add_get(cfg, "gas%pressure", gas_pressure, &
         "The gas pressure (bar)")
    call CFG_add_get(cfg, "gas%temperature", gas_temperature, &
         "The gas temperature (Kelvin)")
    call CFG_add_get(cfg, "gas%molecular_weight", gas_molecular_weight, &
         "Gas mean molecular weight (kg), for gas dynamics")

    ! Ideal gas law
    gas_number_density = 1e5_dp * gas_pressure / &
         (UC_boltzmann_const * gas_temperature)
    gas_inverse_number_density = 1/gas_number_density
    !print *, "Inside gas init, gas number density: ", gas_number_density

    call CFG_add(cfg, "gas%components", ["N2", "O2"], &
         "Gas component names", .true.)
    call CFG_add(cfg, "gas%fractions", [0.8_dp, 0.2_dp], &
         "Gas component fractions", .true.)

    call CFG_get_size(cfg, "gas%components", n)
    allocate(gas_components(n+1))
    allocate(gas_fractions(n+1))
    call CFG_get(cfg, "gas%components", gas_components(1:n))
    call CFG_get(cfg, "gas%fractions", gas_fractions(1:n))

    gas_components(n+1) = "M"
    gas_fractions(n+1)  = 1.0_dp

    if (any(gas_fractions < 0.0_dp)) &
         error stop "gas%fractions has negative value"
    if (abs(sum(gas_fractions(1:n)) - 1.0_dp) > 1e-4_dp) &
         error stop "gas%fractions not normalized"

    gas_densities = gas_fractions * gas_number_density

  end subroutine gas_initialize

 !> Find index of a gas component, return -1 if not found
  elemental integer function gas_index(name)
    character(len=*), intent(in) :: name
    do gas_index = 1, size(gas_components)
       if (gas_components(gas_index) == name) exit
    end do
    if (gas_index == size(gas_components)+1) gas_index = -1
  end function gas_index
end module m_gas
