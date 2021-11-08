program zeroDimPlasmaChem
    use m_config
    use m_chemistry
    use m_types
    implicit none

    type(CFG_t) :: cfg 
    type(ode_sys) :: odes
    integer, parameter :: dt_num_cond = 2
    real(dp) :: dt_array(2)
    integer, parameter :: dt_ix_chem = 1
    integer, parameter :: dt_ix_drt = 2
    real(dp), parameter :: dt_safety_factor = 0.9_dp
    real(dp) :: time, global_dt
    real(dp) :: dt_max = 1e-11_dp
    real(dp) :: dt_min = 1e-14_dp
    character(len=name_len) :: integrator = "heuns_method"
    integer :: time_integrator = -1
    integer, parameter :: integration_advance_steps(3) = [1,2,2]
    integer :: i_electron = -1
    integer :: i_1pos_ion = -1
    integer :: ix_electron = -1
    integer :: i_e_fld = -1
    logical :: field_table_use
    real(dp) ::  field_amplitude
    real(dp), allocatable :: field_table_times(:)
    real(dp), allocatable :: field_table_fields(:)
    real(dp) :: init_specie_density(2)

    print *, "Inside main prog"
    call CFG_update_from_arguments(cfg)
    call init_modules(cfg, odes)
    print *, "Integration method to be used ", integrator !Debug line
    print *, "ODE system stuff: ", odes%n_vars !Debug linr
    print *, "ODE system stuff: ", odes%var_names(1:odes%n_vars) !Debug linr
    print *, "Initial field value: ", field_amplitude
    print *, "First positive ion: ", odes%var_names(i_1pos_ion)
    print *, "Initial densities: ", odes%vars(i_1pos_ion), odes%vars(i_electron)


    time = 0.0_dp

    global_dt = minval(dt_array)
    ! Setting the initial conditions
    call init_cond_initialize(odes, cfg)
    ! Time integration loop here




    print *, "End of simulation"

    contains

    subroutine init_modules(cfg, odes)
        use m_table_data
        use m_transport_data
        use m_gas
        implicit none
        type(CFG_t),intent(inout) :: cfg
        type(ode_sys),intent(inout) :: odes
        !Initialize the time steps here

        !Initialize the tables used for some reaction rates
        call table_data_initialize(cfg)
        call transport_data_initialize(cfg)
        call gas_initialize(odes, cfg)
        !Read the input reactions
        call chemistry_initialize(odes, cfg)
        !Initialize dt values and the type of time integrator to be used
        call dt_initialize(cfg)

        call cfg_init(odes, cfg)
        call field_initialize(odes, cfg)
        ! Initializing the ode variables
        allocate(odes%vars(odes%n_vars))
    
    end subroutine init_modules

    subroutine dt_initialize(cfg)
        use m_config
        implicit none
        type(CFG_t),intent(inout) :: cfg
        call CFG_add_get(cfg, "dt_max", dt_max, &
         "The maximum timestep (s)")
        call CFG_add_get(cfg, "dt_min", dt_min, &
         "The minimum timestep (s)")
        !call CFG_add_get(cfg, "dt_safety_factor", dt_safety_factor, &
         !"Safety factor for the time step")
        call CFG_add_get(cfg, "time_integrator", integrator, &
         "Time integrator (forward_euler, rk2, heuns_method)")
        !> [integrators]
        select case (integrator)
        case ("forward_euler")
           time_integrator = 1
        case ("rk2")
           time_integrator = 2
        case ("heuns_method")
           time_integrator = 3
        case default
           print *, "Time integrator: ", trim(integrator)
           error stop "Invalid time integrator"
        end select
        dt_array = dt_max
    
    end subroutine dt_initialize

    subroutine cfg_init(odes,  cfg)
      use iso_fortran_env, only: int64
      use m_config
      use m_chemistry
      use m_units_constants
      use m_gas
      use m_transport_data
      implicit none
      type(ode_sys),intent(inout) :: odes
      type(CFG_t),intent(inout) ::  cfg
      integer :: n, k, ix_chemistry
      real(dp) :: tmp

      ! Set electron index
      i_electron = find_ode_var(odes, "e")
      ix_electron = species_index("e")
      do n = n_gas_species+1, n_species
         if (species_charge(n) == 1) then
            i_1pos_ion = species_itree(n)
            !ix_1pos_ion = n
            exit
         end if
      end do
  
      if (i_1pos_ion == -1) error stop "No positive ion species (1+) found"
  

      !Used to initialize MANY other input conditions/etc. Afivo has stuff like Bcs, ion motion stuff, etc, which we dont have
      call add_ode_var(odes, "electric_fld", ix=i_e_fld)


    
    end subroutine cfg_init

    subroutine field_initialize(odes,  cfg)
      use m_config
      use m_table_data
      implicit none
      type(ode_sys),intent(inout) :: odes
      type(CFG_t),intent(inout) ::  cfg
      character(len=string_len) :: field_table

      field_table = undefined_str
      call CFG_add_get(cfg, "field_table", field_table, "File containing applied electric field (V/m) versus time")
      if (field_table /= undefined_str) then
         field_table_use = .true.
         call table_from_file(field_table, "field_vs_time", &
              field_table_times, field_table_fields)
      else
         field_table_use = .false.
         call CFG_add_get(cfg, "field_amplitude", field_amplitude, &
         "The initial applied electric field(V/m)")
      end if
    
    end subroutine field_initialize

    subroutine init_cond_initialize(odes,  cfg)
      use m_config
      use m_table_data
      implicit none
      type(ode_sys),intent(inout) :: odes
      type(CFG_t),intent(inout) ::  cfg

      call CFG_add_get(cfg, "init_electron_density", &
      odes%vars(i_electron), "Initial electron density")
      call CFG_add_get(cfg, "init_first_posIon_density", &
      odes%vars(i_1pos_ion), "Initial first positive ion density")
    
    end subroutine init_cond_initialize

    !Write a subroutine that outputs the simulation data as a .txt file
end program zeroDimPlasmaChem

