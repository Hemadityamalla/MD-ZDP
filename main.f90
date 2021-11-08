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
    real(dp) :: dt_max = 1e-11_dp
    real(dp) :: dt_min = 1e-14_dp
    character(len=name_len) :: integrator = "heuns_method"
    integer :: time_integrator = -1




    
    
    print *, "Inside main prog"
    call CFG_update_from_arguments(cfg)
    call init_modules(cfg)
    print *, "Integration method to be used ", integrator !Debug line

    ! Time integration loop here



    print *, "End of simulation"

    contains

    subroutine init_modules(cfg)
        use m_table_data
        use m_transport_data
        use m_gas
        implicit none
        type(CFG_t),intent(inout) :: cfg
        !Initialize the time steps here

        !Initialize the tables used for some reaction rates
        call table_data_initialize(cfg)
        call transport_data_initialize(cfg)
        call gas_initialize(cfg)
        !Read the input reactions
        call chemistry_initialize(cfg)
        !Initialize dt values and the type of time integrator to be used
        call dt_initialize(cfg)
    
    end subroutine init_modules

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

    !Write a subroutine that outputs the simulation data as a .txt file
end program zeroDimPlasmaChem

