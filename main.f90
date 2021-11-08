program zeroDimPlasmaChem
    use m_config
    use m_chemistry
    use m_types
    implicit none

    type(CFG_t) :: cfg 
    type ode_sys
        integer :: n_vars
        character(len=name_len) :: var_names(max_var_lim)
        integer :: var_num_copies(max_var_lim) = 1
        real(dp), allocatable :: vars(:)

    end type ode_sys

    
    
    print *, "Inside main prog"
    call CFG_update_from_arguments(cfg)
    call init_modules(cfg)

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
end program zeroDimPlasmaChem

