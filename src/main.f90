program zeroDimPlasmaChem
    use m_config
    use m_chemistry
    use m_types
    use m_output
    implicit none

    type(CFG_t) :: cfg 
    type(ode_sys) :: odes
    integer, parameter :: n_cell = 1
    integer, parameter :: dt_num_cond = 2
    real(dp) :: dt_array(1:2)
    integer, parameter :: dt_ix_chem = 1
    integer, parameter :: dt_ix_drt = 2
    real(dp), parameter :: dt_safety_factor = 0.9_dp
    real(dp), parameter :: dt_chem_nmin = 1e11_dp
    real(dp) :: time, global_dt
    real(dp) :: end_time
    real(dp) :: dt_max = 1e-11_dp
    real(dp) :: dt_min = 1e-12_dp
    character(len=name_len) :: integrator = "heuns_method"
    integer :: time_integrator = -1
    integer, parameter :: n_integrators = 3
    integer, parameter :: fwd_euler = 1
    integer, parameter :: rk2 = 2
    integer, parameter :: heuns = 3
    integer, parameter :: integration_advance_steps(n_integrators) = [1,2,2]
    integer :: i_electron = -1
    integer :: i_1pos_ion = -1
    integer :: ix_electron = -1
    integer :: i_e_fld = -1
    logical :: field_table_use
    real(dp) ::  field_amplitude
    real(dp), allocatable :: field_table_times(:)
    real(dp), allocatable :: field_table_fields(:)
    real(dp) :: init_specie_density(2)
    integer :: test_iterator
    integer :: output_number = 1
    !character(len=string_len) :: output_name

    print *, "Inside main prog:"
    call CFG_update_from_arguments(cfg)
    call init_modules(cfg, odes)
    
    print *, "Integration method to be used ", integrator !Debug line
    print *, "ODE system number of variables: ", odes%n_vars !Debug linr
    !print *, "ODE system variable names: ", odes%var_names(1:odes%n_vars) !Debug linr
    print *, "Initial field value: ", field_amplitude
    print *, "First positive ion: ", odes%var_names(i_1pos_ion)
    print *, "N_gas_species: ", n_gas_species
    print *, "N_species: ", n_species
    print *, "Species list: ", species_list(n_gas_species+1:n_species)
    print *, "species_itree: ", species_itree(n_gas_species+1:n_species)
    !print *, "Index of electric field, electron: ", find_ode_var(odes, "electric_fld"), find_ode_var(odes, "e")
    ! End debug lines


    time = 0.0_dp
    test_iterator = 0

    global_dt = 0.1_dp*maxval(dt_array)
    ! Setting the initial conditions
    print *, "Initial densities: ", &
       odes%vars(species_itree(n_gas_species+1:n_species))
    ! Time integration loop here
    do while (time < end_time)
      !print *, "Error: ", exp(1.0e-16*2.4143212320551650E+25*time)-odes%vars(i_electron)
      ! Add functionality to compute the wall clock time 
      if (global_dt < dt_min) then
         print *, "dt(", global_dt, ") smaller than dt_min(", dt_min, ")"
         error stop "dt too small"

      endif

      if (mod(test_iterator, output_number) == 0) then
         print *, "Time: ", time
         call output_HDF5_solution(odes, trim(output_name), time, test_iterator)
      endif
      call ode_advance(odes, global_dt, &
         species_itree(n_gas_species+1:n_species), time_integrator)

      !print *, "Electron density: ", odes%vars(i_electron)
      !print *, "Exact: ", exp(1.0e-16*2.4143212320551650E+25*time)
      !test_varnames = pack(odes%var_names, odes%var_matrix==1)
      !print *, "Varnames: ",pack(odes%vars, odes%var_matrix==1) 
      !print *, "Num actual vars: ", sum(odes%var_matrix)
      call dt_constraints(dt_array, odes)
      !print *, "dt array: ", dt_array
      !global_dt = min(dt_safety_factor*minval(dt_array), 2*dt_max)
      time = time + global_dt
      test_iterator = test_iterator + 1



    end do
    !Outputting the last value of the simulation
    call output_HDF5_solution(odes, trim(output_name), time, test_iterator)
    

    print *, "End of simulation at t= ", time

    contains

    subroutine init_modules(cfg, odes)
        use m_table_data
        use m_transport_data
        use m_gas
        use m_output
        implicit none
        type(CFG_t),intent(inout) :: cfg
        type(ode_sys),intent(inout) :: odes
        !Initialize the time steps here

        print *, "Beginning init modules:", odes%n_vars
        !Initialize the tables used for some reaction rates
        call table_data_initialize(cfg)
        call transport_data_initialize(cfg)
        call gas_initialize(odes, cfg)
        call field_initialize(odes, cfg)
        !Initialize dt values and the type of time integrator to be used
        call dt_initialize(cfg)
        !Read the input reactions
        call chemistry_initialize(odes, cfg, & 
         integration_advance_steps(time_integrator))

        call cfg_init(odes, cfg)
        call init_variables(odes)
        call init_cond_initialize(odes, cfg)
        call output_initialize(cfg)
    
    end subroutine init_modules
    
    subroutine init_variables(odes)
      use m_chemistry
      use m_gas
      implicit none
      type(ode_sys),intent(inout) :: odes
      integer :: i

      ! Initializing the ode variables
      allocate(odes%vars(1:odes%n_vars))
      ! Initializing the ode variables rhs -- for the rates
      allocate(odes%vars_rhs(1:odes%n_vars))
      ! Initializing the output variable location array
      allocate(odes%var_matrix(1:odes%n_vars))
      !Initializing the above stuff to zero just to be clear
      odes%vars(1:odes%n_vars) = 0.0_dp
      odes%vars_rhs(1:odes%n_vars) = 0.0_dp
      odes%var_matrix(1:odes%n_vars) = 0

      !The standard gas properties (when they are constant)
      do i=1,4
         odes%var_matrix(i) = 1
      end do
      do i= 1, n_species - n_gas_species
         odes%var_matrix(species_itree(n_gas_species+i)) = 1
      end do
      !print *,"Var matrix: ", odes%var_matrix(:)
      

    end subroutine init_variables


    !subroutine print_var_vals(o_s)
    !  implicit none
    !  type(ode_sys),intent(in) :: o_s
    !  integer, allocatable :: var_matrix(:)
    !  integer :: n_vars
    !  integer :: i
    !  n_vars = o_s%n_vars
    !  !allocate(var_matrix(n_vars))
    !  var_matrix = o_s%var_matrix

    !  
    
    !end subroutine print_var_vals

    subroutine dt_initialize(cfg)
        use m_config
        implicit none
        type(CFG_t),intent(inout) :: cfg
        call CFG_add_get(cfg, "dt_max", dt_max, &
         "The maximum timestep (s)")
        call CFG_add_get(cfg, "dt_min", dt_min, &
         "The minimum timestep (s)")
        call CFG_add_get(cfg, "end_time", end_time, &
         "The end time for the simulation (s)")
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
        dt_array = (/dt_max, dt_min/)
    
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
      !print *, "i_electron", i_electron
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


    
    end subroutine cfg_init

    subroutine field_initialize(odes,  cfg)
      use m_config
      use m_table_data
      implicit none
      type(ode_sys),intent(inout) :: odes
      type(CFG_t),intent(inout) ::  cfg
      character(len=string_len) :: field_table

      call add_ode_var(odes, "electric_fld", ix=i_e_fld)
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
      use m_gas
      use m_chemistry
      implicit none
      type(ode_sys),intent(inout) :: odes
      type(CFG_t),intent(inout) ::  cfg
      real(dp) :: empty_real(0)
      integer :: n_cond, varsize, n, sp_idx, nn
      real(dp), allocatable :: ic_dens(:)
      real(dp), allocatable :: ic_dens2(:)
      character(len=name_len), allocatable :: ic_names2(:)
      character(len=string_len) :: empty_string(0)
      ! The idea is to have an array of size 3 for the initial conditions
      ! Array ic_dens: [electron dens, positive ion dens]
      ! To initialize the densities of other species,
      ! they need to be supplied separately using ic_dens2

      !call CFG_add_get(cfg, "init_electron_density", &
      !odes%vars(i_electron), "Initial electron density")
      !call CFG_add_get(cfg, "init_first_posIon_density", &
      !odes%vars(i_1pos_ion), "Initial first positive ion density")

      call CFG_add(cfg, "init_density", empty_real, &
         "Initial density of electron and first positive ion (1/m3)", .true.)
      call CFG_get_size(cfg, "init_density", n_cond)
      if (n_cond /= 2) &
         stop "init_dens size has incompatible size (must be 2)"
      allocate(ic_dens(1:n_cond))
      call CFG_get(cfg, "init_density", ic_dens)
      odes%vars(i_electron) = ic_dens(1)
      odes%vars(i_1pos_ion) = ic_dens(n_cond)
      call CFG_add(cfg, "init_density2", empty_real, &
         "Initial densities of other species (1/m3)", .true.)
      call CFG_add(cfg, "init_density2_names", empty_string, &
         "Names of other species", .true.)
      call CFG_get_size(cfg, "init_density2", varsize)
      if (varsize > 0) then
         allocate(ic_dens2(1:varsize))
         allocate(ic_names2(1:varsize))
         call CFG_get(cfg, "init_density2", ic_dens2)
         call CFG_get(cfg, "init_density2_names", ic_names2)
         sp_idx = -1
         do n=1,varsize
            do nn=1, n_species
               if (species_list(nn) /= trim(ic_names2(n))) then
                  cycle
               else
                  sp_idx = species_itree(nn)
                  exit
               end if
            end do
            if (sp_idx == -1) &
               stop "Other specie name incorrect or doesnt exist"
            odes%vars(sp_idx) = ic_dens2(n)
         end do

      end if


      !Initializing the field amplitude.
      !TODO: This will change if a field table is supplied
      odes%vars(i_e_fld) = field_amplitude
    
    end subroutine init_cond_initialize

    subroutine compute_rhs(ode_s, specie_idx)
      use m_gas
      use m_chemistry
      use m_transport_data
      use m_lookup_table
      use m_units_constants
      implicit none
      type(ode_sys),intent(inout) :: ode_s
      integer, intent(in) :: specie_idx(:)
      real(dp)                   :: rates(n_cell,n_reactions)
      real(dp)                   :: derivs(n_cell,1:n_species)
      real(dp)                   :: dens(n_cell,n_species)
      real(dp)                   :: field(n_cell), tmp
      real(dp), parameter :: eps = 1e-100_dp

      !print *, "Checking sizes: ", n_species, size(ode_s%vars(specie_idx))
      !print *, "Species list: ", species_list(:)
      !print *, "Specie idx:", specie_idx
      !Obtain the field (E/N) in Townsend units
      tmp = 1 / gas_number_density
      field(n_cell) = SI_to_Townsend * tmp * field_amplitude
      
      !if (gas_constant_density) then
      !   tmp = 1 / gas_number_density
      !   field(n_cell) = SI_to_Townsend * tmp * field_amplitude
      !end if

      ! Get the specie densities for a given time step or intermediate time step
      dens(n_cell,n_gas_species+1:n_species) = ode_s%vars(specie_idx)

      call get_rates(field, rates, n_cell)
      !print *, "Rates: ", rates

      call get_derivatives(dens, rates, derivs, n_cell)
      !print *, "Derivatives: ", derivs
      !print *, "Derivatives specie idx: ", derivs(1, specie_idx)

      ode_s%vars_rhs(specie_idx) = derivs(1,n_gas_species+1:n_species)

    
    end subroutine compute_rhs

    subroutine ode_advance(ode_s,  dt, var_idx, integrator_type)
      implicit none
      type(ode_sys),intent(inout) :: ode_s
      integer, intent(in) :: var_idx(:), integrator_type
      real(dp), intent(in) :: dt
   
      if (integrator_type < 1 .or. integrator_type > n_integrators) &
         error stop "Invalid time integrator"
      
      if (any(ode_s%var_num_copies(var_idx) < &
         integration_advance_steps(integrator_type))) &
         error stop "Not enough copies available"
   
      select case (integrator_type)
      case(fwd_euler)
         !print *, "Inside ode_advance, fwd_euler, var_idx: ", var_idx
         !print *, "Print ode_rhs using var_idx", ode_s%vars_rhs(var_idx)
         call compute_rhs(ode_s, var_idx)
         !print *, "Print ode_rhs using var_idx after update", ode_s%vars_rhs(var_idx)
         ode_s%vars(var_idx) = ode_s%vars(var_idx) + dt*ode_s%vars_rhs(var_idx)
         
      case(rk2)
         call compute_rhs(ode_s, var_idx)
         ode_s%vars(var_idx+1) = ode_s%vars(var_idx) + &
            0.5_dp*dt*ode_s%vars_rhs(var_idx)
         call compute_rhs(ode_s, var_idx+1)
         ode_s%vars(var_idx) = ode_s%vars(var_idx) + &
            dt*ode_s%vars_rhs(var_idx+1)

      case(heuns)
         call compute_rhs(ode_s, var_idx)
         ode_s%vars(var_idx+1) = ode_s%vars(var_idx) + &
            dt*ode_s%vars_rhs(var_idx)
         call compute_rhs(ode_s, var_idx+1)
         ode_s%vars(var_idx) = ode_s%vars(var_idx) + &
            0.5_dp*dt*(ode_s%vars_rhs(var_idx) + ode_s%vars_rhs(var_idx+1))

      end select
    
    
    end subroutine ode_advance

    subroutine dt_constraints(dt, o_s)
      use m_chemistry
      use m_transport_data
      use m_units_constants
      use m_table_data
      use m_lookup_table
      use m_gas
      implicit none
      type(ode_sys),intent(inout) :: o_s
      real(dp),intent(inout) :: dt(2)
      real(dp) ::  ne, E_vm
      real(dp) :: mobility
      real(dp), parameter :: eps = 1e-100_dp
      real(dp) :: ratio = 0.0_dp
      real(dp) :: Td

    !Computing the dielectric relaxation time
    ne = o_s%vars(i_electron)
    E_vm = o_s%vars(4) !TODO:make this better
    Td = E_vm*SI_to_Townsend/gas_number_density
    !print *,"Td: ", Td, gas_number_density
    mobility =  LT_get_col(td_tbl, td_mobility, Td)/gas_number_density
    !dt(dt_ix_drt) = UC_eps0/(UC_elem_charge*max(mobility*ne, eps))
    dt(dt_ix_drt) = 1e+19
    !Computing the chemistry time
    ratio = minval((abs(o_s%vars(species_itree(n_gas_species+1:n_species))) &
      + dt_chem_nmin)/& 
      max(abs(o_s%vars_rhs(species_itree(n_gas_species+1:n_species))),eps))

    dt(dt_ix_chem)  = ratio    
    end subroutine dt_constraints

    !Subroutine that outputs the simulation data as a .txt file
    subroutine output_solution(odes_s, filename, t, idx)
      use m_chemistry
      use m_types
      implicit none
      type(ode_sys), intent(in) :: odes_s
      character(len=*), intent(in) :: filename
      character(len=50), save :: fmt, fmt_header
      character(len=50) :: test_fmt
      integer :: my_unit, n, i, n_vars
      integer, intent(in) :: idx
      real(dp), intent(in) :: t
      logical, save :: first_time = .true.

      n_vars = odes_s%n_vars
      if (first_time) then
         first_time = .false.
         
         open(newunit=my_unit, file=trim(filename), action="write")
         write(my_unit, "(A)", advance="no") "time electron_density"


         write(my_unit, *) ""
         close(my_unit)

        !write(test_fmt, "(A,I0,A)"), "(I6,", 32, "E16.8,I12,1E16.8,I3)"
        !print *, "PRINT TEST FORMAT", test_fmt

      end if

      !write(fmt, "A, I0, A"), 
      fmt = "(E16.8, E16.8)"
      open(newunit=my_unit, file=trim(filename), action="write", &
         position="append")

      write(my_unit, fmt) t, odes_s%vars(idx)
      close(my_unit)
    end subroutine output_solution
    
end program zeroDimPlasmaChem

