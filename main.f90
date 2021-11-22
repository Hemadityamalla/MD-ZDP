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
    real(dp) :: dt_array(2)
    integer, parameter :: dt_ix_chem = 1
    integer, parameter :: dt_ix_drt = 2
    real(dp), parameter :: dt_safety_factor = 0.9_dp
    real(dp) :: time, global_dt
    real(dp) :: end_time
    real(dp) :: dt_max = 1e-11_dp
    real(dp) :: dt_min = 1e-14_dp
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
    !character(len=string_len) :: output_name

    print *, "Inside main prog"
    call CFG_update_from_arguments(cfg)
    call init_modules(cfg, odes)
    
    print *, "Integration method to be used ", integrator !Debug line
    print *, "ODE system number of variables: ", odes%n_vars !Debug linr
    print *, "ODE system variable names: ", odes%var_names(1:odes%n_vars) !Debug linr
    print *, "Initial field value: ", field_amplitude
    print *, "First positive ion: ", odes%var_names(i_1pos_ion)
    print *, "N_gas_species: ", n_gas_species
    print *, "N_species: ", n_species
    print *, "Species list: ", species_list(n_gas_species+1:n_species)
    print *, "species_itree: ", species_itree(n_gas_species+1:n_species)
    print *, "Index of electric field, electron: ", find_ode_var(odes, "electric_fld"), find_ode_var(odes, "e")
    ! End debug lines


    time = 0.0_dp

    global_dt = minval(dt_array)
    ! Setting the initial conditions
    call init_cond_initialize(odes, cfg)
    print *, "Initial densities: ", odes%vars(i_1pos_ion), odes%vars(i_electron)
    ! Time integration loop here
    do while (time < end_time)
      ! Add functionality to compute the wall clock time 

      print *, "Time: ", time
      call output_HDF5_solution(odes, trim(output_name), time, i_electron)
      call ode_advance(odes, global_dt, &
         species_itree(n_gas_species+1:n_species), time_integrator)

      print *, "Electron density: ", odes%vars(i_electron)
      time = time + global_dt



    end do
    

    print *, "End of simulation"

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
        ! Initializing the ode variables rhs -- for the rates
        allocate(odes%vars_rhs(odes%n_vars))

        call output_initialize(cfg)
    
    end subroutine init_modules

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
      print *, "i_electron", i_electron
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
         print *, "Inside ode_advance, fwd_euler, var_idx: ", var_idx
         print *, "Print ode_rhs using var_idx", ode_s%vars_rhs(var_idx)
         call compute_rhs(ode_s, var_idx)
         print *, "Print ode_rhs using var_idx after update", ode_s%vars_rhs(var_idx)
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
    
    !subroutine output_HDF5_solution(odes_s, filename, t, idx)
    !  use m_chemistry
    !  use m_types
    !  use HDF5
    !  implicit none
    !  type(ode_sys), intent(in) :: odes_s
    !  character(len=*), intent(in) :: filename
    !  character(len=50), save :: fmt, fmt_header
    !  character(len=50) :: test_fmt
    !  integer :: my_unit, n, i, n_vars, error, space_rank
    !  integer, intent(in) :: idx
    !  real(dp), intent(in) :: t
    !  logical, save :: first_time = .true.
    !  
    !  integer(HSIZE_T) :: data_dims(1)
    !  integer(HID_T) :: file_id, dspace_id
    !  integer(HID_T), allocatable :: dset_id(:)
    !  integer :: max_rows = 5

    !  

    !  n_vars = odes_s%n_vars
    !  allocate(dset_id(n_vars))
    !  if (first_time) then
    !     first_time = .false.

    !     call h5open_f(error)

    !     call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error)

    !     space_rank = 1
    !     data_dims(1) = max_rows

    !     call h5screate_simple_f(space_rank, data_dims, dspace_id, error)
    !     do i=1,n_vars
    !            call h5dcreate_f(file_id, trim(odes_s%var_names(i)), H5T_NATIVE_DOUBLE, dspace_id, dset_id(i), error)
    !            call h5dclose_f(dset_id(i), error)

    !     end do

    !     




    !     call h5fclose_f(file_id, error)

    !     call h5close_f(error)
    !     

    !  end if

    !  !write(fmt, "A, I0, A"), 
    !  !fmt = "(E16.8, E16.8)"
    !  !open(newunit=my_unit, file=trim(filename), action="write", &
    !  !   position="append")

    !  !write(my_unit, fmt) t, odes_s%vars(idx)
    !  !close(my_unit)
    !end subroutine output_HDF5_solution
end program zeroDimPlasmaChem

