MODULE BOLTZMANN
    !This module calculates the velocity (A/fs) of each atom in a given molecule
    !using the Boltzmann distribtion. 
    !Input file required for initialization: "dft.inp"
    !Run with command line arguments or interactive user input
    !  Command line arguments: <program_executable> <temperature_in_kelvin> <seed>
    USE MOD_IO
    USE MOD_GLOBAL
    USE MARSAGLIA
    implicit none

CONTAINS
    
SUBROUTINE initialize(input_filename)
  character(len=*), intent(in) :: input_filename
  character(256) :: temp_str, seed_str

  call read_input_file(input_filename)

  ! Check command line arguments. If two are not given, then user must input from the terminal
  if (command_argument_count().NE.2) then
      write(*,*) '### USER-INPUT MODE ###'
      write(*,*) 'NOTE: Alternatively, run the program as follows w/ command arguments: '
      write(*,*) '  program_name <temperature> <seed>'
      call input_user_terminal
  else
    ! Get the temperature & seed argument
    call get_command_argument(1, temp_str)
    call get_command_argument(2, seed_str)

    read(temp_str,*) temperature_ions
    read(seed_str,*) ion_velocity_init_seed
  end if

  call compute_atomic_masses

END SUBROUTINE

SUBROUTINE calculate_atomic_velocities()
    !This subroutine gives each atom (ion), which is allowed to move,
    !a random velocity in x,y,z directions according to the Maxwell-Boltzmann
    !distribution at certain temperature
    !Input parameters:
    !  temperature - temperature [in Kelvins]
    !Output data:
    !  v(1:3,1:N_total_atom) - an array where the random velocities are returned [in A/femtosec]
    !,N_total_atom),temperature !NOTE TO SAM: THIS SHOULD ALREADY BE DEFINED AS SHOULD BE ALTERED
    integer :: i,ai,j,elem
    real(8) :: f1,f2,temp,av_kin,tot_kin
    
    !Seed uniform random generator
    call kiss_setseed(ion_velocity_init_seed,43,59,171)
    !Generate velocities

    f1=k_Boltzmann*temperature_ions
    do i=1,N_total_atom
      elem=atom_atomic_number(i)
      if (elem/=0) then
        f2=sqrt(f1/atom_mass(i))
        do j=1,3
          atom_velocity(j,i)=rndnor()*f2
        enddo
      else
        atom_velocity(1:3,i)=0.0d0
      endif
    enddo
    
    !Since the number of ions in not infinite, after generating the
    !x,y,z-velocities with Gaussian distribution we will find the average
    !kinetic energy per ion to be slightly different from k*T (due to
    !statistical error). Here we simply normalize all the velocities in
    !such a way that the relation E_kin = k*T holds exactly.
    call compute_atomic_temperature(temp,av_kin,tot_kin)
    f1=sqrt(temperature_ions/temp)
    do i=1,N_total_atom
        elem=atom_atomic_number(i)
        if (elem/=0) then
          atom_velocity(1:3,i)=f1*atom_velocity(1:3,i)
        endif
    enddo
END SUBROUTINE calculate_atomic_velocities

SUBROUTINE compute_atomic_temperature(temper,av_kin,tot_kin)
    !This subroutine computes the temperature, as well as the total
    !and average kinetic energies of ions. It only counts the ions
    !that are allowed to move.
    !Input parameters:
    !  v(1:3,1:N_total_atoms) - atomic (i.e. ionic) velocities
    !Output data is returned in three variables:
    !  temper -- temperature of ions [in K]
    !  av_kin -- average kinetic energy of ions, which is the basically
    !                        the same as temper [in eV]
    !  tot_kin -- total kinetic energy of ions [in eV]
    real(8) :: temper,av_kin,tot_kin
    integer :: i,j,elem,atom_count
    real(8) :: mass,f,te
    
    atom_count=0
    te=0.0d0
    do i=1,N_total_atom
      elem=atom_atomic_number(i)
      if (elem/=0) then
        atom_count=atom_count+1
        mass=atom_mass(i)
        do j=1,3
          te=te+mass*atom_velocity(j,i)*atom_velocity(j,i)
        enddo
      endif
    enddo
    te=0.5d0*te
    tot_kin=te
    av_kin=tot_kin/atom_count
    temper=(2.0d0/3.0d0)*av_kin/k_Boltzmann
END SUBROUTINE compute_atomic_temperature

SUBROUTINE compute_atomic_masses
  !This subroutine computes the atomic masses of all atoms present in
  !the system
  integer :: i, elem

  do i=1,N_total_atom
    elem=atom_atomic_number(i)
    if (elem/=0) then
      if (use_average_atomic_mass) then
        atom_mass(i)=element_num_nucleons(elem)*mass_convfactor
      else
        atom_mass(i)=element_average_atomic_mass(elem)*mass_convfactor
      endif
    endif
  enddo
END SUBROUTINE compute_atomic_masses  

SUBROUTINE cleanup    
    deallocate(atom_position)
    deallocate(atom_velocity)
    deallocate(atom_atomic_number)
    deallocate(atom_mass)
END SUBROUTINE cleanup


END MODULE BOLTZMANN