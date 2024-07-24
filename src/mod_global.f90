
MODULE MOD_GLOBAL
  !This module contains global variables to be used in the program
  USE PARAMETERS
  implicit none
    
  integer :: N_total_atom
  real*8, dimension(:,:), allocatable :: atom_position
  real*8, dimension(:,:), allocatable :: atom_velocity
  real*8, dimension(:), allocatable   :: atom_mass
  integer, dimension(:), allocatable  :: atom_atomic_number

  logical :: use_average_atomic_mass=.TRUE. 

  ! Command Arguments or user input
  real*8  :: temperature_ions=300     ! The temperature of the ions in K. Default value 300
  integer :: ion_velocity_init_seed=1 ! The seed for the calculation. Default value 1 

END MODULE MOD_GLOBAL
  