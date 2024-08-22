
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

  ! Set in control_boltzmann.inp
  character(len=255) :: dft_input_path = "./"
  character(len=255) :: velocity_output_path = "./"
  logical :: add_proton_velocity = .FALSE.
  real*8 :: proton_x_velocity = 0.0
  real*8 :: proton_y_velocity = 0.0
  real*8 :: proton_z_velocity = 0.0
  logical :: print_info_to_terminal = .FALSE.




  ! Command Arguments or user input
  real*8  :: temperature_ions=300     ! The temperature of the ions in K. Default value 300
  integer :: ion_velocity_init_seed=1 ! The seed for the calculation. Default value 1 

END MODULE MOD_GLOBAL