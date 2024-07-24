! Program that generates the 'velocity.inp' file for a given molecule. 
! Generates random velocities (A/fs) of a molecule based on the Boltzmann distribution 
! given an initial temperature (in Kelvin) and an initial seed
! by Samuel S. Taylor (July 2024)

PROGRAM MAIN
    USE BOLTZMANN
    implicit none
    
    call initialize('dft.inp')
    call calculate_atomic_velocities
    call print_to_terminal
    call write_velocity_to_file('velocity.inp') 
END PROGRAM MAIN