MODULE MOD_IO
    !This module contains subroutines for reading input files and generating output files
    USE MOD_GLOBAL
    implicit none

CONTAINS

SUBROUTINE read_input_file(input_filename)
    character(len=*), intent(in) :: input_filename
    integer :: i, error_code, atom_counter
    character(len=256) :: line
    integer :: unit_num

    ! Open the file
    unit_num = 10
    open(unit=unit_num, file=trim(adjustl(input_filename)), status='old', action='read', iostat=error_code)
    if (error_code /= 0) then
        print *, "Error opening dft.inp"
        print *, "Ensure that dft.inp is in the same directory as the executable"
        stop
    end if

    ! Read lines and skip those starting with '#'
    do
        read(unit_num, '(A)', iostat=error_code) line
        if (error_code /= 0) exit

        ! Ignore anything on the line that appears after a '#'
        i = index(line, '#')
        if (i > 0) then
            line = line(1:i-1)
        end if

        ! Trim leading and trailing spaces
        line = adjustl(line)
        if (len_trim(line) == 0) cycle  ! Skip empty lines

        ! First non-comment line should be the number of atoms
        read(line, *, iostat=error_code) N_total_atom
        if (error_code == 0) exit
        if (error_code /= 0) then
        print *, "Error reading N_total_atom"
        stop
        end if
    end do



    ! Allocate arrays based on N_total_atom
    allocate(atom_position(3, N_total_atom))
    allocate(atom_velocity(3, N_total_atom))
    allocate(atom_atomic_number(N_total_atom))
    allocate(atom_mass(N_total_atom))

    ! Initialize all arrays to zero
    atom_position = 0.0
    atom_velocity = 0.0
    atom_atomic_number = 0
    atom_mass = 0.0

    atom_counter = 0
    ! Read lines and skip those starting with '#'
    do while (atom_counter < N_total_atom)
        read(unit_num, '(A)', end=110,iostat=error_code) line
        if (error_code /= 0) exit

        ! Ignore anything on the line that appears after a '#'
        i = index(line, '#')
        if (i > 0) then
            line = line(1:i-1)
        end if

        ! Trim leading and trailing spaces
        line = adjustl(line)
        if (len_trim(line) == 0) cycle  ! Skip empty lines

        ! First non-comment line should be the number of atoms
        ! Read the atom positions and atomic numbers
        atom_counter = atom_counter + 1        

        read(line, *, iostat=error_code) atom_position(1, atom_counter), atom_position(2, atom_counter), &
                    atom_position(3, atom_counter), atom_atomic_number(atom_counter)
        if (error_code /= 0) then
            print *, "Error reading atom data at atom", atom_counter
            stop
        end if
    end do

    ! Close the file
    110 close(unit_num)    
END SUBROUTINE read_input_file

SUBROUTINE input_user_terminal
    write(*,*) 'Enter the temperature of the molecule in Kelvin:'
    read(*,*) temperature_ions
    write(*,*) 'Enter an integer value for the seed:'
    read(*,*) ion_velocity_init_seed

END SUBROUTINE input_user_terminal


SUBROUTINE print_to_terminal
    integer :: i

    write(*,*) 'Atom positions (A) in x,y,z:'
    do i=1, N_total_atom
        write(*,*) atom_position(:,i)
    end do

    write(*,*) 'Atom atomic number and mass (eV_fs^2/A^2):'
    do i=1, N_total_atom
        write(*,*) atom_atomic_number(i), atom_mass(i)
    end do

    write(*,*) 'Atom velocities (A/fs) in x,y,z:'
    do i=1, N_total_atom
        write(*,*) atom_velocity(:,i)
    end do

END SUBROUTINE print_to_terminal

SUBROUTINE write_velocity_to_file(output_filename)
    character(len=*), intent(in) :: output_filename
    integer :: i, unit_num, error_code

    ! Open the file
    unit_num = 20
    open(unit=unit_num, file=trim(adjustl(output_filename)), status='replace', action='write', iostat=error_code)
    if (error_code /= 0) then
        print *, "Error opening file for writing"
        stop
    end if

    ! Write velocities to the file
    do i = 1, N_total_atom
        write(unit_num, *) atom_velocity(1, i), atom_velocity(2, i), atom_velocity(3, i)
    end do

    ! Close the file
    close(unit_num)

    print *, "Successfully wrote velocities to file: ", output_filename
END SUBROUTINE write_velocity_to_file

END MODULE MOD_IO