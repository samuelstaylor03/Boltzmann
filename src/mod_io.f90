MODULE MOD_IO
  !This module contains subroutines for reading input files and generating output files
  USE MOD_GLOBAL
  implicit none

CONTAINS


SUBROUTINE convert_to_lowercase(string)
  integer          :: i
  character(len=*) :: string
  character        :: the_char

  do i=1,len(string)
    the_char=string(i:i)
    if((iachar(the_char)>=iachar("A")).AND.(iachar(the_char)<=iachar("Z"))) then
      the_char=achar(iachar(the_char)+(iachar('a')-iachar('A')))
    endif
    string(i:i)=the_char
  enddo
END SUBROUTINE convert_to_lowercase


SUBROUTINE read_control_input_file(input_filename)
  character(len=*), intent(in) :: input_filename
  character*256 :: file_line, the_key, value_string
  integer :: i, error_code

  ! Open the file
  open(unit=control_file, file=trim(adjustl(input_filename)), status='old', action='read', iostat=error_code)
  if (error_code /= 0) then
      write(*,*) "***ERROR*** Error opening: ", input_filename
      write(*,*) "***ERROR*** Ensure that ", input_filename, " is in your running directory"
      stop 
  end if

  do
    read(control_file,'(a)',end=110)file_line

    ! Ignore anything on the line that appears after a #
    i=index(file_line,'#')
    if(i>0) then
      file_line=file_line(1:i-1)
    endif

    i=index(file_line,'=')
    if(i==0) then
      cycle
    else
      read(file_line(1:i-1),'(a)')the_key
      read(file_line(i+1:),'(a)')value_string
      call convert_to_lowercase(the_key)

      select case(trim(adjustl(the_key)))
        case("dft_input_path")
          dft_input_path = trim(adjustl(value_string))            
        case("velocity_output_path")
          velocity_output_path = trim(adjustl(value_string))
        case("add_proton_velocity")
          read(value_string,*)add_proton_velocity
        case("proton_x_velocity")
          read(value_string,*)proton_x_velocity
        case("proton_y_velocity")
          read(value_string,*)proton_y_velocity
        case("proton_z_velocity")
          read(value_string,*)proton_z_velocity
        case("print_info_to_terminal")
          read(value_string,*)print_info_to_terminal
        case default
          write(*,*)"ERROR: Invalid variable name: ",trim(adjustl(the_key))
          stop 		   
      end select
    endif
  enddo
  110 close(control_file)

    
END SUBROUTINE read_control_input_file


SUBROUTINE read_input_file(input_filename)
  character(len=*), intent(in) :: input_filename
  integer :: i, error_code, atom_counter
  character(len=256) :: line
  character(len=256) :: input_filename_full
  integer :: unit_num

  input_filename_full = trim(adjustl(dft_input_path)) // trim(adjustl(input_filename))

  ! Open the file
  unit_num = 10
  open(unit=unit_num, file=trim(adjustl(input_filename_full)), status='old', action='read', iostat=error_code)
  if (error_code /= 0) then
      print *, "Error opening ", input_filename_full
      print *, "Ensure that dft.inp is in the directory ", dft_input_path
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
      print *, "Error reading N_total_atom from first line of ", input_filename, " file"
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

  if (add_proton_velocity) then
      write(*,*) proton_x_velocity, proton_y_velocity, proton_z_velocity
  endif

END SUBROUTINE print_to_terminal

SUBROUTINE write_velocity_to_file(output_filename)
  character(len=*), intent(in) :: output_filename
  character(len=255) :: output_filename_full 
  integer :: i, unit_num, error_code, total_number_atoms

  output_filename_full = trim(adjustl(velocity_output_path)) // trim(adjustl(output_filename))

  ! Open the file
  unit_num = 20
  open(unit=unit_num, file=trim(adjustl(output_filename_full)), status='replace', action='write', iostat=error_code)
  if (error_code /= 0) then
      print *, "Error opening velocity output file for writing"
      stop
  end if

  if (add_proton_velocity) then
    total_number_atoms = N_total_atom + 1
  else
    total_number_atoms = N_total_atom
  endif

  write(unit_num, *) total_number_atoms

  ! Write velocities to the file
  do i = 1, N_total_atom
      write(unit_num, *) atom_velocity(1, i), atom_velocity(2, i), atom_velocity(3, i)
  end do

  if (add_proton_velocity) then
      write(unit_num, *) proton_x_velocity, proton_y_velocity, proton_z_velocity
  endif

  ! Close the file
  close(unit_num)

  write(*,'(A,A,A,F8.2,A,I0)') "Finished. Output: ", adjustl(trim(output_filename_full)), &
                               ", temp=", temperature_ions, ", r=", ion_velocity_init_seed
END SUBROUTINE write_velocity_to_file

END MODULE MOD_IO