! Program that generates a random velocities of a molecule based
! on the Boltzmann distribution given an initial temperature
! in Kelvin and a seed.
! Output velocities in A/fs in an output file called velocity.inp

module marsaglia
!This module contains several very efficient and widely used random number 
!generators by G. Marsaglia (Florida State University). 
!All of the subroutines were taken from public sources and modified 
!by Sergiy Bubin (April-June 2009)
implicit none
private
public :: rkiss, ikiss, kiss_setseed, kiss_getseed, ziggurat_set, rndnor, rndexp, rndunitsph

!Seeds for kiss generator
integer(4),save :: ks_x=123456789
integer(4),save :: ks_y=362436069
integer(4),save :: ks_z=521288629
integer(4),save :: ks_w=916191069

!Parameters, tables, and seeds for random normals and random exponentials
real(8),parameter :: m1=2147483648.0D0   
real(8),parameter :: m2=2147483648.0D0
real(8)           :: dn=3.442619855899D0
real(8)           :: tn=3.442619855899D0
real(8),parameter :: vn=0.00991256303526217D0
real(8)           :: de=7.697117470131487D0
real(8)           :: te=7.697117470131487D0
real(8),parameter :: ve=0.003949659822581572D0
integer(4)        :: ip
integer(4)        :: kn(0:127)
integer(4)        :: ke(0:255)
integer(4)        :: hz
real(8)           :: wn(0:127)
real(8)           :: fn(0:127)
real(8)           :: we(0:255)
real(8)           :: fe(0:255)
logical           :: init_done=.false.

contains


!The  KISS (Keep It Simple Stupid) random number generator of G. Marsaglia. 
!It combines:
!(1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
!(2) A 3-shift shift-register generator, period 2^32-1,
!(3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
!Overall period>2^123;  Default seeds ks_x,ks_y,ks_z,ks_w. 
!The real number generator is called rkiss(), the integer generator is called 
!ikiss(). Both ikiss() and rkiss() share the same set of seeds.
!Set your own seeds with 
!call kiss_setseed(ix,iy,iz,iw).
!Get current seeds with
!call kiss_getseed(ix,iy,iz,iw)


function ikiss()
!Function ikiss generates a random integer. 
integer(4)  ikiss
ks_x = 69069 * ks_x + 1327217885
ks_y = m (m (m (ks_y, 13), - 17), 5)
ks_z = 18000 * iand (ks_z, 65535) + ishft (ks_z, - 16)
ks_w = 30903 * iand (ks_w, 65535) + ishft (ks_w, - 16)
ikiss = ks_x + ks_y + ishft (ks_z, 16) + ks_w
contains
  function m(k, n)
  integer m, k, n
  m = ieor (k, ishft (k, n) )
  end function m
end function ikiss


function rkiss()
!Function rkiss generates a random number between 0 and 1 
!(numbers exactly equal to 0.0D0 and 1.0D0 do not happen)
real(8)  rkiss 
integer(4)  ikiss
ks_x = 69069 * ks_x + 1327217885
ks_y = m (m (m (ks_y, 13), - 17), 5)
ks_z = 18000 * iand (ks_z, 65535) + ishft (ks_z, - 16)
ks_w = 30903 * iand (ks_w, 65535) + ishft (ks_w, - 16)
ikiss = ks_x + ks_y + ishft (ks_z, 16) + ks_w
!rkiss = ikiss * 1/(2^32-1) + 2^31/(2^32-1)
rkiss = ikiss*2.32830643D-10 + 0.500000000116415D0
contains
  function m(k, n)
  integer m, k, n
  m = ieor (k, ishft (k, n) )
  end function m
end function rkiss


subroutine kiss_setseed(ix, iy, iz, iw)
!Sets seeds for kiss generator
integer(4) :: ix, iy, iz, iw
ks_x = ix
ks_y = iy
ks_z = iz
ks_w = iw
end subroutine kiss_setseed


subroutine kiss_getseed(ix, iy, iz, iw)
!Gets current seeds of kiss generator 
integer(4) :: ix, iy, iz, iw
ix = ks_x
iy = ks_y
iz = ks_z
iw = ks_w
end subroutine kiss_getseed


!Marsaglia & Tsang generator for random normals & random exponentials,
!which uses the so called ziggurat algorithm.
!Translated from C by Alan Miller. Modified by Sergiy Bubin.
!Reference:
!Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
!random variables', J. Statist. Software, v5(8).


subroutine ziggurat_set()
!This subroutine initializes tables for the generators
!of random normals and random exponentials. It is called once and 
!automatically inside either rndnor or rndexp, when they are called 
!for the first time.
integer(4)             :: i
real(8)                :: q

!Tables for rndnor
q = vn*exp(0.5D0*dn*dn)
kn(0) = (dn/q)*m1
kn(1) = 0
wn(0) = q/m1
wn(127) = dn/m1
fn(0) = 1.0D0
fn(127) = exp( -0.5D0*dn*dn )
do  i = 126, 1, -1
  dn = sqrt( -2.0D0 * log( vn/dn + exp(-0.5D0*dn*dn ) ) )
  kn(i+1) = (dn/tn)*m1
  tn = dn
  fn(i) = exp(-0.5D0*dn*dn)
  wn(i) = dn/m1
end do

!Tables for rndexp
q = ve*exp( de )
ke(0) = (de/q)*m2
ke(1) = 0
we(0) = q/m2
we(255) = de/m2
fe(0) = 1.0D0
fe(255) = exp( -de )
do  i = 254, 1, -1
  de = -log( ve/de + exp( -de ) )
  ke(i+1) = m2 * (de/te)
  te = de
  fe(i) = exp( -de )
  we(i) = de/m2
end do
init_done = .true.
return

end subroutine ziggurat_set


function rndnor() result(fn_val)
!Function  rndnor generates random normals
!Note that rndnor calls ikiss and rkiss. Hence, a call of rndnor changes
!current seeds for those generators.
real(8)             ::  fn_val
real(8), parameter  ::  r = 3.442620D0
real(8)             ::  x, y

if ( .not. init_done ) call ziggurat_set()
hz = ikiss()
ip = iand( hz, 127 )
if ( abs( hz ) < kn(ip) ) then
  fn_val = hz * wn(ip)
else
  do
    if ( ip == 0 ) then
        do
          x = -0.2904764D0* log(rkiss())
          y = -log(rkiss())
          if( y+y >= x*x ) exit
        end do
        fn_val = r+x
        if ( hz <= 0 ) fn_val = -fn_val
        return
    end if
    x = hz * wn(ip)
    if ( fn(ip) + rkiss()*(fn(ip-1)-fn(ip)) < exp(-0.5D0*x*x) ) then
      fn_val = x
      return
    end if
    hz = ikiss()
    ip = iand( hz, 127 )
    if ( abs( hz ) < kn(ip) ) then
      fn_val = hz * wn(ip)
      return
    end if
  end do
end if
return

end function rndnor


function rndexp() result(fn_val)
!Function  rndnor generates random exponentials
!Note that rndexp calls ikiss and rkiss. Hence, a call of rndexp changes
!current seeds for those generators.
real(8)    ::  fn_val
real(8)    ::  x
integer(4) ::  jz

if( .not. init_done ) call ziggurat_set()
jz = ikiss()
ip = iand( jz, 255 )
if( abs( jz ) < ke(ip) ) then
  fn_val = abs(jz) * we(ip)
  return
end if
do
  if( ip == 0 ) then
    fn_val = 7.69711 - log( rkiss() )
    return
  end if
  x = abs( jz ) * we(ip)
  if( fe(ip) + rkiss()*(fe(ip-1) - fe(ip)) < exp( -x ) ) then
    fn_val = x
    return
  end if
  jz = ikiss()
  ip = iand( jz, 255 )
  if( abs( jz ) < ke(ip) ) then
    fn_val = abs( jz ) * we(ip)
    return
  end if
end do
return

end function rndexp


!Marsaglia's method of picking random points on a surface of a unit sphere
!Reference:
!Marsaglia, G. "Choosing a Point from the Surface of a Sphere." 
!Ann. Math. Stat. 43, 645-646, 1972. 

subroutine rndunitsph(x,y,z)
!Subroutine rndunitsph returns in x,y,z a unit vector of random
!direction. Note that rndunitsph calls rkiss. Hence, a call of rndunitsph 
!changes current seeds for that generator.
real(8)  x,y,z
real(8)  x1,x2,t,r

x1=2*(rkiss()-0.5D0)
x2=2*(rkiss()-0.5D0)
t=x1**2+x2**2
do while (t>1.0D0)
  x1=2*(rkiss()-0.5D0)
  x2=2*(rkiss()-0.5D0)
  t=x1**2+x2**2
enddo
r=2.0D0*sqrt(1.0D0-t)
x=x1*r
y=x2*r
z=1.0D0-2.0D0*t
end subroutine rndunitsph


end module marsaglia


MODULE MOD_GLOBAL
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

  ! PARAMETERS:
  integer,parameter       :: N_elements=118
  ! Conversion of the mass from SI units to units corresponding
  ! to Angstrom,eV,femtosec (called here "nano"):
  ! M(nano) = M(kilogram) / ( q * 10^{-10}),
  !where q is the elementary charge in SI units, q=1.602176487d-19
  real*8,parameter      :: mass_convfactor=1.66053886d0/1.602176487d-02
  real*8,parameter      :: k_Boltzmann=8.617343d-05 !in eV/Kelvin

  ! Approximate masses of elements (number of nucleons in the most common isotope)
  real(8),parameter :: element_num_nucleons(N_elements)=(/ &
    1.0d0,   4.0d0, &                                !H,He
    7.0d0,   9.0d0,  11.0d0,  12.0d0,  14.0d0, &     !Li,Be,B,C,N
    16.0d0,  19.0d0,  20.0d0,  23.0d0,  24.0d0, &     !O,F,Ne,Na,Mg
    27.0d0,  28.0d0,  31.0d0,  32.0d0,  35.0d0, &     !Al,Si,P,S,Cl
    40.0d0,  39.0d0,  40.0d0,  45.0d0,  48.0d0, &     !Ar,K,Ca,Sc,Ti
    51.0d0,  52.0d0,  55.0d0,  56.0d0,  59.0d0, &     !V,Cr,Mn,Fe,Co
    58.0d0,  63.0d0,  64.0d0,  69.0d0,  74.0d0, &     !Ni,Cu,Zn,Ga,Ge
    75.0d0,  80.0d0,  79.0d0,  84.0d0,  85.0d0, &     !As,Se,Br,Kr,Rb
    88.0d0,  89.0d0,  90.0d0,  93.0d0,  98.0d0, &     !Sr,Y,Zr,Nb,Mo
    97.0d0, 102.0d0, 103.0d0, 106.0d0, 107.0d0, &     !Tc,Ru,Rh,Pd,Ag
    114.0d0, 115.0d0, 120.0d0, 121.0d0, 130.0d0, &     !Cd,In,Sn,Sb,Te
    127.0d0, 132.0d0, 133.0d0, 138.0d0, 139.0d0, &     !I,Xe,Cs,Ba,La
    140.0d0, 141.0d0, 142.0d0, 145.0d0, 152.0d0, &     !Ce,Pr,Nd,Pm,Sm
    153.0d0, 158.0d0, 159.0d0, 164.0d0, 165.0d0, &     !Eu,Gd,Tb,Dy,Ho
    166.0d0, 169.0d0, 174.0d0, 175.0d0, 180.0d0, &     !Er,Tm,Yb,Lu,Hf
    181.0d0, 184.0d0, 187.0d0, 192.0d0, 193.0d0, &     !Ta,W,Re,Os,Ir
    195.0d0, 197.0d0, 202.0d0, 205.0d0, 208.0d0, &     !Pt,Au,Hg,Tl,Pb
    209.0d0, 209.0d0, 210.0d0, 222.0d0, 223.0d0, &     !Bi,Po,At,Rn,Fr
    226.0d0, 227.0d0, 232.0d0, 231.0d0, 238.0d0, &     !Ra,Ac,Th,Pa,U
    237.0d0, 244.0d0, 243.0d0, 247.0d0, 247.0d0, &     !Np,Pu,Am,Cm,Bk
    251.0d0, 252.0d0, 257.0d0, 258.0d0, 259.0d0, &     !Cf,Es,Fm,Md,No
    262.0d0, 261.0d0, 262.0d0, 266.0d0, 264.0d0, &     !Lr,Rf,Db,Sg,Bh
    277.0d0, 268.0d0, 281.0d0, 272.0d0, 285.0d0, &     !Hs,Mt,Ds,Rg,Cn
    0.0d0, 289.0d0,   0.0d0, 292.0d0,   0.0d0, &     !Uut,Uuq,Uup,Uuh,Uus
    0.0d0 /)   !Uuo

  ! Average masses of elements (in g/mol)
  real(8),parameter :: element_average_atomic_mass(N_elements)=(/ &
    1.00794d0,  4.002602d0,  6.941d0,  9.012182d0,  10.811d0,           & !  H, He, Li, Be, B
    12.0107d0,  14.0067d0,  15.9994d0,  18.9984032d0,  20.1797d0,       & !  C, N, O, F, Ne
    22.98976928d0,  24.3050d0,  26.9815386d0,  28.0855d0,  30.973762d0, & !  Na, Mg, Al, Si, P
    32.065d0,  35.453d0,  39.948d0,  39.0983d0,  40.078d0,              & !  S, Cl, Ar, K, Ca
    44.955912d0,  47.867d0,  50.9415d0,  51.9961d0,  54.938045d0,       & !  Sc, Ti, V, Cr, Mn
    55.845d0,  58.933195d0,  58.6934d0,  63.546d0,  65.38d0,            & !  Fe, Co, Ni, Cu, Zn
    69.723d0,  72.64d0,  74.92160d0,  78.96d0,  79.904d0,               & !  Ga, Ge, As, Se, Br
    83.798d0,  85.4678d0,  87.62d0,  88.90585d0,  91.224d0,             & !  Kr, Rb, Sr, Y, Zr
    92.90638d0,  95.96d0,  98.9063d0,  101.07d0,  102.90550d0,          & !  Nb, Mo, Tc, Ru, Rh
    106.42d0,  107.8682d0,  112.411d0,  114.818d0,  118.710d0,          & !  Pd, Ag, Cd, In, Sn
    121.760d0,  127.60d0,  126.90447d0,  131.293d0,  132.9054519d0,     & !  Sb, Te, I, Xe, Cs
    137.327d0,  138.90547d0,  140.116d0,  140.90765d0,  144.242d0,      & !  Ba, La, Ce, Pr, Nd
    146.9151d0,  150.36d0,  151.964d0,  157.25d0,  158.92535d0,         & !  Pm, Sm, Eu, Gd, Tb
    162.500d0,  164.93032d0,  167.259d0,  168.93421d0,  173.054d0,      & !  Dy, Ho, Er, Tm, Yb
    174.9668d0,  178.49d0,  180.9479d0,  183.84d0,  186.207d0,          & !  Lu, Hf, Ta, W, Re
    190.23d0,  192.217d0,  195.084d0,  196.966569d0,  200.59d0,         & !  Os, Ir, Pt, Au, Hg
    204.3833d0,  207.2d0,  208.98040d0,  208.9824d0,  209.9871d0,       & !  Tl, Pb, Bi, Po, At
    222.0176d0,  223.0197d0,  226.0254d0,  227.0278d0,  232.03806d0,    & !  Rn, Fr, Ra, Ac, Th
    231.03588d0,  238.02891d0,  237.0482d0,  244.0642d0,  243.0614d0,   & !  Pa, U, Np, Pu, Am
    247.0704d0,  247.0703d0,  251.0796d0,  252.0829d0,  257.0951d0,     & !  Cm, Bk, Cf, Es, Fm
    258.0986d0,  259.1009d0,  264d0,  265d0,  268d0,                    & !  Md, No, Lr, Rf, Db
    272d0,  273d0,  276d0,  279d0,  278d0,                              & !  Sg, Bh, Hs, Mt, Ds
    283d0,  285d0,  287d0,  289d0,  291d0,                              & !  Rg, Cn, Uut, Uuq, Uup
    293d0,  295d0,  294d0 /)                                              !  Uuh, Uus, Uuo

END MODULE MOD_GLOBAL


MODULE MOD_IO
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
  
  write(*,*) 'Atom positions:'
  do i=1, N_total_atom
      write(*,*) atom_position(:,i)
  end do
  
  write(*,*) 'Atom atomic number and mass:'
  do i=1, N_total_atom
    write(*,*) atom_atomic_number(i), atom_mass(i)
  end do
  
  write(*,*) 'Atom velocities:'
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


MODULE BOLTZMANN
  USE MOD_GLOBAL
  USE MARSAGLIA
  USE MOD_IO
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
  ! close files and deallocate memory
  
  deallocate(atom_position)
  deallocate(atom_velocity)
  deallocate(atom_atomic_number)
  deallocate(atom_mass)
END SUBROUTINE cleanup

END MODULE BOLTZMANN


PROGRAM MAIN
  USE BOLTZMANN
  implicit none
  
  call initialize('dft.inp')
  call calculate_atomic_velocities
  call print_to_terminal
  call write_velocity_to_file('velocity.inp')
  
END PROGRAM MAIN