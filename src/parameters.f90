MODULE PARAMETERS
    ! contains constant parameters used in various calculations throughout the program
    implicit none

    integer,parameter       :: N_elements=118
    !Conversion of the mass from SI units to units corresponding
    !to Angstrom,eV,femtosec (called here "nano"):
    !M(nano) = M(kilogram) / ( q * 10^{-10}),
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
END MODULE PARAMETERS  