!  FluxPAW - Calculates Plant Available Water (PAW) based on Matric Flux Potential - 2021
!
!  Coded by: Quirijn de Jong van Lier and Marina L A de Melo
!  Revised in Feb 2022
!  Contact: melo.marina@usp.br; marinaluciana94@gmail.com
!
!  FUNCTIONS:
!
!  MatricFlux - Calculates M from Van Genuchten - Mualem parameters using a solution by 
!               De Jong van Lier et al., 2009 (doi:10.1029/2008WR006938)
! 
!  FnFC - Calculates Theta-FC according to Inforsato & De Jong van Lier, 2021 
!         (doi:10.1016/j.geoderma.2021.115308)
!
!***************************************************************************************************************
!***************************************************************************************************************

    program FluxPAW
    implicit none

	integer i 
	real*8 VGa, VGn, VGm, VGl, Ks, RAW, TAW
	real*8 phi, L, p, qfc, zfc
    real*8 FnFC, ThFC, hFC, Thl, hl
    real*8 TH, THw, href, MFP, FNl, h, hm, THm, Km
    real*8 dMdTH, MFP1, MFP2, MFP0
    real*8 crit, dTH, ThR, ThS, Theta, Zsoil, Tp, RLD, Wfactor, Ml, Mw, Hw
    character*20 OutFile, Parm, dum, InpFile
   
   InpFile = 'Input.inp'
   OPEN (unit=1, FILE = InpFile , STATUS='old')  ! Input file
   read (1,*) dum
   read (1,*) dum
   read (1,*) zsoil, ThR, ThS, VGa, VGn, VGl, href, RLD, ks, qfc, Tp, Wfactor
    
    Parm = 'Results'                                           
    OutFile = trim(Parm) // '.out'

    OPEN (unit=2, FILE = OutFile , STATUS='replace')  ! Output file
    write (2, '(8A12)') 'ThFC', 'Thl', 'Thw', 'hFC', 'hl', 'hw', 'RAW', 'TAW'

    ThFC = FnFC (VGa, VGl, VGn, Ks, qfc, zsoil)
    hfc = (ThFC**(VGn/(1-VGn)) - 1)**(1/VGn) / VGa
    ThFC = ThR + (ThS-ThR) * ThFC

    Ml = 1.69 * (Tp/1000.) / Zsoil / (RLD*1.e4)    
    Mw = Ml * Wfactor    
        
    call MatricFLuxInv(Ml, VGa, VGn, VGl, Ks, href, hl)
    call MatricFLuxInv(Mw, VGa, VGn, VGl, Ks, href, hw)
    ThL = ThR + (ThS-ThR) * (1+(VGa*hl)**VGn)**((1-VGn)/VGn)
    ThW = ThR + (ThS-ThR) * (1+(VGa*hw)**VGn)**((1-VGn)/VGn)
    
    if (ThL .gt. ThFC) then
        ThL = ThFC
        hl = hfc
    else
    endif
    
    RAW = ThFC - ThL
    TAW = ThFC - ThW
    
    write (2, '(8F12.5)') ThFC, Thl, Thw, hFC, hl, hw, RAW, TAW 
        
    write (*,*) 'Calculation complete!'

    close(1)
    close(2)
    
    read*
    
    end program FluxPAW
    
 !***************************************************************************************************************
 !***************************************************************************************************************

