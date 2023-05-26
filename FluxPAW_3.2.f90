!  Flux_PAW - Calculates Plant Available Water (PAW) based on Matric Flux Potential - 2022
!
!  Coded by: Quirijn de Jong van Lier and Marina L A de Melo
!  Revised in March 2022
!  Contact: melo.marina@usp.br; marinaluciana94@gmail.com
!  Main reference: paper doi:10.1016/j.geoderma.2022.116253
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
    
    program FluxPAW_32

    implicit none

	integer i, eof, eof2, NSpa, NDat, Spa, Dat, j
	real*8, dimension(:), allocatable :: VGa, VGn, VGm, VGl, Ks, ThR, ThS
	real*8, dimension(:), allocatable :: zsoil, RLD, qfc, Tp, Wfactor
    character*12, dimension(:), allocatable :: IDsoil
    real*8 TAW, RAW
	real*8 phi, L, p, zfc
    real*8 FnFC, ThFC, hFC, Thl, hl
    real*8 TH, THw, href, MFP, FNl, h, hm, THm, Km
    real*8 dMdTH, MFP1, MFP2, MFP0 
    real*8 crit, dTH, Theta, Ml, Mw, Hw
    character*60 OutFile, Parm, dum, InpFile, InpFile2, OutFile2
    character*1 cm

    cm = ','
    NSpa=0
    NDat=0
  
! Fixed parameter
    href = 20000
    
! Identify input file of Soil-Plant-Atmosphere parameters   
    InpFile = 'Input.txt'
    OPEN (unit=1, FILE = InpFile, STATUS='old')
     do j = 1,6  
       read (1,*) dum
     enddo
      eof=0
        do while (eof.ge.0)
           read (1,'(A60)',iostat=eof) dum

           if (trim(dum).ne."") then
             Nspa=Nspa+1
           endif 
        enddo
        NSpa=NSpa-1
        close(1)

    
! Identify soil input file and read header    
    InpFile2 = 'Soil.txt'   
    
        OPEN (unit=2, FILE = InpFile2, STATUS='old')  ! Input file of soil hydraulic parameters
             do j = 1,6  
               read (2,*) dum
             enddo
            eof = 0
            do while (eof.ge.0)
            read (2,'(A70)',iostat=eof) dum
            if (trim(dum).ne."") then
             Ndat=Ndat+1
            endif 
        enddo
        Ndat=Ndat-1
        close(2)

    allocate (IDsoil(NDat), VGa(NDat), VGn(NDat), VGm(NDat), VGl(NDat), Ks(NDat), ThR(NDat), ThS(NDat))
	allocate (zsoil(NSpa), RLD(NSpa), qfc(NSpa), Tp(NSpa), Wfactor(NSpa))


! Create outout file and write header
    Parm = 'RESULTS'                                           
    OutFile = trim(Parm) // '.txt'
    OPEN (unit=3, FILE = OutFile, STATUS='replace')  ! Output file
       
    OPEN (unit=1, FILE = InpFile, STATUS='old')
     do j = 1,6  
      read (1,*) dum
     enddo
      do Spa = 1, NSpa
        read (1,*) zsoil(Spa), RLD(Spa), qfc(Spa), Tp(Spa), Wfactor(Spa)
        ! Header of RESULTS.txt
           if (Spa==1) then
            write (3,'(a109)') '* Calculated water contents (Th), pressure heads (h) and ranges of water availability to plants (RAW and TAW)'  
            write (3,'(a65)') '* FC is field capacity; LP is limiting point; WP is wilting point'
            write (3,'(a1)')'*'   
            write (3,'(a35,i1,a5,f7.2,2(a8,f5.3),a7,f5.3,a12,f5.3)') '*** Results for SPA parameters set ',Spa,': z =',zsoil(Spa),&
                &', RLD = ',RLD(Spa),', qFC = ',qFC(Spa),', Tp = ',Tp(Spa),', Wfactor = ',Wfactor(Spa)
            write (3,'(a109)') '*************************************************************************************************************'
           endif           
      enddo
      close(1)

      write (3,'(9a12)') 'ID_soil     ','ThFC   ', 'ThLP   ', 'ThWP   ', 'hFC(cm)', 'hLP(cm)  ', 'hWP(cm)  ', 'RAW    ', 'TAW    '
      
    OPEN (unit=2, FILE = InpFile2, STATUS='old')  ! Input file of soil hydraulic parameters
     do j = 1,6  
      read (2,*) dum
     enddo
      do Dat = 1, NDat
        read (2,*) IDsoil(Dat), ThR(Dat), ThS(Dat), VGa(Dat), VGn(Dat), VGl(Dat), Ks(Dat)
      enddo
      close(2)
    
        do Spa = 1, NSpa             
          
          do Dat = 1, NDat
                 ThFC = FnFC (VGa(Dat), VGl(Dat), VGn(Dat), Ks(Dat), qfc(Spa), zsoil(Spa))
                   if (ThFC.ge.0.9999) then
                     ThFC = 0.9999 
                   else
                    if (ThFC.le.0.0001) then
                     ThFC = 0.0001   
                     endif    
                   endif
                 hfc = (ThFC**(VGn(Dat)/(1-VGn(Dat))) - 1)**(1/VGn(Dat)) / VGa(Dat)  
                 ThFC = ThR(Dat) + (ThS(Dat)-ThR(Dat)) * ThFC 
               
                  Ml = 1.69 * (Tp(Spa)/10.) / Zsoil(Spa) / RLD(Spa)   ! ***** March 2023: units correction
                  Mw = Ml * Wfactor(Spa)    
        
                 call MatricFLuxInv(Ml, VGa(Dat), VGn(Dat), VGl(Dat), Ks(Dat), href, hl)
                 call MatricFLuxInv(Mw, VGa(Dat), VGn(Dat), VGl(Dat), Ks(Dat), href, hw)
                 ThL = ThR(Dat) + (ThS(Dat)-ThR(Dat)) * (1+(VGa(Dat)*hl)**VGn(Dat))**((1-VGn(Dat))/VGn(Dat))
                 ThW = ThR(Dat) + (ThS(Dat)-ThR(Dat)) * (1+(VGa(Dat)*hw)**VGn(Dat))**((1-VGn(Dat))/VGn(Dat))
    
                 if (ThL .gt. ThFC) then
                     ThL = ThFC
                     hl = hfc
                 else
                 endif
    
                 RAW = ThFC - ThL
                 TAW = ThFC - ThW
                 
                 write (3,'(a12,3f12.5,3f12.3,2f12.5)') IDsoil(Dat), ThFC, Thl, Thw, hFC, hl, hw, RAW, TAW 
          enddo

         if (NSpa .eq. Spa) then
         else
            write (3,'(a109)') '*************************************************************************************************************'
            do i=1,2
             write (3,'(a1)')'*'
            enddo
            write (3,'(a35,i1,a5,f7.2,2(a8,f5.3),a7,f5.3,a12,f5.3)') '*** Results for SPA parameters set ',Spa+1,': z =',zsoil(Spa+1),&
                &', RLD = ',RLD(Spa+1),', qFC = ',qFC(Spa+1),', Tp = ',Tp(Spa+1),', Wfactor = ',Wfactor(Spa+1)
            write (3,'(a109)') '*************************************************************************************************************'
            write (3,'(9a12)') 'ID_soil     ','ThFC   ', 'ThLP   ', 'ThWP   ', 'hFC(cm)', 'hLP(cm)  ', 'hWP(cm)  ', 'RAW    ', 'TAW    '
         endif
       
        enddo

        write (3,'(a17)') '**** end of file!'
        
     write (*,*) 'Calculation complete. RESULTS.txt created.'
     write (*,*) 'Press enter to exit'
     
     read*

    end program FluxPAW_32
