subroutine MatricFLuxInv(MFP0, VGa, VGn, VGl, Ks, href, h)
    
    implicit none

	integer i 
	real*8 VGa, VGn, VGm, VGl, Ks, phi, L 
    real*8 TH, THref, href, MFP, FNl, h, hm, THm, Km
    real*8 dMdTH, MFP1, MFP2, MFP0
    real*8 crit, dTH
    
    VGm = 1-1/VGn  !Mualem
    phi = VGm * (VGl+1)
    
    THref = (1+(VGa*href)**VGn)**(-VGm)
    
    crit = 1e-6
    dTH=.00001

    MFP1=1e8
    
    !First estimate of TH
    TH=.5
    
    do while ((abs((MFP1-MFP0)/MFP0)).gt.crit)
      MFP1 = VGm*(1-VGm)*Ks/2/VGa/(phi+1) * (FnL(TH, VGl, VGm,phi) - FnL(THref, VGl, VGm,phi))  
      MFP2 = VGm*(1-VGm)*Ks/2/VGa/(phi+1) * (FnL(TH+dTH, VGl, VGm,phi) - FnL(THref, VGl, VGm,phi))  
      dMdTH = (MFP2-MFP1)/dTH
      TH = TH - (MFP1-MFP0)/dMdTH
    enddo
    
    h = ((TH**(-1/VGm))-1)**(1/VGn)/VGa
    
    return
    
    end subroutine
    
    ! LAMBDA (returns LAMBDA as function of THETA) 
    real*8 function FnL (TH, VGl, VGm, phi)
	implicit none
	real*8 VGl, VGm, phi, a(3), B(4), TH
    integer i
    
    do i=1, 3
        a(i) = i/VGm + VGl + 1
    enddo
    
    do i=1, 2
        B(i) = (i+phi)*(i+1+VGm)/(i+2)/(i+1+phi)
        B(i+2) = (i+phi)*(i+1-VGm)/(i+2)/(i+1+phi)
    enddo
    
    FnL = 2*VGm*(TH**a(1)) + ((1+VGm)*B(1)-(1-VGm)*B(3))*(TH**a(2)) + ((1+VGm)*B(1)*B(2)-(1-VGm)*B(3)*B(4))*(TH**a(3))
    
    return
    end function
    