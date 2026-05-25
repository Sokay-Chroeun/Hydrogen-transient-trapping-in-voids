	 ! UMATHT : for hydrogen transport and trapping in voids
	 ! UMAT : for linear elastic behaviour for 3D element
	 ! SDV1 : CL
	 ! SDV2 : CT
	 ! SDV3 : CH
	 ! SDV4 : Nb
	 ! SDV5 : P
	 !	  
	  subroutine umatht(u,dudt,dudg,flux,dfdt,dfdg,statev,temp,dtemp,
     1 dtemdx,time,dtime,predef,dpred,cmname,ntgrd,nstatv,props,nprops,
     2 coords,pnewdt,noel,npt,layer,kspt,kstep,kinc)

      include 'aba_param.inc'

      character*80 cmname
      dimension dudg(ntgrd),flux(ntgrd),dfdt(ntgrd),dfdg(ntgrd,ntgrd),
     1 statev(nstatv),dtemdx(ntgrd),time(2),predef(1),dpred(1),
     2 props(nprops),coords(3)

	  real*8 CL, CT, dCT, dudt2, CL0, xT, D, xR, Na, xkb
	  real*8 ks, kr
	  real*8 xIntP, rb, Nb, dNb, Cb, Vb, fTotal, xVm, xfug, xVm_h, xIntP_h, xfug_h
	  real*8 n_old, n_guess, Tol, fn, gn, hh, n_h, fn_h, df_dn, gn_prime, n_new
	  integer Iter
	  
	  ! Normalized parameter
	  CL0 = 1.
	  ! Temperature
	  if (time(2).le.200.) then
		xT = 293.
	  else
		xT = 0.625*time(2) + 168.
	  endif
	  ! Universal constants
	  xR = 8.314          ! J/mol.K
	  Na = 6.022e23       ! H/mol
	  xkb = 8.62e-5		  ! eV/K
	  ! Model parameters
	  ks = 4.4e17*exp(-0.32/xkb/xT)
	  kr = 2e-14*exp(-0.19/xkb/xT)
	  D = 7.2e-2*exp(-0.058/xkb/xT)
	  ! Total porosity & equivalent void radius
	  fTotal = 4e-4
	  rb = 8e-3
	  ! Assign lattice H concentration & transport flux
	  CL = temp + dtemp
	  do i=1,ntgrd
       dudg(i)= 0.
       flux(i)= -D*dtemdx(i)
       dfdt(i)= 0.
       dfdg(i,i)=-D
      end do
	  ! Initialisation
	  Vb = 4./3.*3.1416*rb**3.
	  Cb = fTotal/Vb
      Nb = max(statev(4), 1e-16)
	  xIntP = statev(5)
	  ! Start Newton-raphson to find Nb
	  n_old = Nb
      n_guess = n_old
      Iter = 0
      Tol  = 1.
	  do while (Iter.le.500.and.Tol.gt.1e-4)
		Iter = Iter + 1
		! Fugacity and pressure at n_guess
		xVm = max(1., (8. * 3.1416 * rb **3 * Na) / (3. * n_guess) * 1e-3)  ! cm3/mol
		if (xVm.lt.(8.314*xT/350. + 14.598)) then
			xIntP = P_Tkacz(xT, xVm)
		else
			xIntP = 8.314*xT/(xVm - 14.598)
		endif
		xfug = calcul_fug(xT, xIntP)
		! compute f(n) and g(n)
		fn = (4. * 3.1416 * rb**2) * (kr * (CL*CL0)**2 - kr * ks**2 * xfug)
		gn = n_guess - n_old - fn * dtime
		! compute numerical derivative of fn : df/dn
		hh = n_guess * 1e-5
		n_h = n_guess + hh
		xVm_h = max(1., (8. * 3.1416 * rb **3 * Na) / (3. * n_h) * 1e-3)  ! cm3/mol
		if (xVm_h.lt.(8.314*xT/350. + 14.598)) then
			xIntP_h = P_Tkacz(xT, xVm_h)
		else
			xIntP_h = 8.314*xT/(xVm_h - 14.598)
		endif
		xfug_h = calcul_fug(xT, xIntP_h)
		fn_h = (4. * 3.1416 * rb**2) * (kr * (CL*CL0)**2 - kr * ks**2 * xfug_h)
		df_dn = (fn_h - fn) / hh
		! compute g'(n)
		gn_prime = 1. - df_dn * dtime
		if (gn_prime.eq.0) write(*,*) "cannot find the solution: g'(n) = 0 "
		! Newton update
		n_new = n_guess - gn / gn_prime
		Tol = abs(n_new - n_guess) / (n_guess + 1e-16)
		!
		n_guess = n_new
	  enddo
	  if (Iter.gt.500) write(*,*) ' cannot find the solution: max iteration reaches '
	  ! end Newton-raphson
	  ! update state: Nb, xIntP, CT
	  Nb = n_new
	  dNb = Nb - statev(4)
	  xVm = max(1., (8. * 3.1416 * rb **3 * Na) / (3. * Nb) * 1e-3)  ! cm3/mol
	  if (xVm.lt.(8.314*xT/350. + 14.598)) then
		xIntP = P_Tkacz(xT, xVm)
	  else
		xIntP = 8.314*xT/(xVm - 14.598)
	  endif	
	  xfug = calcul_fug(xT, xIntP)
	  !CT = Nb*Cb/CL0
	  CT = 3.*fTotal/(4.*3.1416*(1. - fTotal)) * Nb/rb**3 /CL0
	  dCT = CT - statev(2)
	  ! Jacobian & internal energy 
	  dudt2 = Cb*(4.*3.1416*rb**2)*kr*( 2.*(CL*CL0) )*dtime
	  dudt = 1. + dudt2
	  u = u + dtemp + dCT
	  ! 
	  statev(1) = CL
	  statev(2) = CT
	  statev(3) = CL + CT
	  
	  statev(4) = Nb
	  statev(5) = xIntP
	  statev(6) = Cb
	  statev(7) = fTotal
	  statev(8) = xT
	  statev(9) = ks
	  statev(10) = kr
	  
      return
      end

! ------------------------------------------------------------------

      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1 drplde,drpldt,stran,dstran,time,dtime,temp2,dtemp,predef,dpred,
     2 cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      include 'aba_param.inc'

      character*8 cmname
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1 ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),
     2 time(2),predef(1),dpred(1),props(nprops),coords(3),drot(3,3),
     3 dfgrd0(3,3),dfgrd1(3,3),jstep(4)
      
	  real*8 E, xnu, eg, elam
	  
      E = 210e3
      xnu = 0.3
	  ddsdde=0.
      
      eg=E/(1.d0+xnu)/2.
      elam=(E/(1.-2.*xnu)-2.*eg)/3.
      
      do i=1,3
       do j=1,3
        ddsdde(j,i)=elam
       end do
       ddsdde(i,i)=2.0*eg+elam
      end do
      do i=4,ntens
       ddsdde(i,i)=eg
      end do

      stress=stress+matmul(ddsdde,dstran)      

      return
      end

! ------------------------------------------------------------------	
  
	  real*8 function P_Tkacz(xT, xV)
	  real*8 xT, xV, xP, Err, Tol, fp, dfdp
	  real*8 A,B,C,D,E,R,bb
	  integer iter, iterMax
	  
	  A = 176.330 
	  B = -633.675
	  C = -304.574
	  D = 731.393
	  E = 8.59805
	  R = 8.314
	  bb = 14.598
	  
	  if (xT.le.0.) then
	    write(*,*) '----------- xT = 0 K in calcul_Pression  -----------', xT
	    stop
	  endif

	  iterMax = 400
	  Tol = 1e-4
	  iter = 0
	  Err = 1.
	  xP = 350.
	  do while (abs(Err).gt.Tol)
	      if(xP.le.0.) write(*,*) ' p <= 0 '
		  iter = iter + 1
		  fp = A*xp**(-0.3333) + B*xp**(-0.6667) + (D+E*xT)/xp + C*xp**(-1.3333)  - xV
		  dfdp = -0.3333*A*xp**(-1.3333) - 0.6667*B*xp**(-1.6667) - 1.3333*C*xp**(-2.3333) - (D+E*xT)*xp**(-2.0)
		  Err = fp/dfdp
		  xP = xP - Err
		  if (iter.ge.iterMax) write(*,*) 'Maximum iteration reached for xV = ', xV, ' and T = ', xT, ' and P = ', xP
		  if (iter.ge.iterMax) stop
	  enddo
	  P_Tkacz = xP
	  
	  return
	  end

! ------------------------------------------------------------------

	  real*8 function calcul_fug(xT, xP)
	  real*8 xT, xV, xP, fug
	  real*8 A,B,C,D,E,R,bb
	  
	  A = 176.330 
	  B = -633.675
	  C = -304.574
	  D = 731.393
	  E = 8.59805
	  R = 8.314
	  bb = 14.598
	  
	  if (xT.le.0.) then
	    write(*,*) '-----------  xT = 0 K in calcul_fug  -----------', xT
	    stop
	  endif
      
	  if (xP.le.0.) xP = 1e-6
	  if (xP.lt.350.) then
		fug = xP * exp(bb * xP / (R*xT))
	  else
	    xV = 1.5*A*xP**(0.6667) + 3.0*B*xP**(0.3333) + (D+E*xT)*log(xP) - 3.0*C*xP**(-0.3333) + 0.22*R*(550. - xT)
		fug = exp(xV / (R * xT) )
	  endif
	  calcul_fug = fug
	  
	  return
	  end