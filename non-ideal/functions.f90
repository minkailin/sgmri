SUBROUTINE F (NEQ, x, Y, YDOT, RPAR, IPAR)
  use global
  implicit none
  integer :: NEQ, IPAR
  real*8 :: x, Y(NEQ), YDOT(NEQ), RPAR
  real*8 :: logrho, phidot
  
  logrho = Y(1)
  phidot = Y(2)
  
  YDOT(1) = -omegaz2*x - phidot
  YDOT(2) = exp(logrho)/bigQ
  
  return
end subroutine F

subroutine jac
  use global
  implicit none
  return
end subroutine jac

subroutine setup_isosg
  use global
  implicit none
  integer, parameter :: NEQ = 2, ITOL = 1 , MF = 10, LIW = 60, IOPT = 0
  integer, parameter :: LRW = 2*(20 + 16*NEQ)
  integer :: i, ITASK, ISTATE, IWORK(LIW), IPAR
  real*8, parameter :: RTOL = 1d-15, ATOL=1d-15
  real*8 :: Y(NEQ), x, xOUT, RWORK(LRW), RPAR
  external :: F, JAC
  
  Y(1) = 0d0 !log rho_hat = 0 at z=0
  Y(2) = 0d0 !dphi/dz =0 at z=0
  
  dnorm(1)   = 1d0
  drho(1)    = 0d0
  
  ITASK = 1
  ISTATE= 1
  
  do i = 1, nz-1
     x    = zaxis(i)
     xOUT = zaxis(i+1)
     
     call DVODE (F, NEQ, Y, x, xOUT, ITOL, RTOL, ATOL, ITASK, &
          ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF, &
          RPAR, IPAR)
     
     if(ISTATE .ne. 2) then
        print*, 'setup unsuccessful'
        stop
     endif
     
     dnorm(i+1) = exp(Y(1))
     drho(i+1)  = -omegaz2*zaxis(i+1) - Y(2)

  enddo
  d2rho= ( drho**2d0 - (omegaz2 + dnorm/bigQ)/fQ**2d0 )
  valf = 1d0/(beta*dnorm)
  csq  = 1d0 
  return
end subroutine setup_isosg

subroutine get_basic
  use global
  implicit none
  integer :: i 
  real*8, external :: density, dlogrho, resistivity, dresistivity, &
       d2resistivity 
  
  if(eos .eq. 'pol') then
     fQ      = sqrt(bigQ)*acos(omegaz2*bigQ/(1d0 + omegaz2*bigQ)) !H*Omega/cs0 
     fQ      = 1d0/fQ 
  else if((eos.eq.'iso').or.(eos.eq.'uns')) then
     fQ = 1d0
  endif
  
  if(eos .eq. 'iso') then !isothermal 
     call setup_isosg
  else if(eos .eq. 'pol') then !n=1 polytrope 
     do i=1, nz 
        dnorm(i)   = density(zaxis(i)) 
        drho(i)    = dlogrho(zaxis(i)) 
     enddo
     d2rho = -(omegaz2 + dnorm/bigQ)/fQ**2d0
     d2rho = d2rho/dnorm 
     valf  = 1d0/(beta*dnorm) 
     csq   = dnorm  
     
  else if(eos .eq. 'uns') then !unstratified 
     dnorm = 1d0
     drho  = 0d0 
     csq   = 1d0
     valf  = 1d0/beta 
  endif
  
  if(beta .lt. 0d0) then !unmagnetized
     valf = 0d0 
  endif
  
   call fill_resistivity 

end subroutine get_basic


real*8 function density(znorm)
  !density normalized to midplane density, as function of normalized height
  !this is also cs^2/cs0^2 for polyn=1
  use global
  implicit none
  real*8 :: ahat, omegaz2_Q, znorm 
  
  omegaz2_Q = omegaz2*bigQ
  ahat      = acos( omegaz2_Q/(omegaz2_Q + 1d0) )
  density   = (1d0 + omegaz2_Q)*cos(znorm*ahat) - omegaz2_Q
  return
end function density
 
real*8 function dlogrho(znorm)
  !dlog(rho)/dz (normalized units) 
  use global
  implicit none
  real*8 :: ahat, omegaz2_Q, znorm 
  real*8, external :: density 
  
  omegaz2_Q = omegaz2*bigQ
  ahat = acos( omegaz2_Q/(omegaz2_Q + 1d0) )
  dlogrho = -ahat*(1d0 + omegaz2_Q)*sin(znorm*ahat)/density(znorm) 
  return
end function dlogrho

subroutine fill_resistivity!based on density field. only iso for now 
  use global
  implicit none
  integer :: i 
  real*8 :: surf, surf0, gplus, gminus, gsum, gdiff
  real*8 :: eta0 , surf_star 
  real*8, external :: density, dlogrho
  
  eta0 = fQ**2d0/(Rm*beta)

  if(eos.eq.'iso') then
     surf0 = -bigQ*( fQ*fQ*drho(nz) + omegaz2*zmax )
  else if(eos.eq.'pol') then
     surf0 = -bigQ*( fQ*fQ*( dnorm(nz)*drho(nz) ) + omegaz2*zmax )
  endif

  if(amp.gt.1d0) then
     zstar = acosh(amp**2d0)
     zstar = 1d0/zstar 
     
     surf_star = zstar*surf0
     
     do i=1, nz
        
        if(eos.eq.'iso') then
           gplus  = -bigQ*( fQ**2d0*(drho(nz) - drho(i)) + omegaz2*(zmax - zaxis(i))  )
           gminus = -bigQ*( fQ**2d0*(drho(i) + drho(nz)) + omegaz2*(zaxis(i) + zmax)  )
        else if(eos.eq.'pol') then
           gplus  = -bigQ*( fQ**2d0*(dnorm(nz)*drho(nz) - dnorm(i)*drho(i) ) + omegaz2*(zmax - zaxis(i))  )
           gminus = -bigQ*( fQ**2d0*(dnorm(i)*drho(i) + dnorm(nz)*drho(nz)) + omegaz2*(zaxis(i) + zmax)  )
        endif

        gplus  = (gplus/surf0 - 1d0 )/zstar 
        gminus = (gminus/surf0 - 1d0)/zstar 
        
        gsum  = exp(-gplus) + exp(-gminus)
        gdiff = exp(-gplus) - exp(-gminus)
        
        eta(i)  = eta0*sqrt(2d0)/sqrt(gsum)
        
        deta(i) = -(eta0/(sqrt(2d0)*surf_star))*dnorm(i)*gsum**(-1.5d0)*gdiff 
        
        d2eta(i) = -(eta0/(sqrt(2d0)*surf_star))*dnorm(i)*gsum**(-1.5d0)
        d2eta(i) = d2eta(i)*( drho(i)*gdiff - 1.5d0*(dnorm(i)/surf_star)&
             &*gdiff**2d0/gsum + (dnorm(i)/surf_star)*gsum) 
        
     enddo
  endif
  if(amp.eq.1d0) then
     eta  = eta0
     deta = 0d0 
     d2eta= 0d0
  endif
  return
end subroutine fill_resistivity

real*8 function resistivity(znorm)
  use global
  implicit none
  real*8 :: znorm

  resistivity = fQ**2d0/(Rm*beta)
  return
end function resistivity

real*8 function dresistivity(znorm)
  use global
  implicit none
  real*8 :: znorm

  dresistivity = 0d0 
  return
end function dresistivity

real*8 function d2resistivity(znorm)
  use global
  implicit none
  real*8 :: znorm

  d2resistivity = 0d0 
  return
end function d2resistivity


subroutine chebyshev_poly(l, zbar, T_l, dT_l, d2T_l)
  implicit none
  real*8, intent(in) :: l, zbar
  real*8, intent(out):: T_l, dT_l, d2T_l
  real*8 :: t, lt, lsq

  lsq = l*l

  t  = acos(zbar)
  lt = l*t

  T_l = cos(lt)

  if(abs(zbar).lt.1d0) then
     dT_l = l*sin(lt)/sin(t)
     d2T_l= -lsq*cos(lt)/sin(t)**2 + l*cos(t)*sin(lt)/sin(t)**3
  else
     dT_l = lsq
     d2t_l =lsq*(lsq-1d0)/3d0
  endif

  return
end subroutine chebyshev_poly

subroutine construct_matrix(big_matrix) 
  use global
  implicit none 
  integer :: i 
  complex*16 :: L11(nz, nz), L12(nz, nz), L13(nz, nz), L14(nz, nz), &
       L15(nz, nz) 
  complex*16 :: L21(nz, nz), L22(nz, nz), L23(nz, nz), L24(nz, nz), &
       L25(nz, nz) 
  complex*16 :: L31(nz, nz), L32(nz, nz), L33(nz, nz), L34(nz, nz), &
       L35(nz, nz) 
  complex*16 :: L41(nz, nz), L42(nz, nz), L43(nz, nz), L44(nz, nz), &
       L45(nz, nz)
  complex*16 :: L51(nz, nz), L52(nz, nz), L53(nz, nz), L54(nz, nz), &
       L55(nz, nz)
  complex*16, intent(out) :: big_matrix(big_nz, big_nz) 
  complex*16 :: sbar , dsbar, isig, ifkx, c2, c1, c0, Dhat(nz), Dhattil(nz)
  complex*16 :: c2til, c1til, c0til
  real*8 :: kxsq, fkx, rho, dden, fvasq, vasq, kappa2
  real*8 :: Diff(nz)
  real*8 :: D2(nz), D1(nz), D0(nz), D2til(nz), D1til(nz), D0til(nz) 
  
  L11 = 0d0 ; L12 = 0d0 ; L13 = 0d0 ; L14 = 0d0 ; L15 = 0d0

  L21 = 0d0 ; L22 = 0d0 ; L23 = 0d0 ; L24 = 0d0 ; L25 = 0d0
  
  L31 = 0d0 ; L32 = 0d0 ; L33 = 0d0 ; L34 = 0d0 ; L35 = 0d0

  L41 = 0d0 ; L42 = 0d0 ; L43 = 0d0 ; L44 = 0d0 ; L45 = 0d0
 
  L51 = 0d0 ; L52 = 0d0 ; L53 = 0d0 ; L54 = 0d0 ; L55 = 0d0

  isig  = ii*sig 
  fkx   = fQ*kx
  ifkx  = ii*fkx 
  kxsq  = kx*kx 
  kappa2= kappa*kappa

  do i=1, nz-1
     rho  = dnorm(i)  
     dden = drho(i) !this is dlogrho/dz   
     vasq = valf(i)
     fvasq = fQ**2d0*vasq
     
     sbar  = sig - ii*eta(i)*kxsq
     dsbar = -ii*deta(i)*kxsq 
     
!     Diff(:) = d2rho(i)*T(i,:) + 2d0*drho(i)*Tp(i,:) + Tpp(i,:)

!EQ1, operators on vx 
!     L11(i,:) = fvasq*Tpp(i,:) + (sbar*sig - fvasq*kxsq)*T(i,:) &
!          +ii*eta(i)*Diff(:)*sig 
!EQ1, operators on vy 
!     L12(i,:) = -2d0*eta(i)*Diff(:) + 2d0*ii*sbar*T(i,:)
!EQ1, operators on w
!     L13(i,:) = ifkx*eta(i)*Diff(:) + fkx*sbar*T(i,:)
!EQ1, operators on phi
!     L14(i,:) = L13(i,:)


!EQ2, operators on vx 
!     L21(i,:) = 0.5d0*eta(i)*kappa2*Diff(:)*sig/shear &
!          +( ii*kxsq*fvasq - kxsq*sig*eta(i) - isig*sig -isig*sbar*kappa2/(2d0*shear) )*T(i,:)
!EQ2, operators on vy
!     L22(i,:) = fvasq*Tpp(i,:)*sig/shear + isig*eta(i)*Diff(:)*sig/shear &
!          +( 2d0*sig + sbar*sig*sig/shear - 2d0*ii*kxsq*eta(i) )*T(i,:)
!EQ2, operators on w
!     L23(i,:) = -ifkx*sbar*T(i,:)
!EQ2, operators on phi
!     L24(i,:) = L23(i,:)                 
     
     D0 = T(i,:)
     D1 = drho(i)*T(i,:) + Tp(i,:)
     D2 = d2rho(i)*T(i,:)+2d0*drho(i)*Tp(i,:)+Tpp(i,:)

     D2til = d2eta(i)*D0 + 2d0*deta(i)*D1 + eta(i)*D2
     D1til = deta(i)*D0 + eta(i)*D1
     D0til = eta(i)*D0 
 
     Dhat = D2til - kxsq*D0til - isig*D0

     L11(i,:) = fvasq*Tpp(i,:) - fvasq*kxsq*T(i,:) &
          +isig*Dhat
     L12(i,:) = -2d0*Dhat
     L13(i,:) = ifkx*Dhat
!     L11(i,:) = fvasq*Tpp(i,:) - fvasq*(sig*dsbar/sbar**2d0 + c1til&
!          &*ii*kxsq/sbar)*Tp(i,:) - fvasq*kxsq*T(i,:) + isig*sig*Dhat/sbar&
!          & + sig*kxsq*Dhattil/sbar 
!     L12(i,:) = -2d0*sig*Dhat/sbar + 2d0*ii*kxsq*Dhattil/sbar 
!     L13(i,:) = ifkx*(sig*Dhat/sbar - ii*kxsq*Dhattil/sbar)
     L14(i,:) = L13(i,:) 

     
     Dhat = D2til -dsbar*D1til/sbar
 
     L21(i,:)= ii*fvasq*Tpp(i,:) + sig*(0.5d0*kappa2/shear - &
          1d0)*Dhat - ii*fvasq*dsbar*Tp(i,:)/sbar &
          -0.5d0*ii*sbar*kappa2*T(i,:)*sig/shear
     L22(i,:)= sig*fvasq*Tpp(i,:)/shear + ii*(sig*sig/shear - &
          2d0)*Dhat - (dsbar/sbar)*fvasq*Tp(i,:)&
          &*sig/shear + sbar*sig*sig*T(i,:)/shear 
     L23(i,:) = -fkx*Dhat
!     L21(i,:)= 0.5d0*kappa2*Dhat*sig/shear + ii*kxsq*fvasq*T(i,:) -&
!          & ii*dsbar*fvasq*Tp(i,:)/sbar + isig*Dhattil 
!     L22(i,:)= fvasq*Tpp(i,:)*sig/shear - fvasq*dsbar*T(i,:)*sig/(sbar*shear) + isig*sig*Dhat/shear &
!               -2d0*Dhattil
!     L23(i,:)= ifkx*Dhattil
     L24(i,:)= L23(i,:) 


!EQ3
     L31(i,:) = (sig*kx/fQ)*T(i,:)
     L33(i,:) = Tpp(i,:) + dden*Tp(i,:) + ( (  sig**2d0/csq(i)   &  
          + rho/(csq(i)*bigQ))/fQ**2d0 )*T(i,:) 
     L34(i,:) = dden*Tp(i,:) + kxsq*T(i,:) 
     
!EQ4
     L43(i,:) = -T(i,:)*rho/(csq(i)*fQ**2d0*bigQ)
     L44(i,:) = Tpp(i,:) - kxsq*T(i,:)
     
  enddo

!upper disk boundary conditions
  i = nz 
  rho  = dnorm(i)
  dden = drho(i)
  vasq = valf(i)
  fvasq = fQ**2d0*vasq

  if((eos.eq.'pol').or.(eos.eq.'iso')) then
!     L11(i,:) = Tp(i,:) !dvx/dz = 0
!     L22(i,:) = Tp(i,:) !dvy/dz = 0 
!bx=0
!     L11(i,:) = (fvasq + isig*eta(i))*Tp(i,:) + isig*drho(i)*eta(i)&
!          &*T(i,:)
!     L12(i,:) = -2d0*eta(i)*( drho(i)*T(i,:) +Tp(i,:) )
!     L13(i,:) = ifkx*eta(i)*(drho(i)*T(i,:) + Tp(i,:) )
!     L14(i,:) = L13(i,:)
!by=0
!     L21(i,:) = 0.5*kappa2*eta(i)*(drho(i)*T(i,:) + Tp(i,:) ) 
!     L22(i,:) = fvasq*Tp(i,:) + isig*eta(i)*(drho(i)*T(i,:) + Tp(i,:) )

     D0 = T(i,:)
     D1 = drho(i)*T(i,:) + Tp(i,:)
     Dhat = eta(i)*D1 + deta(i)*D0 
 
     L11(i,:) = fvasq*Tp(i,:) +isig*Dhat
     L12(i,:) = -2d0*Dhat
     L13(i,:) = ifkx*Dhat
     L14(i,:) = L13(i,:) 
    
     L21(i,:) = 0.5*kappa2*Dhat
     L22(i,:) = fvasq*Tp(i,:) + isig*Dhat 

     
     L33(i,:) = Tp(i,:) !these two imply dw/dz + dphi/dz = 0 => vz=0 
     L34(i,:) = Tp(i,:) 
     
     L44(i,:) = Tp(i,:) + kx*T(i,:) !goldreich and lynden-bell potential bc 
  endif
  
  if(eos.eq.'uns') then !then no gradient at top
     L11(i,:) = Tp(i,:)
     L22(i,:) = Tp(i,:)
     L33(i,:) = Tp(i,:)
     L44(i,:) = Tp(i,:)
  endif
  
  
!put small matrices into big matrix
   big_matrix = 0d0 

   big_matrix(1:nz, 1:nz)         = L11
   big_matrix(1:nz, nz+1 : 2*nz)  = L12
!   big_matrix(1:nz, 2*nz+1:3*nz)  = L13
!   big_matrix(1:nz, 3*nz+1:4*nz)  = L14 

   big_matrix(nz+1:2*nz, 1:nz)         = L21
   big_matrix(nz+1:2*nz, nz+1 : 2*nz)  = L22
!   big_matrix(nz+1:2*nz, 2*nz+1:3*nz)  = L23
!   big_matrix(nz+1:2*nz, 3*nz+1:4*nz)  = L24

!   big_matrix(2*nz+1:3*nz, 1:nz)         = L31
!   big_matrix(2*nz+1:3*nz, 2*nz+1:3*nz)  = L33
!   big_matrix(2*nz+1:3*nz, 3*nz+1:4*nz)  = L34

!   big_matrix(3*nz+1:4*nz, 2*nz+1:3*nz)  = L43
!   big_matrix(3*nz+1:4*nz, 3*nz+1:4*nz)  = L44
  return
end subroutine construct_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








subroutine eigenvalue_problem(output)
  use global
  implicit none 
  real*8, intent(out) :: output(2)
  integer :: info, loc(1) 
  complex*16 :: input(big_nz, big_nz), big_matrix(big_nz, big_nz)
  complex*16 :: fsigma, w(big_nz), vl(big_nz, big_nz), & 
       vr(big_nz,big_nz), work2(4*big_nz)  
  real*8 :: rwork2(2*big_nz), absw(big_nz) 
  complex*16 :: work(3*big_nz)
  complex*16 :: U(1,big_nz), VT(big_nz, big_nz),Mwvec(big_nz)
  real*8 :: singular(big_nz), rwork(5*big_nz), w2
  complex*16, external :: ZDOTC
  
  call construct_matrix(big_matrix)
  input = big_matrix  
  
  if(method .eq. 'eig') then 
     call zgeev('N', 'V', big_nz, input, big_nz, w, vl, big_nz, vr, big_nz, work2, & 
          4*big_nz, rwork2 , info) 
     absw= abs(w)
     loc = minloc(absw)
     fsigma = w(loc(1))/maxval(absw) 
     wvector = vr(:,loc(1)) 
  else if(method .eq. 'svd') then 
     call ZGESVD('N', 'A', big_nz, big_nz, input, big_nz, singular, U, 1, VT, big_nz, &
          WORK, 3*big_nz, RWORK, INFO)
     wvector(:)     = conjg(VT(big_nz, :))
     w2   = ZDOTC(big_nz,wvector,1,wvector,1)
     call  ZGEMV('N',big_nz,big_nz,(1d0,0d0),big_matrix,big_nz,wvector,1,(0d0,0d0),Mwvec,1)
     fsigma = ZDOTC(big_nz,wvector,1,Mwvec,1)/w2 
  endif
  
  output(1) = dble(fsigma)
  output(2) = dimag(fsigma)
  
  return
end subroutine eigenvalue_problem


subroutine svd_final_matrix 
  use global
  implicit none
  integer :: info 
  complex*16 :: input(big_nz, big_nz), big_matrix(big_nz, big_nz)
  complex*16 :: work(3*big_nz) 
  complex*16 :: U(1,big_nz), VT(big_nz, big_nz) 
  real*8 :: singular(big_nz), rwork(5*big_nz)  
  complex*16, external :: ZDOTC

  call construct_matrix(big_matrix)
  input = big_matrix
  call ZGESVD('N', 'N', big_nz, big_nz, input, big_nz, singular, U, 1, VT, big_nz, &
       WORK, 3*big_nz, RWORK, INFO)
  condition = minval(singular)/maxval(singular)
  return
end subroutine svd_final_matrix


SUBROUTINE FCN(N,X,FVEC,IFLAG)
  use global
  integer :: N, IFLAG
  real*8 :: X(N),FVEC(N)
  real*8 :: output(2), err_init(2)
  common /share/ err_init  
  
  output = 0d0 
  sig   = dcmplx(x(1), x(2))
  call eigenvalue_problem(output)
  output  = output/err_init 
  fvec(1) = output(1)
  fvec(2) = output(2) 
  
  RETURN
END SUBROUTINE FCN

subroutine get_eigenvalue
  use global
  implicit none 
  integer, parameter :: N = 2, LWA = N*(3*N+13)
  integer :: INFO 
  real*8, parameter :: TOL = 1d-15
  real*8 :: x(N) , fvec(N), wa(LWA) 
  real*8 :: err_init(2), output(2)
  external :: FCN 
  common /share/ err_init 
 
  output = 0d0 
  !call to get initial error 
  call eigenvalue_problem(output)
  err_init = dabs(output) 
  
  !iterate to get eigenvalue 
  x(1) = dble(sig)
  x(2) = dimag(sig)
  call HYBRD1(FCN,N,X,FVEC,TOL,INFO,WA,LWA)
  sig   = dcmplx(x(1), x(2))
  
  return  
end subroutine get_eigenvalue

subroutine maximize_growth_rate
  use global
  implicit none
  integer, parameter :: N=2 
  integer :: i, loc(1)
  real*8, parameter :: tol=1d-15, rate_lim=1d1
  real*8 :: output(N), a_x, b_x, rates(ntrials), dsig, &
       guess_growth
  complex*16 :: eigen(ntrials)
 
  !trial growth rate read in as table 
  if(refine .eqv. .true.) then
     
     sig = sig_trials(kcount) !dcmplx(0d0, dimag(sig_trials(k)))
     write(6,fmt='(A)') '---------- trial eigen ----------'
     write(6,fmt='(3(e10.3,x))') dble(sig), dimag(sig)
     call get_eigenvalue
     
  else  !trial growth rate divided uniformly between min and max rates
     eigen   = dcmplx(0d0, 0d0) 
     rates   = 0d0 
     
     a_x = -max_growth
     b_x = -min_growth
     dsig = (b_x - a_x)/(ntrials - 1d0) 
     
     do i=1, ntrials 
        guess_growth = a_x + dsig*dble(i) 
        sig = dcmplx(0d0, guess_growth)     

        call get_eigenvalue
        eigen(i) =  sig 
        rates(i) = -dimag(sig)
        
        call svd_final_matrix

        if( (condition .ge. tol) .or. (rates(i) .gt. rate_lim )) then
           rates(i) = 0d0 
           if(prin .eqv. .true.) write(6,fmt='(A,I3,2(e22.15))') 'REJECT  ',i,condition,rates(i)
        else
           if(prin .eqv. .true.) write(6,fmt='(A,I3,2(e22.15))') 'accept  ',i,condition, rates(i)
        endif      

     enddo
      
     loc = maxloc(rates)
     sig = eigen(loc(1))
  endif
  
  !get eigenfunctions
  call eigenvalue_problem(output)
  vx(:)  = wvector(1:nz)
  vy(:)  = wvector(nz+1:2*nz)
!  den(:) = wvector(2*nz+1:3*nz)
!  pot(:) = wvector(3*nz+1:4*nz)
  
  call svd_final_matrix 
  
  write(6,fmt='(A)') '---------- sig_re, sig_im, inv. cond. no. ----------'
  write(6,fmt='(3(e10.3,x))') dble(sig), dimag(sig), condition 
  write(6,fmt='(A)') '----------------------------------------------------' 
end subroutine maximize_growth_rate
