SUBROUTINE F (NEQ, x, Y, YDOT, RPAR, IPAR)
  use global
  implicit none
  integer :: NEQ, IPAR
  real*8 :: x, Y(NEQ), YDOT(NEQ), RPAR
  real*8 :: logrho, phidot
  
  logrho = Y(1)
  phidot = Y(2)
  
  YDOT(1) = -omegaz2*x/fQ**2d0 - phidot
  YDOT(2) = exp(logrho)/(bigQ*fQ**2d0)
  
  return
end subroutine F

subroutine jac
  use global
  implicit none
  return
end subroutine jac

SUBROUTINE FCN_FQ_ISO(N,X,FVEC,IFLAG)
  use global
  implicit none
  integer :: N, IFLAG
  real*8 :: X(N),FVEC(N)
  integer, parameter :: NEQ = 2, ITOL = 1 , MF = 10, LIW = 60, IOPT = 0
  integer, parameter :: LRW = 2*(20 + 16*NEQ)
  integer :: i, ITASK, ISTATE, IWORK(LIW), IPAR
  real*8, parameter :: RTOL = 1d-15, ATOL=1d-15, eps=1d-2
  real*8 :: Y(NEQ), xSTART, xEND, RWORK(LRW), RPAR
  external :: F, JAC

  Y(1) = 0d0 !log rho_hat = 0 at z=0
  Y(2) = 0d0 !dphi/dz =0 at z=0

  ITASK = 1
  ISTATE= 1

  xSTART  = 0d0 
  xEND    = x(1)

  call DVODE (F, NEQ, Y, xSTART, xEND, ITOL, RTOL, ATOL, ITASK, &
          ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF, &
          RPAR, IPAR)

  FVEC(1) = exp(Y(1)) - eps 

  return 
END SUBROUTINE FCN_FQ_ISO


subroutine get_fQ_iso
   use global
   implicit none
   integer, parameter :: N = 1, LWA = N*(3*N+13)
   integer :: INFO
   real*8, parameter :: TOL = 1d-15
   real*8 :: x(N) , fvec(N), wa(LWA)
   external :: FCN_FQ_ISO 

   fQ   = 1d0 

   x(1) = 2d0 !trial value of zmax in units of non-sg scale heights 
100   call HYBRD1(FCN_FQ_ISO, N,X,FVEC,TOL,INFO,WA,LWA) 
    if(x(1).le.0d0) then
       x(1) = 1d0
       goto 100
    endif
  
   fQ = 1d0/x(1)

   return 
end subroutine get_fQ_iso 




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
  
  dnorm(Nzmid)   = 1d0
  drho(Nzmid)    = 0d0
  
  ITASK = 1
  ISTATE= 1
  
  do i = Nzmid, nz-1
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
     drho(i+1)  = -omegaz2*zaxis(i+1)/fQ**2d0 - Y(2)

     dnorm(2*Nzmid - i - 1) = dnorm(i+1)
     drho(2*Nzmid -i - 1)   =-drho(i+1)
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
  real*8, external :: density, dlogrho

  if(eos .eq. 'pol') then
     fQ      = sqrt(bigQ)*acos(omegaz2*bigQ/(1d0 + omegaz2*bigQ)) !H*Omega/cs0 
     fQ      = 1d0/fQ 
  else if((eos.eq.'iso').or.(eos.eq.'uns')) then
     call get_fQ_iso
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
     d2rho = 0d0
     csq   = 1d0
     valf  = 1d0/beta 
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

subroutine fill_resistivity
  use global
  implicit none
  integer :: i 
  real*8 :: surf, surf0, gplus, gminus, gsum, gdiff
  real*8 :: eta0 , surf_star 
  real*8, external :: density, dlogrho
  
  eta0 = fQ**2d0/(Rm*beta)

  if(amp.eq.1d0) then
     eta  = eta0
     deta = 0d0 
     d2eta= 0d0
  else
     print*, 'code set up for uniform resis only'
     stop
  endif

  return
end subroutine fill_resistivity

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
     dT_l =lsq
     d2t_l =lsq*(lsq-1d0)/3d0

     if(zbar.eq.-1d0) then
        dT_l = (-1d0)**(l+1d0)*dT_l
        d2T_l= (-1d0)**(l+2d0)*d2T_l
     endif

  endif

  return
end subroutine chebyshev_poly

subroutine construct_matrix(big_matrix) 
  use global
  implicit none 
  integer :: i 
  complex*16 :: L11(nz, nz), L12(nz, nz), L13(nz, nz), L14(nz, nz)
  complex*16 :: L21(nz, nz), L22(nz, nz), L23(nz, nz), L24(nz, nz)
  complex*16 :: L31(nz, nz), L32(nz, nz), L33(nz, nz), L34(nz, nz)
  complex*16 :: L41(nz, nz), L42(nz, nz), L43(nz, nz), L44(nz, nz)
  complex*16, intent(out) :: big_matrix(big_nz, big_nz) 
  complex*16 :: sbar , dsbar, isig, ifkx, c2, c1, c0, Dhat(nz), Dhattil(nz)
  complex*16 :: c2til, c1til, c0til, U(3), V(2)
  real*8 :: kxsq, fkx, rho, dden, fvasq, vasq, kappa2, d2logrho,&
       & dlogcsq, sgn
  real*8 :: Diff(nz)
  real*8 :: D2(nz), D1(nz), D0(nz), D2til(nz), D1til(nz), D0til(nz) 
  
  L11 = 0d0 ; L12 = 0d0 ; L13 = 0d0 ; L14 = 0d0 

  L21 = 0d0 ; L22 = 0d0 ; L23 = 0d0 ; L24 = 0d0
  
  L31 = 0d0 ; L32 = 0d0 ; L33 = 0d0 ; L34 = 0d0 

  L41 = 0d0 ; L42 = 0d0 ; L43 = 0d0 ; L44 = 0d0 
 
  isig  = ii*sig 
  fkx   = fQ*kx
  ifkx  = ii*fkx 
  kxsq  = kx*kx 
  kappa2= kappa*kappa

  U(1) = isig
  U(2) = -2d0
  U(3) = ifkx

  V(1) = 0.5*kappa2
  V(2) = isig

  !interior
  do i=2, nz-1
     rho  = dnorm(i)  
     dden = drho(i) !this is dlogrho/dz   
     vasq = valf(i)
     fvasq = fQ**2d0*vasq
     d2logrho = d2rho(i) - drho(i)**2d0
     if((eos.eq.'iso').or.(eos.eq.'uns')) then 
        dlogcsq = 0d0
     endif
     if(eos.eq.'pol') then
        dlogcsq = drho(i)
     endif

     sbar  = sig - ii*eta(i)*kxsq
   
     D0 = T(i,:)
     D1 = drho(i)*T(i,:) + Tp(i,:)
     D2 = d2rho(i)*T(i,:)+2d0*drho(i)*Tp(i,:)+Tpp(i,:)
     
     Dhat = eta(i)*(sig*D2 + tilt*kx*shear*D1) - ii*sig*sbar*D0 

     L11(i,:) = fvasq*((1d0+tilt*tilt)*kxsq*T(i,:)*sig - tilt*kx*shear*Tp(i,:) - Tpp(i,:)*sig) &
          -Dhat*U(1) + tilt*tilt*kx*fvasq*tilt*Tp(i,:)*V(1) - ii*tilt**2d0*kxsq*shear*eta(i)*D0*V(1)
     L12(i,:) = ii*tilt*kx*fvasq*sig*Tp(i,:) -Dhat*U(2) + tilt*tilt*kx*fvasq*tilt*Tp(i,:)*V(2) - ii*tilt**2d0*kxsq*shear*eta(i)*D0*V(2)
     L13(i,:) = -Dhat*U(3) + tilt*tilt*kx*fvasq*fQ*Tpp(i,:)
     L14(i,:) = L13(i,:)
 


     Dhat = eta(i)*(sig*D2/shear - tilt*kx*D1) + (ii*tilt**2d0*fvasq&
          &/shear)*(d2logrho*D0  + dden*Tp(i,:)) - ii*sig&
          &*sbar*D0/shear

     L21(i,:) = ii*fvasq*Tpp(i,:) + Dhat*V(1) + ii*eta(i)*D2*U(1)
     L22(i,:) = fvasq*Tpp(i,:)*sig/shear + Dhat*V(2) + ii*eta(i)*D2*U(2)
     L23(i,:) = ii*tilt*fQ*vasq*( (sig**2d0/csq(i))*(Tp(i,:) - &
          & dlogcsq*T(i,:)) + d2logrho*fQ*fQ*Tp(i,:) &
          & + dden*fQ*fQ*Tpp(i,:))/shear + ii*eta(i)*D2*U(3)
     L24(i,:) = ii*tilt*fQ*vasq*( d2logrho*fQ*fQ*Tp(i,:) &
          & + dden*fQ*fQ*Tpp(i,:))/shear + ii*eta(i)*D2*U(3)


 
     L31(i,:) = tilt*fQ*Tp(i,:)*V(1) + dden*tilt*fQ*T(i,:)*V(1) + sig*fQ*kx*T(i,:) 
     L32(i,:) = tilt*fQ*Tp(i,:)*V(2) + dden*tilt*fQ*T(i,:)*V(2)
     L33(i,:) = sig*sig*T(i,:)/csq(i) + fQ*fQ*Tpp(i,:) + dden*fQ*fQ&
          &*Tp(i,:)
     L34(i,:) = fQ*fQ*Tpp(i,:) + dden*fQ*fQ&
          &*Tp(i,:)


     L43(i,:) = -T(i,:)*rho/(csq(i)*fQ**2d0*bigQ)
     L44(i,:) = Tpp(i,:) - kxsq*T(i,:)     
  enddo
  
 
! upper and lower boundary conditions 

  if(eos.eq.'uns') then 
     
     do i = 1, nz, nz-1
        L11(i,:) = Tp(i,:)
        L22(i,:) = Tp(i,:)
        L33(i,:) = Tp(i,:)
        L44(i,:) = Tp(i,:)
     enddo

  else
     do i=1, nz, nz-1

        rho  = dnorm(i)
        dden = drho(i)
        vasq = valf(i)
        fvasq = fQ**2d0*vasq
        
        
        D0 = T(i,:)
        D1 = drho(i)*T(i,:) + Tp(i,:)
        D2 = d2rho(i)*T(i,:)+2d0*drho(i)*Tp(i,:)+Tpp(i,:)
        
        if(vbc.eq.'wall') then
           Dhat = D1til
           !bx = 0
           L11(i,:) = fvasq*Tp(i,:) + eta(i)*(D1*U(1) + ii*tilt*kx*D0*V(1))
           L12(i,:) = eta(i)*(D1*U(2) + ii*tilt*kx*D0*V(2))
           L13(i,:) = eta(i)*D1*U(3)
           L14(i,:) = L13(i,:) 

           !by=0, assuming bx=0
           L21(i,:) = (ii*tilt*tilt*fvasq*dden*D0 + sig*eta(i)*D1)*V(1)
           L22(i,:) = (ii*tilt*tilt*fvasq*dden*D0 + sig*eta(i)*D1)&
                &*V(2) + sig*fvasq*Tp(i,:)
           L23(i,:) = ii*tilt*fQ*vasq*(sig**2d0*T(i,:)/csq(i) + fQ*fQ&
                &*dden*Tp(i,:))
           L24(i,:) = ii*tilt*fQ*fvasq*dden*Tp(i,:) 
           
           !vz = 0
           L31(i,:) = tilt*D0*V(1)
           L32(i,:) = tilt*D0*V(2)
           L33(i,:) = fQ*Tp(i,:)
           L34(i,:) = L33(i,:)
        endif
        
        if(vbc.eq.'nodv') then!      !dv=0 (all components)
           
           L11(i,:) = T(i,:) !vx=0
           
           L22(i,:) = T(i,:) !vy=0
           
           L31(i,:) = tilt*D0*V(1)
           L32(i,:) = tilt*D0*V(2)
           L33(i,:) = fQ*Tp(i,:)
           L34(i,:) = L33(i,:)
           
        else if(vbc.eq.'gb94') then !gammie/balbus hot halo model
           
           !zero divergence condition(identical to zero gas lagrangian pressure pert)
           L11(i,:) = sig*kx*T(i,:) + tilt*Tp(i,:)*V(1)
           L12(i,:) = tilt*Tp(i,:)*V(2)
           L13(i,:) = fQ*Tpp(i,:)
           L14(i,:) = L13(i,:)
           
           !by=0 condition

           L21(i,:) = ii*fvasq*Tp(i,:) + ii*eta(i)*D1*U(1) + eta(i)&
                &*(sig*D1/shear - tilt*kx*D0)*V(1)
           L22(i,:) = (sig/shear)*fvasq*Tp(i,:) + ii*eta(i)*D1*U(2) + eta(i)&
                &*(sig*D1/shear - tilt*kx*D0)*V(2)
           L23(i,:) = ii*eta(i)*D1*U(3)
           L24(i,:) = L23(i,:)
           
           !bz - ii*bx = 0 relation  (upper)
           !bz + ii*bx = 0 relation  (lower)

           if(i.eq.nz) sgn =  1d0
           if(i.eq.1)  sgn = -1d0

              L31(i,:) = sgn*sbar*(fvasq*Tp(i,:) + eta(i)*(D1*U(1) +&
                   & ii*tilt*kx*D0*V(1))) + sig*kx*fvasq*T(i,:) &
                   & - ii*eta(i)*kx*(fvasq*Tpp(i,:) + eta(i)*(D2*U(1)&
                   & + ii*tilt*kx*D1*V(1)))
              L32(i,:) =  sgn*sbar*(eta(i)*(D1*U(2) +&
                   & ii*tilt*kx*D0*V(2))) &
                   & - ii*eta(i)*kx*(eta(i)*(D2*U(2)&
                   & + ii*tilt*kx*D1*V(2)))
              L33(i,:) =  sgn*sbar*eta(i)*D1*U(3) &
                   & - ii*eta(i)*kx*eta(i)*D2*U(3) 
              L34(i,:) = L33(i,:) 
              
           
        endif
     
     !GLB65 potential bc. allowing for vertical motion
     L41(i,:) = rho*tilt*D0*V(1)
     L42(i,:) = rho*tilt*D0*V(2)
     L43(i,:) = rho*fQ*Tp(i,:)

     if(i.eq.1)  L44(i,:) = bigQ*sig**2d0*fQ*(Tp(i,:) - kx*T(i,:)) + rho*fQ*Tp(i,:)
     if(i.eq.nz) L44(i,:) = bigQ*sig**2d0*fQ*(Tp(i,:) + kx*T(i,:)) + rho*fQ*Tp(i,:)
  enddo
endif
  
!put small matrices into big matrix
   big_matrix = 0d0 

   big_matrix(1:nz, 1:nz)         = L11
   big_matrix(1:nz, nz+1 : 2*nz)  = L12
   big_matrix(1:nz, 2*nz+1:3*nz)  = L13
!   big_matrix(1:nz, 3*nz+1:4*nz)  = L14 

   big_matrix(nz+1:2*nz, 1:nz)         = L21
   big_matrix(nz+1:2*nz, nz+1 : 2*nz)  = L22
   big_matrix(nz+1:2*nz, 2*nz+1:3*nz)  = L23
!   big_matrix(nz+1:2*nz, 3*nz+1:4*nz)  = L24

   big_matrix(2*nz+1:3*nz, 1:nz)         = L31
   big_matrix(2*nz+1:3*nz, nz+ 1 :2*nz)  = L32 
   big_matrix(2*nz+1:3*nz, 2*nz+1:3*nz)  = L33
!   big_matrix(2*nz+1:3*nz, 3*nz+1:4*nz)  = L34

!   big_matrix(3*nz+1:4*nz, 1:nz)         = L41
!   big_matrix(3*nz+1:4*nz, nz+ 1 :2*nz)  = L42
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
  integer :: info, loc(1), k
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
  output  = output/err_init!sqrt(abs(err_init(1))**2d0 +  abs(err_init(2))**2d0)
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
     
  
     if (use_old.eq..false.) then
        do i=1, ntrials
           
           guess_growth = a_x + dsig*dble(i)
           sig = dcmplx(0d0, guess_growth)
           
           call get_eigenvalue
           eigen(i) =  sig
           rates(i) = -dimag(sig)
           
           call svd_final_matrix
           
           if( (condition .ge. tol) .or. (rates(i) .gt. rate_lim )) then
              rates(i) = 0d0
!              if(prin .eqv. .true.) write(6,fmt='(A,I3,2(e22.15))') 'REJECT ',i,condition,rates(i)
           else
              if((prin .eqv. .true.).and.(rates(i).ge.1d-2)) write(6,fmt='(A,I3,2(e22.15))') 'accept ',i,condition, rates(i)
           endif
           
        enddo
        loc = maxloc(rates)
        sig = eigen(loc(1))
        
     else
        sig = sig_old
        write(6,fmt='(A)') '---------- trial eigen ----------'
        write(6,fmt='(3(e10.3,x))') dble(sig), dimag(sig)
        call get_eigenvalue

        call svd_final_matrix
        rates(1) = -dimag(sig)    
!        if( (condition .ge. tol) .or. (rates(1) .gt. rate_lim )) then
!        print*, 'bad eigenvlaue, do trials'
!        do i=1, ntrials
!           guess_growth = a_x + dsig*dble(i)
!           sig = dcmplx(0d0, guess_growth)
!           call get_eigenvalue
!           eigen(i) =  sig
!           rates(i) = -dimag(sig)
!           call svd_final_matrix
!           if( (condition .ge. tol) .or. (rates(i) .gt. rate_lim ))  rates(i) = 0d0
!        enddo
!        loc = maxloc(rates)
!        sig = eigen(loc(1))
!        endif        
        endif
     
  endif
  
  !get eigenfunctions
  call eigenvalue_problem(output)
  vx(:)  = wvector(1:nz)
  vy(:)  = wvector(nz+1:2*nz)
  den(:) = wvector(2*nz+1:3*nz)
  pot(:) = 0d0 !wvector(3*nz+1:4*nz)
  
  call svd_final_matrix 
  
  write(6,fmt='(A)') '---------- sig_re, sig_im, inv. cond. no. ----------'
  write(6,fmt='(3(e10.3,x))') dble(sig), dimag(sig), condition 
  write(6,fmt='(A)') '----------------------------------------------------' 
end subroutine maximize_growth_rate
