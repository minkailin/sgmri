module global
  logical :: refine, prin, use_old
  character*3 :: eos, method
  character*4 :: var, vbc 
  integer, parameter :: nvar=3
  integer :: nz, nk, big_nz, ntrials , kcount, nuse_old 
  real*8, parameter :: omega=1.0, bigG = 1.0, mstar = 1.0, r0=1.0
  real*8, parameter :: pi = 2d0*acos(0d0) 
  real*8 :: bigQ, beta, shear, kmin, kmax, kx, dk, qmin, qmax, bmin, &
       bmax, var_min, var_max 
  real*8 :: zmax, dz    
  real*8 :: kappa, omegak, omegaz2, fQ, Rm, amp, zstar, Rm_max, Rm_min
  real*8 :: condition, min_growth, max_growth, sig_re, sig_im    
  real*8, allocatable :: dnorm(:), drho(:), d2rho(:), valf(:), zaxis(:),&
       kaxis(:), csq(:), eta(:), deta(:), d2eta(:) 
  real*8, allocatable :: T(:,:), Tp(:,:), Tpp(:,:)  
  complex*16, parameter :: ii = (0d0, 1d0)      
  complex*16 :: sig, sig_old
  complex*16, allocatable :: vx(:), vy(:), bz(:)
  complex*16, allocatable :: den(:), pot(:), wvector(:), sig_trials(:)
end module global

program sgmri
  use global 
  implicit none 
  integer :: i, j, k  
  integer :: haf_lmax, lmax  
  character*3 :: str_mode 
  real*8  :: zbar, T_l, dT_l, d2T_l, lmode, dummy
  real*8, external :: density, dlogrho, alfven
  namelist /params/ bigQ, beta, kx, var, shear, kmin, kmax, qmin,&
       & qmax, bmin, bmax, eos 
  namelist /resis/ Rm, amp, Rm_min, Rm_max
  namelist /grid/ zmax, nz, nk, vbc, nuse_old   
  namelist /eigen/ method, min_growth, max_growth, ntrials , refine, prin
  
  write(6,fmt='(A)') '--------------- welcome ---------------'
  
  !read input parameters.
  open(7, file="input")
  read(7, nml=params)
  read(7, nml=resis)
  read(7, nml=grid)
  read(7, nml=eigen) 
  close(7)

  if( (eos.eq.'pol').and.(zmax .ge. 1d0)) then
     print*, 'polytropic disk does not allow zmax>1'
     stop
  endif

  if(nk .eq. 1) then 
     write(6,fmt='(A)') ' fixed parameters'
     write(6,fmt='(A,f5.2)') 'self-gravity Q =', bigQ
     write(6,fmt='(A,f6.0)') 'plasma beta    =', beta
     write(6,fmt='(A,f4.1)') 'horizontal k   =', kx
     write(6,fmt='(A,f6.2)') 'elsasser num   =', Rm
  else
     if(var.eq.'grav') then
        write(6,fmt='(A)') ' varying gravity'
        write(6,fmt='(A,f5.2,A,f5.2)') 'self-gravity Q =', qmin, ' to ' ,qmax 
        write(6,fmt='(A,f6.0)')    'plasma beta    =', beta
        write(6,fmt='(A,f4.1)')    'horizontal k   =', kx
        write(6,fmt='(A,f6.2)') 'elsasser num   =', Rm
        var_min = qmin
        var_max = qmax 
     endif
     if(var.eq.'beta') then
        write(6,fmt='(A)') ' varying B field'
        write(6,fmt='(A,f5.2)')   'self-gravity Q =', bigQ
        write(6,fmt='(A,f6.0,A,f6.0)') 'plasma beta    =',bmin, ' to ' ,bmax 
        write(6,fmt='(A,f4.1)')    'horizontal k   =', kx
        write(6,fmt='(A,f6.2)') 'elsasser num   =', Rm
        var_min = log10(bmin)
        var_max = log10(bmax)
     endif
     if(var.eq.'wave') then
        write(6,fmt='(A)') ' varying horizontal wave number'
        write(6,fmt='(A,f5.2)')    'self-gravity Q =', bigQ
        write(6,fmt='(A,f6.0)')    'plasma beta    =', beta
        write(6,fmt='(A,f4.1,A,f4.1)') 'horizontal k   =', kmin, ' to ' ,kmax 
        write(6,fmt='(A,f6.2)') 'elsasser num   =', Rm
        var_min = kmin
        var_max = kmax
     endif
     if(var.eq.'resi') then
        write(6,fmt='(A)') ' varying resistivity wave number'
        write(6,fmt='(A,f5.2)')    'self-gravity Q =', bigQ
        write(6,fmt='(A,f6.0)')    'plasma beta    =', beta
        write(6,fmt='(A,f4.1)') 'horizontal k   =', kx
        write(6,fmt='(A,f6.2,A,f6.2)') 'elsasser num   =', Rm_min, ' to ' , Rm_max
        var_min = log10(Rm_min)
        var_max = log10(Rm_max) 
     endif
  endif
  write(6,fmt='(A,A)')    'e.o.s.         = ', eos
  write(6,fmt='(A,f4.0)')    'conduct. boost = ', amp
  write(6,fmt='(A,A)')    'method         = ', method
  write(6,fmt='(A,A)')    'vert. b.c.     = ', vbc
  !setup parameters (dimensionless units)
  omegak  = sqrt(bigG*mstar/r0**3d0)/omega 
  kappa   = sqrt(2.0*(2.0 - shear))
  omegaz2 = omegak**2d0 
 
  big_nz = nvar*nz 

  !allocate grids
  allocate(zaxis(nz))
  allocate(kaxis(nk))
  allocate(dnorm(nz))!rho/rho0, rho0=rho(z=0)
  allocate(drho(nz))!dlogrho/dz
  allocate(d2rho(nz))!d^2 (rho/rho0)/dz^2 / rho0
  allocate(valf(nz))!1/beta*rho
  allocate(csq(nz))
  allocate(eta(nz))
  allocate(deta(nz))
  allocate(d2eta(nz))
  allocate(T(nz, nz))
  allocate(Tp(nz,nz))  
  allocate(Tpp(nz,nz))

  allocate(vx(nz)) !solve for
  allocate(vy(nz)) !solve for
  allocate(bz(nz)) !solve for 
  allocate(den(nz))!solve for
  allocate(pot(nz))!solve for
  allocate(wvector(big_nz))   

  if(refine .eqv. .true.) allocate(sig_trials(nk))

  !setup Z axis. these are extrema of T_lmax(Z/zmax) plus end point
  haf_lmax = nz-1
  lmax = 2*haf_lmax
  zaxis(1)  = 0d0
  zaxis(nz) = zmax
  do j = haf_lmax+1, lmax-1
     zaxis(j-haf_lmax+1) = -zmax*cos(dble(j)*pi/lmax)
  enddo

  !setup the chebyshev matrices
  do j=1, nz !jth physical grid
     zbar = zaxis(j)/zmax
     do k=1, nz !kth basis
        lmode = 2d0*(k-1d0)!lth chebyshev mode
        call chebyshev_poly(lmode, zbar, T_l, dT_l, d2T_l)
        T(j,k)  = T_l
        Tp(j,k) = dT_l/zmax     !division needed to get d/dz
        Tpp(j,k)= d2T_l/zmax**2
     enddo
  enddo

  !setup the variable axis if we are varying. could be Q, beta, or kx, or Rm
  if(nk .gt. 1) then 
     dk = (var_max - var_min)/(nk - 1d0) 
     open(10, file='var_axis.dat')
     do i=1, nk
        kaxis(i)   = var_min + dk*(i - 1d0)
        write(10, fmt = '(2(e22.15,x))') dble(i), kaxis(i)
     enddo
     close(10)
  endif
  
  !if sig_table = true then read in table of trial values 
  if(refine .eqv. .true.) then
     open(20, file='eigenvalues.dat')
     do i=1, nk 
        read(20, fmt='(6(e22.15,x))') dummy, dummy, dummy, sig_re, sig_im, dummy
        sig_trials(i) = dcmplx(sig_re, sig_im)
     enddo
     close(20) 
  endif
  
  !solve the eigenvalue problem for each k (= bigQ, beta or kx)
  open(99, file='basic.dat')
  open(100, file='params.dat')
  
  if(nk .eq. 1) then

     call get_basic 
     
     write(100, fmt='(5(e22.15,x))') zmax, dble(nz), fQ, shear, Rm
     
     do j=1, nz 
        write(99, fmt='(9(e22.15,x))') zaxis(j), dnorm(j), valf(j),&
             & csq(j), eta(j), drho(j), deta(j), d2rho(j), d2eta(j)
     enddo
     
     
  endif
  
  open(20, file='eigenvalues.dat')  
  do i=1, nk 
     kcount = i      
     if(nk .gt. 1) then
        
        if(var.eq.'grav') then
           bigQ = kaxis(i)
           write(6,fmt='(A,A,A,f5.2,A)') '---------------', var, '=',bigQ,'---------------'
           call get_basic
        endif
        if(var.eq.'beta') then
           beta = 10d0**kaxis(i)
           write(6,fmt='(A,A,A,f7.1,A)') '---------------', var, '= ',beta,'---------------'
           call get_basic 
        endif
        if(var.eq.'wave') then
           kx   = kaxis(i)
           write(6,fmt='(A,A,a,f5.2,A)') '---------------', var, '= ',kx,'----------------'
           call get_basic 
        endif
        if(var.eq.'resi') then
           Rm   = 10d0**kaxis(i)
           write(6,fmt='(A,A,a,f5.2,A)') '---------------', var, '=', Rm,'----------------'
           call get_basic
        endif 

        write(100, fmt='(5(e22.15,x))') zmax, dble(nz), fQ, shear, Rm
        do j=1, nz 
           write(99, fmt='(9(e22.15,x))') zaxis(j), dnorm(j), valf(j)&
                &, csq(j), eta(j), drho(j), deta(j), d2rho(j), d2eta(j)
        enddo

        
     endif
    
     write(str_mode, fmt='(I3)') i 

     if(i.ge.nuse_old) then
     use_old = .true.
     else
     use_old = .false.
     endif

     !get largest growth mode 

     call maximize_growth_rate
     sig_old = sig
 
     !output
     write(20, fmt='(6(e22.15,x))') bigQ, beta, kx, dble(sig), dimag(sig), condition 
     
     open(10, file='eigenfunctions_'//trim(adjustl(str_mode))//'.dat')      
     do j=1, big_nz 
        write(10, fmt='(2(e22.15,x))') dble(wvector(j)), dimag(wvector(j))
     enddo
     close(10)
     
  enddo
  close(20)!eigenvalue.dat
  close(99)!basic.dat 
  close(100)!params.dat
end program sgmri
