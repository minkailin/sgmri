function cheby, lmode, zbar
  
  lsq = lmode*lmode
  t  = acos(zbar)
  tl = t*lmode
  
  T_l = cos(tl)
  
  if(abs(zbar) lt 1d0) then begin
     dT_l = lmode*sin(tl)/sin(t)
     d2T_l= -lsq*cos(tl)/sin(t)^2 + lmode*cos(t)*sin(tl)/sin(t)^3
  endif else begin
     dT_l = lsq
     d2t_l =lsq*(lsq-1d0)/3d0
  endelse
 
  if(zbar eq -1d0) then begin
        dT_l = (-1d0)^(lmode+1d0)*dT_l
        d2T_l= (-1d0)^(lmode+2d0)*d2T_l
  endif

 
  return, [t_l, dt_l, d2t_l]
end

pro error, loc=loc, nk=nk, nvar=nvar, nz=nz, chebaxis=chebaxis, nderiv=nderiv
  
  if not keyword_set(nvar) then nvar=4

  location =strcompress(loc,/remove_all)

  ii = dcomplex(0d0, 1d0)
  
  kappa2 = 1.0
  
  for k=1, nk do begin 

  kmode = k

  file = filepath('params.dat',root_dir='.',subdir=[location(0)])
  nlines = file_lines(file)
  param = dblarr(6,nlines)
  openr,1,file
  readf,1, param, format='(6(e22.15,x))'
  close,1
  zmax   = param(0,k-1) 
  lmodes = param(1,k-1)
  fq     = param(2,k-1)
  shear  = param(3,k-1) 
  rm     = param(4,k-1)     ;this is actually midplane elsasser
  tilt   = param(5,k-1)
  
 if not keyword_set(chebaxis) then begin
  zaxis = -zmax + 2d0*zmax*dindgen(nz)/(nz-1d0)
endif else begin
  nz = lmodes 
  zaxis = dblarr(nz)
  lmax = nz-1
  for j=0, nz-1 do begin
     zaxis(j) = -zmax*cos(!dpi*j/lmax)
  endfor
endelse





  
;basic state for this mode 
  file = filepath('basic.dat',root_dir='.',subdir=[location])
  nlines = file_lines(file)
  array = dblarr(7, nlines)
  openr,1,file
  readf,1,array,format='(7(e22.15,x))'
  close,1

  z_small = array(0,lmodes*(k-1): lmodes*kmode(0)- 1)
  rho     = array(1,lmodes*(k-1): lmodes*kmode(0)- 1)
  vasq    = array(2,lmodes*(kmode(0)-1): lmodes*kmode(0)- 1)
  csq     = array(3,lmodes*(kmode(0)-1): lmodes*kmode(0)- 1)
  eta     = array(4,lmodes*(kmode(0)-1): lmodes*kmode(0)- 1)
  drho    = array(5,lmodes*(kmode(0)-1): lmodes*kmode(0)- 1)
  deta    = array(6,lmodes*(kmode(0)-1): lmodes*kmode(0)- 1)

;interpolate to fine grid
  rho  = interpol(rho, z_small, zaxis,  /spline)
  drho = interpol(drho, z_small, zaxis,  /spline)
  vasq = interpol(vasq, z_small, zaxis, /spline)
  csq  = interpol(csq, z_small, zaxis, /spline) 
  eta  = interpol(eta, z_small, zaxis, /spline) 
  deta = interpol(deta, z_small, zaxis, /spline)

;additional basic state variables 
  dlogcsq  = deriv(zaxis, alog(csq))
  d2logrho = deriv(zaxis, drho) ;(d/dz)(dlogrho/dz)
  d2rho    = d2logrho + drho^2  ; (d^2rho/dz^2)/rho 

;eigenvalue information 
  file = filepath('eigenvalues.dat',root_dir='.',subdir=[location])
  nlines = file_lines(file)
  eigen = dblarr(7, nlines)
  openr,1,file
  readf,1, eigen, format='(7(e22.15,x))'
  close,1
  
  bigQ = eigen(0, kmode(0)-1)
  bta  = eigen(1, kmode(0)-1)
  tilt = eigen(2, kmode(0)-1)
  kx   = eigen(3, kmode(0)-1)
  sig  = dcomplex( eigen(4, kmode(0)-1), eigen(5, kmode(0)-1) )
  sbar = sig -ii*eta*kx^2

;eigenfunction information 
  array = dblarr(2,nvar*lmodes)
  filename =  strcompress('eigenfunctions_'+string(kmode(0), format='(I3)')+'.dat', /remove_all)
  file = filepath(filename,root_dir='.',subdir=[location])
  openr, 1, file 
  readf,1,array,format='(2(e22.15,x))'
  close, 1
 
 
  vx = dcomplexarr(nz) 
  vy = dcomplexarr(nz)
  vz = dcomplexarr(nz)
  w  = dcomplexarr(nz);this is actually w=dp/rho
  pot= dcomplexarr(nz)
  
  dvx = vx
  dvy = vy
  dvz = vz
  dw  = w
  dpot= pot 

  d2vx = vx
  d2vy = vy
  d2w  = w
  d2pot= pot

  bx  = dcomplexarr(nz)
  by  = dcomplexarr(nz)
  bz  = dcomplexarr(nz) 

  dbx = bx
  dby = by
  dbz = bz
 
  d2bx = bx
  d2by = by
  d2bz = bz  
  
  for j=0, nz-1 do begin
     zbar = zaxis(j)/zmax 
     vx(j)  = 0d0 
     vy(j)  = 0d0 
     w(j) = 0d0 
     pot(j) = 0d0 
     
     dvx(j) = 0d0 
     dvy(j) = 0d0 
     dw(j)  = 0d0 
     dpot(j)= 0d0

     d2vx(j) = 0d0
     d2vy(j) = 0d0
     d2w(j)  = 0d0
     d2pot(j)= 0d0

     for i=1, lmodes do begin
        arr = cheby(i-1d0, zbar)
        
        vx(j) +=dcomplex( array(0, i-1),  $
                          array(1, i-1)) $
                *arr(0)         ;this is T_l
        dvx(j) +=dcomplex( array(0, i-1),  $
                           array(1, i-1)) $
                 *arr(1)
        d2vx(j) +=dcomplex( array(0, i-1),  $
                           array(1, i-1)) $
                 *arr(2)       

        vy(j) +=dcomplex( array(0, lmodes + i-1),  $
                          array(1, lmodes + i-1 )) $
                *arr(0)         
        
        dvy(j) +=dcomplex( array(0, lmodes + i-1),  $
                           array(1, lmodes + i-1 )) $
                 *arr(1)
        d2vy(j) +=dcomplex( array(0, lmodes + i-1),  $
                           array(1, lmodes + i-1 )) $
                 *arr(2)

        if(nvar ge 3) then begin
        w(j) +=dcomplex( array(0, 2*lmodes + i-1  ),$
                           array(1, 2*lmodes + i-1) ) $
                 *arr(0)
        dw(j) +=dcomplex( array(0, 2*lmodes + i-1  ),$
                            array(1, 2*lmodes + i-1) ) $
                  *arr(1)
        d2w(j) +=dcomplex( array(0, 2*lmodes + i-1  ),$
                            array(1, 2*lmodes + i-1) ) $
                  *arr(2)        
        endif
        if(nvar eq 4 ) then begin
        pot(j) +=dcomplex( array(0, 3*lmodes + i-1),  $
                           array(1, 3*lmodes + i-1)) $
                 *arr(0)
        
        dpot(j) +=dcomplex( array(0, 3*lmodes + i-1),  $
                            array(1, 3*lmodes + i-1)) $
                  *arr(1)
        d2pot(j) +=dcomplex( array(0, 3*lmodes + i-1),  $
                            array(1, 3*lmodes + i-1)) $
                  *arr(2)
        endif

     endfor
  endfor 
  
  dvx  /= zmax
  dvy  /= zmax
  dw   /= zmax
  dpot /= zmax
  
  d2vx  /= zmax^2
  d2vy  /= zmax^2
  d2w /= zmax^2
  d2pot /= zmax^2

  wtilde  = pot + w
  dwtilde = dpot + dw
  d2wtilde = d2pot + d2w

  V  = 0.5*kappa2*vx + ii*sig*vy
  dV = 0.5*kappa2*dvx + ii*sig*dvy
  ddV= 0.5*kappa2*d2vx+ ii*sig*d2vy

  U  = ii*sig*vx - 2d0*vy + ii*kx*fQ*wtilde
  dU = ii*sig*dvx - 2d0*dvy + ii*kx*fQ*dwtilde
  ddU= ii*sig*d2vx -2d0*d2vy+ ii*kx*fQ*d2wtilde

  D1U = drho*U + dU
  D2U = d2rho*U + 2d0*drho*dU + ddU 

  D1V = drho*V + dV
  D2V = d2rho*V + 2d0*drho*dV + ddV

  vz = (ii/sig)*(fQ*dwtilde + tilt*V)
  dvz= (ii/sig)*(fQ*d2wtilde + tilt*dV)

;calculate  mag field  
  bx = fQ*dvx + eta*(D1U + ii*tilt*kx*V)/(fQ*vasq)
  bx/= ii*sig
 
  by = tilt*(ii*sig*w/csq + fQ*drho*vz) + fQ*dvy - shear*bx + eta*D1V/(fQ*vasq)
  by /= ii*sbar 

  bz = fQ*kx*vx + eta*kx*( U/(fQ*vasq) + ii*tilt*kx*by )
  bz/= -sig 

if keyword_set(nderiv) then begin
  dbx = dcomplex(dydx(zaxis, real_part(bx)), dydx(zaxis, imaginary(bx)))
  dby = dcomplex(dydx(zaxis, real_part(by)), dydx(zaxis, imaginary(by)))
  dbz = dcomplex(dydx(zaxis, real_part(bz)), dydx(zaxis, imaginary(bz)))

  d2bx = dcomplex(dydx(zaxis, real_part(dbx)), dydx(zaxis, imaginary(dbx)))
  d2by = dcomplex(dydx(zaxis, real_part(dby)), dydx(zaxis, imaginary(dby)))
  d2bz = dcomplex(dydx(zaxis, real_part(dbz)), dydx(zaxis, imaginary(dbz))) 
endif else begin
 
  dbx = fQ*d2vx + (eta/(fQ*vasq))*(D2U + ii*tilt*kx*D1V)
  dbx/= ii*sig 

  dby = tilt*( (ii*sig/csq)*(dw - w*dlogcsq) + fq*(d2logrho*vz + drho*dvz)) $
        + fq*d2vy - shear*dbx + eta*D2V/(fQ*vasq)
  dby/= ii*sbar

  dbz = fQ*kx*dvx + eta*kx*( D1U/(fQ*vasq) + ii*tilt*kx*dby )
  dbz/= -sig

  
  d2bx = D1U/(fQ*vasq) + ii*kx*(dbz + tilt*dby)
  d2by = D1V/(fQ*vasq)
  d2bz = -ii*kx*dbx
 
endelse


norm = max(abs(vx))

test_pot = abs(d2pot - kx*kx*pot - rho*w/csq/fq^2/bigQ)
test_den = abs(ii*sig*w/fq/csq + ii*kx*vx + drho*vz + dvz)

test_vx = abs(ii*sig*vx - 2*vy + ii*fQ*kx*wtilde - fQ*vasq*(dbx - ii*kx*(bz + tilt*by)))
test_vy = abs(ii*sig*vy + 0.5*kappa2*vx - fq*vasq*dby)
test_vz = abs(ii*sig*vz + fQ*dwtilde + tilt*fQ*vasq*dby)

test_bx = abs(ii*sbar*bx - fQ*dvx - eta*d2bx)
div = ii*kx*vx + dvz
test_by = abs(ii*sbar*by - fq*dvy + tilt*fq*div + shear*bx - eta*d2by)
test_bz = abs(ii*sbar*bz + ii*fq*kx*vx - eta*d2bz)

test_divb = abs(ii*kx*bx + dbz)

test_pot /= norm
test_den /= norm
test_vx  /= norm
test_vy  /= norm
test_vz  /= norm
test_bx  /= norm
test_by  /= norm
test_bz  /= norm
test_divb/= norm

m1 = max(test_pot)
m2 = max(test_den)
m3 = max(test_vx)
m4 = max(test_vy)
m5 = max(test_vz)
m6 = max(test_bx)
m7 = max(test_by)
m8 = max(test_bz)
m9 = max(test_divb)

if(nvar eq 3) then m1=0d0
err = max([m1, m2, m3, m4, m5, m6, m7, m8, m9])

m1 = max(test_pot(1:nz-2))
m2 = max(test_den(1:nz-2))
m3 = max(test_vx(1:nz-2))
m4 = max(test_vy(1:nz-2))
m5 = max(test_vz(1:nz-2))
m6 = max(test_bx(1:nz-2))
m7 = max(test_by(1:nz-2))
m8 = max(test_bz(1:nz-2))
m9 = max(test_divb(1:nz-2))

if(nvar eq 3) then m1=0d0
err_int = max([m1, m2, m3, m4, m5, m6, m7, m8, m9])


bctest_div = abs(div)

bctest_by = abs(by)

bctest_bxmin = abs(bz + ii*bx)

bctest_bxmax = abs(bz - ii*bx)

bctest_potxmax = abs(dpot +kx*pot - ii*rho*vz/bigQ/sig/fQ)

bctest_potxmin = abs(dpot  -kx*pot - ii*rho*vz/bigQ/sig/fQ)

m1 = max([ bctest_div(0), bctest_div(nz-1)])/norm 
m2 = max([ bctest_by(0), bctest_by(nz-1)])/norm
m3 = max([ bctest_bxmin(0), bctest_bxmax(nz-1)])/norm
m4 = max([ bctest_potxmin(0), bctest_potxmax(nz-1)])/norm

if(nvar eq 3) then m4=0d0
err_bc = max([m1, m2, m3, m4])

print, "k, max error (all), max error (int), bc error", k, err, err_int, err_bc

endfor

end
