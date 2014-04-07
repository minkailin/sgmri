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
  
  return, [t_l, dt_l, d2t_l]
end

pro error, loc=loc, nk=nk, nz=nz, nvar=nvar, chebaxis=chebaxis
  
  if not keyword_set(nvar) then nvar=4


  location =strcompress(loc,/remove_all)

  ii = dcomplex(0d0, 1d0)
  
  kappa2 = 1.0

  for k=1, nk do begin
  
  kmode = k 
  
  file = filepath('params.dat',root_dir='.',subdir=[location(0)])
  nlines = file_lines(file)
  param = dblarr(5,nlines)
  openr,1,file
  readf,1, param, format='(5(e22.15,x))'
  close,1
  zmax   = param(0,kmode(0)-1) 
  lmodes = param(1,kmode(0)-1)
  fq     = param(2,kmode(0)-1)
  shear  = param(3,kmode(0)-1) 
  rm     = param(4,kmode(0)-1)     ;this is actually midplane elsasser
 
if not keyword_set(chebaxis) then begin 
  zaxis = zmax*dindgen(nz)/(nz-1d0)
endif else begin
;   !setup Z axis. these are extrema of T_lmax(Z/zmax) plus end point

  nz = lmodes   

  zaxis = dblarr(nz) 
  haf_lmax = nz-1
  lmax = 2*haf_lmax
  zaxis(0)  = 0d0
  zaxis(nz-1) = zmax
  for j = haf_lmax+1, lmax-1 do begin
     zaxis(j-haf_lmax) = -zmax*cos(j*!dpi/lmax)
     ;print, zaxis(j-haf_lmax)
  endfor
endelse

  
;basic state for this mode 
  file = filepath('basic.dat',root_dir='.',subdir=[location])
  nlines = file_lines(file)
  array = dblarr(7, nlines)
  openr,1,file
  readf,1,array,format='(7(e22.15,x))'
  close,1

  z_small = array(0,lmodes*(kmode(0)-1): lmodes*kmode(0)- 1)
  rho     = array(1,lmodes*(kmode(0)-1): lmodes*kmode(0)- 1)
  vasq    = array(2,lmodes*(kmode(0)-1): lmodes*kmode(0)- 1)
  csq     = array(3,lmodes*(kmode(0)-1): lmodes*kmode(0)- 1)
  eta   = array(4,lmodes*(kmode(0)-1): lmodes*kmode(0)- 1)
  drho    = array(5,lmodes*(kmode(0)-1): lmodes*kmode(0)- 1)
  deta    = array(6,lmodes*(kmode(0)-1): lmodes*kmode(0)- 1)


;interpolate to fine grid
  rho  = interpol(rho, z_small, zaxis,  /spline)
  drho = interpol(drho, z_small, zaxis,  /spline)
  vasq = interpol(vasq, z_small, zaxis, /spline)
  csq  = interpol(csq, z_small, zaxis, /spline) 
  eta  = interpol(eta, z_small, zaxis, /spline) 
  deta = interpol(deta, z_small, zaxis, /spline)


;plot elsasser number; omega=1
  elsas = alog10(fQ^2*vasq/eta)
  temp = min(abs(elsas),grid)
  zcrit = zaxis(grid)
  set_plot, 'ps'
  file = filepath('elsas.ps',root_dir='.',subdir=[location])
  device, filename=file, $
          bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, elsas,xmargin=[6,2],ymargin=[3.,1.75], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=textoidl('z/H'), $
        linestyle = 0, xtickinterval=xtickinterval, ytitle=textoidl('Log_{10}\Lambda'), charthick=6
  oplot, [zcrit, zcrit], [-1,1]*1d10, thick=4, linestyle=2
  device,/close
  
;eigenvalue information 
  file = filepath('eigenvalues.dat',root_dir='.',subdir=[location])
  nlines = file_lines(file)
  eigen = dblarr(6, nlines)
  openr,1,file
  readf,1, eigen, format='(6(e22.15,x))'
  close,1
;write(20, fmt='(6(e22.15,x))') bigQ, beta, kx, dble(sig), dimag(sig), condition
  bigQ = eigen(0, kmode(0)-1)
  bta  = eigen(1, kmode(0)-1)
  kx   = eigen(2, kmode(0)-1)
  sig  = dcomplex( eigen(3, kmode(0)-1), eigen(4, kmode(0)-1) )
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
        arr = cheby(2d0*(i-1d0), zbar)
        
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

  vz = ii*fq*(dw + dpot)/sig 
  dvz= ii*fq*(d2w + d2pot)/sig



;calculate  mag field  
  g  = ii*sig*vx - 2d0*vy + ii*fQ*kx*(w+pot)
  dg = ii*sig*dvx- 2d0*dvy+ ii*fQ*kx*(dw+dpot)

  bx = fq*dvx + (bta*rho/fQ)*(eta*(drho*g + dg) + deta*g)
  bx *= -ii/sig
 
  g  = ii*sig*vy + 0.5*kappa2*vx
  dg = ii*sig*dvy+ 0.5*kappa2*dvx

  by = fQ*dvy - shear*bx + (bta*rho/fQ)*(eta*(drho*g + dg) + deta*g)
  by /= ii*sbar

  bz = -fQ*kx*vx - (kx*eta*bta*rho/fQ)*( ii*sig*vx - 2d0*vy + ii*fQ*kx*(w+pot) )
  bz /= sig 
;calculate first deriv of magn field
  dbx = (bta*rho/fq)*(ii*sig*vx - 2d0*vy + ii*fq*kx*(w+pot)) + ii*kx*bz
  dby = (bta*rho/fq)*(ii*sig*vy + 0.5*kappa2*vx)

  g  = ii*sig*vx - 2d0*vy + ii*fQ*kx*(w+pot)
  dg = ii*sig*dvx- 2d0*dvy+ ii*fQ*kx*(dw+dpot)

  dbz = fq*kx*dvx + (kx*bta*rho/fq)*( (deta + eta*drho)*g + eta*dg) 
  dbz /= -sig

;calculate second deriv of mag field 
  g  = ii*sig*vx - 2d0*vy + ii*fQ*kx*(w+pot)
  dg = ii*sig*dvx- 2d0*dvy+ ii*fQ*kx*(dw+dpot)
  d2bx= (bta*rho/fq)*(drho*g + dg) + ii*kx*dbz

  g  = ii*sig*vy + 0.5*kappa2*vx
  dg = ii*sig*dvy+ 0.5*kappa2*dvx
  d2by= (bta*rho/fq)*(drho*g + dg)

  g  = ii*sig*vx - 2d0*vy + ii*fQ*kx*(w+pot)
  dg = ii*sig*dvx- 2d0*dvy+ ii*fQ*kx*(dw+dpot)
  d2g= ii*sig*d2vx- 2d0*d2vy+ ii*fQ*kx*(d2w+d2pot)


  d2bz = -ii*kx*dbx

;tests   
  wtilde = pot + w
  dwtilde = dpot + dw
 
  norm = max(abs(vx))

  test_pot = abs(d2pot - kx*kx*pot - rho*w/csq/fq^2/bigQ)
  test_den = abs(ii*sig*w/fq/csq + ii*kx*vx + drho*vz + dvz)
  
  test_vx = abs(ii*sig*vx - 2*vy + ii*fQ*kx*wtilde - fQ*vasq*(dbx - ii*kx*(bz)))
  test_vy = abs(ii*sig*vy + 0.5*kappa2*vx - fq*vasq*dby)
  test_vz = abs(ii*sig*vz + fQ*dwtilde)

  test_bx = abs(ii*sbar*bx - fQ*dvx - eta*d2bx -  deta*dbx + ii*kx*deta*bz )
  test_by = abs(ii*sbar*by - fq*dvy + shear*bx - eta*d2by - deta*dby)
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

  bctest_bx = abs(bx)
  
  bctest_by = abs(by)

;wall bc below  
  bctest_vz = abs(vz)
  
  bctest_potxmax = abs(dpot + kx*pot)

m1 =  bctest_bx(nz-1)/norm   
m2 =  bctest_by(nz-1)/norm
m3 =  bctest_vz(nz-1)/norm
m4 =  bctest_potxmax(nz-1)/norm

if(nvar eq 3) then m4=0d0
err_bc = max([m1, m2, m3, m4])

print, "k, max error (all), max error (int), bc error", k, err, err_int, err_bc


endfor
end
