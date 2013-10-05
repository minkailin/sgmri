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

pro result, loc=loc, kmode=kmode, nz=nz, xrange=xrange, legend=legend
  
  location =strcompress(loc,/remove_all)

  ii = dcomplex(0d0, 1d0)
  
  kappa2 = 1.0
  
  param = dblarr(5,1)
  file = filepath('params.dat',root_dir='.',subdir=[location])
  openr,1,file
  readf,1, param, format='(5(e22.15,x))'
  close,1
;write(10, fmt='(4(e22.15,x))') zmax, dble(nz), fQ, shear
  zmax   = param(0) 
  lmodes = param(1)
  fq     = param(2)
  shear  = param(3) 
  rm     = param(4)

  zaxis = zmax*dindgen(nz)/(nz-1d0)
  
;focus on a specific mode (among diff kx, bigQ or beta) 
  
;basic state for this mode 
  file = filepath('basic.dat',root_dir='.',subdir=[location])
  nlines = file_lines(file)
  array = dblarr(7, nlines)
  openr,1,file
  readf,1,array,format='(7(e22.15,x))'
  close,1
;write(99, fmt='(4(e22.15,x))') zaxis(j), dnorm(j), valf(j), csq(j)
  z_small = array(0,lmodes*(kmode-1): lmodes*kmode- 1)
  rho     = array(1,lmodes*(kmode-1): lmodes*kmode- 1)
  vasq    = array(2,lmodes*(kmode-1): lmodes*kmode- 1)
  csq     = array(3,lmodes*(kmode-1): lmodes*kmode- 1)
  eta   = array(4,lmodes*(kmode-1): lmodes*kmode- 1)
  drho    = array(5,lmodes*(kmode-1): lmodes*kmode- 1)
  deta    = array(6,lmodes*(kmode-1): lmodes*kmode- 1)
  
;interpolate to fine grid
  rho  = interpol(rho, z_small, zaxis,  /spline)
  drho = interpol(drho, z_small, zaxis,  /spline)
  vasq = interpol(vasq, z_small, zaxis, /spline)
  csq  = interpol(csq, z_small, zaxis, /spline) 
  eta  = interpol(eta, z_small, zaxis, /spline) 
  deta = interpol(deta, z_small, zaxis, /spline)
;  deta = deriv(zaxis, eta) 

;plot elsasser number; omega=1
  elsas = alog10(vasq/eta)
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
  bigQ = eigen(0, kmode-1)
  bta  = eigen(1, kmode-1)
  kx   = eigen(2, kmode-1)
  sig  = dcomplex( eigen(3, kmode-1), eigen(4, kmode-1) )
  sbar = sig -ii*eta*kx^2

  if(bta le 0d0) then bta =1d0 
  
;eigenfunction information 
  array = dblarr(2,4*lmodes)
  filename =  strcompress('eigenfunctions_'+string(kmode, format='(I2)')+'.dat', /remove_all)
  file = filepath(filename,root_dir='.',subdir=[location])
  openr, 1, file 
  readf,1,array,format='(2(e22.15,x))'
  close, 1
  
  vx = dcomplexarr(nz) 
  vy = dcomplexarr(nz)
  vz = dcomplexarr(nz)
  den= dcomplexarr(nz);this is actually w=dp/rho
  pot= dcomplexarr(nz)
  
  dvx = vx
  dvy = vy
  dden= den
  dpot= pot 
  bx  = dcomplexarr(nz)
  by  = dcomplexarr(nz)
  bz  = dcomplexarr(nz) 
  
  
  for j=0, nz-1 do begin
     zbar = zaxis(j)/zmax 
     vx(j)  = 0d0 
     vy(j)  = 0d0 
     den(j) = 0d0 
     pot(j) = 0d0 
     
     dvx(j) = 0d0 
     dvy(j) = 0d0 
     
     for i=1, lmodes do begin
        arr = cheby(2d0*(i-1d0), zbar)
        
        vx(j) +=dcomplex( array(0, i-1),  $
                          array(1, i-1)) $
                *arr(0)         ;this is T_l
        dvx(j) +=dcomplex( array(0, i-1),  $
                           array(1, i-1)) $
                 *arr(1)
        
        vy(j) +=dcomplex( array(0, lmodes + i-1),  $
                          array(1, lmodes + i-1 )) $
                *arr(0)         
        
        dvy(j) +=dcomplex( array(0, lmodes + i-1),  $
                           array(1, lmodes + i-1 )) $
                 *arr(1)
        
        den(j) +=dcomplex( array(0, 2*lmodes + i-1  ),$
                           array(1, 2*lmodes + i-1) ) $
                 *arr(0)
        dden(j) +=dcomplex( array(0, 2*lmodes + i-1  ),$
                            array(1, 2*lmodes + i-1) ) $
                  *arr(1)
        
        pot(j) +=dcomplex( array(0, 3*lmodes + i-1),  $
                           array(1, 3*lmodes + i-1)) $
                 *arr(0)
        
        dpot(j) +=dcomplex( array(0, 3*lmodes + i-1),  $
                            array(1, 3*lmodes + i-1)) $
                  *arr(1)
     endfor
  endfor 
  
  dvx  /= zmax
  dvy  /= zmax
  dden /= zmax 
  dpot /= zmax 
  
  vz = ii*fq*(dden + dpot)/sig 
  
;  bx = -ii*fq*dvx/(sig*zmax)
;  by = -ii*fq*(dvy + ii*shear*dvx/sig)/(sig*zmax)
;  bz = -kx*fq*vx/sig 
  
  g  = ii*sig*vx - 2d0*vy + ii*fQ*kx*(den+pot)
  dg = ii*sig*dvx- 2d0*dvy+ ii*fQ*kx*(dden+dpot)

  bx = fq*dvx + (bta*rho/fQ)*(eta*(drho*g + dg) + deta*g)
  bx *= -ii/sig

 
  g  = ii*sig*vy + 0.5*kappa2*vx
  dg = ii*sig*dvy+ 0.5*kappa2*dvx

  by = fQ*dvy - shear*bx + (bta*rho/fQ)*(eta*(drho*g + dg) + deta*g)
  by /= ii*sbar

  bz = -fQ*kx*vx - (kx*eta*bta*rho/fQ)*( ii*sig*vx - 2d0*vy + ii*fQ*kx*(den+pot) )
  bz /= sig 

  v2  = abs(vx)^2 + abs(vy)^2 + abs(vz)^2
  b2  = abs(bx)^2 + abs(by)^2 + abs(bz)^2 
  gpe = real_part(den*conj(pot))
  
  en1 = b2/(2d0*bta)
  en2 = rho*abs(gpe)/csq 
  en3 = rho*v2/2d0 
  en4 = rho*abs(den)^2/(2d0*csq)

  en1 /= max(en1)
  en2 /= max(en1)
  en3 /= max(en1)

  
  set_plot,'x'
  plot, zaxis, en1, xrange=xrange, ystyle=1, charsize=1.5, thick=2
  oplot, zaxis, en2, linestyle=1, thick=2
  oplot, zaxis, en3, linestyle=2, thick=2
  oplot, [zcrit, zcrit], [-1,1]*1d10, thick=4, linestyle=2
  
  set_plot, 'ps'
  rate_string = strcompress(textoidl('\gamma=')+ string(-imaginary(sig), format='(f4.2)')+textoidl('\Omega'), /remove_all)
  kx_string   = strcompress(textoidl('k_xH=')+ string(kx, format='(f4.2)'), /remove_all)
  q3d_string  = strcompress(textoidl('Q_{3D}=')+ string(bigQ, format='(f5.2)'), /remove_all)
  beta_string = strcompress(textoidl('\beta=')+ string(bta, format='(f5.1)'), /remove_all)
  rm_string = strcompress(textoidl('log\Lambda=')+ string(alog10(rm), format='(f4.1)'), /remove_all)
  
  title = q3d_string +', '+ beta_string  +', '+ rm_string+ ', '   +  kx_string +', '+ rate_string
  
  file = filepath('result.ps',root_dir='.',subdir=[location])
  device, filename=file, $
          bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, en1,xmargin=[4,2],ymargin=[3.,1.75], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=[0,1], xtitle=textoidl('z/H'), $
        linestyle = 0, xtickinterval=xtickinterval, title=title, charthick=6
  oplot, zaxis, en2, thick=4, linestyle=1
  oplot, zaxis, en3, thick=4, linestyle=2
  oplot, [zcrit, zcrit], [-1,1]*1d10, thick=1, linestyle=0
  if keyword_set(legend) then begin
     label   = strarr(3)
     label(0)='magnetic';+textoidl('|\deltaB^2|/2\mu_0')
     label(1)='gravity';+textoidl('|Re(\delta\rho\delta\Phi_d^*)|')
     label(2)='kinetic';+textoidl('|\deltav^2|/2')

     x0=legend(0)
     x1=legend(1)
     y0=legend(2)
     dy=legend(3)
     for j=0, n_elements(label)-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5,charthick=6
     endfor
  endif
  device,/close
end
