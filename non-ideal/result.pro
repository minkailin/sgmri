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

pro result, loc=loc, kmode=kmode, nz=nz, xrange=xrange, legend=legend, yrange=yrange, nvar=nvar, nob=nob, hold=hold, rescale=rescale, chebaxis=chebaxis
  
  if not keyword_set(nvar) then nvar=4
  if not keyword_set(rescale) then rescale=[1,1,1,1]


  location =strcompress(loc,/remove_all)

  ii = dcomplex(0d0, 1d0)
  
  kappa2 = 1.0
  
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


  Q2d = 2d0*fq*bigQ*int_tabulated(zaxis, sqrt(csq)*rho)/int_tabulated(zaxis, rho)^2
  print, Q2d

  lam2= bta*bigQ/Q2d^2

  kx_opt = int_tabulated(zaxis, rho)^2/int_tabulated(zaxis, csq*rho)
  kx_opt/= 2d0*bigQ*fq^2

 ; kx_opt*= 1d0 - 1d0/lam2
 ; kx_opt/= 1d0 + 1d0/bta

  ;print,'kx, kx_opt, growth', kx, kx_opt, -imaginary(sig)
   print,'Q, kx, beta, growth', bigQ, kx, bta, -imaginary(sig)

  if(bta le 0d0) then bta =1d0 
  
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

  if keyword_set(nob) then begin
  bx(*) = 0d0
  by(*) = 0d0
  bz(*) = 0d0
  endif

;stop

; calculate energies 

  v2  = abs(vx)^2 + abs(vy)^2 + abs(vz)^2
  b2  = abs(bx)^2 + abs(by)^2 + abs(bz)^2 
  gpe = real_part(w*conj(pot))
  
  en1 = b2/(2d0*bta)
  en2 = (rho/csq)*abs(gpe)/2d0
  en3 = rho*v2/2d0  
  en4 = rho*abs(w)^2/(2d0*csq)

  etot = en1+en2+en3+en4

;  m1 = mean(etot)
   
   m1 = int_tabulated(zaxis, etot)/zmax 

  en1 /= m1
  en2 /= m1
  en3 /= m1
  en4 /= m1

  en1 *= rescale[0]
  en2 *= rescale[1]
  en3 *= rescale[2]
  en4 *= rescale[3]

  set_plot,'x'
  window,0, retain=2
  plot, zaxis, en1, xrange=xrange, ystyle=1, charsize=1.5, thick=2, yrange=yrange;[0,max([en1,en2,en3])]
  oplot, zaxis, en2, linestyle=1, thick=2
  oplot, zaxis, en3, linestyle=2, thick=2
  oplot, zaxis, en4, linestyle=3, thick=2
  oplot, [zcrit, zcrit], [-1,1]*1d10, thick=4, linestyle=2


  
;  div = ii*kx*vx + dvz  
 
;stresses  
  maxwell  = - real_part(bx*conj(by))*vasq  
  reynolds =   real_part(vx*conj(vy))

  stress = reynolds + maxwell 


;vorticity
 vort2 = abs(dvy)^2 +abs(ii*kx*vz - dvx)^2 + kx*kx*abs(vy)^2
 vort2/=max(vort2)

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
  
  bctest_bx = abs(bx)
  
  bctest_by = abs(by)

;wall bc below  
  bctest_vz = abs(vz)
  
  bctest_potxmax = abs(dpot + kx*pot)

if keyword_set(hold) then stop

  set_plot, 'ps'
  rate_string = strcompress(textoidl('\gamma=')+ string(-imaginary(sig), format='(f4.2)')+textoidl('\Omega'), /remove_all)
  kx_string   = strcompress(textoidl('k_xH=')+ string(kx, format='(f4.2)'), /remove_all)
  q3d_string  = strcompress(textoidl('Q=')+ string(bigQ, format='(f5.2)'), /remove_all)
  beta_string = strcompress(textoidl('\beta=')+ string(bta, format='(f5.1)'), /remove_all)
  rm_string = strcompress(textoidl('log\Lambda_0=')+ string(alog10(rm), format='(f4.1)'), /remove_all)
  
  title = q3d_string +', '+ beta_string  +', '+ rm_string+ ', '   +  kx_string +', '+ rate_string
  
  file = filepath('result.ps',root_dir='.',subdir=[location])
  device, filename=file, $
          bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, en1,xmargin=[6,2],ymargin=[3.,1.75], ystyle=1  $
        ,charsize=1.5, thick=4, xrange=xrange,yrange=yrange, xtitle=textoidl('z/H'), $
        linestyle = 0, xtickinterval=xtickinterval, title=title, charthick=6, ytitle=textoidl('E_i/<E>') 
  oplot, zaxis, en2, thick=4, linestyle=1
  oplot, zaxis, en3, thick=4, linestyle=2
  oplot, zaxis, en4, thick=4, linestyle=3
  oplot, [zcrit, zcrit], [-1,1]*1d10, thick=1, linestyle=0
  if keyword_set(legend) then begin
     label   = strarr(4)
     label(0)='magnetic';+textoidl('|\deltaB^2|/2\mu_0')
     label(1)='gravity';+textoidl('|Re(\delta\rho\delta\Phi_d^*)|')
     label(2)='kinetic';+textoidl('|\deltav^2|/2')
     label(3)='thermal'

     for i=0, 3 do begin
        if (rescale(i) ne 1d0) then begin
        label(i) += textoidl('\times10^{'+string(alog10(rescale[i]),format='(I2)')+'}')
        endif
     endfor


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
