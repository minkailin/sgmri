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

pro get_stress, loc=loc, modes=modes, nz=nz, nvar=nvar, kmode= kmode, chebgrid=chebgrid
  common stress, maxwell, reynolds
  
  ii = dcomplex(0d0, 1d0)
  kappa2 = 1.0
  
     
  file = filepath('params.dat',root_dir='.',subdir=[loc])
  nlines = file_lines(file)
  param = dblarr(5,nlines)
  openr,1,file
  readf,1, param, format='(5(e22.15,x))'
  close,1
  zmax   = param(0,kmode-1) 
  lmodes = param(1,kmode-1)
  fq     = param(2,kmode-1)
  shear  = param(3,kmode-1) 
  rm     = param(4,kmode-1)
  
  if (chebgrid eq 0)  then begin
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

  vx = dcomplexarr(nz)
  vy = dcomplexarr(nz)
  vz = dcomplexarr(nz)
  w  = dcomplexarr(nz)
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


;basic state for this mode 
  file = filepath('basic.dat',root_dir='.',subdir=[loc])
  nlines = file_lines(file)
  array = dblarr(7, nlines)
  openr,1,file
  readf,1,array,format='(7(e22.15,x))'
  close,1
     
  z_small = array(0,lmodes*(kmode-1): lmodes*kmode- 1)
  rho     = array(1,lmodes*(kmode-1): lmodes*kmode- 1)
  vasq    = array(2,lmodes*(kmode-1): lmodes*kmode- 1)
  csq     = array(3,lmodes*(kmode-1): lmodes*kmode- 1)
  eta     = array(4,lmodes*(kmode-1): lmodes*kmode- 1)
  drho    = array(5,lmodes*(kmode-1): lmodes*kmode- 1)
  deta    = array(6,lmodes*(kmode-1): lmodes*kmode- 1)

;interpolate to fine grid
  rho  = interpol(rho, z_small, zaxis,  /spline)
  drho = interpol(drho, z_small, zaxis,  /spline)
  vasq = interpol(vasq, z_small, zaxis, /spline)
  csq  = interpol(csq, z_small, zaxis, /spline) 
  eta  = interpol(eta, z_small, zaxis, /spline) 
  deta = interpol(deta, z_small, zaxis, /spline)
  
;eigenvalue information 
  file = filepath('eigenvalues.dat',root_dir='.',subdir=[loc])
  nlines = file_lines(file)
  eigen = dblarr(6, nlines)
  openr,1,file
  readf,1, eigen, format='(6(e22.15,x))'
  close,1
  
  bigQ = eigen(0, kmode-1)
  bta  = eigen(1, kmode-1)
  kx   = eigen(2, kmode-1)
  sig  = dcomplex( eigen(3, kmode-1), eigen(4, kmode-1))
  sbar = sig -ii*eta*kx^2
  
;eigenfunction information 
  array = dblarr(2,nvar*lmodes)
  filename =  strcompress('eigenfunctions_'+string(kmode, format='(I3)')+'.dat', /remove_all)
  file = filepath(filename,root_dir='.',subdir=[loc])
  openr, 1, file 
  readf,1,array,format='(2(e22.15,x))'
  close, 1
  
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
  
; calculate stress



  temp     = real_part(vx*conj(vy))*rho
  reynolds = int_tabulated(zaxis, temp) 
  
  temp    =  real_part(bx*conj(by))
  maxwell = -int_tabulated(zaxis, temp)/bta 

end

pro compare_stress, loc=loc, nz=nz, xrange=xrange, legend=legend, label=label, title=title, yrange=yrange, var=var, nvar=nvar, chebaxis=chebaxis
  common stress, maxwell, reynolds

  location = strcompress(loc,/remove_all)  
  ncases   = n_elements(loc)
  
  if keyword_set(chebaxis) then begin
  chebgrid = 1
  endif else chebgrid = 0 


  file = filepath('eigenvalues.dat',root_dir='.',subdir=[location(0)])
  nlines = file_lines(file)
  
  eigen = dblarr(6,nlines)
  
  openr,1,file
  readf,1, eigen, format='(6(e22.15,x))'
  close,1
  
  case var of
     
     'toomre':  begin
        xaxis = eigen(0,*)      ; toomre q
        xtitle = 'Q'
     end
     'kx':      begin
        xaxis = eigen(2,*)      ; kx
        xtitle=textoidl('k_xH')
     end
     'beta':   begin
        xaxis = eigen(1,*)      ; kx
        xtitle=textoidl('\beta')
     end
     
  endcase
  
  stress = dblarr(nlines)

  for i=0, nlines-1 do begin
     ip1 = i+1
     
     get_stress, loc=location(0), nz=nz, nvar=nvar, kmode=ip1, chebgrid=chebgrid 

     stress(i) = maxwell/reynolds

  endfor

; ytitle = '-<Re(\deltav_{Ax}\deltav_{Ay}^{*})>/<Re(\deltav_x\deltav_y^{*})>'

  set_plot, 'ps'
  file = filepath('compare_stress.ps',root_dir='.',subdir=[location(0)])
  device, filename=file, $
          bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  
  plot, xaxis, stress, xmargin=[8,2],ymargin=[3.5,0.5]  $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, $
        linestyle = 0, xtickinterval=xtickinterval, ytitle=textoidl('M_{xy}/R_{xy}'), charthick=6, yrange=yrange, ystyle=1, xstyle=1 
  for i=1, ncases-1 do begin
     file = filepath('eigenvalues.dat',root_dir='.',subdir=[location(i)])
     nlines = file_lines(file)
     eigen = dblarr(6,nlines)
     
     openr,1,file
     readf,1, eigen, format='(6(e22.15,x))'
     close,1
     case var of
        'toomre':  begin
           xaxis = eigen(0,*)   ; toomre q
           xtitle = 'Q'
        end
        'kx':      begin
           xaxis = eigen(2,*)   ; kx
           xtitle=textoidl('k_xH')
        end
        'beta':   begin
           xaxis = eigen(1,*)   ; kx
           xtitle=textoidl('\beta')
        end
     endcase

     stress = dblarr(nlines)

     for j=0, nlines-1 do begin
        ip1 = j+1
        
        get_stress, loc=location(i), nz=nz, nvar=nvar, kmode=ip1, chebgrid=chebgrid

       stress(j) = maxwell/reynolds 

     endfor
     
     oplot, xaxis , stress , thick=4, linestyle=i

  endfor
  if keyword_set(legend) then begin
     x0=legend(0)
     x1=legend(1)
     y0=legend(2)
     dy=legend(3)
     for j=0, n_elements(label)-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5,charthick=6
     endfor
     if keyword_set(title) then xyouts, x1, y0+dy, textoidl(title),charsize=1.5,charthick=6
  endif
  device,/close

end

