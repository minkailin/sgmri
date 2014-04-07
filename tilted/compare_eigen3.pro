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

pro get_energies, loc=loc, nz=nz, nvar=nvar, kmode= kmode, chebaxis=chebaxis
  common energy, emag, epot, ekin, etherm
  
  ii = dcomplex(0d0, 1d0)
  kappa2 = 1.0
       
  file = filepath('params.dat',root_dir='.',subdir=[loc])
  nlines = file_lines(file)
  param = dblarr(6,nlines)
  openr,1,file
  readf,1, param, format='(6(e22.15,x))'
  close,1
  zmax   = param(0,kmode(0)-1) 
  lmodes = param(1,kmode(0)-1)
  fq     = param(2,kmode(0)-1)
  shear  = param(3,kmode(0)-1) 
  rm     = param(4,kmode(0)-1)  ;this is actually midplane elsasser
  tilt   = param(5,kmode(0)-1)
  
  if (chebaxis eq 0) then begin
     zaxis = -zmax + 2d0*zmax*dindgen(nz)/(nz-1d0)
  endif else begin
     nz = lmodes 
     zaxis = dblarr(nz)
     lmax = nz-1 
     for j=0, nz-1 do begin
        zaxis(j) = -zmax*cos(!dpi*j/lmax)
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
  
  dlogcsq = deriv(zaxis, alog(csq))

;eigenvalue information 
  file = filepath('eigenvalues.dat',root_dir='.',subdir=[loc])
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
        arr = cheby((i-1d0), zbar)
        
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

  D1V = drho*V + dV
  
  vz = (ii/sig)*(fQ*dwtilde + tilt*V)
  dvz= (ii/sig)*(fQ*d2wtilde + tilt*dV)
  
;calculate  mag field  
  bx = fQ*dvx + eta*(D1U + ii*tilt*kx*V)/(fQ*vasq)
  bx/= ii*sig
  
  by = tilt*(ii*sig*w/csq + fQ*drho*vz) + fQ*dvy - shear*bx + eta*D1V/(fQ*vasq)
  by /= ii*sbar
  
  bz = fQ*kx*vx + eta*kx*( U/(fQ*vasq) + ii*tilt*kx*by )
  bz/= -sig

;energies  
  v2  = abs(vx)^2 + abs(vy)^2 + abs(vz)^2
  b2  = abs(bx)^2 + abs(by)^2 + abs(bz)^2 
  gpe = real_part(w*conj(pot))
  
  en1 = b2/(2d0*bta)
  en2 = rho*abs(gpe)/csq/2d0 
  en3 = rho*v2/2d0 
  en4 = rho*abs(w)^2/(2d0*csq)
  etot = en1+en2+en3+en4
 
  
  emag = mean(en1)
  epot = mean(en2)
  ekin = mean(en3)
  etherm=mean(en4) 

end

pro compare_eigen3, loc=loc, nz=nz, xrange=xrange, legend=legend, label=label, title=title, yrange=yrange, var=var, sample=sample, psym=psym, color=color, nvar=nvar, solid=solid, plotrange=plotrange, chebaxis=chebaxis 
  common energy, emag, epot, ekin, etherm                          

  location = strcompress(loc,/remove_all)  
  ncases   = n_elements(loc)
  samp_arr = dblarr(ncases)
  
  if keyword_set(chebaxis) then begin
     chebgrid = 1
  endif else chebgrid = 0 

  if not keyword_set(nvar) then begin
    nvar = indgen(ncases) 
    nvar(*) = 4
  endif

  if not keyword_set(solid) then begin
    solid = intarr(ncases)
  endif

  if not keyword_set(plotrange) then plotrange=[0,1]

  if not keyword_set(sample) then begin
     samp_arr(*) = 1
  endif else begin
     if (n_elements(sample) eq 1) then begin
        samp_arr(*) = sample
     endif else samp_arr(*) = sample(*)
  endelse

  xstart = dblarr(ncases)
  ystart = dblarr(ncases)

  xfinal = dblarr(ncases)
  yfinal = dblarr(ncases)

  

  i = 0 
  
  file   = filepath('eigenvalues.dat',root_dir='.',subdir=[location(i)])
  nlines = file_lines(file)
  modenum   = dindgen(nlines)+1   ;store mode number 

  eigen = dblarr(7,nlines)  
  openr,1,file
  readf,1, eigen, format='(7(e22.15,x))'
  close,1
  
  case var of     
     'toomre':  begin
        grid = 0
        xtitle = 'Q'
     end
     'kx':      begin
        grid = 3 
        xtitle=textoidl('k_xH')
     end
     'beta':   begin
        grid = 1
        xtitle=textoidl('\beta')
     end
     'tilt':   begin
        grid = 2      ; kx
        xtitle=textoidl('B_y/B_z')
     end
     
     'rm':     
  endcase  
  if(var ne 'rm') then begin
     xaxis =  eigen(grid,*)
  endif else begin
     file = filepath('params.dat',root_dir='.',subdir=[location(i)])
     params = dblarr(6,nlines)
     openr,1,file
     readf,1, params, format='(6(e22.15,x))'
     close,1
     xaxis = params(4,*)
     xtitle=textoidl('\Lambda_0')
  endelse

  data  = -eigen(5,*)*100d0     ;growth rate
  
  nsamp  = samp_arr(i)*nlines

  xaxis    = congrid(xaxis(*),   nsamp)
  data     = congrid(data(*),    nsamp)
  modenum  = congrid(modenum(*), nsamp)
  
  

  if not keyword_set(color) then begin
     xmargin=[6,2]
  endif else begin
     x = cos(2d0*!dpi*dindgen(16)/15d0)
     y = sin(2d0*!dpi*dindgen(16)/15d0)
     usersym, x, y, /fill 
     xmargin=[6,6]
  endelse 
  
  cut = 2.1
;  cut = 2

  set_plot, 'ps'
  file = filepath('compare_growth3.ps',root_dir='.',subdir=[location(0)])
  device, filename=file, $
          bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, xaxis, data,xmargin=xmargin,ymargin=[4.,1]  $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, $
        linestyle = 0, xtickinterval=xtickinterval, ytitle=textoidl('10^2\gamma/\Omega'), charthick=6, yrange=yrange, ystyle=1, xstyle=1,/nodata 
  
  
;  temp = min(abs(xaxis - cut), nsamp)
  xfinal(0) = xaxis(nsamp-1)
  yfinal(0) = data(nsamp-1)

  xstart(0) = xaxis(0)
  ystart(0) = data(0)

  if(solid(0) eq 0 ) then begin
  if keyword_set(color) then loadct, color,/silent
  for j=0, nsamp-1 do begin
     get_energies, loc=location(i), nz=nz, nvar=nvar(0), kmode=modenum(j), chebaxis=chebgrid 

     if not keyword_set(color) then begin
        if (emag gt epot) then begin
           temp=psym(0)
        endif else temp=psym(1)
        oplot, [1,1]*xaxis(j), [1,1]*data(j), psym=temp,thick=2
     endif else begin 
        col = epot/(epot+emag)*256/max(plotrange)
        oplot, [1,1]*xaxis(j), [1,1]*data(j), psym=8,thick=4, color = col
     endelse
  endfor
  endif else begin
  loadct, 0, /silent
  oplot, xaxis, data, thick=4 
  endelse   


  for i=1, ncases-1 do begin
     file = filepath('eigenvalues.dat',root_dir='.',subdir=[location(i)])
     nlines = file_lines(file)
     eigen = dblarr(7,nlines)
     modenum   = dindgen(nlines)+1   ;store mode number 
     
     openr,1,file
     readf,1, eigen, format='(7(e22.15,x))'
     close,1

     if(var ne 'rm') then begin
        xaxis =  eigen(grid,*)
     endif else begin
        file = filepath('params.dat',root_dir='.',subdir=[location(i)])
        params = dblarr(6,nlines)
        openr,1,file
        readf,1, params, format='(6(e22.15,x))'
        close,1
        xaxis = params(4,*)
     endelse
     
     data  = -eigen(5,*)*100d0  ;growth rate
     
     nsamp  = samp_arr(i)*nlines
     
     xaxis    = congrid(xaxis(*),   nsamp)
     data     = congrid(data(*),    nsamp)
     modenum  = congrid(modenum(*), nsamp)

;     temp = min(abs(xaxis - cut), nsamp)
     xfinal(i) = xaxis(nsamp-1)
     yfinal(i) = data(nsamp-1)
     
     xstart(i) = xaxis(0)
     ystart(i) = data(0)

     if(solid(i) eq 0 ) then begin
        if keyword_set(color) then loadct, color,/silent
     for j=0, nsamp-1 do begin
        get_energies, loc=location(i), nz=nz, nvar=nvar(i), kmode=modenum(j), chebaxis=chebgrid 
        if not keyword_set(color) then begin
        if (emag gt epot) then begin
           temp=psym(0)
        endif else temp=psym(1)
        oplot, [1,1]*xaxis(j), [1,1]*data(j), psym=temp,thick=2
        endif else begin
        col = epot/(epot+emag)*256/max(plotrange)
        oplot, [1,1]*xaxis(j), [1,1]*data(j), psym=8,thick=4, color = col
        endelse
     endfor     
     endif else begin
     loadct, 0, /silent
     oplot, xaxis, data, thick=4                          
     endelse

  endfor


  if keyword_set(color) then begin
     loadct, color,/silent
     colorbar, position=[0.91, 0.192, 0.93, 0.95],/vertical,/right,format='(f4.2)', thick=4, charsize=1.5, range=plotrange
  endif
  
  if keyword_set(legend) then begin
     x0=legend(0)
     x1=legend(1)
     y0=legend(2)
     dy=legend(3)

     nlabel=n_elements(label)

     if not keyword_set(color) then begin
        for j=0, nlabel-2 do begin
           oplot, [x0,x0], [y0,y0]-dy*j, thick=4, psym=psym(j)
           xyouts, x0, y0-dy*j,textoidl(label(j)),charsize=1.5,charthick=6
        endfor
     endif else begin
        loadct, color,/silent
        for j=0, nlabel-2 do oplot, [x1,x1], [y0,y0]-dy*j, thick=4, psym=8, color=j*256
        loadct, 0, /silent
        for j=0, nlabel-2 do  xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5,charthick=6
     endelse

      ;; if keyword_set(solid) then begin
      ;;    j=nlabel-1
      ;;    oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=0
      ;;    xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5,charthick=6
      ;; endif

;     x0=1.2
     if keyword_set(title) then xyouts, x0, y0+dy, textoidl(title),charsize=1.5,charthick=6


      ;; names ='B_y/B_z=' + ['0','1','2','3']
      ;; names = textoidl(names)
      ;; loadct, 0, /silent

      ;; xfinal *= 1.03

      ;; i=0
      ;; xyouts, xfinal(i), yfinal(i)*.97, names(i), charsize=1.5,charthick=6
     
      ;; i=1
      ;; xyouts, xfinal(i), yfinal(i)*.98, names(i), charsize=1.5,charthick=6
     
      ;; i=2
      ;; xyouts, xfinal(i), yfinal(i), names(i), charsize=1.5,charthick=6
     
      ;;  i=3
      ;; xyouts, xfinal(i), yfinal(i), names(i), charsize=1.5,charthick=6


     

;      names = ['odd','even','overstable']
     names = 'E_m(0): '+['local max.','local min.']
     names = textoidl(names)
     loadct, 0, /silent
     xstart*=1.5
;      for i=0, 1 do begin
     
     i=0
     xyouts, xstart(i), ystart(i)*1.01, names(i), charsize=1.5,charthick=6
     i=1
     xyouts, xstart(i), ystart(i), names(i), charsize=1.5,charthick=6

;      endfor
;      xstart *= 1.4
;      ystart *= 0.9
;      xyouts, xstart(3), ystart(3), names(2), charsize=1.5,charthick=6

  endif
  device,/close
  
  loadct,0,/silent
end
