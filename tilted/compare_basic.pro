pro compare_basic, loc=loc, cases=cases, nz=nz, xrange=xrange, legend=legend, label=label, title=title, yrange=yrange
  
  location = strcompress(loc,/remove_all)  
  ncases   = n_elements(cases)
  
  density = dblarr(ncases, nz)
  elsasser= dblarr(ncases, nz)
  resis = dblarr(ncases, nz)
  lambda = dblarr(ncases, nz)
 
  for i=0, ncases-1 do begin
     kmode = cases(i)
     
     file = filepath('params.dat',root_dir='.',subdir=[location])
     nlines = file_lines(file)
     param = dblarr(5,nlines)
     openr,1,file
     readf,1, param, format='(5(e22.15,x))'
     close,1
     zmax   = param(0,kmode-1) 
     lmodes = param(1,kmode-1)
     fq     = param(2,kmode-1)
     shear  = param(3,kmode-1) 
     rm     = param(4,kmode-1)          ;this is actually midplane elsasser
     
     zaxis = zmax*dindgen(nz)/(nz-1d0)

;basic state for this mode 
     file = filepath('basic.dat',root_dir='.',subdir=[location])
     nlines = file_lines(file)
     array = dblarr(9, nlines)
     openr,1,file
     readf,1,array,format='(9(e22.15,x))'
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

     density(i,*)  =  rho
     elsasser(i,*) = alog10(fq*fq*vasq/eta)
     

;ideal and resistive mri wavelengths. FOR POLYTROPIC DISK ONLY right now
     for j=0, nz-1 do begin
        lambda(i,j) =  max([ fq*sqrt(vasq(j))*4d0*!dpi/sqrt(15d0), (2d0*!dpi/sqrt(3d0))*eta(j)/(fq*sqrt(vasq(j)))])
     endfor
     
  endfor

  set_plot, 'ps'
  file = filepath('lambda.ps',root_dir='.',subdir=[location])
  device, filename=file, $
          bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, lambda(0,*),xmargin=[6,2],ymargin=[3.,1.75], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=textoidl('z/H'), $
        linestyle = 0, xtickinterval=xtickinterval, ytitle=textoidl('\lambda/2H'), charthick=6, yrange=yrange
;  temp = min(abs(lambda(0,*)-1d0), grid)
;  if((grid lt nz-1) and (grid gt 0)) then oplot,[zaxis(grid),zaxis(grid)],[-1,1]*1d10,thick=1
  for i=1, ncases-1 do begin
     oplot, zaxis, lambda(i,*), thick=4, linestyle=i
;     temp = min(abs(lambda(i,*)-1d0), grid)
;  if((grid lt nz-1) and (grid gt 0)) then   oplot,[zaxis(grid),zaxis(grid)],[-1,1]*1d10,linestyle=i,thick=1
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
  oplot, [0,10], [1,1], linestyle=1, thick=1
  device,/close



;plot elsasser number; omega=1
  set_plot, 'ps'
  file = filepath('elsasser.ps',root_dir='.',subdir=[location])
  device, filename=file, $
          bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, elsasser(0,*),xmargin=[6,2],ymargin=[3.,1.75], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=textoidl('z/H'), $
        linestyle = 0, xtickinterval=xtickinterval, ytitle=textoidl('Log_{10}\Lambda'), charthick=6, yrange=yrange
  temp = min(abs(elsasser(0,*)), grid)
  oplot,[zaxis(grid),zaxis(grid)],[-1,1]*1d10,linestyle=0,thick=1
  for i=1, ncases-1 do begin
     oplot, zaxis, elsasser(i,*), thick=4, linestyle=i
     temp = min(abs(elsasser(i,*)), grid)
     oplot,[zaxis(grid),zaxis(grid)],[-1,1]*1d10,linestyle=i,thick=1
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
 
  set_plot, 'ps'
  file = filepath('density.ps',root_dir='.',subdir=[location])
  device, filename=file, $
          bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, density(0,*),xmargin=[6,2],ymargin=[3.,1.75], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=textoidl('z/H'), $
        linestyle = 0, xtickinterval=xtickinterval, ytitle=textoidl('\rho/\rho_0'), charthick=6
  for i=1, ncases-1 do begin
     oplot, zaxis, density(i,*), thick=4, linestyle=i
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
