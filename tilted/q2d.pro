function get_q2d, location

  file = filepath('params.dat',root_dir='.',subdir=[location])
  nlines = file_lines(file)
  param = dblarr(5,nlines)
  openr,1,file
  readf,1, param, format='(5(e22.15,x))'
  close,1
  
  file = filepath('eigenvalues.dat',root_dir='.',subdir=[location])
  nk   = file_lines(file)
  
  eigen = dblarr(6,nk)  
  q3d2d = dblarr(2,nk)
  
  openr,1,file
  readf,1, eigen, format='(6(e22.15,x))'
  close,1
  q3d2d(0,*)= eigen(0,*)
  
                                ;basic state for this mode 
  file = filepath('basic.dat',root_dir='.',subdir=[location])
  nlines = file_lines(file)
  basic = dblarr(7, nlines)
  openr,1,file
  readf,1,basic,format='(7(e22.15,x))'
  close,1
  
  for k=0, nk-1 do begin

     kmode = k+1
     
     lmodes = param(1,kmode-1)
     fq     = param(2,kmode-1) 
     
     z_small = basic(0,lmodes*(kmode-1): lmodes*kmode- 1)
     rho     = basic(1,lmodes*(kmode-1): lmodes*kmode- 1)
     csq     = basic(3,lmodes*(kmode-1): lmodes*kmode- 1) 

     q3d2d(1,k) = 2d0*fq*q3d2d(0,k)*int_tabulated(z_small, sqrt(csq)*rho)/int_tabulated(z_small, rho)^2
  endfor

  return, q3d2d
end

pro q2d, loc=loc, xrange=xrange, title=title, yrange=yrange, legend=legend, label=label
  
  location = strcompress(loc,/remove_all)  
  ncases = n_elements(location)
;assume kappa=1
  
  
;plot Q2d and Q 
  set_plot, 'ps'
  file = filepath('q2d.ps',root_dir='.',subdir=[location(0)])
  device, filename=file, $
          bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  q3d2d = get_q2d(location(0))
  plot, q3d2d(0,*), q3d2d(1,*),xmargin=[5,2],ymargin=[3,1]  $
        ,charsize=1.5, thick=4, xrange=xrange, $
        linestyle = 0, xtickinterval=xtickinterval, ytitle=textoidl('Q_{2D}'), charthick=6, yrange=yrange, ystyle=1, xstyle=1 $
        , xtitle='Q'
  for j=1, ncases-1 do begin
     q3d2d = get_q2d(location(j))
     oplot, q3d2d(0,*), q3d2d(1,*), thick=4, linestyle=j
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
  endif

  device,/close

;stop
end
