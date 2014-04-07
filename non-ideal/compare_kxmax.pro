pro compare_kxmax, loc=loc, files=files, nz=nz, xrange=xrange, legend=legend, label=label, title=title, yrange=yrange 
; plot optimal kx (max growth rate) as function of  beta 

  location = strcompress(loc,/remove_all);+strcompress('eigenvalues.dat.'+string(files),/remove_all) 
  ncases   = n_elements(files)

  array = dblarr(2, ncases)
 
  for i = 0, ncases-1 do begin 
     
      
     file = filepath(strcompress('eigenvalues.dat.'+string(files(i)),/remove_all),root_dir='.',subdir=[location])
     nlines = file_lines(file)

     eigen = dblarr(6,nlines)
 
     openr,1,file
     readf,1, eigen, format='(6(e22.15,x))'
     close,1

     rates   = -eigen(4,*)*100d0 ;growth rate
     kxaxis = eigen(2,*) 
     btaaxis= eigen(1,*) 

     temp  = max(rates, grid)
;     kxopt = kxaxis(grid)

     result = poly_fit(kxaxis, rates, 2)
     kxopt = -result(1)/(2d0*result(2))
     bta   = btaaxis(grid)  

     array(0, i) = kxopt
     array(1, i) = bta 
  endfor

;plot growth rates
  set_plot, 'ps'
  file = filepath('compare_kxmax.ps',root_dir='.',subdir=[location(0)])
  device, filename=file, $
          bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, array(1,*), array(0,*), xmargin=[6,2],ymargin=[4.,1]  $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=textoidl('\beta'), $
        linestyle = 0, xtickinterval=xtickinterval, ytitle=textoidl('(k_xH)_{opt}'), charthick=6, yrange=yrange, ystyle=1, xstyle=1

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
