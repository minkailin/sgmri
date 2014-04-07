pro compare_eigen, loc=loc, nz=nz, xrange=xrange, legend=legend, label=label, title=title, yrange=yrange, var=var
  
  location = strcompress(loc,/remove_all)  
  ncases   = n_elements(loc)
 
  i = 0 
 
     file = filepath('eigenvalues.dat',root_dir='.',subdir=[location(i)])
     nlines = file_lines(file)

     eigen = dblarr(7,nlines)
 
     openr,1,file
     readf,1, eigen, format='(7(e22.15,x))'
     close,1

     case var of

'toomre':  begin
           xaxis = eigen(0,*) ; toomre q
           xtitle = 'Q'
           end
'kx':      begin
           xaxis = eigen(2,*); kx
           xtitle=textoidl('k_xH')
        end
'beta':   begin
   xaxis = eigen(1,*)           ; kx
   xtitle=textoidl('\beta')
end
'tilt':   begin
   xaxis = eigen(2,*)           ; kx
   xtitle=textoidl('B_y/B_z')
end

      endcase

     data  = -eigen(5,*)*100d0 ;growth rate
     

;plot growth rates
  set_plot, 'ps'
  file = filepath('compare_growth.ps',root_dir='.',subdir=[location(0)])
  device, filename=file, $
          bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, xaxis, data,xmargin=[6,2],ymargin=[4.,1]  $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, $
        linestyle = 0, xtickinterval=xtickinterval, ytitle=textoidl('10^2\gamma/\Omega'), charthick=6, yrange=yrange, ystyle=1, xstyle=1
  for i=1, ncases-1 do begin
     file = filepath('eigenvalues.dat',root_dir='.',subdir=[location(i)])
     nlines = file_lines(file)
     eigen = dblarr(7,nlines)

     openr,1,file
     readf,1, eigen, format='(7(e22.15,x))'
     close,1
          case var of

'toomre':  begin
           xaxis = eigen(0,*) ; toomre q
           xtitle = 'Q'
           end
'kx':      begin
           xaxis = eigen(2,*); kx
           xtitle=textoidl('k_xH')
        end
'beta':   begin
   xaxis = eigen(1,*)           ; kx
   xtitle=textoidl('\beta')
end
'tilt':   begin
   xaxis = eigen(2,*)           ; kx
   xtitle=textoidl('B_y/B_z')
end
endcase


     data  = -eigen(5,*)*100d0 ;growth rate

     oplot, xaxis , data , thick=4, linestyle=i
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
