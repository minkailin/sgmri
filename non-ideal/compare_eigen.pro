pro compare_eigen, loc=loc, nz=nz, xrange=xrange, legend=legend, label=label, title=title, yrange=yrange, var=var, solid=solid, psym=psym
                             
  location = strcompress(loc,/remove_all)  
  ncases   = n_elements(loc)
 
  i = 0 
 
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

        'rm': begin
           file = filepath('params.dat',root_dir='.',subdir=[location(i)])
           params = dblarr(5,nlines)
           
           openr,1,file
           readf,1, params, format='(5(e22.15,x))'
           close,1
           
           xaxis = params(4,*)
           xtitle=textoidl('\Lambda_0')
        end
     endcase
     
     data  = -eigen(4,*)*100d0 ;growth rate
     
     if keyword_set(psym) then begin
        xaxis = congrid(xaxis(*), 32)
        data  = congrid(data(*), 32)
     endif

;stop
;plot growth rates
  set_plot, 'ps'
  file = filepath('compare_growth.ps',root_dir='.',subdir=[location(0)])
  device, filename=file, $
          bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, xaxis, data,xmargin=[6,2],ymargin=[4.,1]  $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, psym=psym, $
        linestyle = 0, xtickinterval=xtickinterval, ytitle=textoidl('10^2\gamma/\Omega'), charthick=6, yrange=yrange, ystyle=1, xstyle=1
  for i=1, ncases-1 do begin
     file = filepath('eigenvalues.dat',root_dir='.',subdir=[location(i)])
     nlines = file_lines(file)
     eigen = dblarr(6,nlines)

     openr,1,file
     readf,1, eigen, format='(6(e22.15,x))'
     close,1
     
     case var of
        
        'toomre':  begin
           xaxis = eigen(0,*) 
        end
        'kx':      begin
           xaxis = eigen(2,*)
        end
        'beta':   begin
           xaxis = eigen(1,*)
        end
        'rm': begin
           file = filepath('params.dat',root_dir='.',subdir=[location(i)])
           params = dblarr(5,nlines)
           
           openr,1,file
           readf,1, params, format='(5(e22.15,x))'
           close,1
           
           xaxis = params(4,*)
        end
        
     endcase
     
     data  = -eigen(4,*)*100d0 ;growth rate

     if keyword_set(psym) then begin
        xaxis = congrid(xaxis(*), 32)
        data  = congrid(data(*), 32)
     endif


     if keyword_set(solid) then begin
        temp = 0
     endif else temp = i

     oplot, xaxis , data , thick=4, linestyle=temp, psym=psym
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
