pro compare_eigen2, loc=loc, panels=panels, xrange=xrange, legend=legend, label=label, title=title, yrange=yrange, var=var, $
                    custom_margin=custom_margin
  
  npanels  = n_elements(panels)
  ncases   = n_elements(loc)
  location = strcompress(loc,/remove_all)  
  
  i = 0 
  j = 0
  
  file = filepath('eigenvalues.dat',root_dir='.',subdir=[location(i+j)])
  nlines = file_lines(file)
  
  eigen = dblarr(6,nlines)
  
  openr,1,file
  readf,1, eigen, format='(6(e22.15,x))'
  close,1
  
  case var of
     'toomre':  begin
        grid = 0                ; toomre q
        xtitle1= 'Q'
     end
     'kx':      begin
        grid = 2                ; kx
        xtitle1=textoidl('k_xH')
     end
     'beta':   begin
        grid = 1                ; beta
        xtitle1=textoidl('\beta')
     end
  endcase
  xaxis = eigen(grid,*)
  data  = -eigen(4,*)*100d0     ;growth rate
  
  set_plot, 'ps'
  file = filepath('compare_growth2.ps',root_dir='.',subdir=[location(0)])
  device, filename=file, $
          bits_per_pixel=8,xsize=8, ysize=npanels*2,xoffset=0,yoffset=0,/inches
  multiplot, [1,npanels],custom_margin=custom_margin
  
  count=0
  for i=0, npanels-1 do begin     
     xtitle = ''
     ytitle = ''
     if(i eq npanels-1) then xtitle=xtitle1
     if(i eq npanels/2) then ytitle=textoidl('10^2\gamma/\Omega') 
     
     plot, xaxis, data, xmargin=[6,2],ymargin=[4.,1]  $
           ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, $
           linestyle = 0, xtickinterval=xtickinterval, ytitle=ytitle, charthick=6, yrange=yrange, xstyle=1,ystyle=1, /nodata  
     
     if keyword_set(legend) then begin
        x0=legend(0)
        y0=legend(1)
        ;; y0=legend(2)
        ;; dy=legend(3)
        ;; for j=0, n_elements(label)-1 do begin
        ;;    oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
        ;;    xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5,charthick=6
        ;; endfor
        if keyword_set(title) then xyouts, x0, y0, textoidl(title(i)),charsize=1.5,charthick=6
     endif
     
     for j=0, panels(i)-1 do begin
        file = filepath('eigenvalues.dat',root_dir='.',subdir=[location(count)])
        nlines = file_lines(file)
        
        eigen = dblarr(6,nlines)
        
        openr,1,file
        readf,1, eigen, format='(6(e22.15,x))'
        close,1
        
        xaxis = eigen(grid,*)
        data  = -eigen(4,*)*100d0 ;growth rate
        
        oplot, xaxis, data, thick=4
        count += 1
     endfor
     multiplot
  endfor
  
  ;; if keyword_set(legend) then begin
  ;;    x0=legend(0)
  ;;    x1=legend(1)
  ;;    y0=legend(2)
  ;;    dy=legend(3)
  ;;    for j=0, n_elements(label)-1 do begin
  ;;       oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
  ;;       xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5,charthick=6
  ;;    endfor
  ;;    if keyword_set(title) then xyouts, x1, y0+dy, textoidl(title),charsize=1.5,charthick=6
  ;; endif
  device,/close
  multiplot,/reset
end
