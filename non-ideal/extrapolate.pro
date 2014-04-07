pro extrapolate, loc=loc, nk=nk, var=var, polyn=polyn, xrange=xrange, xsample=xsample 

 location = strcompress(loc,/remove_all)

 file = filepath('eigenvalues.dat',root_dir='.',subdir=[location])
 nlines = file_lines(file)

 eigen = dblarr(6,nlines)
 
 openr,1,file
 readf,1, eigen, format='(6(e22.15,x))'
 close,1


 eigen1 = dblarr(6, nk)
 xaxis1 = dblarr(nk)
        
  case var of

        'toomre':  begin
           xaxis = eigen(0,*)
           
           if not keyword_set(xrange) then begin 
           xmin  = xaxis(0)
           xmax  = xaxis(nlines-1)
           endif else begin
           xmin = xrange(0)
           xmax = xrange(1) 
           endelse         


           dx    = (xmax - xmin)/(nk-1)
           xaxis1 = xmin + dx*dindgen(nk)
           grid = 0
        end
        'kx':      begin
           xaxis = eigen(2,*)

 	   if not keyword_set(xrange) then begin
           xmin  = xaxis(0)
           xmax  = xaxis(nlines-1)
           endif else begin
           xmin = xrange(0)
           xmax = xrange(1)
           endelse

           dx    = (xmax - xmin)/(nk-1)
           xaxis1 = xmin + dx*dindgen(nk)
           grid = 2
        end
        'beta':   begin
           xaxis = eigen(1,*)
           
           if not keyword_set(xrange) then begin
           xmin  = alog10(xaxis(0))
           xmax  = alog10(xaxis(nlines-1))
           endif else begin
           xmin = xrange(0)
           xmax = xrange(1)
           endelse

           dx    = (xmax - xmin)/(nk-1)
           xaxis1 =10^( xmin + dx*dindgen(nk) )
           grid = 1
        end
        'rm': begin
           file = filepath('params.dat',root_dir='.',subdir=[location(i)])
           params = dblarr(5,nlines)

           openr,1,file
           readf,1, params, format='(5(e22.15,x))'
           close,1

           xaxis = params(4,*)


           if not keyword_set(xrange) then begin
           xmin  = alog10(xaxis(0))
           xmax  = alog10(xaxis(nlines-1))
           endif else begin
           xmin = xrange(0)
           xmax = xrange(1)
           endelse
           
           dx    = (xmax - xmin)/(nk-1)
           xaxis1 =10^( xmin + dx*dindgen(nk) )
        end

 endcase

 eigen1(grid,*) = xaxis1

 temp = min(abs(xaxis - xsample(0)),x1)
 temp = min(abs(xaxis - xsample(1)),x2)

 result = poly_fit(xaxis(x1:x2), eigen(4,x1:x2), polyn, /double) 
 for i = 0, nk-1 do begin
 eigen1(4,i) = 0d0 
 for j = 0, polyn do eigen1(4,i) += result(j)*xaxis1(i)^j
 endfor 
  
 result = poly_fit(xaxis(x1:x2), eigen(3,x1:x2), polyn, /double)
 for i = 0, nk-1 do begin
 eigen1(3,i) = 0d0
 for j = 0, polyn do eigen1(3,i) += result(j)*xaxis1(i)^j
 endfor

  

       file = filepath('eigenvalues.dat',root_dir='.',subdir=[location])
       openw,1,file
       for j=0, nk-1 do begin
       printf,1, eigen1(*,j), format='(6(e22.15,x))'
       endfor
       close,1

end 
