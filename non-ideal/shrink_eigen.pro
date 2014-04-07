pro shrink_eigen, loc=loc, nk=nk, var=var

 location = strcompress(loc,/remove_all)
 ncases   = n_elements(loc)

 for  i = 0, ncases-1 do begin 

 file = filepath('eigenvalues.dat',root_dir='.',subdir=[location(i)])
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
           xmin  = xaxis(0)
           xmax  = xaxis(nlines-1)
           dx    = (xmax - xmin)/(nk-1) 
           xaxis1 = xmin + dx*dindgen(nk)            
           grid = 0 
        end
        'kx':      begin
           xaxis = eigen(2,*)
           xmin  = xaxis(0)
           xmax  = xaxis(nlines-1)
           dx    = (xmax - xmin)/(nk-1)
           xaxis1 = xmin + dx*dindgen(nk)
           grid = 2 
        end
        'beta':   begin
           xaxis = eigen(1,*)
           xmin  = alog10(xaxis(0))
           xmax  = alog10(xaxis(nlines-1))
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
           xmin  = alog10(xaxis(0))
           xmax  = alog10(xaxis(nlines-1))
           dx    = (xmax - xmin)/(nk-1)
           xaxis1 =10^( xmin + dx*dindgen(nk) )
        end
        
 endcase
          sig_im = interpol(eigen(4,*), xaxis, xaxis1,  /spline)
          sig_re = interpol(eigen(3,*), xaxis, xaxis1,  /spline)

          eigen1(3,*) = sig_re
          eigen1(4,*) = sig_im 
          eigen1(grid,*) = xaxis1     

          

       file = filepath('eigenvalues.dat',root_dir='.',subdir=[location(i)])
       openw,1,file
       for j=0, nk-1 do begin
       printf,1, eigen1(*,j), format='(6(e22.15,x))'
       endfor
       close,1
endfor

end 
