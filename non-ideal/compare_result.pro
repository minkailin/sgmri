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

pro compare_result, loc=loc, kmode=kmode, nz=nz, xrange=xrange, legend=legend, label=label, en=en, cowling=cowling 
  
  location =strcompress(loc,/remove_all)

  ii = dcomplex(0d0, 1d0)  
  kappa2 = 1.0
  
  ncases = n_elements(kmode)
  en1 = dblarr(ncases, nz)
  en2 = en1
  en3 = en1
  en4 = en1
  elsasser = en1 
  data = en1
  if not keyword_set(cowling) then cowling = intarr(ncases) 

  vort2 = en1

  vx = dcomplexarr(nz) 
  vy = dcomplexarr(nz)
  vz = dcomplexarr(nz)
  den= dcomplexarr(nz)          ;this is actually w=dp/rho
  pot= dcomplexarr(nz)
  
  dvx = vx
  dvy = vy
  dden= den
  dpot= pot 
  bx  = dcomplexarr(nz)
  by  = dcomplexarr(nz)
  bz  = dcomplexarr(nz) 
  



  for k=0, ncases-1 do begin
     
     file = filepath('params.dat',root_dir='.',subdir=[location(k)])
     nlines = file_lines(file)
     param = dblarr(5,nlines)
     openr,1,file
     readf,1, param, format='(5(e22.15,x))'
     close,1
     zmax   = param(0,kmode(k)-1) 
     lmodes = param(1,kmode(k)-1)
     fq     = param(2,kmode(k)-1)
     shear  = param(3,kmode(k)-1) 
     rm     = param(4,kmode(k)-1) ;this is actually midplane elsasser
     
     zaxis = zmax*dindgen(nz)/(nz-1d0)
     
;focus on a specific mode (among diff kx, bigQ or beta) 
     
     file = filepath('var_axis.dat',root_dir='.',subdir=[location(k)])
     nlines = file_lines(file)
     array = dblarr(2, nlines)
     openr,1,file
     readf,1,array,format='(2(e22.15,x))'
     close,1
     ;print, 'cases, var, growth', kmode(k), array(1,kmode(k)-1) 
     
;basic state for this mode 
     file = filepath('basic.dat',root_dir='.',subdir=[location(k)])
     nlines = file_lines(file)
     array = dblarr(7, nlines)
     openr,1,file
     readf,1,array,format='(7(e22.15,x))'
     close,1
;write(99, fmt='(4(e22.15,x))') zaxis(j), dnorm(j), valf(j), csq(j)
     z_small = array(0,lmodes*(kmode(k)-1): lmodes*kmode(k)- 1)
     rho     = array(1,lmodes*(kmode(k)-1): lmodes*kmode(k)- 1)
     vasq    = array(2,lmodes*(kmode(k)-1): lmodes*kmode(k)- 1)
     csq     = array(3,lmodes*(kmode(k)-1): lmodes*kmode(k)- 1)
     eta   = array(4,lmodes*(kmode(k)-1): lmodes*kmode(k)- 1)
     drho    = array(5,lmodes*(kmode(k)-1): lmodes*kmode(k)- 1)
     deta    = array(6,lmodes*(kmode(k)-1): lmodes*kmode(k)- 1)
     
;interpolate to fine grid
     rho  = interpol(rho, z_small, zaxis,  /spline)
     drho = interpol(drho, z_small, zaxis,  /spline)
     vasq = interpol(vasq, z_small, zaxis, /spline)
     csq  = interpol(csq, z_small, zaxis, /spline) 
     eta  = interpol(eta, z_small, zaxis, /spline) 
     deta = interpol(deta, z_small, zaxis, /spline)
     
     
     elsasser(k,*) = alog10(fq*fq*vasq/eta)

;eigenvalue information 
     file = filepath('eigenvalues.dat',root_dir='.',subdir=[location(k)])
     nlines = file_lines(file)
     eigen = dblarr(6, nlines)
     openr,1,file
     readf,1, eigen, format='(6(e22.15,x))'
     close,1
;write(20, fmt='(6(e22.15,x))') bigQ, beta, kx, dble(sig), dimag(sig), condition
     bigQ = eigen(0, kmode(k)-1)
     bta  = eigen(1, kmode(k)-1)
     kx   = eigen(2, kmode(k)-1)
     sig  = dcomplex( eigen(3, kmode(k)-1), eigen(4, kmode(k)-1) )
     sbar = sig -ii*eta*kx^2
     print, 'growth, kx, kmode(k)', -imaginary(sig), kx, kmode(k)
;eigenfunction information 

if(cowling(k) eq 0) then begin
    array = dblarr(2,4*lmodes)
endif else begin
    array = dblarr(2,3*lmodes)
endelse

     filename =  strcompress('eigenfunctions_'+string(kmode(k), format='(I2)')+'.dat', /remove_all)
     file = filepath(filename,root_dir='.',subdir=[location(k)])
     openr, 1, file 
     readf,1,array,format='(2(e22.15,x))'
     close, 1
     
     for j=0, nz-1 do begin
        zbar = zaxis(j)/zmax 
        vx(j)  = 0d0 
        vy(j)  = 0d0 
        den(j) = 0d0 
        pot(j) = 0d0 
        
        dvx(j) = 0d0 
        dvy(j) = 0d0 
        dden(j)= 0d0
        dpot(j)= 0d0 
        for i=1, lmodes do begin
           arr = cheby(2d0*(i-1d0), zbar)
           
           vx(j) +=dcomplex( array(0, i-1),  $
                             array(1, i-1)) $
                *arr(0)         ;this is T_l
           dvx(j) +=dcomplex( array(0, i-1),  $
                              array(1, i-1)) $
                    *arr(1)
           
           vy(j) +=dcomplex( array(0, lmodes + i-1),  $
                             array(1, lmodes + i-1 )) $
                   *arr(0)         
           
           dvy(j) +=dcomplex( array(0, lmodes + i-1),  $
                              array(1, lmodes + i-1 )) $
                    *arr(1)
           
           den(j) +=dcomplex( array(0, 2*lmodes + i-1  ),$
                              array(1, 2*lmodes + i-1) ) $
                    *arr(0)
           dden(j) +=dcomplex( array(0, 2*lmodes + i-1  ),$
                               array(1, 2*lmodes + i-1) ) $
                     *arr(1)

           if( cowling(k) eq 0) then begin
           pot(j) +=dcomplex( array(0, 3*lmodes + i-1),  $
                              array(1, 3*lmodes + i-1)) $
                    *arr(0)
           
           dpot(j) +=dcomplex( array(0, 3*lmodes + i-1),  $
                               array(1, 3*lmodes + i-1)) $
                     *arr(1)
           endif 

        endfor
     endfor 
     
     dvx  /= zmax
     dvy  /= zmax
     dden /= zmax 
     dpot /= zmax 
     
     vz = ii*fq*(dden + dpot)/sig 
  
     g  = ii*sig*vx - 2d0*vy + ii*fQ*kx*(den+pot)
     dg = ii*sig*dvx- 2d0*dvy+ ii*fQ*kx*(dden+dpot)
     
     bx = fq*dvx + (bta*rho/fQ)*(eta*(drho*g + dg) + deta*g)
     bx *= -ii/sig
     
     
     g  = ii*sig*vy + 0.5*kappa2*vx
     dg = ii*sig*dvy+ 0.5*kappa2*dvx
     
     by = fQ*dvy - shear*bx + (bta*rho/fQ)*(eta*(drho*g + dg) + deta*g)
     by /= ii*sbar
     
     bz = -fQ*kx*vx - (kx*eta*bta*rho/fQ)*( ii*sig*vx - 2d0*vy + ii*fQ*kx*(den+pot) )
     bz /= sig 
     
     v2  = abs(vx)^2 + abs(vy)^2 + abs(vz)^2
     b2  = abs(bx)^2 + abs(by)^2 + abs(bz)^2 
     gpe = real_part(den*conj(pot))
     
     en1(k,*) = b2/(2d0*bta)
     en2(k,*) = rho*abs(gpe)/csq/2d0 
     en3(k,*) = rho*v2/2d0 
     en4(k,*) = rho*abs(den)^2/(2d0*csq)
     
     m1 = max(en1(k,*))
     
     en1(k,*) /=  max(en1(k,*))
if(cowling(k) eq 0) then begin
     en2(k,*) /=  max(en2(k,*))
endif else en2(k,*) = 0d0  
     en3(k,*) /=  max(en3(k,*))
     en4(k,*) /=  max(en4(k,*))
        

     data(k,*) = b2/max(b2)


;vorticity
; vort2(k,*) = abs(dvy)^2 +abs(ii*kx*vz - dvx)^2 + kx*kx*abs(vy)^2

     domega_z = ii*kx*fQ*vy
     dvorten = domega_z
     vort2(k,*)= abs(dvorten)

;     vort2(k,*) = abs(vy)

     vort2(k,*)/= max(vort2(k,*))
  endfor

  case en of
     1 : begin
;        data = en1
        ytitle = textoidl('E_m/max(E_m)')
     end
     2 : begin
        data = en2
        ytitle = textoidl('E_g/max(E_g)')
     end
     3 : begin
        data = en3
        ytitle = textoidl('E_k/max(E_k)')
     end
     4 : begin
        data = en4
        ytitle = textoidl('E_t/max(E_t)')
     end
     'vort' :begin
      data = vort2 
      ytitle = textoidl('|\omega|^2/max(|\omega|^2)')
     end
  endcase

  set_plot, 'ps'

  kx_string   = strcompress(textoidl('k_xH=')+ string(kx, format='(f4.2)'), /remove_all)
  beta_string = strcompress(textoidl('\beta=')+ string(bta, format='(f5.1)'), /remove_all)
  rm_string = strcompress(textoidl('log\Lambda_0=')+ string(alog10(rm), format='(f4.1)'), /remove_all)
  
  title = beta_string  +', '+ rm_string+ ', '   +  kx_string
  
  file = filepath('compare_result.ps',root_dir='.',subdir=[location(0)])
  device, filename=file, $
          bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  
  plot, zaxis, data(0,*), xmargin=[6.5,1.5],ymargin=[3.,1.75], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=[0,1], xtitle=textoidl('z/H'), $
        linestyle = 0, xtickinterval=xtickinterval, title=title, charthick=6, ytitle=ytitle

;  temp = min(abs(elsasser(0,*)),z1)
;  zcrit = zaxis(z1)
;  oplot, [zcrit, zcrit], [-1,1]*1d10, thick=1, linestyle=0

  for k=1, ncases-1 do begin
     oplot, zaxis, data(k,*), thick=4, linestyle=k
;     temp = min(abs(elsasser(k,*)),z1)
;     zcrit = zaxis(z1)
;     oplot, [zcrit, zcrit], [-1,1]*1d10, thick=1, linestyle=k
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
end
