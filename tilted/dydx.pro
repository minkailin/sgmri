function dydx,x,y,lin_end=lin_end
;+
; ROUTINE:    dydx
;
; PURPOSE:    compute numerical derivative accurate to second order. 
;             In most cases this should be more accurate than the 
;             IDL userlib procedure DERIV.
;
; USEAGE:     dydx(x,y)
;
; INPUT:
;    x        vector of independent variable (monitonic)
;    y        vector of dependent variable
;
; KEYWORD INPUT:
;   lin_end   if set, use a one sided linear approximation for end points
;
; OUTPUT:     numerical derivative dy/dx
;
; EXAMPLE:    
;     !p.multi=[0,1,2]
;     x=[0.,1.,3.,8.,10.,13.,17.,19.,22.,25.]
;     y=x*(25-x) & yp=25-2*x & labels=['analytic','dydx','deriv']
;     plot,x,y,title='quadratic [dydx(x,y) is exact in this case]'
;     plot,x,yp & oplot,x,dydx(x,y),psym=2 & oplot,x,deriv(x,y),psym=4
;     legend,labels,li=[0,0,0],psym=[0,2,4],pos=[0.5,0.7,0.75,0.98]
;
;     x=[0.0,0.5,1.5,2.,3.,3.6,4.0,5.,6.,6.3]
;     y=sin(x) & yp=cos(x) 
;     plot,x,y,title='sin(x)'
;     plot,x,yp & oplot,x,dydx(x,y),psym=2 & oplot,x,deriv(x,y),psym=4
;     legend,labels,li=[0,0,0],psym=[0,2,4],pos=[0.25,0.7,0.50,0.98]
;
;     x=[-2.,-1.,-.7,-.3,-.1,0.,.1,.3,.7,1.,2]
;     y=exp(-x^2) & yp=-2*x*y
;     plot,x,y,title='gausian'
;     plot,x,yp & oplot,x,dydx(x,y),psym=2 & oplot,x,deriv(x,y),psym=4
;     legend,labels,li=[0,0,0],psym=[0,2,4],pos=[0.5,0.7,0.75,0.98]
;             
;
; REVISIONS:
;
;  author:  Paul Ricchiazzi                            jan94
;           Institute for Computational Earth System Science
;           University of California, Santa Barbara
;-
;
n=n_elements(y)

if n eq 0 then begin
  n=n_elements(x)
  y=x
  x=findgen(n)
endif

if n_elements(x) ne n then $
     message,'number of elements of x and y vectors must match'

ym=shift(y,1)
yp=shift(y,-1)
xm=shift(x,1)
xp=shift(x,-1)

denom=(x-xm)*(xp-xm)*(xp-x)
denom(0)=denom(1)
denom(n-1)=denom(n-2)

dfdx=yp*(x-xm)^2+y*(xp-xm)*(xp+xm-2*x)-ym*(xp-x)^2

if keyword_set(lin_end) then begin
  
  dfdx=dfdx/denom
  
  dfdx(0)=(y(1)-y(0))/(x(1)-x(0))

  dfdx(n-1)=(y(n-1)-y(n-2))/(x(n-1)-x(n-2))


endif else begin

  dfdx(0)=-y(2)*(x(1)-x(0))^2+y(1)*(x(2)-x(0))^2 $
          -y(0)*(x(2)-x(1))*(x(1)+x(2)-2*x(0))

  dfdx(n-1)= y(n-1)*(x(n-2)-x(n-3))*(2*x(n-1)-x(n-2)-x(n-3)) $
            -y(n-2)*(x(n-1)-x(n-3))^2 $
            +y(n-3)*(x(n-1)-x(n-2))^2

  dfdx=dfdx/denom

endelse


return,dfdx

end
