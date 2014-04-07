pro test, var=var

case var of 

1: plot, abs(d2pot - kx*kx*pot - rho*w/csq/fq^2/bigQ)/abs(kx*kx*pot(0)), charsize=2

2: plot, abs(ii*sig*w/fq/csq + ii*kx*vx + drho*vz + dvz)/abs(ii*kx*vx(0)), charsize=2

3: plot, abs(ii*sig*vx - 2d0*vy - fq*( vasq*dbx - ii*kx*(w + pot + vasq*bz) ))/abs(vy(0)),charsize=2

4: plot, abs(ii*sig*vy + 0.5*kappa2*vx - fq*vasq*dby)/abs(0.5*kappa2*vx(0)),charsize=2

5: plot, abs(ii*sig*vz + fq*(dw+dpot))/abs(vz(nz/2)),charsize=2

6: plot, abs(ii*sbar*bx - fq*dvx - eta*d2bx - deta*dbx + ii*kx*deta*bz)/abs(dvx(nz/2)), charsize=2

7: plot, abs(ii*sbar*by - fq*dvy +shear*bx - eta*d2by - deta*dby)/abs(dvy(nz/2)),charsize=2

8: plot, abs(ii*sbar*bz + ii*fq*kx*vx - eta*d2bz)/abs(ii*fq*kx*vx(0)),charsize=2

9: plot, abs(ii*sbar*bz + ii*fq*kx*vx + eta*ii*kx*dbx)/abs(ii*fq*kx*vx(0)),charsize=2

10: plot, abs(ii*kx*bx+dbz)/abs(ii*kx*bx(nz/2)),charsize=2

endcase


end
