pro plot_fq, qrange=qrange

nx = 128
xaxis = qrange(0) + (qrange(1)-qrange(0))*dindgen(nx)/(nx-1d0)
yaxis = sqrt(xaxis)*acos(xaxis/(1d0+xaxis))
yaxis = 1d0/yaxis

set_plot, 'ps'
file = 'plot_fq.ps'
device, filename=file, $
        bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
plot, xaxis, yaxis,xmargin=[6,2],ymargin=[3.,1.75], ystyle=1  $
      ,charsize=1.5, thick=4, yrange=[0.5,4], xtitle=textoidl('Q'), $
      linestyle = 0, xtickinterval=xtickinterval, ytitle=textoidl('polytropic f (solid)'), charthick=6
oplot, [0,100], [1,1]*1d0/sqrt(2d0), linestyle=2, thick=4

device,/close
end
