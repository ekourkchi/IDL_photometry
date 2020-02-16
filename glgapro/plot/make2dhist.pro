pro make2dhist, x, y, xlimit, dx, ylimit, dy, f, xreverse=xreverse, $
                yreverse=yreverse
;;;Compute a 2-D histogram.

nx=round(abs(xlimit[1]-xlimit[0])/dx)
ny=round(abs(ylimit[1]-ylimit[0])/dy)
if keyword_set(xreverse) then xpts=xlimit[0]-dx*(findgen(nx)-0.5) $
else xpts=xlimit[0]+dx*(findgen(nx)+0.5)
if keyword_set(yreverse) then ypts=ylimit[0]-dy*(findgen(ny)-0.5) $
else ypts=ylimit[0]+dy*(findgen(ny)+0.5)

f=create_struct('f',intarr(nx,ny),'x',xpts,'y',ypts,'nx',nx,'ny',ny,$
                'xlimit',xlimit,'ylimit',ylimit,'dx',dx,'dy',dy)
for i=0,nx-1 do for j=0,ny-1 do begin
    x1=xpts[i]-0.5*dx
    x2=xpts[i]+0.5*dx
    y1=ypts[j]-0.5*dy
    y2=ypts[j]+0.5*dy
    w=where((x ge x1) and (x lt x2) and (y ge y1) and (y lt y2),n)
    f.f[i,j]=n
endfor

return
end

