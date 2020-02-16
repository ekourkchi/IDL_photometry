PRO tvcircle, x0, y0, maj, nx, ny, linestyle=linestyle, color=color

;----------------------------------
; plot circle at x0,y0 with size
; maj 
;
; scott t. brown - june 1996
;----------------------------------

npoints  = 30.
interval = (2. * !pi) / npoints

t = findgen(npoints + 1) * interval

x = x0 + maj * cos(t) 
y = y0 + maj * sin(t)

c = 255

if keyword_set(color) then c = color

plot,x,y,/noerase,/nodata,xran=[0,nx-1],yran=[0,ny-1], $
	pos=[0,0,1,1],xsty=5,ysty=5

if keyword_set(linestyle) then $
        oplot, x, y, color=c, linestyle=linestyle $
else    oplot, x, y, color=c,thick=3

END
