PRO ellipse, x0, y0, maj, min, t0, linestyle=linestyle, color=color, $
	init=init, device=device

;----------------------------------
; plot ellipse at x0,y0 with axes
; maj and min, inclination angle t0
; t0 in degrees east of north or
; degrees left of up
;
; scott t. brown - june 1996
;----------------------------------

tz = -(t0 + 90.0)/!radeg

npoints  = 30.
interval = (2. * !pi) / npoints

t = findgen(npoints + 1) * interval

x = x0 + ((maj * cos(t) * cos(tz)) - (min * sin(t) * sin(tz)))
y = y0 + ((maj * cos(t) * sin(tz)) + (min * sin(t) * cos(tz)))

c = 255

if keyword_set(color) then c = color

if keyword_set(init) then begin
	if keyword_set(linestyle) then $
		plot, x, y, color=c, linestyle=linestyle,device=device $
	else	plot, x, y, color=c,device=device

endif else begin
	if keyword_set(linestyle) then $
		oplot, x, y, color=c, linestyle=linestyle $
	else	oplot, x, y, color=c
endelse

END
