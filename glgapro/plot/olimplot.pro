pro olimplot,x,y,yerr,thick=thick,lower=lower,upper=upper,color=color, $
	scale=scale,xlimit=xlimit
;
;
xx = [-1,1,0, 0,-1, 0, 1]
yy = [ 0,0,0,-5,-3,-5,-3]
if keyword_set(upper) then begin
	yy = -yy
endif
if keyword_set(scale) then begin
	yy = yy * scale
	xx = xx * scale
endif
if keyword_set(xlimit) then begin
	tmp = xx
	xx = yy
	yy = tmp
endif
;
usersym,xx,yy,thick=thick
;
if keyword_set(color) then begin
	for i=0,n_elements(x)-1 do $
		oplot,[x(i),x(i)],[y(i),y(i)],thick=thick,psym=8,color=color
	if n_params(0) ge 3 then errplot,x,y-yerr,y+yerr,color=color
endif else begin
	for i=0,n_elements(x)-1 do $
		oplot,[x(i),x(i)],[y(i),y(i)],thick=thick,psym=8
	if n_params(0) ge 3 then errplot,x,y-yerr,y+yerr
endelse
;
return
end
