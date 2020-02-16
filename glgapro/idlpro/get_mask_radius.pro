function get_mask_radius,im,xcen,ycen,dtype
;+
; get_mask_radius - get the default radius for masking point sources
;
; return value in pixels
;-
radius = 6.	; default
if strpos(dtype,'galex') ge 0 then $
	psffile = !GLGA_ROOT+'psf/CDFS_00-nd-prof.dat' $
else if strpos(dtype,'sdss') ge 0 then $
	psffile = !GLGA_ROOT+'psf/sdss_prof.dat' $
else if strpos(dtype,'panstarrs') ge 0 then $
	psffile = !GLGA_ROOT+'psf/sdss_prof.dat' $	
else if strpos(dtype,'2mass') ge 0 then $
	psffile = !GLGA_ROOT+'psf/2mass_prof.dat' $
else if strpos(dtype,'wise') ge 0 then $
	psffile = !GLGA_ROOT+'psf/wise_prof.dat' $
else if strpos(dtype,'irac') ge 0 then $
	psffile = !GLGA_ROOT+'psf/irac_prof.dat' $
else if strpos(dtype,'csp') ge 0 then $
	psffile = !GLGA_ROOT+'psf/csp_prof.dat' $
else begin
	print,'No profile for data type, returning default'
	return,radius
endelse

readcol,psffile,r,psf,format='f,f',/silent

if strpos(dtype,'galex') ge 0 then begin
	meanclip, im, nsky1, nskysigma, clipsig = 3
	a = where(im le nsky1 + 3.*nskysigma)
	nsky2 = median(im[a])

	if nsky1 gt 0 and nsky2 gt 0 then nsky = nsky1 < nsky2
	if nsky1 eq 0 then nsky = nsky2 
	if nsky2 eq 0 then nsky = nsky1
	if nsky1 eq 0 and nsky2 eq 0 then nsky=0.002

	val = im[xcen>0, ycen>0]
	radius = ceil(r[value_to_index(psf*val, nsky)])+2.
	;what peak value should one do max mask?????
	;yes 3.8 (UGC 24), no 2.0 (UGC 34), no 3.3 (UGC 15)
	if val ge 3.5 then radius = 30.
endif else if strpos(dtype,'sdss') ge 0 then begin
	MMM,im,nsky,nsky_sig,/silent
	aper,im,xcen,ycen,mag,merr,sky,skyerr,1.,2.5,[5.,10.],[0.,0.],/silent
	if not finite(sky) or sky lt 0. then sky = 0.
	if nsky_sig lt 0. then nsky_sig = skyerr
	val = ((im[xcen>0, ycen>0]-sky[0])>nsky_sig)/10.0
	radius = ceil(r[value_to_index(psf*val, nsky_sig)])+2.
endif else if strpos(dtype,'panstarrs') ge 0 then begin
	MMM,im,nsky,nsky_sig,/silent
	aper,im,xcen,ycen,mag,merr,sky,skyerr,1.,2.5,[5.,10.],[0.,0.],/silent
	if not finite(sky) or sky lt 0. then sky = 0.
	if nsky_sig lt 0. then nsky_sig = skyerr
	val = ((im[xcen>0, ycen>0]-sky[0])>nsky_sig)/10.0
	radius = ceil(r[value_to_index(psf*val, nsky_sig)])+2.
	endif else if strpos(dtype,'2mass') ge 0 then begin
	MMM,im,nsky,nsky_sig,/silent
	aper,im,xcen,ycen,mag,merr,sky,skyerr,1.,2.5,[5.,10.],[0.,0.],/silent
	if not finite(sky) or sky lt 0. then sky = 0.
	if nsky_sig lt 0. then nsky_sig = skyerr
	val = ((im[xcen>0, ycen>0]-sky[0])>nsky_sig)/10.0
	radius = ceil(r[value_to_index(psf*val, nsky_sig)])+2.
endif else if strpos(dtype,'wise') ge 0 then begin
	MMM,im,nsky,nsky_sig,/silent
	aper,im,xcen,ycen,mag,merr,sky,skyerr,1.,2.5,[5.,10.],[0.,0.],/silent
	if not finite(sky) or sky lt 0. then sky = 0.
	if nsky_sig lt 0. then nsky_sig = skyerr
	val = ((im[xcen>0, ycen>0]-sky[0])>nsky_sig)/10.0
	radius = ceil(r[value_to_index(psf*val, nsky_sig)])+2.
endif else if strpos(dtype,'irac') ge 0 then begin
	MMM,im,nsky,nsky_sig,/silent
	aper,im,xcen,ycen,mag,merr,sky,skyerr,1.,2.5,[5.,10.],[0.,0.],/silent
	if not finite(sky) or sky lt 0. then sky = 0.
	if nsky_sig lt 0. then nsky_sig = skyerr
	val = ((im[xcen>0, ycen>0]-sky[0])>nsky_sig);/10.0
	radius = ceil(r[value_to_index(psf*val, nsky_sig)])+2.
endif else if strpos(dtype,'csp') ge 0 then begin
	MMM,im,nsky,nsky_sig,/silent
	aper,im,xcen,ycen,mag,merr,sky,skyerr,1.,2.5,[5.,10.],[0.,0.],/silent
	if not finite(sky) or sky lt 0. then sky = 0.
	if nsky_sig lt 0. then nsky_sig = skyerr
	val = ((im[xcen, ycen]-sky[0])>nsky_sig)/10.0
	radius = ceil(r[value_to_index(psf*val, nsky_sig)])+2.
endif else begin
	print,'No profile for data type: ',dtype
	print,'Returning default: 6.'
endelse

return,radius

end
