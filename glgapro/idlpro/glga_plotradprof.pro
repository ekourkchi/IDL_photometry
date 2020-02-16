pro glga_plotradprof, id, band, survey=survey, type=type, $
    intfile=intfile, maskimgfile=maskimgfile, bandmaskimgfile=bandmaskimgfile, $
    nonlinearity=nonlinearity, brite=brite, origin=origin, $
    quality=quality, outpath=outpath, verbose=verbose, $
    pathtoprofile=pathtoprofile, yuan13=yuan13, update=update, rotatee=rotatee
;
; note: set yuan13 keyword to use Yuan et al. 2013 empirical extintion coeff's
;	Yuan coeff's assume SFD98 over-estimates E(B-V) by 14% (see table 2)
;
; Warning this procedure spawns convert and ps2pdf14
; (without any checking) to convert postscript output 
; to pdf and jpeg files. Can just comment out if postscript
; is OK. 
;

if not keyword_set(type) then typ = '-' else typ = strtrim(type,2)
if not keyword_set(survey) then srv = '' else srv = strtrim(survey,2)
if not keyword_set(outpath) then outpath='./'
if not keyword_set(pathtoprofile) then pathtoprofile='./'
if not keyword_set(quality) then quality = 75
if not keyword_set(brite) then brite = 1.
if not keyword_set(nonlinearity) then nonlinearity = 2.5
if not keyword_set(origin) then origin = [0,0,0]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in profile data

read_radprof,id,band, tot_a, tot_mag, tot_mag_e, ann_mu, ann_mu_e, $
	ra_cen, dec_cen, semimajor, semiminor, pa, as_pix, $
	tf_mag, tf_mag_e, af_mag, af_mag_e, mu_bg, mu_bg_e, $
	skyradius_in, skyradius_out, tf_a, asymag, asymag_e, asyma, $
	tot_int=tot_int, nbg_flx=bg, nbg_e=bg_e, phot_ts=phot_ts, $
	annuli_size=annuli_size, pathtoprofile=pathtoprofile, $
	int_to_mjy=mjc
; check inputs
if tot_a[0] eq -1 then begin
	if keyword_set(verbose) then print,'Error reading profile data for: ',$
		id,' ',band
	return
endif
;
; check plot grade
gflist = file_search(outpath+id+'_'+band+'_grade*',count=ngf)
if ngf gt 0 then begin
	grade = fix(strmid(gflist[0],strpos(gflist[0],'grade')+5,1))
	gfl_fi = file_info(gflist[0])
	grade_ts = gfl_fi.ctime
endif else begin
	grade = 1
	grade_ts = 0.
endelse
;
; check for output existence
ofl_fi = file_info(outpath+id+'_'+band+'.jpg')
;
; are we updating?
if ofl_fi.exists and keyword_set(update) then begin
	if phot_ts lt ofl_fi.ctime and grade_ts lt ofl_fi.ctime then begin
		if keyword_set(verbose) then $
			print,'No update required for: ',band,' ',id
		return
	endif
endif
;
; get surface brightness range
	sbdel = 17.
	sbran = [13.+sbdel, 13.]
	if band eq 'FUV' or band eq 'NUV' then begin
		sbran = [18.+sbdel, 18.]
	endif
	if band eq 'u' or band eq 'g' or band eq 'r' or $
		band eq 'i' or band eq 'z' then begin
		sbran = [15.+sbdel, 15.]
	endif
	if band eq 'j' or band eq 'h' or band eq 'k' then begin
		sbran = [13.+sbdel, 13.]
	endif
	if band eq 'w1' or band eq 'w2' then begin
		sbran = [13.+sbdel, 13.]
	endif
	if band eq 'w3' or band eq 'w4' then begin
		sbran = [13.+sbdel, 13.]
	endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set some useful values to
; set up plots regardless of band

a=tot_a			; semi-major axes of apertures
ra_cen=ra_cen[0]
dec_cen=dec_cen[0]
semimajor=semimajor[0]
semiminor=semiminor[0]
pa=pa[0]>0.
asymag = asymag[0]
;
; get index of aperture magnitude semi-major axis
idxap=where(a ge semimajor, nrap)
if nrap gt 0 then $
	idxap=idxap[0] $
else	idxap=-1
;
; Asymptotic values
asyma=asyma[0]
;
; central sufrace brightness
mucent = ann_mu[0]
mucente = ann_mu_e[0]
;
; get indices based on asymtotic magnitude semi-major axis
idxsy=where(a eq asyma, nat)
if nat eq 1 and asyma gt 0. and finite(asyma) and $
	asymag gt 0. and grade le 2 then begin
	;
	; get magnitudes at fractional fluxes
	amg20 = asymag + 1.747	; -2.5 * alog10(0.20) [20%]
	amg50 = asymag + 0.753	; -2.5 * alog10(0.50) [50%]
	amg80 = asymag + 0.242	; -2.5 * alog10(0.80) [80%]
	amg90 = asymag + 0.114	; -2.5 * alog10(0.90) [90%]
	;
	; find monotonic section of growth curve
	del = tot_int - shift(tot_int,1)
	del[0] = 0.
	;
	; get radii from growth curve
	;
	; handle innermost radius
	; 20%
	ia = where(tot_mag lt amg20 and tot_mag gt 0. and del gt 0.)
	ia = ia[0]
	if ia le 0 then $
		ar20 = a[0] $
	else	linterp,tot_mag[ia-1:ia],a[ia-1:ia],amg20,ar20
	; 50%
	ia = where(tot_mag lt amg50 and tot_mag gt 0. and del gt 0.)
	ia = ia[0]
	if ia le -1 then $
		ar50 = -1 $
	else	linterp,tot_mag[ia-1:ia],a[ia-1:ia],amg50,ar50
	; 80%
	ia = where(tot_mag lt amg80 and tot_mag gt 0. and del gt 0.)
	ia = ia[0]
	if ia le -1 then $
		ar80 = -1. $
	else	linterp,tot_mag[ia-1:ia],a[ia-1:ia],amg80,ar80
	;
	; 90%
	ia = where(tot_mag lt amg90 and tot_mag gt 0. and del gt 0.)
	ia = ia[0]
	if ia le -1 then $
		ar90 = -1 $
	else	linterp,tot_mag[ia-1:ia],a[ia-1:ia],amg90,ar90
	;
	; get surface brightnesses
	ia = where(a gt ar20) & ia = ia[0]
	linterp,a[ia-1:ia],ann_mu[ia-1:ia],ar20,mu20
	linterp,a[ia-1:ia],ann_mu_e[ia-1:ia],ar20,mu20e
	if ar50 gt 0. then begin
		ia = where(a gt ar50) & ia = ia[0]
		linterp,a[ia-1:ia],ann_mu[ia-1:ia],ar50,mu50
		linterp,a[ia-1:ia],ann_mu_e[ia-1:ia],ar50,mu50e
	endif else begin
		mu50 = -1.
		mu50e = -1.
	endelse
	if ar80 gt 0. then begin
		ia = where(a gt ar80) & ia = ia[0]
		linterp,a[ia-1:ia],ann_mu[ia-1:ia],ar80,mu80
		linterp,a[ia-1:ia],ann_mu_e[ia-1:ia],ar80,mu80e
	endif else begin
		mu80 = -1.
		mu80e = -1.
	endelse
	if ar90 gt 0. then begin
		ia = where(a gt ar90) & ia = ia[0]
		linterp,a[ia-1:ia],ann_mu[ia-1:ia],ar90,mu90
		linterp,a[ia-1:ia],ann_mu_e[ia-1:ia],ar90,mu90e
	endif else begin
		mu90 = -1.
		mu90e = -1.
	endelse
	;
	; check for noisy surface brightness
	if mu90e gt 1.0 then begin
		while ann_mu_e[ia] gt 1.0 and ia ge 1 do ia = ia - 1
		if ia gt 0 then begin
			mu90 = ann_mu[ia]
			mu90e = 1.0
		endif
	endif
endif else begin
	print,'Unable to set values at asymptotic semi-major axis = ',asyma
	idxsy = -1
	mu20 = -1
	mu20e = -1
	mu50 = -1
	mu50e = -1
	mu80 = -1
	mu80e = -1
	mu90 = -1
	mu90e = -1
	asyma = -1
	ar20 = -1
	ar50 = -1
	ar80 = -1
	ar90 = -1
endelse

;;;;;;;;;;;;;;;;;;;;;;
; start profile plots
if keyword_set(verbose) then $
	print,band+' plot for '+id+' a = ',semimajor/60.,"' ...", $
		format='($,a,f6.2,a)'

;
; output file
filename=outpath+'/'+id+'_'+band+'.ps'
ms_ps_start, filename=filename, xsize=7,ysize=10,/inch,$
	/color,/true,bits=8, xoffset=0.75, yoffset = 0.75

!p.multi=[0,2,3]
!p.charsize=1.2

psize=0.5

;              20%                 50%              80% 
cs = [cgcolor('magenta'), cgcolor('blue'), cgcolor('brown'), $
;              90%               ASY                AP
      cgcolor('forest green'), cgcolor('orange'), cgcolor('red')]
;
ls = [2,1,4,3,5,0]

;growth curve
bad=where(tot_mag lt 0., nbad)
if nbad gt 0 then tot_mag[bad] = !values.f_nan
if grade gt 2 then tot_mag[*] = -9.99

if asymag gt 0. then $
	nminmax=minmax([tot_mag,asymag], /nan) $
else	nminmax=minmax(tot_mag, /nan)
min=nminmax[0]
max=nminmax[1]
del = max-min

xrng = [min(a/60.)-(max(a/60.)*.05), max(a/60.)*1.025]
yrng = [max+del*0.05, min-del*0.1]

plot,[a/60],[tot_mag],yr=yrng,/ys,/nodata,xr = xrng, /xs,$
	xtitle='R!da!n [arcmin]',ytitle=band+' Growth Curve [AB mag]'

;
; make sure we are including the sky annulus in the plot
if a[n_elements(a)-1] gt skyradius_in[0] then begin
	oplot,[skyradius_in[0],skyradius_out[0]]/60., $
		[min-del*0.03,min-del*0.03], color=cgcolor('sky blue')
	xyouts,skyradius_in[0]/60.,min-del*0.04,'BG ANN', $
		color=cgcolor('sky blue'), charsize=0.5
endif

if ar20 ge 0 then begin
	oplot,[-100,ar20/60],[amg20,amg20],color=cs[0],linesty=2
	oplot,[ar20/60,ar20/60],[100,amg20],color=cs[0],linesty=2
endif
if ar50 ge 0 then begin
	oplot,[-100,ar50/60],[amg50,amg50],color=cs[1],linesty=0
	oplot,[ar50/60,ar50/60],[100,amg50],color=cs[1],linesty=0
endif
if ar80 ge 0 then begin
	oplot,[-100,ar80/60],[amg80,amg80],color=cs[2],linesty=4
	oplot,[ar80/60,ar80/60],[100,amg80],color=cs[2],linesty=4
endif
if ar90 ge 0 then begin
	oplot,[-100,ar90/60],[amg90,amg90],color=cs[3],linesty=0
	oplot,[ar90/60,ar90/60],[100,amg90],color=cs[3],linesty=0
endif
;
; plot semimajor axis
sbeps = (yrng[0] - yrng[1]) * (3./sbdel)
oplot,[semimajor/60,semimajor/60],[yrng[0],yrng[0]-sbeps],color=cs[5],linesty=0
oplot,[semimajor/60,semimajor/60],[yrng[1],yrng[1]+sbeps],color=cs[5],linesty=0
;
; asymptotic mag
oplot,[-100,100],[asymag,asymag],color=cs[4],linesty=5
;
; plot total profile
if grade le 2 then begin
	oplot,[a/60.],[tot_mag+tot_mag_e],color=cgcolor('grey')
	oplot,[a/60.],[tot_mag-tot_mag_e],color=cgcolor('grey')
	oplot,[a/60.],[tot_mag]
endif

legend,['20%','50%','80%','90%','ASY','APR'],$
	textcolor=cs,/right,/bot,box=0,spac=1.5,charsize=0.75
;
;annular surface brightness
bad=where(ann_mu lt 0., nbad)
if nbad gt 0 then ann_mu[bad] = !values.f_nan
if grade gt 2 then ann_mu = -9.99

; get range from within ap
good = where(a lt a[idxap], ngood)
if ngood le 0 then $
	good = lindgen(n_elements(a))
yrng = sbran


plot,[a/60],[ann_mu],yr=yrng,/ys,/nodata,$
	xr = xrng, /xs,xtit='R!da!n [arcmin]', $ ;/xlog,$
	ytit=greek('mu',/append)+'!D'+band+'!N [AB mag/arcsec!u2!n]'

if ar50 ge 0 then begin
	oplot,[xrng[0],ar50/60],[mu50,mu50],color=cs[1],linesty=0
	oplot,[ar50/60,ar50/60],[100,mu50],color=cs[1],linesty=0
endif
if ar90 ge 0 then begin
	oplot,[xrng[0],ar90/60],[mu90,mu90],color=cs[3],linesty=0
	oplot,[ar90/60,ar90/60],[100,mu90],color=cs[3],linesty=0
endif
;
; plot annular profile
if grade le 2 then begin
	oplot,[a/60.],[ann_mu],color=cgcolor('grey')
	oploterror, [a/60.],[ann_mu],[ann_mu_e*0],[ann_mu_e],$
		psym=sym(1, psize = psize), /nohat
endif
;
; plot semimajor axis
oplot,[semimajor/60,semimajor/60],[yrng[0],yrng[0]-3.],color=cs[5],linesty=0
oplot,[semimajor/60,semimajor/60],[yrng[1],yrng[1]+3.],color=cs[5],linesty=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in image

if keyword_set(intfile) and file_exist(intfile) then begin
	img=mrdfits(intfile,0,hdr,/silent,/fscale)
	extast,hdr,astr
	AD2XY, ra_cen, dec_cen, astr, x0, y0
	getrot,hdr,rot,cdelt
	pixscl=abs(cdelt[0]*3600.)
	diff = abs(pixscl - as_pix[0])
	if diff gt 1.e-2 then $
		print,'Warning - image/photometry pixel scale mismatch: ',diff
	x0=x0[0]
	y0=y0[0]
	sz=size(img,/dim)
	sqasec = sz[0]*pixscl * sz[1]*pixscl
;
; get band-specific infor and perform band-specific processing
	if band eq 'FUV' or band eq 'NUV' then begin
		nimgs = sxpar(hdr,'nadded')
		exptim = sxpar(hdr,'medrr')
		if band eq 'FUV' then $
			img = smooth(img,2,/nan)>0
	endif
	if band eq 'u' or band eq 'g' or band eq 'r' or $
		band eq 'i' or band eq 'z' then begin
		nimgs = 1
		exptim = 53.9
	endif
	if band eq 'j' or band eq 'h' or band eq 'k' then begin
		nimgs = 1
		exptim = 1.3
	endif
	if band eq 'w1' or band eq 'w2' then begin
		nimgs = sxpar(hdr,'numfrms') / ( (sqasec / 7806805.1d0) > 1.)
		exptim = 7.7 * nimgs
	endif
	if band eq 'w3' or band eq 'w4' then begin
		nimgs = sxpar(hdr,'numfrms') / ( (sqasec / 7806805.1d0) > 1.)
		exptim = 8.8 * nimgs
	endif
;
; read in mask
	if keyword_set(verbose) then $
		print,' mask',format='($,a)'
	maskimg = glga_getmask(maskimgfile,sz,astr,pixscl)
	mhdr=hdr
	
	if keyword_set(bandmaskimgfile) then begin
	    if file_test(bandmaskimgfile) eq 1 then begin
	        maskimg += glga_getmask(bandmaskimgfile,sz,astr,pixscl)
	    endif
	endif
	
	
	
;
; subtract sky
	if mjc le 0. then mjc = 1.0		; conversion from DN to mJy
	skypp = (bg[0]*as_pix[0]^2)/mjc		; sky per pixel
	skyepp= (bg_e[0]*as_pix[0]^2)/mjc	; sky err per pixel
	img = (img - skypp) + skyepp		; leave a small floor
	if keyword_set(verbose) then $
		print,' sky,sig:',bg,bg_e,format='($,a,2g11.5)'
;
; create false three color for scaling
	image=fltarr(sz[0],sz[1],3)
	image[*,*,0]=img
	image[*,*,1]=img
	image[*,*,2]=img

	scale = glga_getimplotscl(band)
	image = nw_scale_rgb(image,scales=([scale,scale,scale])*brite*5.)
	image = nw_arcsinh_fit(image,nonlinearity=nonlinearity)
	image = nw_fit_to_box(image,origin=origin)
	image = nw_float_to_byte(image)
	image = reform(image[*,*,0])
;
; get plotting params

	delta = skyradius_out[0]/60. * 1.01
	imscale = pixscl[0]/60.0 ; arcmins
	ximgrng = ([x0,x0-sz[0]]) * imscale
	yimgrng = ([0-y0, sz[1]-y0]) * imscale
	xpltrng = [delta,0-delta]
	ypltrng = [0-delta,delta]
	delta = semimajor/60. * 1.01
	xpltrng2 = [delta,0-delta]
	ypltrng2 = [0-delta,delta]

;;;;;;;;;;;;;;;;;;;;;;;;
; scale with mask

	img2=bytarr(sz[0],sz[1],3)
	mindx=where(maskimg ge 1,nmask)
	nandx = where(finite(img) ne 1,nnans)

	;
	; reverse scale iamges
	red=not bytscl(image)
	green=not bytscl(image)
	blue=not bytscl(image)

	; are there pixels to mask?
	if nmask gt 0 then begin
		c=cgcolor('Dodger Blue', /Triple)
		red[mindx]=c[0,0]
		green[mindx]=c[0,1]
		blue[mindx]=c[0,2]
	endif

	img2[*,*,0]=red
	img2[*,*,1]=green
	img2[*,*,2]=blue

	loadct,0,/silent
	plotimage, img2, /preserve, color=cgcolor('black'), $
		xran=xpltrng, yran=ypltrng, $
		imgxrange=ximgrng, imgyrange=yimgrng, $
 		xtitle='arcmin', ytitle='arcmin'

	ratio=semiminor/semimajor

	ela=[ar50/60.0,ratio*ar50/60.0,0,0,float(90.-pa,0)]
	tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
		linestyle=0,color=cgcolor("blue"), thick=4,noclip=0

	ela=[ar90/60.0,ratio*ar90/60.0,0,0,float(90.-pa,0)]
	tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
		linestyle=0,color=cgcolor("forest green"), thick=4,noclip=0

	ela=[semimajor/60.0,semiminor/60.0,0,0,float(90.-pa,0)]
	tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
	   linestyle=0,color=cgcolor("red"), thick=4,noclip=0

	eli=[skyradius_out/60.0,ratio*skyradius_out/60.0,0,0,float(90.-pa,0)]
	tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
	  linestyle=2,color=cgcolor('black'), thick=4,noclip=0
 
	elo=[skyradius_in/60.0,ratio*skyradius_in/60.0,0,0,float(90.-pa,0)]
	tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
	  linestyle=2,color=cgcolor('black'), thick=4,noclip=0
;
; create false three color for scaling
	image=fltarr(sz[0],sz[1],3)
	image[*,*,0]=img
	image[*,*,1]=img
	image[*,*,2]=img

	image = nw_scale_rgb(image,scales=([scale,scale,scale])*brite)
	image = nw_arcsinh_fit(image,nonlinearity=nonlinearity)
	image = nw_fit_to_box(image,origin=origin)
	image = nw_float_to_byte(image)
	image = reform(image[*,*,0])
	loadct,0,/silent
	;
	; reverse scale images
	red=not bytscl(image)
	green=not bytscl(image)
	blue=not bytscl(image)

	img2[*,*,0]=red
	img2[*,*,1]=green
	img2[*,*,2]=blue

	plotimage, img2, /preserve, color=cgcolor('black'), $
		xran=xpltrng2, yran=ypltrng2, $
		imgxrange=ximgrng, imgyrange=yimgrng, $
 		xtitle='arcmin', ytitle='arcmin'

;;;;;;;;;;;;;;;
;plot ellipses

	ela=[ar50/60.0,ratio*ar50/60.0,0,0,float(90.-pa,0)]
	tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
		linestyle=0,color=cgcolor("blue"), thick=4,noclip=0

	ela=[ar90/60.0,ratio*ar90/60.0,0,0,float(90.-pa,0)]
	tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
		linestyle=0,color=cgcolor("forest green"), thick=4,noclip=0

	ela=[semimajor/60.0,semiminor/60.0,0,0,float(90.-pa,0)]
	tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
		linestyle=0,color=cgcolor("red"), thick=4,noclip=0

	eli=[skyradius_out/60.0,ratio*skyradius_out/60.0,0,0,float(90.-pa,0)]
	tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
	  linestyle=2,color=cgcolor('black'), thick=4,noclip=0
 
	elo=[skyradius_in/60.0,ratio*skyradius_in/60.0,0,0,float(90.-pa,0)]
	tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
	  linestyle=2,color=cgcolor('black'), thick=4,noclip=0
;
endif else begin

	 plot,[0,1],[0,1],/nodata,xs=4,ys=4
	 xyouts,[0.32],[0.5],'No '+band+' Image',charsize=1, $
		 color=cgcolor("black")

	 plot,[0,1],[0,1],/nodata,xs=4,ys=4
	 xyouts,[0.3],[0.5],'No Masked Image',charsize=1, $
		 color=cgcolor("black")

endelse

;;;;;;;;;;;;;;;;;;
; text

!p.charsize=0.8

plot,[0,100],[0,100],xs=4,ys=4,/nodata, pos = [0.085, 0.01, 0.98, 0.33]

; galactic coords to index dust maps
euler, ra_cen, dec_cen, l, b, 1
ebv=dust_getval(l,b,/interp)			; SFD98 E(B-V)
a_gal=glga_getextin(ebv,band,yuan13=yuan13)	; get MW extinction for band
if a_gal lt 1.e-3 then a_gal = 0.

chw = 1.24
xyouts,0,85,id,charsi=1.5
if strlen(srv) gt 0 then $
	xyouts,0,75,srv+' '+band,charsi=1.25 $
else	xyouts,0,75,band,charsi=1.25
rastr = strtrim(strn(ra_cen,format='(f12.6)'),2)
decstr = strtrim(strn(dec_cen,format='(f12.6)'),2)
nrachr = strlen(rastr)
ndechr = strlen(decstr)
coolab = 'R.A., Dec. [J2K, deg]:'
xyouts,30,85,coolab
off = chw*((10-nrachr)>0)
xyouts,55+off,85,rastr
off = chw*((10-ndechr)>0)
xyouts,71+off,85,decstr
lonstr = strtrim(strn(l,format='(f12.6)'),2)
latstr = strtrim(strn(b,format='(f12.6)'),2)
nlochr = strlen(lonstr)
nlachr = strlen(latstr)
gallab = 'Lon., Lat. [Gal., deg]:'
xyouts,30,75,gallab
off = chw*((10-nlochr)>0)
xyouts,55+off,75,lonstr
off = chw*((10-nlachr)>0)
xyouts,71+off,75,latstr

y0 = 60
del = 9
xyouts,0,y0,'GAL E(B-V):  '+strn(ebv,format='(f5.3)')
xyouts,0,y0-del,'GAL A!d'+band+'!n:  '+strn(a_gal,format='(f5.3)')
xyouts,0,y0-del*2.,"SEMIMAJOR:  "+strn(semimajor/60,format='(f6.2)')+"'"
 
xyouts,0,y0-del*3.,'RATIO (a/b):  '+strn(semimajor/semiminor,format='(f6.2)')
xyouts,0,y0-del*4.,'P.A.:  '+strn(pa,format='(f6.2)')+greek("degrees",/append)
xyouts,0,y0-del*5.,'Annuli:  '+strn(annuli_size,format='(f6.2)')+'"'
xyouts,0,y0-del*6.,'TYPE:  '+typ


; time stamp
xyouts,30,y0-del*6.,'Phot Date: '+systime(0,phot_ts)

if tf_mag lt 0. or grade gt 2 then  $
	outstr = 'APR:  m!d0!n = -9.99   m!ddr!n = -9.99 ' + $
		greek("plus_minus", /append) + ' -9.99' $
else	outstr = 'APR:  m!do!n = ' + $
 		strn(tf_mag>0.,format='(f5.2)')+'   m!ddr!n = ' + $
 		strn((tf_mag-a_gal)>0.,format='(f5.2)') + ' ' + $
		greek("plus_minus", /append) + $
 		' '+strn(tf_mag_e>0.01,format='(F4.2)')
xyouts,30,y0,outstr

if asymag lt 0. or grade gt 2 then $
	outstr = 'ASY:  m!do!n = -9.99   m!ddr!n = -9.99 ' + $
 		greek("plus_minus", /append) + ' -9.99' $
else	outstr = 'ASY:  m!do!n = ' + $
		strn(asymag>0.,format='(f5.2)')+'   m!ddr!n = ' + $
		strn((asymag-a_gal)>0.,format='(f5.2)') +' ' + $
		greek("plus_minus", /append) + $
		' '+string(asymag_e>0.01,format='(F4.2)')
xyouts,30,y0-del,outstr

if mu90 lt 0 then $
	outstr = '90%:  '+greek("mu",/append)+'!do!n = -9.99' + $
 		'   '+greek("mu",/append)+'!ddr!n = -9.99 ' + $
 		greek("plus_minus", /append) + ' -9.99' $
else	outstr = '90%:  '+greek("mu",/append)+'!do!n = ' + $
		strn(mu90>0.,format='(f5.2)')+'   '+ $
		greek("mu",/append)+'!ddr!n = ' + $
		strn((mu90-a_gal)>0.,format='(f5.2)') +' '+ $
		greek("plus_minus", /append) + $
 		' '+string(mu90e>0.01,format='(F4.2)')
xyouts,30,y0-del*2.,outstr

if mu50 lt 0 then $
	outstr = '50%:  '+greek("mu",/append)+'!do!n = -9.99' + $
		'   '+greek("mu",/append)+'!ddr!n = -9.99 ' + $
		greek("plus_minus", /append) + ' -9.99' $
else	outstr = '50%:  '+greek("mu",/append)+'!do!n = ' + $
		strn(mu50>0.,format='(f5.2)')+'   '+ $
		greek("mu",/append)+'!ddr!n = ' + $
		strn((mu50-a_gal)>0.,format='(f5.2)') +' '+ $
		greek("plus_minus", /append) + $
		' '+string(mu50e>0.01,format='(F4.2)')
xyouts,30,y0-del*3.,outstr

if ar20 le 0 or ar80 le 0 or grade gt 2 then $
	outstr = 'Concentration (C82) = -9.99' $
else	outstr = 'Concentration (C82) = '+ $
		strn((ar80/ar20)>0.,format='(f5.2)')
xyouts,30,y0-del*4.,outstr

xyouts,30,y0-del*5.,'EXPTIME: '+strn(exptim,format='(F8.1)') + ' s' + $
	'  # imgs: '+strn(fix(nimgs))

if mu_bg lt 0. or grade gt 2 then $
	outstr = greek("mu", /append)+'BG = -9.99 ' + $
		greek("plus_minus", /append)+ ' -9.99' $
else	outstr = greek("mu", /append)+'BG = ' + $
		strn(mu_bg>0.,format='(f5.2)')+' '+ $
		greek("plus_minus", /append)+$
		' '+strn(mu_bg_e>0.01<9.99,format='(F5.2)')
xyouts,78,y0,outstr

xyouts,78,y0-del, 'R!dASY!n = '+string(asyma/60,format='(F6.2)')+"'"
xyouts,78,y0-del*2., 'R!d90!n = '+string(ar90/60,format='(F6.2)')+"'"
xyouts,78,y0-del*3., 'R!d80!n = '+string(ar80/60,format='(F6.2)')+"'"
xyouts,78,y0-del*4., 'R!d50!n = '+string(ar50/60,format='(F6.2)')+"'"
xyouts,78,y0-del*5., 'R!d20!n = '+string(ar20/60,format='(F6.2)')+"'"

xyouts,78,y0-del*6., 'Grade = '+string(grade,format='(I3)')

plots,[0.06,.98,.98,0.06,0.06],[.32,.32,.01,0.01,0.32],/norm, $
	color=cgcolor('Charcoal')

plots,[0.06,0.98],[0.235,0.235],/norm,color=cgcolor('Charcoal')

xyouts,0.06, 0.002,systime(),charsize=.5,/norm

!p.multi=0

ms_ps_end

;;;;;;;;;;;;;;;;;;;;;;
; convert to pdf & jpg

p=outpath+'/'
name=id+'_'+band

if keyword_set(rotatee) then $
    spawn,'convert +matte -rotate -90 -density 196 -resample 72 '+p+name+'.ps '+p+name+'.jpg' $
else $
    spawn,'convert +matte -density 196 -resample 72 '+p+name+'.ps '+p+name+'.jpg'
   
spawn,'ps2pdf '+p+name+'.ps '+p+name+'.pdf'
spawn,'rm '+p+name+'.ps'

;;;;;;;;;;;;
; all done
if keyword_set(verbose) then print,' Done.'

return
end
