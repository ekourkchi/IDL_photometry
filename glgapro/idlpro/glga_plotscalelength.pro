pro glga_plotscalelength, id, band, survey=survey, $
    outpath=outpath, verbose=verbose, $
    pathtoprofile=pathtoprofile, yuan13=yuan13


if not keyword_set(survey) then srv = '' else srv = strtrim(survey,2)
if not keyword_set(outpath) then outpath='./'
if not keyword_set(pathtoprofile) then pathtoprofile='./'


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
	if band eq 'V' or band eq 'I' then begin
		sbran = [17.+sbdel, 20.]
	endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
;
;;;;;;;;;;;;;;;;;;;;;;
; start profile plots
if keyword_set(verbose) then $
	print,band+' plot for '+id+' a = ',semimajor/60.,"' ...", $
		format='($,a,f6.2,a)'


; output file
filename=outpath+'/'+id+'_'+band+'.scales.ps'
ms_ps_start, filename=filename, xsize=7,ysize=10,/inch,$
	/color,/true,bits=8, xoffset=0.85, yoffset = 0.75
!p.multi=[0,1,2]	
!p.charsize=1.2
psize=0.5	
	
;              20%                 50%              80% 
cs = [cgcolor('magenta'), cgcolor('blue'), cgcolor('brown'), $
;              90%               ASY                AP
      cgcolor('forest green'), cgcolor('orange'), cgcolor('red')]
;
ls = [2,1,4,3,5,0]	


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
xrng = [min([min(a/60.)-(max(a/60.)*.05),0.]), max(a/60.)*1.025]

band_txt = band
if band eq 'V' then band_txt = 'F606W'
if band eq 'I' then band_txt = 'F814W'

plot,[a/60],[ann_mu],yr=yrng,/ys,/nodata,$
	xr = xrng, /xs,xtit='R!da!n [arcmin]', $ ;/xlog,$
	ytit=greek('mu',/append)+'!D'+band_txt+'!N [AB mag/arcsec!u2!n]'

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

; Ra_min = ar50
; Ra_max = ar90
; ind = where(a ge Ra_min and a le Ra_max)
; Ra   = a[ind]
; Mu   = ann_mu[ind]
; Mu_e = ann_mu_e[ind]
; 
; line_u = LINFIT(Ra, Mu, MEASURE_ERRORS=Mu_e)
; 
; Mu_0 = line_u[0]
; beta = line_u[1]
; Xu = [0, 1000]
; Yu = beta*Xu+Mu_0
; 
; oplot, Xu/60.,Yu, linestyle = 2, color=cgcolor('brown')


Ra_min = 29.0 ;ar50
ind = where(a ge Ra_min and ann_mu le 30.0)
Ra   = a[ind]
Mu   = ann_mu[ind]
Mu_e = ann_mu_e[ind]

m_pivot = Mu[0]
Ra_   = []
Mu_   = []
Mu_e_ = []
n_i = 0
for i=0, n_elements(Mu)-1 do begin
   if Mu[i] ge m_pivot then begin
       m_pivot = Mu[i] 
       Ra_ = [Ra_, Ra[i]]
       Mu_ = [Mu_ , Mu[i]]
       Mu_e_ = [Mu_e_, Mu_e[i]]
       n_i++
   endif
endfor


if n_i ge 2 then begin
    line_u = LINFIT(Ra_/60., Mu_)

    Mu_0 = line_u[0]
    beta = line_u[1]
    Xu = [0, 1000]
    Yu = beta*Xu+Mu_0
    oplot, Xu,Yu, linestyle = 2, color=cgcolor('brown')

    c = 2.5*alog10(exp(1.))
    h = c/beta
endif else begin
    c = -9.99
    h = -9.99
    Mu_0 = -9.99
endelse



mu_i = []
a_i = []
n_i = 0
mag_i = []

if n_elements(ann_mu) ge 2 and n_elements(ann_mu) eq n_elements(a) and n_elements(ann_mu) eq n_elements(tot_mag) then begin
    for i=0, n_elements(ann_mu)-2 do begin
       if ann_mu[i] le 29.0 and ann_mu[i+1] gt 29.0 then begin 
            mu_i=[ann_mu[i],ann_mu[i+1]]
            a_i = [a[i],a[i+1]]
            mag_i = [tot_mag[i], tot_mag[i+1]]
            n_i++
            break
        endif
    endfor  
endif


if n_i gt 0 then begin 
    linterp,mu_i,a_i,28.,ar255   ; 29.0
    linterp,a_i,mag_i,ar255, mag255
endif else begin 
    ar255=-9.99
    mag255=-9.99
endelse



mu_i = []
mu_ei = []
a_i = []
n_i = 0
if n_elements(ann_mu) ge 2 and n_elements(ann_mu) eq n_elements(a) then begin
    for i=0, n_elements(ann_mu)-1 do begin
       if ann_mu[i] gt 0. then begin 
            mu_i=[mu_i,ann_mu[i]]
            a_i = [a_i,a[i]]
            n_i++
        endif else break
    endfor  
endif

if n_i ge 2 then begin 
    l = LINFIT(a_i[0:1],mu_i[0:1])
    mu_00 = l[0]
endif else begin 
    mu_00=-9.99
endelse



if Mu_0 gt 0 then begin
    delta_n = (29.0-Mu_0)/c
    delta_m_ext = 2.5*alog10(1.-(1.+delta_n)*exp(-1.*delta_n))
endif else begin
    delta_n  = -9.99
    delta_m_ext = -9.99
endelse



x0_legend = xrng[0]+0.8*(xrng[1]-xrng[0])
y0_legend = yrng[0]+0.85*(yrng[1]-yrng[0])


legend,['50%','90%','Aperture', 'Best Fit'],$
	textcolor=[cs[1],cs[3],cs[5],cs[2]],box=0,spac=1.5,charsize=0.85, position=[x0_legend,y0_legend]

;;;;;;;;;;;;;;;;;;
; text in the plot

!p.charsize=0.8

plot,[0,100],[0,100],xs=4,ys=4,/nodata, pos = [0.085, 0.11, 0.98, 0.43]

; galactic coords to index dust maps
euler, ra_cen, dec_cen, l, b, 1
ebv=dust_getval(l,b,/interp)			; SFD98 E(B-V)
a_gal=glga_getextin(ebv,band,yuan13=yuan13)	; get MW extinction for band
if a_gal lt 1.e-3 then a_gal = 0.

chw = 1.24
xyouts,0,95,id,charsi=1.5
if strlen(srv) gt 0 then $
	xyouts,0,85,srv+' '+band,charsi=1.25 $
else	xyouts,0,85,band,charsi=1.25
rastr = strtrim(strn(ra_cen,format='(f12.6)'),2)
decstr = strtrim(strn(dec_cen,format='(f12.6)'),2)
nrachr = strlen(rastr)
ndechr = strlen(decstr)
coolab = 'R.A., Dec. [J2K, deg]:'
xyouts,30,95,coolab
off = chw*((10-nrachr)>0)
xyouts,55+off,95,rastr
off = chw*((10-ndechr)>0)
xyouts,71+off,95,decstr
lonstr = strtrim(strn(l,format='(f12.6)'),2)
latstr = strtrim(strn(b,format='(f12.6)'),2)
nlochr = strlen(lonstr)
nlachr = strlen(latstr)
gallab = 'Lon., Lat. [Gal., deg]:'
xyouts,30,85,gallab
off = chw*((10-nlochr)>0)
xyouts,55+off,85,lonstr
off = chw*((10-nlachr)>0)
xyouts,71+off,85,latstr

y0 = 70
del = 9
xyouts,0,y0,'GAL E(B-V):  '+strn(ebv,format='(f5.3)')
xyouts,0,y0-del,'GAL A!d'+band+'!n:  '+strn(a_gal,format='(f5.3)')
xyouts,0,y0-del*2.,"SEMIMAJOR:  "+strn(semimajor/60.,format='(f6.2)')+"'"
 
xyouts,0,y0-del*3.,'RATIO (a/b):  '+strn(semimajor/semiminor,format='(f6.2)')
xyouts,0,y0-del*4.,'P.A.:  '+strn(pa,format='(f6.2)')+greek("degrees",/append)
xyouts,0,y0-del*5.,'Annuli:  '+strn(annuli_size,format='(f6.2)')+'"'



if asymag lt 0. or grade gt 2 then $
	outstr = 'ASY:  m = -9.99   ' $
else	outstr = 'ASY:  m = ' + strn(asymag>0.,format='(f5.2)')+' [mag]'
xyouts,30,y0-0.*del,outstr



if asymag lt 0. or grade gt 2 then $
	outstr = 'Central:  '+greek("mu",/append)+'!do!n = -9.99   ' $
else	outstr = 'Central:  '+greek("mu",/append)+'!do!n = ' + strn(mu_00>0.,format='(f5.2)')+' [mag/arcsec!u2!n]'
xyouts,30,y0-1.*del,outstr


if mu90 lt 0 then $
	outstr = '90%:  '+greek("mu",/append)+'= -9.99'  $
else	outstr = '90%:  '+greek("mu",/append)+'= ' + $
		strn(mu90>0.,format='(f5.2)')+' [mag/arcsec!u2!n]'
xyouts,30,y0-del*2.,outstr

if mu50 lt 0 then $
	outstr = '50%:  '+greek("mu",/append)+'= -9.99' $
else	outstr = '50%:  '+greek("mu",/append)+'= ' + $
		strn(mu50>0.,format='(f5.2)')+' [mag/arcsec!u2!n]'
xyouts,30,y0-del*3.,outstr

if ar20 le 0 or ar80 le 0 or grade gt 2 then $
	outstr = 'Concentration (C82) = -9.99' $
else	outstr = 'Concentration (C82) = '+ $
		strn((ar80/ar20)>0.,format='(f5.2)')
xyouts,30,y0-del*4.,outstr

if mag255 lt 0. or grade gt 2 then $
	outstr = 'm!d29.0!n = -9.99   ' $
else	outstr = 'm!d29.0!n = ' + strn(mag255>0.,format='(f5.2)')+' [mag]'
xyouts,30,y0-5.*del,outstr


xyouts,30,y0-6.*del,'Disk '+greek("mu",/append)+'!do!n = '+string(Mu_0,format='(F5.2)')+' [mag/arcsec!u2!n]'
xyouts,30,y0-7.*del,'Disk Scale Length h = '+string(h,format='(F6.2)')+"'"
xyouts,30,y0-8.*del,greek("Delta",/append)+'m!dext!n = '+string(delta_m_ext,format='(F10.4)')+' [mag]'

if (mag255 gt 0. and delta_m_ext gt -2) then $ 
       xyouts,30,y0-9.*del,'m!dtotal!n = m!d29.0!n + '+ greek("Delta",/append)+'m!dext!n = '+string(delta_m_ext+mag255,format='(F10.2)')+' [mag]'

xyouts,78,y0-del*5.,'R!d29.0!n = '+string(ar255/60.,format='(F6.2)')+"'"

xyouts,78,y0, 'R!dASY!n = '+string(asyma/60.,format='(F6.2)')+"'"
xyouts,78,y0-del*1., 'R!d90!n = '+string(ar90/60.,format='(F6.2)')+"'"
xyouts,78,y0-del*2., 'R!d80!n = '+string(ar80/60.,format='(F6.2)')+"'"
xyouts,78,y0-del*3., 'R!d50!n = '+string(ar50/60.,format='(F6.2)')+"'"
xyouts,78,y0-del*4., 'R!d20!n = '+string(ar20/60.,format='(F6.2)')+"'"

xyouts,78,y0-del*6., 'Grade = '+string(grade,format='(I3)')

plots,[0.06,.98,.98,0.06,0.06],[.44,.44,.05,0.05,0.44],/norm, $
	color=cgcolor('Charcoal')

plots,[0.06,0.98],[0.36,0.36],/norm,color=cgcolor('Charcoal')

xyouts,0.07, 0.038,'All magnitudes are raw and should be de-reddened.',charsize=.7,/norm
xyouts,0.07, 0.055,systime(),charsize=.6,/norm


;;;;;;;;;;;;;;;;;;

!p.multi=0

ms_ps_end
;;;;;;;;;;;;;;;;;;;;;;
; convert to pdf & jpg

p=outpath+'/'
name=id+'_'+band

if keyword_set(rotatee) then $
    spawn,'convert +matte -rotate -90 -density 196 -resample 72 '+p+name+'.scales.ps '+p+name+'.scales.jpg' $
else $
    spawn,'convert +matte -density 196 -resample 72 '+p+name+'.scales.ps '+p+name+'.scales.jpg'
   
spawn,'ps2pdf '+p+name+'.scales.ps '+p+name+'.scales.pdf'
spawn,'rm '+p+name+'.scales.ps'

;;;;;;;;;;;;
;; Generating the output text file

if not keyword_set(pathtoprofile) then pathtoprofile="./"
; output file
outfile=pathtoprofile+'/'+id+'_'+band+'.scales.dat'
openw,lun,outfile,/get_lun
printf,lun,'# ID   : '+id
printf,lun,'# Band : '+band_txt
printf,lun,'# RA [deg]: '+rastr
printf,lun,'# Dec [deg]: '+decstr
printf,lun,'# Lon. [Galactic deg]: '+lonstr
printf,lun,'# Lat. [Galactic deg]: '+latstr
printf,lun,'# GAL E(B-V): '+strn(ebv,format='(f5.3)')
; printf,lun,'# Extinction GAL A_'+band+': '+strn(a_gal,format='(f5.3)')
printf,lun,'# Semi-Major [arcmin]: '+strn(semimajor/60.,format='(f6.2)')
printf,lun,'# Ratio (a/b): '+strn(semimajor/semiminor,format='(f6.2)')
printf,lun,'# P.A. [deg]: '+strn(pa,format='(f6.2)')
printf,lun,'# Annuli [arcsec]: '+strn(annuli_size,format='(f6.2)')
printf,lun,'# ASY m_total: '+strn(asymag>0.,format='(f6.3)')
printf,lun,'# Central Surface Brightness mu_0: '+strn(mu_00>0.,format='(f6.3)')
printf,lun,'# 90% mu_90: '+strn(mu90>0.,format='(f6.3)')
printf,lun,'# 50% mu_50: '+strn(mu50>0.,format='(f6.3)')
printf,lun,'# Concentration (C82): '+strn((ar80/ar20)>0.,format='(f6.2)')
printf,lun,'# Magnitude within R_29.0 (m_29.0): '+strn(mag255>0.,format='(f6.3)')
printf,lun,'# Disk mu_0: '+string(Mu_0,format='(F6.3)')
printf,lun,'# Disk Scale Length h [arcmin]: '+string(h,format='(F6.2)')
printf,lun,'# Magnitude correction based on disk extrapolation (delta_m_ext): '+string(delta_m_ext,format='(F10.4)')
printf,lun,'# Asymptotic Radius (R_asy) [arcmin]: '+string(asyma/60.,format='(F6.2)')
printf,lun,'# R_90%  [arcmin]: '+string(ar90/60.,format='(F6.2)')
printf,lun,'# R_80%  [arcmin]: '+string(ar80/60.,format='(F6.2)')
printf,lun,'# R_50%  [arcmin]: '+string(ar50/60.,format='(F6.2)')
printf,lun,'# R_20%  [arcmin]: '+string(ar20/60.,format='(F6.2)')
printf,lun,'# R_29.0 [arcmin]: '+string(ar255/60.,format='(F6.2)')
printf,lun,'# Note: All magnitudes are raw and should be de-reddened.'
printf,lun,'# Date : '+systime()
printf,lun,'# '
printf,lun,'ID,ra,dec,lon,lat,ebv,a_gal,semimajor,a_b,pa,annuli,m_asy,central_mu,mu_90,mu_50,concentration,m_255,disc_mu0,scale_length_h,R_asy,R_90,R_80,R_50,R_20,R_255, d_m_ext'
printf, lun,id+','+rastr+','+decstr+','+lonstr+','+latstr+','+strn(ebv,format='(f5.3)')+','+strn(a_gal,format='(f5.3)')+','+strn(semimajor/60.,format='(f6.2)')+','+strn(semimajor/semiminor,format='(f6.2)')+','+strn(pa,format='(f6.2)')+','+strn(annuli_size,format='(f6.2)')+','+strn(asymag>0.,format='(f6.3)')+','+strn(mu_00>0.,format='(f6.3)')+','+strn(mu90>0.,format='(f6.3)')+','+strn(mu50>0.,format='(f6.3)')+','+strn((ar80/ar20)>0.,format='(f6.2)')+','+strn(mag255>0.,format='(f6.3)')+','+string(Mu_0,format='(F6.3)')+','+string(h,format='(F6.2)')+','+string(asyma/60.,format='(F6.2)')+','+string(ar90/60.,format='(F6.2)')+','+string(ar80/60.,format='(F6.2)')+','+string(ar50/60.,format='(F6.2)')+','+string(ar20/60.,format='(F6.2)')+','+string(ar255/60.,format='(F6.2)')+','+string(delta_m_ext,format='(F10.4)')
 free_lun,lun









; all done
if keyword_set(verbose) then print,' Done.'

return



END




