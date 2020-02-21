


Pro scalelength




pathtoprofile = '/home/ehsan/db_esn/NGC6684/raw/data/282D/photometry'
id = 'NGC6684'
fpath = '/home/ehsan/db_esn/NGC6684/raw/data/282D/acs/fits/'
intfile=fpath+id+'_I.fits'
outpath = '/home/ehsan/db_esn/NGC6684/raw/data/282D/plots/'
filename=outpath+id+'_acs_profile.esn.ps'
jpgpath='/home/ehsan/db_esn/NGC6684/raw/data/282D/acs/jpg'

bands = ['F606W','F814W']
nband = n_elements(bands)

if not keyword_set(type) then typ = '-' else typ = strtrim(type,2)
if not keyword_set(pathtoprofile) then pathtoprofile='./'
if not keyword_set(jpgpath) then jpgpath='./'
b0 = 0
b1 = 1


read_radprof,id,bands[b0], utot_a, utot_mag, utot_mag_e, uann_mu, uann_mu_e, $
	ura_cen, udec_cen, usemimajor, usemiminor, upa, uscale, $
	utf_mag, utf_mag_e, uaf_mag, uaf_mag_e, umu_bg, umu_bg_e, $
	uskyradius_in, uskyradius_out, utf_a, usymag, usymag_e, usyma, $
	pathtoprofile=pathtoprofile
read_radprof,id,bands[b1], gtot_a, gtot_mag, gtot_mag_e, gann_mu, gann_mu_e, $
	gra_cen, gdec_cen, gsemimajor, gsemiminor, gpa, gscale, $
	gtf_mag, gtf_mag_e, gaf_mag, gaf_mag_e, gmu_bg, gmu_bg_e, $
	gskyradius_in, gskyradius_out, gtf_a, gsymag, gsymag_e, gsyma, $
	pathtoprofile=pathtoprofile, annuli_size=annuli_size


;
; check inputs
if utot_a[0] eq -1 then nouband = (1 eq 1) else nouband = (1 eq 0)
if gtot_a[0] eq -1 then nogband = (1 eq 1) else nogband = (1 eq 0)

if nogband then begin
	print,'No primary image photometry for '+id+', returning.'
	return
endif
if nouband and nogband then return
; fix missing values
if nouband then begin
; 	utf_mag=!values.f_nan
	utf_mag_e=!values.f_nan
	uaf_mag=!values.f_nan
	uaf_mag_e=!values.f_nan
	usyma=!values.f_nan
	usymag=!values.f_nan
	usymag_e=!values.f_nan
endif
if nogband then begin
	gtf_mag=!values.f_nan
	gtf_mag_e=!values.f_nan
	gaf_mag=!values.f_nan
	gaf_mag_e=!values.f_nan
	gsyma=!values.f_nan
	gsymag=!values.f_nan
	gsymag_e=!values.f_nan
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set some useful values to
; set up plots regardless of band

 a=gtot_a
 tot_mag=gtot_mag
 ann_mu=gann_mu 
 ra_cen=gra_cen[0]
 dec_cen=gdec_cen[0]
 semimajor=gsemimajor[0]
 semiminor=gsemiminor[0]
 pa=gpa[0]>0.
 r_ap=semimajor
 r_ski=gskyradius_in
 r_sko=gskyradius_out

;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ms_ps_start, filename=filename, xsize=7.,ysize=8.,/inch,$
          /color,bit=8, xoffset=0.25, yoffset = 10.5, /landscape
; ps_on , filename, noask='y', /color, cc=cc



!p.multi=[0,1,2]
!p.charsize=1.5

setplotcolors
psize=0.75

cs = [!blue, !red]

; --------------------------------------------------------

;growth curve
if not nouband then begin
	bad=where(utot_mag lt 0., nbad)
	if nbad gt 0 then utot_mag[bad] = !values.f_nan
endif else utot_mag = !values.f_nan
if not nogband then begin
	bad=where(gtot_mag lt 0., nbad)
	if nbad gt 0 then gtot_mag[bad] = !values.f_nan
endif else gtot_mag = !values.f_nan

nminmax=minmax([utot_mag,gtot_mag], /nan)
min=nminmax[0]
max=nminmax[1]

ydel = max-min
yr0 = max+ydel*0.1
yr1 = min-ydel*0.2

xdel = r_sko[0]/60.
xr0 = min(a/60.) - xdel*0.02
xlr0 = (0.8*(annuli_size/60.))>0.01
xr1 = r_sko[0]/60. + xdel*0.02

; growth curve linear x-axis
plot,[a/60],[tot_mag],yr=[yr0, yr1],/ys,/nodata,$
    xr = [0, xr1], /xs,$
    xtit='R!da!n [arcmin]',ytit='Growth Curve [mag]', charsize=1.5
;
; make sure we are including the sky annulus in the plot
if a[n_elements(a)-1] gt r_ski then begin
	oplot,[r_ski,r_sko]/60.,[min-ydel*0.1,min-ydel*0.1], color=!blue
	xyouts,r_ski/60.,min-ydel*0.12,'SKY ANN', $
		color=!blue, charsize=0.5
endif

if not nouband then begin
 oploterror, [a/60.],[utot_mag],[utot_mag_e*0],[utot_mag_e],$
  color=cs[b0],errcolor=cs[b0],psym=0, /nohat  ; psym=-sym(1, psize = psize)
 oplot,!x.crange,[usymag,usymag],linesty=1,color=cs[b0]
 oplot,[usyma,usyma]/60.,[yr0,usymag],linesty=2,color=cs[b0],thick=4
endif

if not nogband then begin
 oploterror, [a/60.],[gtot_mag],[gtot_mag_e*0],[gtot_mag_e],$
  color=cs[b1],errcolor=cs[b1],psym=0, /nohat ; psym=-sym(1, psize = psize)
 oplot,!x.crange,[gsymag,gsymag],linesty=1,color=cs[b1]
 oplot,[gsyma,gsyma]/60.,[yr0,gsymag],linesty=2,color=cs[b1],thick=4
endif


oplot,[r_ap,r_ap]/60.,[yr0, yr1],thick=4,color=!brown, linestyle=5
xx = r_ap/60. - abs(!x.crange[1]-!x.crange[0])*0.09
yy = !y.crange[0] - abs(!y.crange[0]-!y.crange[1])*0.05
xyouts,0.45,21,'APERTURE',color=!brown,charsiz=0.8, orientation=90

; --------------------------------------------------------
; --------------------------------------------------------
; annular surface brightness

if not nouband then begin
	bad=where(uann_mu lt 0., nbad)
	if nbad gt 0 then uann_mu[bad] = !values.f_nan
endif else uann_mu = !values.f_nan
if not nogband then begin
	bad=where(gann_mu lt 0., nbad)
	if nbad gt 0 then gann_mu[bad] = !values.f_nan
endif else gann_mu = !values.f_nan

nminmax=minmax([uann_mu,gann_mu], /nan)

min=nminmax[0]
max=nminmax[1]

; annular surface brightness linear x-axis
plot,[a/60],[ann_mu],yr=[max+2.0, min-2.],/ys,/nodata,$
    xr = [0, xr1], /xs,$
    xtit='R!da!n [arcmin]',$
    ytit=greek('mu',/append)+' [mag/arcsec!u2!n]'

; --------------------------------------------------------
;
; --------------------------------------------------------

Ra = a/60  ;; --> arcmin

; disclude all bad values
ind = where(uann_mu gt -90.)
Ra = Ra[ind]
Mu = uann_mu[ind]
Mu_e = uann_mu_e[ind]

; Specifying Ra range 
Ra_min = 0.7
Ra_max = 0.9

ind = where(Ra gt Ra_min and Ra lt Ra_max)
Ra = Ra[ind]
Mu = Mu[ind]
Mu_e = Mu_e[ind]

line_u = LINFIT(Ra, Mu, MEASURE_ERRORS=Mu_e)

Mu_0 = line_u[0]
beta = line_u[1]
Xu = [0.,2]
Yu = beta*Xu+Mu_0



; --------------------------------------------------------
; --------------------------------------------------------
;
Ra = a/60

; disclude all bad values
ind = where(uann_mu gt -90.)
Ra = Ra[ind]
Mu = gann_mu[ind]
Mu_e = gann_mu_e[ind]

; Specifying Ra range 
Ra_min = 0.7
Ra_max = 0.9

ind = where(Ra gt Ra_min and Ra lt Ra_max)
Ra = Ra[ind]
Mu = Mu[ind]
Mu_e = Mu_e[ind]

line_g = LINFIT(Ra, Mu, MEASURE_ERRORS=Mu_e)

Mu_0 = line_g[0]
beta = line_g[1]
Xg = [0.,2]
Yg = beta*Xg+Mu_0


; --------------------------------------------------------

if not nouband then begin
 
 Ra = [a[0]/60.]
 Y = [uann_mu[0]]
 Ye = [uann_mu_e[0]]
 step = 0.015
 min_lim = a[0]/60. + step
 for i=1, n_elements(a)-1 do begin
    if a[i]/60. gt min_lim then begin
      Ra = [Ra, a[i]/60.]
      Y  = [Y, uann_mu[i]]
      Ye = [Ye, uann_mu_e[i]]
      min_lim+=step
    endif
  endfor
 
 
 oploterror, Ra,Y,Ye*0,Ye,$
  color=cs[b0],errcolor=cs[b0],psym=sym(1, psize = psize), /nohat
endif
  
  
if not nogband then begin

 Ra = [a[0]/60.]
 Y = [gann_mu[0]]
 Ye = [gann_mu_e[0]]
 step = 0.015
 min_lim = a[0]/60. + step
 for i=1, n_elements(a)-1 do begin
    if a[i]/60. gt min_lim then begin
      Ra = [Ra, a[i]/60.]
      Y  = [Y, gann_mu[i]]
      Ye = [Ye, gann_mu_e[i]]
      min_lim+=step
    endif
  endfor
 
  oploterror, Ra,Y,Ye*0,Ye,$
  color=cs[b1],errcolor=cs[b1],psym=sym(1, psize = psize), /nohat
endif

 oplot,[r_ap,r_ap]/60.,[max+2.,min-2.],thick=4,color=!brown, linestyle=5
 xx = r_ap/60. - abs(!x.crange[1]-!x.crange[0])*0.09
 yy = !y.crange[1] + abs(!y.crange[0]-!y.crange[1])*0.1
 xyouts,1.0,25,'APERTURE',color=!brown,charsiz=0.8, orientation=90

legend,['F606W','F814W'],textcolors=cs[[b0,b1]],/bot,/left, $
	box=0, charsize=1.0


oplot, Xu,Yu, linestyle = 1
oplot, Xg,Yg, linestyle = 1
oplot, [Ra_min, Ra_min],[27,26.5],color=50
oplot, [Ra_max, Ra_max],[27,26.5],color=50
oplot, [Ra_min, Ra_max],[26.75,26.75],color=50
xyouts,0.70,26,'Fitting Range',color=50,charsiz=0.7

c = 2.5*alog10(exp(1.))
Mu_0_u = line_u[0]
beta_u = line_u[1]
scale_u = c/beta_u
Mu_0_g = line_g[0]
beta_g = line_g[1]
scale_g = c/beta_g
print, 'Mu_0_u: ', Mu_0_u
print, 'scale_u: ', scale_u
print, 'Mu_0_g: ', Mu_0_g
print, 'scale_g: ', scale_g


ms_ps_end
;;;;;;;;;;;;;;;;;;;;;;
; convert to pdf & jpg

p=outpath+'/'
name=id+'_acs_profile.esn'
spawn,'convert +matte -density 196 -resample 72 -rotate -90 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'


;;;;;;;;;;;;;;;;;;;;;;

END
