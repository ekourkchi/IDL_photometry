pro wise_plotradprof, id, pathtoprofile=pathtoprofile, type=type, $
    intfile=intfile, maskimgfile=maskimgfile, auxbase=auxbase, $
    jpgpath=jpgpath, outpath=outpath, verbose=verbose, yuan13=yuan13
;
; note: set yuan13 keyword to use Yuan et al. 2013 empirical extintion coeff's
;	Yuan coeff's assume SFD98 over-estimates E(B-V) by 14% (see table 2)
;
; Warning this procedure spawns convert and ps2pdf14
; (without any checking) to convert postscript output 
; to pdf and jpeg files. Can just comment out if postscript
; is OK. 
;

bands = ['w1','w2','w3','w4']
nband = n_elements(bands)

if not keyword_set(type) then typ = '-' else typ = strtrim(type,2)
if not keyword_set(pathtoprofile) then pathtoprofile='./'
if not keyword_set(jpgpath) then jpgpath='./'
b0 = 0
b1 = 1
b2 = 2
b3 = 3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in profile data

read_radprof,id,bands[b0], w1tot_a, w1tot_mag, w1tot_mag_e, $
	w1ann_mu, w1ann_mu_e, $
	w1ra_cen, w1dec_cen, w1semimajor, w1semiminor, w1pa, w1scale, $
	w1tf_mag, w1tf_mag_e, w1af_mag, w1af_mag_e, w1mu_bg, w1mu_bg_e, $
	w1skyradius_in, w1skyradius_out, w1tf_a, w1symag, w1symag_e, w1syma, $
	pathtoprofile=pathtoprofile
read_radprof,id,bands[b1], w2tot_a, w2tot_mag, w2tot_mag_e, $
	w2ann_mu, w2ann_mu_e, $
	w2ra_cen, w2dec_cen, w2semimajor, w2semiminor, w2pa, w2scale, $
	w2tf_mag, w2tf_mag_e, w2af_mag, w2af_mag_e, w2mu_bg, w2mu_bg_e, $
	w2skyradius_in, w2skyradius_out, w2tf_a, w2symag, w2symag_e, w2syma, $
	pathtoprofile=pathtoprofile
read_radprof,id,bands[b2], w3tot_a, w3tot_mag, w3tot_mag_e, $
	w3ann_mu, w3ann_mu_e, $
	w3ra_cen, w3dec_cen, w3semimajor, w3semiminor, w3pa, w3scale, $
	w3tf_mag, w3tf_mag_e, w3af_mag, w3af_mag_e, w3mu_bg, w3mu_bg_e, $
	w3skyradius_in, w3skyradius_out, w3tf_a, w3symag, w3symag_e, w3syma, $
	pathtoprofile=pathtoprofile
read_radprof,id,bands[b3], w4tot_a, w4tot_mag, w4tot_mag_e, $
	w4ann_mu, w4ann_mu_e, $
	w4ra_cen, w4dec_cen, w4semimajor, w4semiminor, w4pa, w4scale, $
	w4tf_mag, w4tf_mag_e, w4af_mag, w4af_mag_e, w4mu_bg, w4mu_bg_e, $
	w4skyradius_in, w4skyradius_out, w4tf_a, w4symag, w4symag_e, w4syma, $
	pathtoprofile=pathtoprofile
; check inputs
if w1tot_a[0] eq -1 then now1band = (1 eq 1) else now1band = (1 eq 0)
if w2tot_a[0] eq -1 then now2band = (1 eq 1) else now2band = (1 eq 0)
if w3tot_a[0] eq -1 then now3band = (1 eq 1) else now3band = (1 eq 0)
if w4tot_a[0] eq -1 then now4band = (1 eq 1) else now4band = (1 eq 0)
if now1band then begin
	print,'No primary image photometry for '+id+', returning.'
	return
endif
if now1band and now2band and now3band and now4band then return

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set some useful values to
; set up plots regardless of band

 a=w1tot_a
 tot_mag=w1tot_mag
 ann_mu=w1ann_mu 
 ra_cen=w1ra_cen[0]
 dec_cen=w1dec_cen[0]
 semimajor=w1semimajor[0]
 semiminor=w1semiminor[0]
 pa=w1pa[0]>0.
 r_ap=semimajor
 r_ski=w1skyradius_in
 r_sko=w1skyradius_out

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in headers

if keyword_set(intfile) and file_exist(intfile) then begin
 w1hdr=headfits(intfile,ext=0)
 extast,w1hdr,w1astr
 AD2XY, w1ra_cen[0] ,w1dec_cen[0], w1astr, nx0, ny0
 sz = [sxpar(w1hdr,'NAXIS1'),sxpar(w1hdr,'NAXIS2')]
 getrot,w1hdr,w1rot,w1cdelt
 as_pix = abs(w1cdelt[0])*3600.
endif else begin
	print,'No image data for '+id+', returning.'
	return
endelse

 if keyword_set(verbose) then $
	print,'wise ',form='($,a)'
;;;;;;;;;;;;;;;;;;;;;;;;;;;

 w1_bmap_avail = (1 eq 0)
 w2_bmap_avail = (1 eq 0)
 w3_bmap_avail = (1 eq 0)
 w4_bmap_avail = (1 eq 0)
 
 w1_mask_avail = (1 eq 0)
 w2_mask_avail = (1 eq 0)
 w3_mask_avail = (1 eq 0)
 w4_mask_avail = (1 eq 0)
 
 if keyword_set(jpgpath) and keyword_set(auxbase) then begin
   rute = strmid(jpgpath,0,strpos(jpgpath,'/jpg'))

   bmap_w1 = auxbase+'_w1_badmap.fits.gz'
   bmap_w2 = auxbase+'_w2_badmap.fits.gz'
   bmap_w3 = auxbase+'_w3_badmap.fits.gz'
   bmap_w4 = auxbase+'_w4_badmap.fits.gz'
   
   if file_exist(bmap_w1) then w1_bmap_avail = (1 eq 1)
   if file_exist(bmap_w2) then w2_bmap_avail = (1 eq 1)
   if file_exist(bmap_w3) then w3_bmap_avail = (1 eq 1)
   if file_exist(bmap_w4) then w4_bmap_avail = (1 eq 1)
   
   rute = strmid(rute,0,strpos(rute,'/wise'))
   mask_w1 = rute+'/aux/'+id+'_wise_mask_w1.fits.gz'
   mask_w2 = rute+'/aux/'+id+'_wise_mask_w2.fits.gz'
   mask_w3 = rute+'/aux/'+id+'_wise_mask_w3.fits.gz'
   mask_w4 = rute+'/aux/'+id+'_wise_mask_w4.fits.gz'
   
   if file_exist(mask_w1) then w1_mask_avail = (1 eq 1)
   if file_exist(mask_w2) then w2_mask_avail = (1 eq 1)
   if file_exist(mask_w3) then w3_mask_avail = (1 eq 1)
   if file_exist(mask_w4) then w4_mask_avail = (1 eq 1)
   
 endif
 
 vertices = ''
 background = ''
 if keyword_set(auxbase) then begin
   vertices = auxbase+'_stv_background_vertices.dat'
   background = auxbase+'_stv_background_output.dat'
 endif
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;

 ; read in jpg images

xdjpgfile =jpgpath+'/'+id+'_w123.jpg'

if file_test(xdjpgfile) then begin
 read_jpeg,xdjpgfile,xdjpg,/true
 xdpic=1

 ;mask the 3color image

 if keyword_set(maskimgfile) then begin
  if keyword_set(verbose) then $
	 print,'masking ... ',form='($,a)'
  mskimg = glga_getmask(maskimgfile,sz,w1astr,as_pix)
  maskidx = where(mskimg ge 1, nmaskidx)
  if nmaskidx gt 0 then begin
	xdjpgr = xdjpg[0, *, *]
	xdjpgg = xdjpg[1, *, *]
	xdjpgb = xdjpg[2, *, *]
	xdjpgr[maskidx] = 244
	xdjpgg[maskidx] = 164
	xdjpgb[maskidx] = 96
	xdjpg[0, *, *] = xdjpgr
	xdjpg[1, *, *] = xdjpgg 
	xdjpg[2, *, *] = xdjpgb
  endif
  delvarx,mskimg,maskidx,/free
 endif

endif else xdpic=0

w1djpgfile =jpgpath+'/'+id+'_w1.jpg'  ;################

if file_test(w1djpgfile) then begin
 read_jpeg,w1djpgfile,w1djpg,/true
 sz=size(w1djpg,/dim)
  w1dimg=bytarr(sz[1],sz[2])
  w1dimg[*,*]=w1djpg[0,*,*]

 if w1_bmap_avail then begin 
      bmap_image =  mrdfits(bmap_w1,0,/silent)
      b_maskidx = where(bmap_image ge 1, nmaskidx)
      if nmaskidx gt 0 then begin
      
	xdjpgr = xdjpg[0, *, *]
	xdjpgg = xdjpg[1, *, *]
	xdjpgb = xdjpg[2, *, *]
  	xdjpgr[b_maskidx] = 244
  	xdjpgg[b_maskidx] = 164
  	xdjpgb[b_maskidx] = 96
	xdjpg[0, *, *] = xdjpgr
	xdjpg[1, *, *] = xdjpgg 
	xdjpg[2, *, *] = xdjpgb        
      
	s = size(bmap_image)
	ncol = s(1)
	w1_col = b_maskidx mod ncol
	w1_row = b_maskidx / ncol
      endif
 endif

 w1dpic=1
endif else w1dpic=0


w2djpgfile =jpgpath+'/'+id+'_w2.jpg'  ;################

if file_test(w2djpgfile) then begin
 read_jpeg,w2djpgfile,w2djpg,/true
 sz=size(w2djpg,/dim)
  w2dimg=bytarr(sz[1],sz[2])
  w2dimg[*,*]=w2djpg[0,*,*]

 if w2_bmap_avail then begin 
      bmap_image =  mrdfits(bmap_w2,0,/silent)
      b_maskidx = where(bmap_image ge 1, nmaskidx)
      if nmaskidx gt 0 then begin
      
	xdjpgr = xdjpg[0, *, *]
	xdjpgg = xdjpg[1, *, *]
	xdjpgb = xdjpg[2, *, *]
  	xdjpgr[b_maskidx] = 244
  	xdjpgg[b_maskidx] = 164
  	xdjpgb[b_maskidx] = 96
	xdjpg[0, *, *] = xdjpgr
	xdjpg[1, *, *] = xdjpgg 
	xdjpg[2, *, *] = xdjpgb        
      
	s = size(bmap_image)
	ncol = s(1)
	w2_col = b_maskidx mod ncol
	w2_row = b_maskidx / ncol
      endif
 endif

 w2dpic=1
endif else w2dpic=0

w3djpgfile =jpgpath+'/'+id+'_w3.jpg'  ;################

if file_test(w3djpgfile) then begin
 read_jpeg,w3djpgfile,w3djpg,/true
 sz=size(w3djpg,/dim)
  w3dimg=bytarr(sz[1],sz[2])
  w3dimg[*,*]=w3djpg[0,*,*]

 if w3_bmap_avail then begin 
      bmap_image =  mrdfits(bmap_w3,0,/silent)
      b_maskidx = where(bmap_image ge 1, nmaskidx)
      if nmaskidx gt 0 then begin
      
	xdjpgr = xdjpg[0, *, *]
	xdjpgg = xdjpg[1, *, *]
	xdjpgb = xdjpg[2, *, *]
  	xdjpgr[b_maskidx] = 244
  	xdjpgg[b_maskidx] = 164
  	xdjpgb[b_maskidx] = 96
	xdjpg[0, *, *] = xdjpgr
	xdjpg[1, *, *] = xdjpgg 
	xdjpg[2, *, *] = xdjpgb        
      
	s = size(bmap_image)
	ncol = s(1)
	w3_col = b_maskidx mod ncol
	w3_row = b_maskidx / ncol
      endif
 endif

 w3dpic=1
endif else w3dpic=0

;;;;;;;;;;;;;;;;;;;;;;
; start profile plots

if keyword_set(verbose) then $
	print,'profiles ',form='($,a)'

filename=outpath+'/'+id+'_wise_profile.ps'
ms_ps_start, filename=filename, xsize=10.5,ysize=8.,/inch,$
          /color,bit=8, xoffset=0.15, yoffset = 10.75,/landscape

!p.multi=[0,2,3]
!p.charsize=2

setplotcolors
psize=0.75

cs = [!green, !orange, !red, !brown]

;growth curve
if not now1band then begin
	bad=where(w1tot_mag lt 0., nbad)
	if nbad gt 0 then w1tot_mag[bad] = !values.f_nan
endif else w1tot_mag = !values.f_nan
if not now2band then begin
	bad=where(w2tot_mag lt 0., nbad)
	if nbad gt 0 then w2tot_mag[bad] = !values.f_nan
endif else w2tot_mag = !values.f_nan
if not now3band then begin
	bad=where(w3tot_mag lt 0., nbad)
	if nbad gt 0 then w3tot_mag[bad] = !values.f_nan
endif else w3tot_mag = !values.f_nan
if not now4band then begin
	bad=where(w4tot_mag lt 0., nbad)
	if nbad gt 0 then w4tot_mag[bad] = !values.f_nan
endif else w4tot_mag = !values.f_nan

; nminmax=minmax([w1tot_mag,w2tot_mag,w3tot_mag,w4tot_mag], /nan)
nminmax=minmax([w1tot_mag,w2tot_mag], /nan)
min=nminmax[0]
max=nminmax[1]

ydel = max-min
yr0 = max+ydel*0.1
yr1 = min-ydel*0.2

xdel = r_sko[0]/60.
xr0 = min(a/60.) - xdel*0.02
xr1 = r_sko[0]/60. + xdel*0.02

; growth curve linear x-axis
plot,[a/60],[tot_mag],yr=[yr0, yr1],/ys,/nodata,$
    xr = [xr0, xr1], /xs,$
    xtitle='R!da!n [arcmin]',ytitle='Growth Curve [AB mag]'
;
; make sure we are including the sky annulus in the plot
if a[n_elements(a)-1] gt r_ski then begin
	oplot,[r_ski,r_sko]/60.,[min-ydel*0.1,min-ydel*0.1], color=!blue
	xyouts,r_ski/60.,min-ydel*0.12,'SKY ANN', $
		color=!blue, charsize=0.5
endif

if not now1band then begin
 oploterror, [a/60.],[w1tot_mag],[w1tot_mag_e*0],[w1tot_mag_e],$
  color=cs[b0],errcolor=cs[b0],psym=-sym(1, psize = psize), /nohat
 oplot,!x.crange,[w1symag,w1symag],linesty=1,color=cs[b0]
 oplot,[w1syma,w1syma]/60.,[yr0,w1symag],linesty=2,color=cs[b0],thick=4
endif

if not now2band then begin
 oploterror, [a/60.],[w2tot_mag],[w2tot_mag_e*0],[w2tot_mag_e],$
  color=cs[b1],errcolor=cs[b1],psym=-sym(1, psize = psize), /nohat
 oplot,!x.crange,[w2symag,w2symag],linesty=1,color=cs[b1]
 oplot,[w2syma,w2syma]/60.,[yr0,w2symag],linesty=2,color=cs[b1],thick=4
endif

if not now3band then begin
 oploterror, [a/60.],[w3tot_mag],[w3tot_mag_e*0],[w3tot_mag_e],$
  color=cs[b2],errcolor=cs[b2],psym=-sym(1, psize = psize), /nohat
 oplot,!x.crange,[w3symag,w3symag],linesty=1,color=cs[b2]
 oplot,[w3syma,w3syma]/60.,[yr0,w3symag],linesty=2,color=cs[b2],thick=4
endif

; if not now4band then begin
;  oploterror, [a/60.],[w4tot_mag],[w4tot_mag_e*0],[w4tot_mag_e],$
;   color=cs[b3],errcolor=cs[b3],psym=-sym(1, psize = psize), /nohat
;  oplot,!x.crange,[w4symag,w4symag],linesty=1,color=cs[b3]
;  oplot,[w4syma,w4syma]/60.,[yr0,w4symag],linesty=2,color=cs[b3],thick=4
; endif

 oplot,[r_ap,r_ap]/60.,[yr0,yr1],thick=4,color=!red
 xx = r_ap/60. - abs(!x.crange[1]-!x.crange[0])*0.09
 yy = !y.crange[0] - abs(!y.crange[0]-!y.crange[1])*0.05
 xyouts,xx,yy,'APR',color=!red,charsiz=0.8

; growth curve log x-axis
plot,[a/60],[tot_mag],yr=[yr0, yr1],/ys,/nodata,$
    xr = [0.08, xr1], /xs, /xlog, $
    xtitle='R!da!n [arcmin]',ytitle='Growth Curve [AB mag]'

xyouts,0.1,yr1+ydel*0.14,'ASY',color=!blue,charsiz=0.8

if not now1band then begin
 oploterror, [a/60.],[w1tot_mag],[w1tot_mag_e*0],[w1tot_mag_e],$
  color=cs[b0],errcolor=cs[b0],psym=-sym(1, psize = psize), /nohat
 oplot,[0.08,1.e9],[w1symag,w1symag],linesty=1,color=cs[b0]
endif

if not now2band then begin
 oploterror, [a/60.],[w2tot_mag],[w2tot_mag_e*0],[w2tot_mag_e],$
  color=cs[b1],errcolor=cs[b1],psym=-sym(1, psize = psize), /nohat
 oplot,[0.08,1.e9],[w2symag,w2symag],linesty=1,color=cs[b1]
endif

if not now3band then begin
 oploterror, [a/60.],[w3tot_mag],[w3tot_mag_e*0],[w3tot_mag_e],$
  color=cs[b2],errcolor=cs[b2],psym=-sym(1, psize = psize), /nohat
 oplot,[0.08,1.e9],[w3symag,w3symag],linesty=1,color=cs[b2]
endif

; if not now4band then begin
;  oploterror, [a/60.],[w4tot_mag],[w4tot_mag_e*0],[w4tot_mag_e],$
;   color=cs[b3],errcolor=cs[b3],psym=-sym(1, psize = psize), /nohat
;  oplot,[0.08,1.e9],[w4symag,w4symag],linesty=1,color=cs[b3]
; endif

 oplot,[r_ap,r_ap]/60.,[yr0,yr1],thick=4,color=!red

;annular surface brightness
if not now1band then begin
	bad=where(w1ann_mu lt 0., nbad)
	if nbad gt 0 then w1ann_mu[bad] = !values.f_nan
endif else w1ann_mu = !values.f_nan
if not now2band then begin
	bad=where(w2ann_mu lt 0., nbad)
	if nbad gt 0 then w2ann_mu[bad] = !values.f_nan
endif else w2ann_mu = !values.f_nan
if not now3band then begin
	bad=where(w3ann_mu lt 0., nbad)
	if nbad gt 0 then w3ann_mu[bad] = !values.f_nan
endif else w3ann_mu = !values.f_nan
if not now4band then begin
	bad=where(w4ann_mu lt 0., nbad)
	if nbad gt 0 then w4ann_mu[bad] = !values.f_nan
endif else w4ann_mu = !values.f_nan

nminmax=minmax([w1ann_mu,w2ann_mu,w3ann_mu,w4ann_mu], /nan)
min=nminmax[0]
max=nminmax[1]
; annular surface brightness linear x-axis
plot,[a/60],[ann_mu],yr=[max+0.5, min-0.5],/ys,/nodata,$
    xr = [xr0, xr1], /xs,$
    xtitle='R!da!n [arcmin]',$
    ytitle=greek('mu',/append)+' [AB mag/arcsec!u2!n]'

if not now1band then $
 oploterror, [a/60.],[w1ann_mu],[w1ann_mu_e*0],[w1ann_mu_e],$
  color=cs[b0],errcolor=cs[b0],psym=-sym(1, psize = psize), /nohat

if not now2band then $
 oploterror, [a/60.],[w2ann_mu],[w2ann_mu_e*0],[w2ann_mu_e],$
  color=cs[b1],errcolor=cs[b1],psym=-sym(1, psize = psize), /nohat

if not now3band then $
 oploterror, [a/60.],[w3ann_mu],[w3ann_mu_e*0],[w3ann_mu_e],$
  color=cs[b2],errcolor=cs[b2],psym=-sym(1, psize = psize), /nohat

; if not now4band then $
;  oploterror, [a/60.],[w4ann_mu],[w4ann_mu_e*0],[w4ann_mu_e],$
;   color=cs[b3],errcolor=cs[b3],psym=-sym(1, psize = psize), /nohat

 oplot,[r_ap,r_ap]/60.,[max+0.5,min-0.5],thick=4,color=!red
 xx = r_ap/60. - abs(!x.crange[1]-!x.crange[0])*0.09
 yy = !y.crange[1] + abs(!y.crange[0]-!y.crange[1])*0.1
 xyouts,xx,yy,'APR',color=!red,charsiz=0.8

legend,bands[[b0,b1,b2,b3]],textcolors=cs[[b0,b1,b2,b3]],/bot,/left,box=0, $
	charsize=1.0

; annular surface brightness log x-axis
plot,[a/60],[ann_mu],yr=[max+0.5, min-0.5],/ys,/nodata,$
    xr = [0.08, xr1], /xs, /xlog, $
    xtitle='R!da!n [arcmin]',$
    ytitle=greek('mu',/append)+' [AB mag/arcsec!u2!n]'

if not now1band then $
 oploterror, [a/60.],[w1ann_mu],[w1ann_mu_e*0],[w1ann_mu_e],$
  color=cs[b0],errcolor=cs[b0],psym=-sym(1, psize = psize), /nohat

if not now2band then $
 oploterror, [a/60.],[w2ann_mu],[w2ann_mu_e*0],[w2ann_mu_e],$
  color=cs[b1],errcolor=cs[b1],psym=-sym(1, psize = psize), /nohat

if not now3band then $
 oploterror, [a/60.],[w3ann_mu],[w3ann_mu_e*0],[w3ann_mu_e],$
  color=cs[b2],errcolor=cs[b2],psym=-sym(1, psize = psize), /nohat

; if not now4band then $
;  oploterror, [a/60.],[w4ann_mu],[w4ann_mu_e*0],[w4ann_mu_e],$
;   color=cs[b3],errcolor=cs[b3],psym=-sym(1, psize = psize), /nohat

 oplot,[r_ap,r_ap]/60.,[max+0.5,min-0.5],thick=4,color=!red

legend,bands[[b0,b1,b2,b3]],textcolors=cs[[b0,b1,b2,b3]],/bot,/left,box=0, $
	charsize=1.0

; text

!p.charsize=1.0

plot,[0,100],[0,100],xs=4,ys=4,/nodata, pos = [0.085, 0.01, 0.98, 0.33]

xyouts,0,80,id,charsize=1.5
xyouts,0,60,'R.A.  [J2K]:  '+strn(ra_cen[0],format='(f12.6)')
xyouts,0,50,'DEC [J2K]:  '+strn(dec_cen[0],format='(f12.6)')
xyouts,0,40,'SEMIMAJOR [arcmin]:  '+strn(semimajor/60,format='(f6.2)')
xyouts,0,30,'RATIO (a/b):  '+strn(semimajor/semiminor,format='(f6.2)')
xyouts,0,20,'P.A. [deg]:  '+strn(pa[0],format='(f6.2)')
xyouts,0,10,'TYP:  '+typ

euler, ra_cen, dec_cen, l, b, 1
ebv=dust_getval(l,b,/interp)
a_w1  = glga_getextin(ebv,'w1',yuan13=yuan13)
a_w2  = glga_getextin(ebv,'w2',yuan13=yuan13)
a_w3  = glga_getextin(ebv,'w3',yuan13=yuan13)
a_w4  = glga_getextin(ebv,'w4',yuan13=yuan13)

a_0 = a_w1
a_1 = a_w2
a_2 = a_w3
a_3 = a_w4

xyouts,28,80,'APR:  m('+bands[b0]+'!do!n) = ' + $
 strn(w1af_mag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(w1af_mag_e>0.01,format='(F4.2)')

xyouts,28,70,'APR:  m('+bands[b1]+'!do!n) = ' + $
 strn(w2af_mag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(w2af_mag_e>0.01,format='(F4.2)')

xyouts,28,60,'APR:  m('+bands[b2]+'!do!n) = ' + $
 strn(w3af_mag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(w3af_mag_e>0.01,format='(F4.2)')

xyouts,28,50,'APR:  m('+bands[b3]+'!do!n) = ' + $
 strn(w4af_mag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(w4af_mag_e>0.01,format='(F4.2)')

xyouts,28,40,'ASY:  m('+bands[b0]+'!do!n) = ' + $
 strn(w1symag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(w1symag_e>0.01,format='(F4.2)') + $
 ' R!da!n = '+strn(w1syma/60.,format='(f8.2)')+"'"

xyouts,28,30,'ASY:  m('+bands[b1]+'!do!n) = ' + $
 strn(w2symag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(w2symag_e>0.01,format='(F4.2)') + $
 ' R!da!n = '+strn(w2syma/60.,format='(f8.2)')+"'"

xyouts,28,20,'ASY:  m('+bands[b2]+'!do!n) = ' + $
 strn(w3symag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(w3symag_e>0.01,format='(F4.2)') + $
 ' R!da!n = '+strn(w3syma/60.,format='(f8.2)')+"'"

xyouts,28,10,'ASY:  m('+bands[b3]+'!do!n) = ' + $
 strn(w4symag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(w4symag_e>0.01,format='(F4.2)') + $
 ' R!da!n = '+strn(w4syma/60.,format='(f8.2)')+"'"

uvcolor=w1af_mag-w3af_mag>(-99.)<99.
err=sqrt(w1af_mag_e^2 + w3af_mag_e^2)>0.01<9.
xyouts,75,80,'APR:  '+bands[b0]+'-'+bands[b2]+' = '+$
 strn(uvcolor,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(err,format='(F4.2)')

uvcolor=w1symag-w3symag>(-99.)<99.
err=sqrt(w1symag_e^2 + w3symag_e^2)>0.01<9.
xyouts,75,70,'ASY:  '+bands[b0]+'-'+bands[b2]+' = '+$
 strn(uvcolor,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(err,format='(F4.2)')

uvcolor=w2af_mag-w4af_mag>(-99.)<99.
err=sqrt(w2af_mag_e^2 + w4af_mag_e^2)>0.01<9.
xyouts,75,60,'APR:  '+bands[b1]+'-'+bands[b3]+' = '+$
 strn(uvcolor,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(err,format='(F4.2)')

uvcolor=w2symag-w4symag>(-99.)<99.
err=sqrt(w2symag_e^2 + w4symag_e^2)>0.01<9.
xyouts,75,50,'ASY:  '+bands[b1]+'-'+bands[b3]+' = '+$
 strn(uvcolor,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(err,format='(F4.2)')

;xyouts,75,40,greek('mu',/append)+' BG '+bands[b1]+':  ' $
; +strn(nmu_bg,format='(f5.2)')$
; +'  ('+greek("plus_minus", /append)+strn(nmu_bg_e,format='(F4.2)')+')'$
; ,charsize=chrsz

xyouts,75,40,'GAL E(B-V):  '+strn(ebv,format='(f5.3)')
xyouts,75,30,'GAL A!d'+bands[b0]+'!n:  '+strn(a_w1,format='(f5.3)'), $
	charsize=chrsz
xyouts,75,20,'GAL A!d'+bands[b1]+'!n:  '+strn(a_w2,format='(f5.3)'), $
	charsize=chrsz

xyouts,0.97,0.077,systime(),charsize=.7,/norm,ori=90.

; read in qa
qafile = strmid(intfile,0,strpos(intfile,'_')+1) +'qa.txt'
qa = glga_read_qa_stat(qafile,stat=qastat)

if qastat eq 1 then begin
	complete = qa.complete
	user = strtrim(qa.user_name,2)
	machine = strtrim(qa.machine_name,2)
	note = strtrim(qa.note,2)
endif else begin
	linfo = get_login_info()
	complete = 0
	user = strtrim(linfo.user_name)
	machine = strtrim(linfo.machine_name)
	note = 'MEAS'
endelse
qastr = 'QA: '+strn(complete)+' '+user+'@'+machine+' '+note
xq0 = 75 - fix(strlen(qastr)*0.3)
xyouts,xq0,10,qastr

ms_ps_end

;;;;;;;;;;;;;;;;;;;;;;
; convert to pdf & jpg

p=outpath+'/'
name=id+'_wise_profile'
spawn,'convert +matte -density 196 -resample 72 -rotate -90 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'rm '+p+name+'.ps'

;;;;;;;;;;;;;;;;;;;;;;
; start image plots

if keyword_set(verbose) then $
	print,'images ',form='($,a)'

filename=outpath+'/'+id+'_wise_images.ps'
ms_ps_start, filename=filename, xsize=10.5,ysize=8.,/inch,$
          /color,bit=8, xoffset=0.15, yoffset = 11,/landscape

!p.charsize=1.0
chrsz=0.75

!p.multi=[0,2,2]

minimsz = 30. / w1scale	; 30 arcsec is minimum image size

; w1

if w1dpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(w1dimg,/dim)
 delta = w1skyradius_out[0]/60. * 1.01 > 0.5
 scale = w1scale[0]/60.0 ; arcmins
 ximgrng = ([nx0,nx0-sz[0]]) * scale
 yimgrng = ([0-ny0, sz[1]-ny0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, w1dimg, /preserve, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='w1'

 setplotcolors
 
 if w1_bmap_avail then begin 
    col = (nx0-w1_col)*scale
    row = (w1_row-ny0)*scale
    a_bol = (1 eq 0)
    b_bol = (1 eq 0)
    indx_tmp = where(col gt xpltrng[1] and col lt xpltrng[0], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    a_bol = (1 eq 1)
    endif
    indx_tmp = where(row gt ypltrng[0] and row lt ypltrng[1], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    b_bol = (1 eq 1)
    endif
    
    if a_bol and b_bol then plots, col, row, color=!red, psym=3
 endif
 
 ela=[w1semimajor/60.0,w1semiminor/60.0,0,0,float(90.-w1pa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=w1semiminor/w1semimajor

 eli=[w1skyradius_out/60.0,ratio*w1skyradius_out/60.0,0,0,float(90.-w1pa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[w1skyradius_in/60.0,ratio*w1skyradius_in/60.0,0,0,float(90.-w1pa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1

;###########


 if file_test(vertices) then begin
 
   polygons = read_back_roi(vertices)
   xvertices = polygons[*,0]
   yvertices = polygons[*,1] 
   split = where(xvertices eq -1)

   for i = 0, n_elements(split)-2 do begin
     XV = xvertices[split[i]+1: split[i+1]-1] 
     YV = yvertices[split[i]+1: split[i+1]-1]
     plots, (nx0-XV)*scale,(YV-ny0)*scale, color=!red, linestyle=2, thick=3, /data
   endfor 
 endif
;###########     
   
   
endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No w1 Image',charsize=1

endelse


; w2

if w2dpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(w2dimg,/dim)
 delta = w2skyradius_out[0]/60. * 1.01 > 0.5
 scale = w2scale[0]/60.0 ; arcmins
 ximgrng = ([nx0,nx0-sz[0]]) * scale
 yimgrng = ([0-ny0, sz[1]-ny0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, w2dimg, /preserve, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='w2'

 setplotcolors

 if w2_bmap_avail then begin 
    col = (nx0-w2_col)*scale
    row = (w2_row-ny0)*scale
    a_bol = (1 eq 0)
    b_bol = (1 eq 0)
    indx_tmp = where(col gt xpltrng[1] and col lt xpltrng[0], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    a_bol = (1 eq 1)
    endif
    indx_tmp = where(row gt ypltrng[0] and row lt ypltrng[1], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    b_bol = (1 eq 1)
    endif
    
    if a_bol and b_bol then plots, col, row, color=!red, psym=3
 endif
 


 ela=[w2semimajor/60.0,w2semiminor/60.0,0,0,float(90.-w2pa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=w2semiminor/w2semimajor

 eli=[w2skyradius_out/60.0,ratio*w2skyradius_out/60.0,0,0,float(90.-w2pa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[w2skyradius_in/60.0,ratio*w2skyradius_in/60.0,0,0,float(90.-w2pa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1
   
;###########


 if file_test(vertices) then begin
 
   polygons = read_back_roi(vertices)
   xvertices = polygons[*,0]
   yvertices = polygons[*,1] 
   split = where(xvertices eq -1)

   for i = 0, n_elements(split)-2 do begin
     XV = xvertices[split[i]+1: split[i+1]-1] 
     YV = yvertices[split[i]+1: split[i+1]-1]
     plots, (nx0-XV)*scale,(YV-ny0)*scale, color=!red, linestyle=2, thick=3, /data
   endfor 
 endif
;###########   
   
endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No w2 Image',charsize=1

endelse


; w3

if w3dpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(w3dimg,/dim)
 delta = w3skyradius_out[0]/60. * 1.01 > 0.5
 scale = w3scale[0]/60.0 ; arcmins
 ximgrng = ([nx0,nx0-sz[0]]) * scale
 yimgrng = ([0-ny0, sz[1]-ny0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, w3dimg, /preserve, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='w3'

 setplotcolors

 if w3_bmap_avail then begin 
    col = (nx0-w3_col)*scale
    row = (w3_row-ny0)*scale
    a_bol = (1 eq 0)
    b_bol = (1 eq 0)
    indx_tmp = where(col gt xpltrng[1] and col lt xpltrng[0], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    a_bol = (1 eq 1)
    endif
    indx_tmp = where(row gt ypltrng[0] and row lt ypltrng[1], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    b_bol = (1 eq 1)
    endif
    
    if a_bol and b_bol then plots, col, row, color=!red, psym=3
 endif
 
 ela=[w1semimajor/60.0,w1semiminor/60.0,0,0,float(90.-w1pa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=w1semiminor/w1semimajor

 eli=[w1skyradius_out/60.0,ratio*w1skyradius_out/60.0,0,0,float(90.-w1pa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[w1skyradius_in/60.0,ratio*w1skyradius_in/60.0,0,0,float(90.-w1pa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1
   
;###########


 if file_test(vertices) then begin
 
   polygons = read_back_roi(vertices)
   xvertices = polygons[*,0]
   yvertices = polygons[*,1] 
   split = where(xvertices eq -1)

   for i = 0, n_elements(split)-2 do begin
     XV = xvertices[split[i]+1: split[i+1]-1] 
     YV = yvertices[split[i]+1: split[i+1]-1]
     plots, (nx0-XV)*scale,(YV-ny0)*scale, color=!red, linestyle=2, thick=3, /data
   endfor 
 endif
;###########   
   
endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No w3 Image',charsize=1

endelse


; COLOR (make use of last derived values)

if xdpic then begin

 loadct,0,/silent

 plotimage, xdjpg, /preserve, $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='Masked Composite'

 setplotcolors

 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3
 plots,ela[2],ela[3],psym=4,symsi=2.,thick=0.5,color=!green,/data

 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!gray, thick=1

 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!gray, thick=1

;###########


 if file_test(vertices) then begin
 
   polygons = read_back_roi(vertices)
   xvertices = polygons[*,0]
   yvertices = polygons[*,1] 
   split = where(xvertices eq -1)

   for i = 0, n_elements(split)-2 do begin
     XV = xvertices[split[i]+1: split[i+1]-1] 
     YV = yvertices[split[i]+1: split[i+1]-1]
     plots, (nx0-XV)*scale,(YV-ny0)*scale, color=!cyan, linestyle=2, thick=3, /data
   endfor 
 endif
;###########
   
 
endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.3],[0.5],'No Composite Image',charsize=1

endelse

if n_elements(xpltrng) eq 2 and n_elements(ypltrng) eq 2 then begin
	xdel = xpltrng[0] - xpltrng[1]
	xx0 = xpltrng[0] + xdel*0.55
	ydel = ypltrng[1] - ypltrng[0]
	yy0 = ypltrng[1] + ydel*0.065
endif else begin
	xx0 = 0.55
	yy0 = 0.065
endelse

xyouts,xx0,yy0,id,charsi=1.5
xyouts,0.96,0.06,systime(),charsize=.7,/norm,ori=90.

!p.multi=0

ms_ps_end

;;;;;;;;;;;;;;;;;;;;;;
; convert to pdf & jpg

p=outpath+'/'
name=id+'_wise_images'
spawn,'convert +matte -density 196 -resample 72 -rotate -90 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'rm '+p+name+'.ps'


;;;;;;;;;;;;;;;;;;;;;;
; start image plots2

if w1_mask_avail or w2_mask_avail or w3_mask_avail or w4_mask_avail then begin

filename=outpath+'/'+id+'_wise_images_masks.ps'
ms_ps_start, filename=filename, xsize=10.5,ysize=8.,/inch,$
          /color,bit=8, xoffset=0.15, yoffset = 11,/landscape

!p.charsize=1.0
chrsz=0.75

!p.multi=[0,2,2]

minimsz = 30. / w1scale	; 30 arcsec is minimum image size

; w1

if w1dpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(w1dimg,/dim)
 delta = w1skyradius_out[0]/60. * 1.01 > 0.5
 scale = w1scale[0]/60.0 ; arcmins
 ximgrng = ([nx0,nx0-sz[0]]) * scale
 yimgrng = ([0-ny0, sz[1]-ny0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, w1dimg, /preserve, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='w1'

 setplotcolors
 
 if w1_bmap_avail then begin 
    col = (nx0-w1_col)*scale
    row = (w1_row-ny0)*scale
    a_bol = (1 eq 0)
    b_bol = (1 eq 0)
    indx_tmp = where(col gt xpltrng[1] and col lt xpltrng[0], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    a_bol = (1 eq 1)
    endif
    indx_tmp = where(row gt ypltrng[0] and row lt ypltrng[1], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    b_bol = (1 eq 1)
    endif
    
    if a_bol and b_bol then plots, col, row, color=!red, psym=3
 endif

 if w1_mask_avail  then begin 
 
      mask_image =  mrdfits(mask_w1,0,/silent)
      maskidx = where(mask_image ge 1, nmaskidx)
      if nmaskidx gt 0  then begin
      
	s = size(mask_image)
	ncol = s(1)
	w1_col_mask = maskidx mod ncol
	w1_row_mask = maskidx / ncol
      endif    
      
    col = (nx0-w1_col_mask)*scale
    row = (w1_row_mask-ny0)*scale
    a_bol = (1 eq 0)
    b_bol = (1 eq 0)
    indx_tmp = where(col gt xpltrng[1] and col lt xpltrng[0], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    a_bol = (1 eq 1)
    endif
    indx_tmp = where(row gt ypltrng[0] and row lt ypltrng[1], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    b_bol = (1 eq 1)
    endif
    
    if a_bol and b_bol then plots, col, row, color=!green, psym=3
 endif 
  

 ela=[w1semimajor/60.0,w1semiminor/60.0,0,0,float(90.-w1pa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=w1semiminor/w1semimajor

 eli=[w1skyradius_out/60.0,ratio*w1skyradius_out/60.0,0,0,float(90.-w1pa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[w1skyradius_in/60.0,ratio*w1skyradius_in/60.0,0,0,float(90.-w1pa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1

   
endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No w1 Image',charsize=1

endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; w2

if w2dpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(w2dimg,/dim)
 delta = w2skyradius_out[0]/60. * 1.01 > 0.5
 scale = w2scale[0]/60.0 ; arcmins
 ximgrng = ([nx0,nx0-sz[0]]) * scale
 yimgrng = ([0-ny0, sz[1]-ny0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, w2dimg, /preserve, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='w2'

 setplotcolors

 if w2_bmap_avail then begin 
    col = (nx0-w2_col)*scale
    row = (w2_row-ny0)*scale
    a_bol = (1 eq 0)
    b_bol = (1 eq 0)
    indx_tmp = where(col gt xpltrng[1] and col lt xpltrng[0], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    a_bol = (1 eq 1)
    endif
    indx_tmp = where(row gt ypltrng[0] and row lt ypltrng[1], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    b_bol = (1 eq 1)
    endif
    
    if a_bol and b_bol then plots, col, row, color=!red, psym=3
 endif
 

 if w2_mask_avail  then begin 
 
      mask_image =  mrdfits(mask_w2,0,/silent)
      maskidx = where(mask_image ge 1, nmaskidx)
      if nmaskidx gt 0  then begin
      
	s = size(mask_image)
	ncol = s(1)
	w2_col_mask = maskidx mod ncol
	w2_row_mask = maskidx / ncol
      endif    
      
    col = (nx0-w2_col_mask)*scale
    row = (w2_row_mask-ny0)*scale
    a_bol = (1 eq 0)
    b_bol = (1 eq 0)
    indx_tmp = where(col gt xpltrng[1] and col lt xpltrng[0], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    a_bol = (1 eq 1)
    endif
    indx_tmp = where(row gt ypltrng[0] and row lt ypltrng[1], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    b_bol = (1 eq 1)
    endif
    
    if a_bol and b_bol then plots, col, row, color=!green, psym=3
 endif 
  

 ela=[w2semimajor/60.0,w2semiminor/60.0,0,0,float(90.-w2pa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=w2semiminor/w2semimajor

 eli=[w2skyradius_out/60.0,ratio*w2skyradius_out/60.0,0,0,float(90.-w2pa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[w2skyradius_in/60.0,ratio*w2skyradius_in/60.0,0,0,float(90.-w2pa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No w2 Image',charsize=1

endelse


; w3

if w3dpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(w3dimg,/dim)
 delta = w3skyradius_out[0]/60. * 1.01 > 0.5
 scale = w3scale[0]/60.0 ; arcmins
 ximgrng = ([nx0,nx0-sz[0]]) * scale
 yimgrng = ([0-ny0, sz[1]-ny0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, w3dimg, /preserve, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='w3'

 setplotcolors

 if w3_bmap_avail then begin 
    col = (nx0-w3_col)*scale
    row = (w3_row-ny0)*scale
    a_bol = (1 eq 0)
    b_bol = (1 eq 0)
    indx_tmp = where(col gt xpltrng[1] and col lt xpltrng[0], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    a_bol = (1 eq 1)
    endif
    indx_tmp = where(row gt ypltrng[0] and row lt ypltrng[1], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    b_bol = (1 eq 1)
    endif
    
    if a_bol and b_bol then plots, col, row, color=!red, psym=3
 endif

 if w3_mask_avail  then begin 
 
      mask_image =  mrdfits(mask_w3,0,/silent)
      maskidx = where(mask_image ge 1, nmaskidx)
      if nmaskidx gt 0  then begin
      
	s = size(mask_image)
	ncol = s(1)
	w3_col_mask = maskidx mod ncol
	w3_row_mask = maskidx / ncol
      endif    
      
    col = (nx0-w3_col_mask)*scale
    row = (w3_row_mask-ny0)*scale
    a_bol = (1 eq 0)
    b_bol = (1 eq 0)
    indx_tmp = where(col gt xpltrng[1] and col lt xpltrng[0], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    a_bol = (1 eq 1)
    endif
    indx_tmp = where(row gt ypltrng[0] and row lt ypltrng[1], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    b_bol = (1 eq 1)
    endif
    
    if a_bol and b_bol then plots, col, row, color=!green, psym=3
 endif 
   
 ela=[w1semimajor/60.0,w1semiminor/60.0,0,0,float(90.-w1pa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=w1semiminor/w1semimajor

 eli=[w1skyradius_out/60.0,ratio*w1skyradius_out/60.0,0,0,float(90.-w1pa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[w1skyradius_in/60.0,ratio*w1skyradius_in/60.0,0,0,float(90.-w1pa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1
   

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No w3 Image',charsize=1

endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


w4djpgfile =jpgpath+'/'+id+'_w4.jpg'

if file_test(w4djpgfile) then begin
 read_jpeg,w4djpgfile,w4djpg,/true
 sz=size(w4djpg,/dim)
  w4dimg=bytarr(sz[1],sz[2])
  w4dimg[*,*]=w4djpg[0,*,*]

 if w4_bmap_avail then begin 
      bmap_image =  mrdfits(bmap_w4,0,/silent)
      b_maskidx = where(bmap_image ge 1, nmaskidx)
      if nmaskidx gt 0 then begin
      
	s = size(bmap_image)
	ncol = s(1)
	w4_col = b_maskidx mod ncol
	w4_row = b_maskidx / ncol
      endif
 endif

w4dpic=1
endif else w4dpic=0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; w4

if w4dpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(w4dimg,/dim)
 delta = w4skyradius_out[0]/60. * 1.01 > 0.5
 scale = w4scale[0]/60.0 ; arcmins
 ximgrng = ([nx0,nx0-sz[0]]) * scale
 yimgrng = ([0-ny0, sz[1]-ny0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, w4dimg, /preserve, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='w4'

 setplotcolors

 if w4_bmap_avail then begin 
    col = (nx0-w4_col)*scale
    row = (w4_row-ny0)*scale
    a_bol = (1 eq 0)
    b_bol = (1 eq 0)
    indx_tmp = where(col gt xpltrng[1] and col lt xpltrng[0], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    a_bol = (1 eq 1)
    endif
    indx_tmp = where(row gt ypltrng[0] and row lt ypltrng[1], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    b_bol = (1 eq 1)
    endif
    
    if a_bol and b_bol then plots, col, row, color=!red, psym=3
 endif

 if w4_mask_avail  then begin 
 
      mask_image =  mrdfits(mask_w4,0,/silent)
      maskidx = where(mask_image ge 1, nmaskidx)
      if nmaskidx gt 0  then begin
      
	s = size(mask_image)
	ncol = s(1)
	w4_col_mask = maskidx mod ncol
	w4_row_mask = maskidx / ncol
      endif    
      
    col = (nx0-w4_col_mask)*scale
    row = (w4_row_mask-ny0)*scale
    a_bol = (1 eq 0)
    b_bol = (1 eq 0)
    indx_tmp = where(col gt xpltrng[1] and col lt xpltrng[0], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    a_bol = (1 eq 1)
    endif
    indx_tmp = where(row gt ypltrng[0] and row lt ypltrng[1], n_tmp)
    if n_tmp gt 0 then begin
    col = col[indx_tmp]
    row = row[indx_tmp]
    b_bol = (1 eq 1)
    endif
    
    if a_bol and b_bol then plots, col, row, color=!green, psym=3
 endif 
   
 ela=[w1semimajor/60.0,w1semiminor/60.0,0,0,float(90.-w1pa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=w1semiminor/w1semimajor

 eli=[w1skyradius_out/60.0,ratio*w1skyradius_out/60.0,0,0,float(90.-w1pa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[w1skyradius_in/60.0,ratio*w1skyradius_in/60.0,0,0,float(90.-w1pa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1
   

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No w4 Image',charsize=1

endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if n_elements(xpltrng) eq 2 and n_elements(ypltrng) eq 2 then begin
	xdel = xpltrng[0] - xpltrng[1]
	xx0 = xpltrng[0] + xdel*0.55
	ydel = ypltrng[1] - ypltrng[0]
	yy0 = ypltrng[1] + ydel*0.065
endif else begin
	xx0 = 0.55
	yy0 = 0.065
endelse

xyouts,xx0,yy0,id,charsi=1.5
xyouts,0.96,0.06,systime(),charsize=.7,/norm,ori=90.

!p.multi=0

ms_ps_end

;;;;;;;;;;;;;;;;;;;;;;
; convert to pdf & jpg

p=outpath+'/'
name=id+'_wise_images_masks'
spawn,'convert +matte -density 196 -resample 72 -rotate -90 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'rm '+p+name+'.ps'

endif else begin
    p=outpath+'/'
    name=id+'_wise_images_masks'
    spawn,  'rm ' + p+name+'.*'
    endelse

if keyword_set(verbose) then $
	print,'Done.'

;;;;;;;;;;;;
; all done

return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function read_back_roi, filename

array = ''
line = ''


OPENR, lun, filename, /GET_LUN

WHILE NOT EOF(lun) DO BEGIN & $
  READF, lun, line & $
  array = [array, line] & $
ENDWHILE


xvertices = [-1]
yvertices = [-1]
i = 0 
while i lt n_elements(array) do begin
  

  first = strmid(array[i], 0, 1)
  while strcmp(first , '#') and i lt n_elements(array) do begin
      i++
      if i eq n_elements(array) then break
      first = strmid(array[i], 0, 1)
      
  endwhile
  
  S = STRSPLIT(array[i],' ', /extract)
  x0 = long(S[0])
  y0 = long(S[-1])

  
  while first ne '#' and i lt n_elements(array) do begin
       S = STRSPLIT(array[i],' ', /extract)
       xvertices = [xvertices, long(S[0])]
       yvertices = [yvertices, long(S[-1])]
       i++
       if i eq n_elements(array) then break
       first = strmid(array[i], 0, 1)
  endwhile
  
  ; to have a closed polygon
  xvertices = [xvertices, x0, -1]
  yvertices = [yvertices, y0, -1]
  

endwhile 

i = 0
while xvertices[i] le 0 do i++

xvertices =  xvertices[i-1:n_elements(xvertices)-1]
yvertices =  yvertices[i-1:n_elements(yvertices)-1]

free_lun,lun

help,  [[xvertices], [yvertices]]
return, [[xvertices], [yvertices]]

end ; program ends here


;;;;;;;;;;;;;;;;
