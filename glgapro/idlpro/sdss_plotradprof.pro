pro sdss_plotradprof, id, pathtoprofile=pathtoprofile, type=type, $
    intfile=intfile, maskimgfile=maskimgfile, auxbase=auxbase,$
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

bands = ['u','g','r','i','z']
nband = n_elements(bands)

if not keyword_set(type) then typ = '-' else typ = strtrim(type,2)
if not keyword_set(pathtoprofile) then pathtoprofile='./'
if not keyword_set(jpgpath) then jpgpath='./'
b0 = 0
b1 = 1
b2 = 2
b3 = 3
b4 = 4

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in profile data

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
read_radprof,id,bands[b2], rtot_a, rtot_mag, rtot_mag_e, rann_mu, rann_mu_e, $
	rra_cen, rdec_cen, rsemimajor, rsemiminor, rpa, rscale, $
	rtf_mag, rtf_mag_e, raf_mag, raf_mag_e, rmu_bg, rmu_bg_e, $
	rskyradius_in, rskyradius_out, rtf_a, rsymag, rsymag_e, rsyma, $
	pathtoprofile=pathtoprofile
read_radprof,id,bands[b3], itot_a, itot_mag, itot_mag_e, iann_mu, iann_mu_e, $
	ira_cen, idec_cen, isemimajor, isemiminor, ipa, iscale, $
	itf_mag, itf_mag_e, iaf_mag, iaf_mag_e, imu_bg, imu_bg_e, $
	iskyradius_in, iskyradius_out, itf_a, isymag, isymag_e, isyma, $
	pathtoprofile=pathtoprofile
read_radprof,id,bands[b4], ztot_a, ztot_mag, ztot_mag_e, zann_mu, zann_mu_e, $
	zra_cen, zdec_cen, zsemimajor, zsemiminor, zpa, zscale, $
	ztf_mag, ztf_mag_e, zaf_mag, zaf_mag_e, zmu_bg, zmu_bg_e, $
	zskyradius_in, zskyradius_out, ztf_a, zsymag, zsymag_e, zsyma, $
	pathtoprofile=pathtoprofile
; check inputs
if utot_a[0] eq -1 then nouband = (1 eq 1) else nouband = (1 eq 0)
if gtot_a[0] eq -1 then nogband = (1 eq 1) else nogband = (1 eq 0)
if rtot_a[0] eq -1 then norband = (1 eq 1) else norband = (1 eq 0)
if itot_a[0] eq -1 then noiband = (1 eq 1) else noiband = (1 eq 0)
if ztot_a[0] eq -1 then nozband = (1 eq 1) else nozband = (1 eq 0)
if norband then begin
	print,'No primary image photometry for '+id+', returning.'
	return
endif
if nogband and norband and noiband then return
; fix missing values
if nouband then begin
	utf_mag=!values.f_nan
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
if norband then begin
	rtf_mag=!values.f_nan
	rtf_mag_e=!values.f_nan
	raf_mag=!values.f_nan
	raf_mag_e=!values.f_nan
	rsyma=!values.f_nan
	rsymag=!values.f_nan
	rsymag_e=!values.f_nan
endif
if noiband then begin
	itf_mag=!values.f_nan
	itf_mag_e=!values.f_nan
	iaf_mag=!values.f_nan
	iaf_mag_e=!values.f_nan
	isyma=!values.f_nan
	isymag=!values.f_nan
	isymag_e=!values.f_nan
endif
if nozband then begin
	ztf_mag=!values.f_nan
	ztf_mag_e=!values.f_nan
	zaf_mag=!values.f_nan
	zaf_mag_e=!values.f_nan
	zsyma=!values.f_nan
	zsymag=!values.f_nan
	zsymag_e=!values.f_nan
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set some useful values to
; set up plots regardless of band

 a=rtot_a
 tot_mag=rtot_mag
 ann_mu=rann_mu 
 ra_cen=rra_cen[0]
 dec_cen=rdec_cen[0]
 semimajor=rsemimajor[0]
 semiminor=rsemiminor[0]
 pa=rpa[0]>0.
 r_ap=semimajor
 r_ski=rskyradius_in
 r_sko=rskyradius_out

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in headers

if keyword_set(intfile) and file_exist(intfile) then begin
 rhdr=headfits(intfile,ext=0)
 extast,rhdr,rastr
 AD2XY, rra_cen[0] ,rdec_cen[0], rastr, rx0, ry0
 sz = [sxpar(rhdr,'NAXIS1'),sxpar(rhdr,'NAXIS2')]
 getrot,rhdr,rrot,rcdelt
 as_pix = abs(rcdelt[0])*3600.
endif else begin
	print,'No image data for '+id+', returning.'
	return
endelse

if keyword_set(verbose) then $
	print,'sdss ',form='($,a)'
;;;;;;;;;;;;;;;;;;;;;;;;;;;
 g_bmap_avail = (1 eq 0)
 r_bmap_avail = (1 eq 0)
 i_bmap_avail = (1 eq 0)
 z_bmap_avail = (1 eq 0)
 
 g_mask_avail = (1 eq 0)
 r_mask_avail = (1 eq 0)
 i_mask_avail = (1 eq 0)
 z_mask_avail = (1 eq 0)
 
 if keyword_set(jpgpath) and keyword_set(auxbase) then begin
   rute = strmid(jpgpath,0,strpos(jpgpath,'/jpg'))
   
   bmap_g = auxbase+'_g_badmap.fits.gz'
   bmap_r = auxbase+'_r_badmap.fits.gz'
   bmap_i = auxbase+'_i_badmap.fits.gz'
   bmap_z = auxbase+'_z_badmap.fits.gz'
   
   if file_exist(bmap_g) then g_bmap_avail = (1 eq 1)
   if file_exist(bmap_r) then r_bmap_avail = (1 eq 1)
   if file_exist(bmap_i) then i_bmap_avail = (1 eq 1)
   if file_exist(bmap_z) then z_bmap_avail = (1 eq 1)
   
   rute = strmid(rute,0,strpos(rute,'/sdss'))
   mask_g = rute+'/aux/'+id+'_sdss_mask_g.fits.gz'
   mask_r = rute+'/aux/'+id+'_sdss_mask_r.fits.gz'
   mask_i = rute+'/aux/'+id+'_sdss_mask_i.fits.gz'
   mask_z = rute+'/aux/'+id+'_sdss_mask_z.fits.gz'
   
   if file_exist(mask_g) then g_mask_avail = (1 eq 1)
   if file_exist(mask_r) then r_mask_avail = (1 eq 1)
   if file_exist(mask_i) then i_mask_avail = (1 eq 1)
   if file_exist(mask_z) then z_mask_avail = (1 eq 1)
   
 endif

 vertices = ''
 background = ''
 if keyword_set(auxbase) then begin
   vertices = auxbase+'_stv_background_vertices.dat'
   background = auxbase+'_stv_background_output.dat'
 endif
  
 
 ; read in jpg images

xdjpgfile =jpgpath+'/'+id+'_gri.jpg'

if file_test(xdjpgfile) then begin
 read_jpeg,xdjpgfile,xdjpg,/true
 xdpic=1

 ;mask the 3color image

 if keyword_set(maskimgfile) then begin
  if keyword_set(verbose) then $
	 print,'masking ... ',form='($,a)'
  mskimg = glga_getmask(maskimgfile,sz,rastr,as_pix)
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

gdjpgfile =jpgpath+'/'+id+'_g.jpg'

if file_test(gdjpgfile) then begin
 read_jpeg,gdjpgfile,gdjpg,/true
 sz=size(gdjpg,/dim)
  gdimg=bytarr(sz[1],sz[2])
  gdimg[*,*]=gdjpg[0,*,*]

 if g_bmap_avail then begin 
      bmap_image =  mrdfits(bmap_g,0,/silent)
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
	g_col = b_maskidx mod ncol
	g_row = b_maskidx / ncol
      endif
 endif

 gdpic=1
endif else gdpic=0

rdjpgfile =jpgpath+'/'+id+'_r.jpg'

if file_test(rdjpgfile) then begin
 read_jpeg,rdjpgfile,rdjpg,/true
 sz=size(rdjpg,/dim)
 rdimg=bytarr(sz[1],sz[2])
 rdimg[*,*]=rdjpg[0,*,*]
 
 if r_bmap_avail then begin 
      bmap_image =  mrdfits(bmap_r,0,/silent)
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
	r_col = b_maskidx mod ncol
	r_row = b_maskidx / ncol
      endif
 endif

 
 if r_mask_avail then begin 
      mask_image =  mrdfits(mask_r,0,/silent)
      maskidx = where(mask_image ge 1, nmaskidx)
      if nmaskidx gt 0  then begin
      
	s = size(mask_image)
	ncol = s(1)
	r_col_mask = maskidx mod ncol
	r_row_mask = maskidx / ncol
      endif
 endif 
 
rdpic=1
endif else rdpic=0

idjpgfile =jpgpath+'/'+id+'_i.jpg'

if file_test(idjpgfile) then begin
 read_jpeg,idjpgfile,idjpg,/true
 sz=size(idjpg,/dim)
 idimg=bytarr(sz[1],sz[2])
 idimg[*,*]=idjpg[0,*,*]
 
  if i_bmap_avail then begin 
      bmap_image =  mrdfits(bmap_i,0,/silent)
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
	i_col = b_maskidx mod ncol
	i_row = b_maskidx / ncol
      endif
 endif

 
 
 idpic=1
endif else idpic=0

;;;;;;;;;;;;;;;;;;;;;;
; start profile plots

if keyword_set(verbose) then $
	print,'profiles ',form='($,a)'

filename=outpath+'/'+id+'_sdss_profile.ps'
ms_ps_start, filename=filename, xsize=10.5,ysize=8.,/inch,$
          /color,bit=8, xoffset=0.15, yoffset = 10.75,/landscape

!p.multi=[0,2,3]
!p.charsize=2

setplotcolors
psize=0.75

cs = [!blue, !green, !orange, !red, !brown]

;growth curve
if not nouband then begin
	bad=where(utot_mag lt 0., nbad)
	if nbad gt 0 then utot_mag[bad] = !values.f_nan
endif else utot_mag = !values.f_nan
if not nogband then begin
	bad=where(gtot_mag lt 0., nbad)
	if nbad gt 0 then gtot_mag[bad] = !values.f_nan
endif else gtot_mag = !values.f_nan
if not norband then begin
	bad=where(rtot_mag lt 0., nbad)
	if nbad gt 0 then rtot_mag[bad] = !values.f_nan
endif else rtot_mag = !values.f_nan
if not noiband then begin
	bad=where(itot_mag lt 0., nbad)
	if nbad gt 0 then itot_mag[bad] = !values.f_nan
endif else itot_mag = !values.f_nan
if not nozband then begin
	bad=where(ztot_mag lt 0., nbad)
	if nbad gt 0 then ztot_mag[bad] = !values.f_nan
endif else ztot_mag = !values.f_nan

nminmax=minmax([utot_mag,gtot_mag,rtot_mag,itot_mag,ztot_mag], /nan)
; nminmax=minmax([gtot_mag,rtot_mag,itot_mag], /nan)
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
    xr = [xr0, xr1], /xs,$
    xtit='R!da!n [arcmin]',ytit='Growth Curve [AB mag]'
;
; make sure we are including the sky annulus in the plot
if a[n_elements(a)-1] gt r_ski then begin
	oplot,[r_ski,r_sko]/60.,[min-ydel*0.1,min-ydel*0.1], color=!blue
	xyouts,r_ski/60.,min-ydel*0.12,'SKY ANN', $
		color=!blue, charsize=0.5
endif

if not nouband then begin
 oploterror, [a/60.],[utot_mag],[utot_mag_e*0],[utot_mag_e],$
  color=cs[b0],errcolor=cs[b0],psym=-sym(1, psize = psize), /nohat
 oplot,!x.crange,[usymag,usymag],linesty=1,color=cs[b0]
 oplot,[usyma,usyma]/60.,[yr0,usymag],linesty=2,color=cs[b0],thick=4
endif

if not nogband then begin
 oploterror, [a/60.],[gtot_mag],[gtot_mag_e*0],[gtot_mag_e],$
  color=cs[b1],errcolor=cs[b1],psym=-sym(1, psize = psize), /nohat
 oplot,!x.crange,[gsymag,gsymag],linesty=1,color=cs[b1]
 oplot,[gsyma,gsyma]/60.,[yr0,gsymag],linesty=2,color=cs[b1],thick=4
endif

if not norband then begin
 oploterror, [a/60.],[rtot_mag],[rtot_mag_e*0],[rtot_mag_e],$
  color=cs[b2],errcolor=cs[b2],psym=-sym(1, psize = psize), /nohat
 oplot,!x.crange,[rsymag,rsymag],linesty=1,color=cs[b2]
 oplot,[rsyma,rsyma]/60.,[yr0,rsymag],linesty=2,color=cs[b2],thick=4
endif

if not noiband then begin
 oploterror, [a/60.],[itot_mag],[itot_mag_e*0],[itot_mag_e],$
  color=cs[b3],errcolor=cs[b3],psym=-sym(1, psize = psize), /nohat
 oplot,!x.crange,[isymag,isymag],linesty=1,color=cs[b3]
 oplot,[isyma,isyma]/60.,[yr0,isymag],linesty=2,color=cs[b3],thick=4
endif

if not nozband then begin
 oploterror, [a/60.],[ztot_mag],[ztot_mag_e*0],[ztot_mag_e],$
  color=cs[b4],errcolor=cs[b4],psym=-sym(1, psize = psize), /nohat
 oplot,!x.crange,[zsymag,zsymag],linesty=1,color=cs[b4]
 oplot,[zsyma,zsyma]/60.,[yr0,zsymag],linesty=2,color=cs[b4],thick=4
endif

oplot,[r_ap,r_ap]/60.,[yr0, yr1],thick=4,color=!red
xx = r_ap/60. - abs(!x.crange[1]-!x.crange[0])*0.09
yy = !y.crange[0] - abs(!y.crange[0]-!y.crange[1])*0.05
xyouts,xx,yy,'APR',color=!red,charsiz=0.8

; growth curve log x-axis
plot,[a/60],[tot_mag],yr=[yr0, yr1],/ys,/nodata,$
    xr = [xlr0, xr1], /xs, /xlog, $
    xtit='R!da!n [arcmin]',ytit='Growth Curve [AB mag]'

xyouts,0.1,yr1+ydel*0.14,'ASY',color=!blue,charsiz=0.8

if not nouband then begin
 oploterror, [a/60.],[utot_mag],[utot_mag_e*0],[utot_mag_e],$
  color=cs[b0],errcolor=cs[b0],psym=-sym(1, psize = psize), /nohat
 oplot,[xlr0,1.e9],[usymag,usymag],linesty=1,color=cs[b0]
endif

if not nogband then begin
 oploterror, [a/60.],[gtot_mag],[gtot_mag_e*0],[gtot_mag_e],$
  color=cs[b1],errcolor=cs[b1],psym=-sym(1, psize = psize), /nohat
 oplot,[xlr0,1.e9],[gsymag,gsymag],linesty=1,color=cs[b1]
endif

if not norband then begin
 oploterror, [a/60.],[rtot_mag],[rtot_mag_e*0],[rtot_mag_e],$
  color=cs[b2],errcolor=cs[b2],psym=-sym(1, psize = psize), /nohat
 oplot,[xlr0,1.e9],[rsymag,rsymag],linesty=1,color=cs[b2]
endif

if not noiband then begin
 oploterror, [a/60.],[itot_mag],[itot_mag_e*0],[itot_mag_e],$
  color=cs[b3],errcolor=cs[b3],psym=-sym(1, psize = psize), /nohat
 oplot,[xlr0,1.e9],[isymag,isymag],linesty=1,color=cs[b3]
endif

if not nozband then begin
 oploterror, [a/60.],[ztot_mag],[ztot_mag_e*0],[ztot_mag_e],$
  color=cs[b4],errcolor=cs[b4],psym=-sym(1, psize = psize), /nohat
 oplot,[xlr0,1.e9],[zsymag,zsymag],linesty=1,color=cs[b4]
endif

 oplot,[r_ap,r_ap]/60.,[yr0, yr1],thick=4,color=!red

;annular surface brightness
if not nouband then begin
	bad=where(uann_mu lt 0., nbad)
	if nbad gt 0 then uann_mu[bad] = !values.f_nan
endif else uann_mu = !values.f_nan
if not nogband then begin
	bad=where(gann_mu lt 0., nbad)
	if nbad gt 0 then gann_mu[bad] = !values.f_nan
endif else gann_mu = !values.f_nan
if not norband then begin
	bad=where(rann_mu lt 0., nbad)
	if nbad gt 0 then rann_mu[bad] = !values.f_nan
endif else rann_mu = !values.f_nan
if not noiband then begin
	bad=where(iann_mu lt 0., nbad)
	if nbad gt 0 then iann_mu[bad] = !values.f_nan
endif else iann_mu = !values.f_nan
if not nozband then begin
	bad=where(zann_mu lt 0., nbad)
	if nbad gt 0 then zann_mu[bad] = !values.f_nan
endif else zann_mu = !values.f_nan

nminmax=minmax([uann_mu,gann_mu,rann_mu,iann_mu,zann_mu], /nan)
min=nminmax[0]
max=nminmax[1]
; annular surface brightness linear x-axis
plot,[a/60],[ann_mu],yr=[max+0.5, min-0.5],/ys,/nodata,$
    xr = [xr0, xr1], /xs,$
    xtit='R!da!n [arcmin]',$
    ytit=greek('mu',/append)+' [AB mag/arcsec!u2!n]'

if not nouband then $
 oploterror, [a/60.],[uann_mu],[uann_mu_e*0],[uann_mu_e],$
  color=cs[b0],errcolor=cs[b0],psym=-sym(1, psize = psize), /nohat

if not nogband then $
 oploterror, [a/60.],[gann_mu],[gann_mu_e*0],[gann_mu_e],$
  color=cs[b1],errcolor=cs[b1],psym=-sym(1, psize = psize), /nohat

if not norband then $
 oploterror, [a/60.],[rann_mu],[rann_mu_e*0],[rann_mu_e],$
  color=cs[b2],errcolor=cs[b2],psym=-sym(1, psize = psize), /nohat

if not noiband then $
 oploterror, [a/60.],[iann_mu],[iann_mu_e*0],[iann_mu_e],$
  color=cs[b3],errcolor=cs[b3],psym=-sym(1, psize = psize), /nohat

if not nozband then $
 oploterror, [a/60.],[zann_mu],[zann_mu_e*0],[zann_mu_e],$
  color=cs[b4],errcolor=cs[b4],psym=-sym(1, psize = psize), /nohat

 oplot,[r_ap,r_ap]/60.,[max+0.5,min-0.5],thick=4,color=!red
 xx = r_ap/60. - abs(!x.crange[1]-!x.crange[0])*0.09
 yy = !y.crange[1] + abs(!y.crange[0]-!y.crange[1])*0.1
 xyouts,xx,yy,'APR',color=!red,charsiz=0.8

legend,bands[[b0,b1,b2,b3,b4]],textcolors=cs[[b0,b1,b2,b3,b4]],/bot,/left, $
	box=0, charsize=1.0

; annular surface brightness log x-axis
plot,[a/60],[ann_mu],yr=[max+0.5, min-0.5],/ys,/nodata,$
    xr = [xlr0, xr1], /xs, /xlog, $
    xtit='R!da!n [arcmin]',$
    ytit=greek('mu',/append)+' [AB mag/arcsec!u2!n]'

if not nouband then $
 oploterror, [a/60.],[uann_mu],[uann_mu_e*0],[uann_mu_e],$
  color=cs[b0],errcolor=cs[b0],psym=-sym(1, psize = psize), /nohat

if not nogband then $
 oploterror, [a/60.],[gann_mu],[gann_mu_e*0],[gann_mu_e],$
  color=cs[b1],errcolor=cs[b1],psym=-sym(1, psize = psize), /nohat

if not norband then $
 oploterror, [a/60.],[rann_mu],[rann_mu_e*0],[rann_mu_e],$
  color=cs[b2],errcolor=cs[b2],psym=-sym(1, psize = psize), /nohat

if not noiband then $
 oploterror, [a/60.],[iann_mu],[iann_mu_e*0],[iann_mu_e],$
  color=cs[b3],errcolor=cs[b3],psym=-sym(1, psize = psize), /nohat

if not nozband then $
 oploterror, [a/60.],[zann_mu],[zann_mu_e*0],[zann_mu_e],$
  color=cs[b4],errcolor=cs[b4],psym=-sym(1, psize = psize), /nohat

 oplot,[r_ap/60,r_ap/60],[max+0.5,min-0.5],thick=4,color=!red

legend,bands[[b0,b1,b2,b3,b4]],textcolors=cs[[b0,b1,b2,b3,b4]],/bot,/left, $
	box=0, charsize=1.0

; text

!p.charsize=1.0

plot,[0,100],[0,100],xs=4,ys=4,/nodata, pos = [0.085, 0.01, 0.98, 0.33]

xyouts,0,80,id,charsize=1.5
xyouts,0,70,'TYP:  '+typ
xyouts,0,60,'R.A.  [J2K]:  '+strn(ra_cen[0],format='(f12.6)')
xyouts,0,50,'DEC [J2K]:  '+strn(dec_cen[0],format='(f12.6)')
xyouts,0,40,'SEMIMAJOR [arcmin]:  '+strn(semimajor/60,format='(f6.2)')
xyouts,0,30,'RATIO (a/b):  '+strn(semimajor/semiminor,format='(f6.2)')
xyouts,0,20,'P.A. [deg]:  '+strn(pa[0],format='(f6.2)')
xyouts,0,10,'Annuli [asec]:  '+strn(annuli_size,format='(f6.2)')

euler, ra_cen, dec_cen, l, b, 1
ebv=dust_getval(l,b,/interp)
a_u  = glga_getextin(ebv,'u',yuan13=yuan13)
a_g  = glga_getextin(ebv,'g',yuan13=yuan13)
a_r  = glga_getextin(ebv,'r',yuan13=yuan13)
a_i  = glga_getextin(ebv,'i',yuan13=yuan13)
a_z  = glga_getextin(ebv,'z',yuan13=yuan13)

xyouts,25,80,'APR:  m('+bands[b0]+') = ' + $
 strn(uaf_mag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(uaf_mag_e>0.01,format='(F4.2)')

xyouts,25,70,'APR:  m('+bands[b1]+') = ' + $
 strn(gaf_mag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(gaf_mag_e>0.01,format='(F4.2)')

xyouts,25,60,'APR:  m('+bands[b2]+') = ' + $
 strn(raf_mag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(raf_mag_e>0.01,format='(F4.2)')

xyouts,25,50,'APR:  m('+bands[b3]+') = ' + $
 strn(iaf_mag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(iaf_mag_e>0.01,format='(F4.2)')

xyouts,25,40,'APR:  m('+bands[b4]+') = ' + $
 strn(zaf_mag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(zaf_mag_e>0.01,format='(F4.2)')

xyouts,50,80,'ASY:  m('+bands[b0]+') = ' + $
 strn(usymag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(usymag_e>0.01,format='(F4.2)') + $
 '  R!da!n = '+strn(usyma/60.,format='(f8.2)')+"'"

xyouts,50,70,'ASY:  m('+bands[b1]+') = ' + $
 strn(gsymag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(gsymag_e>0.01,format='(F4.2)') + $
 '  R!da!n = '+strn(gsyma/60.,format='(f8.2)')+"'"

xyouts,50,60,'ASY:  m('+bands[b2]+') = ' + $
 strn(rsymag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(rsymag_e>0.01,format='(F4.2)') + $
 '  R!da!n = '+strn(rsyma/60.,format='(f8.2)')+"'"

xyouts,50,50,'ASY:  m('+bands[b3]+') = ' + $
 strn(isymag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(isymag_e>0.01,format='(F4.2)') + $
 '  R!da!n = '+strn(isyma/60.,format='(f8.2)')+"'"

xyouts,50,40,'ASY:  m('+bands[b4]+') = ' + $
 strn(zsymag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(zsymag_e>0.01,format='(F4.2)') + $
 '  R!da!n = '+strn(zsyma/60.,format='(f8.2)')+"'"

opcolor=((uaf_mag-a_u)-(gaf_mag-a_g))>(-99.)<99.
err=sqrt(uaf_mag_e^2 + gaf_mag_e^2)>0.01<9.
xyouts,25,30,'APR:  ('+bands[b0]+'-'+bands[b1]+')!d0!n = '+$
 strn(opcolor,format='(f5.2)') + $
 ' '+greek("plus_minus", /append)+' '+strn(err,$
 format='(F4.2)')

opcolor=((gaf_mag-a_g)-(raf_mag-a_r))>(-99.)<99.
err=sqrt(gaf_mag_e^2 + raf_mag_e^2)>0.01<9.
xyouts,25,20,'APR:  ('+bands[b1]+'-'+bands[b2]+')!d0!n = '+$
 strn(opcolor,format='(f5.2)') + $
 ' '+greek("plus_minus", /append)+' '+strn(err,$
 format='(F4.2)')

opcolor=((raf_mag-a_r)-(iaf_mag-a_i))>(-99.)<99.
err=sqrt(raf_mag_e^2 + iaf_mag_e^2)>0.01<9.
xyouts,25,10,'APR:  ('+bands[b2]+'-'+bands[b3]+')!d0!n = '+$
 strn(opcolor,format='(f5.2)') + $
 ' '+greek("plus_minus", /append)+' '+strn(err,$
 format='(F4.2)')

opcolor=((usymag-a_u)-(gsymag-a_g))>(-99.)<99.
err=sqrt(usymag_e^2 + gsymag_e^2)>0.01<9.
xyouts,50,30,'ASY:  ('+bands[b0]+'-'+bands[b1]+')!d0!n = '+$
 strn(opcolor,format='(f5.2)') + $
 ' '+greek("plus_minus", /append)+' '+strn(err,$
 format='(F4.2)')

opcolor=((gsymag-a_g)-(rsymag-a_r))>(-99.)<99.
err=sqrt(gsymag_e^2 + rsymag_e^2)>0.01<9.
xyouts,50,20,'ASY:  ('+bands[b1]+'-'+bands[b2]+')!d0!n = '+$
 strn(opcolor,format='(f5.2)') + $
 ' '+greek("plus_minus", /append)+' '+strn(err,$
 format='(F4.2)')

opcolor=((rsymag-a_r)-(isymag-a_i))>(-99.)<99.
err=sqrt(rsymag_e^2 + isymag_e^2)>0.01<9.
xyouts,50,10,'ASY:  ('+bands[b2]+'-'+bands[b3]+')!d0!n = '+$
 strn(opcolor,format='(f5.2)') + $
 ' '+greek("plus_minus", /append)+' '+strn(err,$
 format='(F4.2)')

xyouts,83,80,'GAL E(B-V):  '+strn(ebv,format='(f5.3)')
xyouts,83,70,'GAL A!d'+bands[b0]+'!n:  '+strn(A_u,format='(f5.3)')
xyouts,83,60,'GAL A!d'+bands[b1]+'!n:  '+strn(A_g,format='(f5.3)')
xyouts,83,50,'GAL A!d'+bands[b2]+'!n:  '+strn(A_r,format='(f5.3)')
xyouts,83,40,'GAL A!d'+bands[b3]+'!n:  '+strn(A_i,format='(f5.3)')
xyouts,83,30,'GAL A!d'+bands[b4]+'!n:  '+strn(A_z,format='(f5.3)')

xyouts,83,20,greek('mu',/append)+'!DBG!N '+bands[b1]+':  ' + $
 strn(gmu_bg>0.<99.,format='(f8.2)') + $
 ' '+greek("plus_minus", /append)+' '+strn(gmu_bg_e>0.<99.,format='(F8.2)')
 

xyouts,83,3,systime(),charsize=.7

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
xyouts,0,3,qastr,charsize=.7

ms_ps_end

;;;;;;;;;;;;;;;;;;;;;;
; convert to pdf & jpg

p=outpath+'/'
name=id+'_sdss_profile'
spawn,'convert +matte -density 196 -resample 72 -rotate -90 '+p+name+'.ps '+p+name+'.jpg'  ; 
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'rm '+p+name+'.ps'

;;;;;;;;;;;;;;;;;;;;;;
; start image plots

if keyword_set(verbose) then $
	print,'images ',form='($,a)'

filename=outpath+'/'+id+'_sdss_images.ps'
ms_ps_start, filename=filename, xsize=10.5,ysize=8.,/inch,$
          /color,bit=8, xoffset=0.15, yoffset = 11,/landscape

!p.charsize=1.0
chrsz=0.75

!p.multi=[0,2,2]

minimsz = 30. / rscale	; 30 arcsec is minimum image size

; g

if gdpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(gdimg,/dim)
 delta = rskyradius_out[0]/60. * 1.01 > 0.5
 scale=rscale[0]/60.0
 ximgrng = ([rx0,rx0-sz[0]]) * scale
 yimgrng = ([0-ry0,sz[1]-ry0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, gdimg, /preserve, charsize=chrsz, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='g'

 setplotcolors

 
 
 if g_bmap_avail then begin 
    col = (rx0-g_col)*scale
    row = (g_row-ry0)*scale
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
 
 
 
 ela=[rsemimajor/60.0,rsemiminor/60.0,0,0,float(90.-rpa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=rsemiminor/rsemimajor

 eli=[rskyradius_out/60.0,ratio*rskyradius_out/60.0,0,0,float(90.-rpa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[rskyradius_in/60.0,ratio*rskyradius_in/60.0,0,0,float(90.-rpa,0)]
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
     plots, (rx0-XV)*scale,(YV-ry0)*scale, color=!red, linestyle=2, thick=3, /data
   endfor 
 endif
;###########   
   
endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No g Image',charsize=1

endelse

; r

if rdpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(rdimg,/dim)
 delta = rskyradius_out[0]/60. * 1.01 > 0.5
 scale=rscale[0]/60.0
 ximgrng = ([rx0,rx0-sz[0]]) * scale
 yimgrng = ([0-ry0,sz[1]-ry0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, rdimg, /preserve, charsize=chrsz, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='r'

 setplotcolors

 if r_bmap_avail then begin 
    col = (rx0-r_col)*scale
    row = (r_row-ry0)*scale
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
 

 
 ela=[rsemimajor/60.0,rsemiminor/60.0,0,0,float(90.-rpa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=rsemiminor/rsemimajor

 eli=[rskyradius_out/60.0,ratio*rskyradius_out/60.0,0,0,float(90.-rpa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[rskyradius_in/60.0,ratio*rskyradius_in/60.0,0,0,float(90.-rpa,0)]
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
     plots, (rx0-XV)*scale,(YV-ry0)*scale, color=!red, linestyle=2, thick=3, /data
   endfor 
 endif
;###########   
   
endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No r Image',charsize=1

endelse

; i

if idpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(idimg,/dim)
 delta = rskyradius_out[0]/60. * 1.01 > 0.5
 scale = rscale[0]/60.0
 ximgrng = ([rx0,rx0-sz[0]]) * scale
 yimgrng = ([0-ry0,sz[1]-ry0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, idimg, /preserve, charsize=chrsz, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='i'

 setplotcolors

 if i_bmap_avail then begin 
    col = (rx0-i_col)*scale
    row = (i_row-ry0)*scale
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
 
 
 ela=[rsemimajor/60.0,rsemiminor/60.0,0,0,float(90.-rpa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=rsemiminor/rsemimajor

 eli=[rskyradius_out/60.0,ratio*rskyradius_out/60.0,0,0,float(90.-rpa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[rskyradius_in/60.0,ratio*rskyradius_in/60.0,0,0,float(90.-rpa,0)]
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
     plots, (rx0-XV)*scale,(YV-ry0)*scale, color=!red, linestyle=2, thick=3, /data
   endfor 
 endif
;###########
   
   
endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No i Image',charsize=1

endelse

; COLOR (make use of last derived values)

if xdpic then begin

 loadct,0,/silent

 plotimage, xdjpg, /preserve, charsize=chrsz, $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='Masked Composite'

 setplotcolors

 ela=[rsemimajor/60.0,rsemiminor/60.0,0,0,float(90.-rpa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3
 plots,ela[2],ela[3],psym=4,symsi=2.,thick=0.5,color=!green,/data
 
 ratio=rsemiminor/rsemimajor
 
 eli=[rskyradius_out/60.0,ratio*rskyradius_out/60.0,0,0,float(90.-rpa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!gray, thick=1
   
 elo=[rskyradius_in/60.0,ratio*rskyradius_in/60.0,0,0,float(90.-rpa,0)]
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
     plots, (rx0-XV)*scale,(YV-ry0)*scale, color=!cyan, linestyle=2, thick=3, /data
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
name=id+'_sdss_images'
spawn,'convert +matte -density 196 -resample 72 -rotate -90 '+p+name+'.ps '+p+name+'.jpg'  ; 
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'rm '+p+name+'.ps'


;;;;;;;;;;;;;;;;;;;;;;
; start image plots2



if g_mask_avail or r_mask_avail or i_mask_avail or z_mask_avail then begin

filename=outpath+'/'+id+'_sdss_images_masks.ps'
ms_ps_start, filename=filename, xsize=10.5,ysize=8.,/inch,$
          /color,bit=8, xoffset=0.15, yoffset = 11,/landscape

!p.charsize=1.0
chrsz=0.75

!p.multi=[0,2,2]

minimsz = 30. / rscale	; 30 arcsec is minimum image size


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if gdpic then begin
 
 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(gdimg,/dim)
 delta = rskyradius_out[0]/60. * 1.01 > 0.5
 scale=rscale[0]/60.0
 ximgrng = ([rx0,rx0-sz[0]]) * scale
 yimgrng = ([0-ry0,sz[1]-ry0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, gdimg, /preserve, charsize=chrsz, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='g'

 setplotcolors


 if g_bmap_avail then begin 
    col = (rx0-g_col)*scale
    row = (g_row-ry0)*scale
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
 
 if g_mask_avail  then begin 
 
      mask_image =  mrdfits(mask_g,0,/silent)
      maskidx = where(mask_image ge 1, nmaskidx)
      if nmaskidx gt 0  then begin
      
	s = size(mask_image)
	ncol = s(1)
	g_col_mask = maskidx mod ncol
	g_row_mask = maskidx / ncol
      endif    
      
    col = (rx0-g_col_mask)*scale
    row = (g_row_mask-ry0)*scale
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
 
 ela=[rsemimajor/60.0,rsemiminor/60.0,0,0,float(90.-rpa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=rsemiminor/rsemimajor

 eli=[rskyradius_out/60.0,ratio*rskyradius_out/60.0,0,0,float(90.-rpa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[rskyradius_in/60.0,ratio*rskyradius_in/60.0,0,0,float(90.-rpa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1
   
endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No g Image',charsize=1

endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if rdpic then begin
 
 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(rdimg,/dim)
 delta = rskyradius_out[0]/60. * 1.01 > 0.5
 scale=rscale[0]/60.0
 ximgrng = ([rx0,rx0-sz[0]]) * scale
 yimgrng = ([0-ry0,sz[1]-ry0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, rdimg, /preserve, charsize=chrsz, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='r'

 setplotcolors


 if r_bmap_avail then begin 
    col = (rx0-r_col)*scale
    row = (r_row-ry0)*scale
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
 
 if r_mask_avail  then begin 
 
      mask_image =  mrdfits(mask_r,0,/silent)
      maskidx = where(mask_image ge 1, nmaskidx)
      if nmaskidx gt 0  then begin
      
	s = size(mask_image)
	ncol = s(1)
	r_col_mask = maskidx mod ncol
	r_row_mask = maskidx / ncol
      endif    
      
    col = (rx0-r_col_mask)*scale
    row = (r_row_mask-ry0)*scale
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
 
 ela=[rsemimajor/60.0,rsemiminor/60.0,0,0,float(90.-rpa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=rsemiminor/rsemimajor

 eli=[rskyradius_out/60.0,ratio*rskyradius_out/60.0,0,0,float(90.-rpa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[rskyradius_in/60.0,ratio*rskyradius_in/60.0,0,0,float(90.-rpa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1
   
endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No r Image',charsize=1

endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if idpic then begin
 
 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(idimg,/dim)
 delta = rskyradius_out[0]/60. * 1.01 > 0.5
 scale=rscale[0]/60.0
 ximgrng = ([rx0,rx0-sz[0]]) * scale
 yimgrng = ([0-ry0,sz[1]-ry0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, idimg, /preserve, charsize=chrsz, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='i'

 setplotcolors


 if i_bmap_avail then begin 
    col = (rx0-i_col)*scale
    row = (i_row-ry0)*scale
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
 
 if i_mask_avail  then begin 
 
      mask_image =  mrdfits(mask_i,0,/silent)
      maskidx = where(mask_image ge 1, nmaskidx)
      if nmaskidx gt 0  then begin
      
	s = size(mask_image)
	ncol = s(1)
	i_col_mask = maskidx mod ncol
	i_row_mask = maskidx / ncol
      endif    
      
    col = (rx0-i_col_mask)*scale
    row = (i_row_mask-ry0)*scale
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
 
 ela=[rsemimajor/60.0,rsemiminor/60.0,0,0,float(90.-rpa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=rsemiminor/rsemimajor

 eli=[rskyradius_out/60.0,ratio*rskyradius_out/60.0,0,0,float(90.-rpa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[rskyradius_in/60.0,ratio*rskyradius_in/60.0,0,0,float(90.-rpa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1
   
endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No i Image',charsize=1

endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


zdjpgfile =jpgpath+'/'+id+'_z.jpg'

if file_test(zdjpgfile) then begin
 read_jpeg,zdjpgfile,zdjpg,/true
 sz=size(zdjpg,/dim)
  zdimg=bytarr(sz[1],sz[2])
  zdimg[*,*]=zdjpg[0,*,*]

 if z_bmap_avail then begin 
      bmap_image =  mrdfits(bmap_z,0,/silent)
      b_maskidx = where(bmap_image ge 1, nmaskidx)
      if nmaskidx gt 0 then begin
      
	s = size(bmap_image)
	ncol = s(1)
	z_col = b_maskidx mod ncol
	z_row = b_maskidx / ncol
      endif
 endif

zdpic=1
endif else zdpic=0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if zdpic then begin
 
 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(zdimg,/dim)
 delta = rskyradius_out[0]/60. * 1.01 > 0.5
 scale=rscale[0]/60.0
 ximgrng = ([rx0,rx0-sz[0]]) * scale
 yimgrng = ([0-ry0,sz[1]-ry0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, zdimg, /preserve, charsize=chrsz, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='z'

 setplotcolors


 if z_bmap_avail then begin 
    col = (rx0-z_col)*scale
    row = (z_row-ry0)*scale
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
 
 if z_mask_avail  then begin 
 
      mask_image =  mrdfits(mask_z,0,/silent)
      maskidx = where(mask_image ge 1, nmaskidx)
      if nmaskidx gt 0  then begin
      
	s = size(mask_image)
	ncol = s(1)
	z_col_mask = maskidx mod ncol
	z_row_mask = maskidx / ncol
      endif    
      
    col = (rx0-z_col_mask)*scale
    row = (z_row_mask-ry0)*scale
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
 
 ela=[rsemimajor/60.0,rsemiminor/60.0,0,0,float(90.-rpa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=rsemiminor/rsemimajor

 eli=[rskyradius_out/60.0,ratio*rskyradius_out/60.0,0,0,float(90.-rpa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[rskyradius_in/60.0,ratio*rskyradius_in/60.0,0,0,float(90.-rpa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1
   
endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No z Image',charsize=1

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
name=id+'_sdss_images_masks'

; ; thanks! I modify the policy.xml, comment this line" <policy domain="coder" rights="none" pattern="LAEBL"> ". It does work ! Thanks ! –
; /etc/ImageMagick/policy.xml – jianwei May 25 '17 at 13:44
; 

spawn,'convert +matte -density 196 -resample 72 -rotate -90 '+p+name+'.ps '+p+name+'.jpg'  ; 
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'rm '+p+name+'.ps'

endif else begin
    p=outpath+'/'
    name=id+'_sdss_images_masks'
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
; file = 'stv_background_vertices.dat'

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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;
