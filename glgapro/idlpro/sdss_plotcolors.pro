pro sdss_plotcolors, id, pathtoprofile=pathtoprofile, type=type, $
    intfile=intfile, maskimgfile=maskimgfile, $
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
	pathtoprofile=pathtoprofile
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
endif
if nogband then begin
	gtf_mag=!values.f_nan
	gtf_mag_e=!values.f_nan
	gaf_mag=!values.f_nan
	gaf_mag_e=!values.f_nan
endif
if norband then begin
	rtf_mag=!values.f_nan
	rtf_mag_e=!values.f_nan
	raf_mag=!values.f_nan
	raf_mag_e=!values.f_nan
endif
if noiband then begin
	itf_mag=!values.f_nan
	itf_mag_e=!values.f_nan
	iaf_mag=!values.f_nan
	iaf_mag_e=!values.f_nan
endif
if nozband then begin
	ztf_mag=!values.f_nan
	ztf_mag_e=!values.f_nan
	zaf_mag=!values.f_nan
	zaf_mag_e=!values.f_nan
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

;;;;;;;;;;;;;;;;;;;;;;
; start profile plots

if keyword_set(verbose) then $
	print,'plot sdss colors ... ',form='($,a)'

filename=outpath+'/'+id+'_sdss_colors.ps'
ms_ps_start, filename=filename, xsize=7,ysize=10,/inch,$
          /color,/true,bits=8, xoffset=0.75, yoffset = 0.75

!p.multi=[0,2,3]
!p.charsize=1.2

psize=0.5
;
; clean total mags
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
;
; clean annular mags
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
;
; galactic coords to index dust maps
euler, ra_cen, dec_cen, l, b, 1
ebv=dust_getval(l,b,/interp)			; SFD98 E(B-V)
a_u  = glga_getextin(ebv,'u',yuan13=yuan13)	; get extinction for u
a_g  = glga_getextin(ebv,'g',yuan13=yuan13)	; get extinction for g
a_r  = glga_getextin(ebv,'r',yuan13=yuan13)	; get extinction for r
a_i  = glga_getextin(ebv,'i',yuan13=yuan13)	; get extinction for i
a_z  = glga_getextin(ebv,'z',yuan13=yuan13)	; get extinction for z
;
; color colors
cs = [cgcolor('forest green'), cgcolor('black'), $
      cgcolor('org5'), cgcolor('brown')]

; plot total color
nn=where(finite(utot_mag) and finite(gtot_mag), cnt1finite)
if cnt1finite gt 0 then begin
	op1color = (utot_mag-a_u) - (gtot_mag-a_g)
	op1cerr  = sqrt(utot_mag_e^2 + gtot_mag_e^2)
endif
nn=where(finite(gtot_mag) and finite(rtot_mag), cnt2finite)
if cnt2finite gt 0 then begin
	op2color = (gtot_mag-a_g) - (rtot_mag-a_r)
	op2cerr  = sqrt(gtot_mag_e^2 + rtot_mag_e^2)
endif
nn=where(finite(rtot_mag) and finite(itot_mag), cnt3finite)
if cnt3finite gt 0 then begin
	op3color = (rtot_mag-a_r) - (itot_mag-a_i)
	op3cerr  = sqrt(rtot_mag_e^2 + itot_mag_e^2)
endif
nn=where(finite(itot_mag) and finite(ztot_mag), cnt4finite)
if cnt4finite gt 0 then begin
	op4color = (itot_mag-a_i) - (ztot_mag-a_z)
	op4cerr  = sqrt(itot_mag_e^2 + ztot_mag_e^2)
endif

apidx = where(a le r_ap)
if total([cnt1finite,cnt2finite,cnt3finite,cnt4finite]) gt 0 then begin
	min=min([op1color[apidx],op2color[apidx], $
		 op3color[apidx],op4color[apidx]],/nan)
	max=max([op1color[apidx],op2color[apidx], $
		 op3color[apidx],op4color[apidx]],/nan)
	yttl = 'Total Color'
	plot, [a/60.], [op1color], psym=-sym(1, psize = psize),$
		charsize=chrsz, ytit = yttl, xtit='R!da!n [arc minutes]',$
		xr = [0, (r_ap/60.)+0.2], /xs, $
		yr = [min-0.05, max+0.05], /ys, /nodata
	if cnt1finite gt 0 then $
	    oploterror, [a/60.],[op1color],[op1cerr*0],[op1cerr],$
		psym=-sym(1, psize = psize), /nohat, $
		color=cs[0], errcolor=cs[0]
	if cnt2finite gt 0 then $
	    oploterror, [a/60.],[op2color],[op2cerr*0],[op2cerr],$
		psym=-sym(1, psize = psize), /nohat, $
		color=cs[1], errcolor=cs[1]
	if cnt3finite gt 0 then $
	    oploterror, [a/60.],[op3color],[op3cerr*0],[op3cerr],$
		psym=-sym(1, psize = psize), /nohat, $
		color=cs[2], errcolor=cs[2]
	if cnt4finite gt 0 then $
	    oploterror, [a/60.],[op4color],[op4cerr*0],[op4cerr],$
		psym=-sym(1, psize = psize), /nohat, $
		color=cs[3], errcolor=cs[3]
endif 
oplot, [0,1e4],[0,0], color=cgcolor('black')
oplot,[r_ap/60,r_ap/60],[-100,100], color=cgcolor('red')
xl = !x.crange[1] - (!x.crange[1]-!x.crange[0]) * 0.035
yl = !y.crange[1]
legend,['(u-g)!Ddr!N','(g-r)!Ddr!N','(r-i)!Ddr!N','(i-z)!Ddr!N'], $
	box=0, charsi=0.9, textcolors=cs, pos=[xl,yl]
;
; plot annular color
nn=where(finite(uann_mu) and finite(gann_mu), cnt1finite)
if cnt1finite gt 0 then begin
	op1color = (uann_mu-a_u) - (gann_mu-a_g)
	op1cerr  = sqrt(uann_mu_e^2 + gann_mu_e^2)
endif
nn=where(finite(gann_mu) and finite(rann_mu), cnt2finite)
if cnt2finite gt 0 then begin
	op2color = (gann_mu-a_g) - (rann_mu-a_r)
	op2cerr  = sqrt(gann_mu_e^2 + rann_mu_e^2)
endif
nn=where(finite(rann_mu) and finite(iann_mu), cnt3finite)
if cnt3finite gt 0 then begin
	op3color = (rann_mu-a_r) - (iann_mu-a_i)
	op3cerr  = sqrt(rann_mu_e^2 + iann_mu_e^2)
endif
nn=where(finite(iann_mu) and finite(zann_mu), cnt4finite)
if cnt4finite gt 0 then begin
	op4color = (iann_mu-a_i) - (zann_mu-a_z)
	op4cerr  = sqrt(iann_mu_e^2 + zann_mu_e^2)
endif

if total([cnt1finite,cnt2finite,cnt3finite,cnt4finite]) gt 0 then begin
	min=min([op1color[apidx],op2color[apidx], $
		 op3color[apidx],op4color[apidx]],/nan)
	max=max([op1color[apidx],op2color[apidx], $
		 op3color[apidx],op4color[apidx]],/nan)
	yttl = 'Annular Color'
	plot, [a/60.], [op1color], psym=-sym(1, psize = psize),$
		charsize=chrsz, ytit = yttl, xtit='R!da!n [arc minutes]',$
		xr = [0, (r_ap/60.)+0.2], /xs, $
		yr = [min-0.05, max+0.05], /ys, /nodata
	if cnt1finite gt 0 then $
	    oploterror, [a/60.],[op1color],[op1cerr*0],[op1cerr],$
		psym=-sym(1, psize = psize), /nohat, $
		color=cs[0], errcolor=cs[0]
	if cnt2finite gt 0 then $
	    oploterror, [a/60.],[op2color],[op2cerr*0],[op2cerr],$
		psym=-sym(1, psize = psize), /nohat, $
		color=cs[1], errcolor=cs[1]
	if cnt3finite gt 0 then $
	    oploterror, [a/60.],[op3color],[op3cerr*0],[op3cerr],$
		psym=-sym(1, psize = psize), /nohat, $
		color=cs[2], errcolor=cs[2]
	if cnt4finite gt 0 then $
	    oploterror, [a/60.],[op4color],[op4cerr*0],[op4cerr],$
		psym=-sym(1, psize = psize), /nohat, $
		color=cs[3], errcolor=cs[3]
endif 
oplot, [0,1e4],[0,0], color=cgcolor('black')
oplot,[r_ap/60,r_ap/60],[-100,100], color=cgcolor('red')

;;;;;;;;;;;;;;;;;;;;;;
; start image plots

;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in gri jpg images

xdjpgfile =jpgpath+'/'+id+'_x.jpg'

xdpic=0
if file_test(xdjpgfile) then begin
	read_jpeg,xdjpgfile,xdjpg,/true
	xdpic=1
endif


; COLOR (make use of last derived values)

if xdpic then begin
	msklab = ''
;
;mask the 3color image if requested
	if keyword_set(maskimgfile) then begin
		if keyword_set(verbose) then $
			print,'mask1 ',form='($,a)'
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
			msklab = 'Masked '
		endif
 	endif

	minimsz = 30. / rscale	; 30 arcsec is minimum image size
	delta=(semimajor/rscale * 1.01) > minimsz 
	sz=size(xdjpg,/dim)
	x1=rx0-delta > 0
	x2=rx0+delta < sz[1]-1
	y1=ry0-delta > 0
	y2=ry0+delta < sz[2]-1
	scale=rscale[0]/60.0
	xdjpg=xdjpg[*,x1:x2,y1:y2]
	sz=size(xdjpg,/dim)
	xrng = ([0+sz[1]/2,0-sz[1]/2]) * scale
	yrng = ([0-sz[2]/2,0+sz[2]/2]) * scale

	loadct,0,/silent

	plotimage,xdjpg,/preserve,color=max(!d.n_colors), $
		imgxrange=xrng, imgyrange=yrng, $
		xtitle='arcmin', ytitle='arcmin', $
		title='gri '+msklab+'Composite'

	ela=[semimajor/60.0,semiminor/60.0,0,0,float(90.-pa,0)]
	tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
		linestyle=0,color=cgcolor('red'), thick=3, noclip=0
	plots,ela[2],ela[3],psym=4,symsi=2.,color=cgcolor('green'),/data, $
		thick=0.5

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.3],[0.5],'No gri Composite Image',charsize=1

endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in urz jpg images

xdjpgfile =jpgpath+'/'+id+'_w.jpg'

xdpic=0
if file_test(xdjpgfile) then begin
	read_jpeg,xdjpgfile,xdjpg,/true
	xdpic=1
endif


; COLOR (make use of last derived values)

if xdpic then begin
	msklab = ''
;
;mask the 3color image if requested
	if keyword_set(maskimgfile) then begin
		if keyword_set(verbose) then $
			print,'mask2 ',form='($,a)'
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
			msklab = 'Masked '
		endif
 	endif

	minimsz = 30. / rscale	; 30 arcsec is minimum image size
	delta=(semimajor/rscale * 1.01) > minimsz 
	sz=size(xdjpg,/dim)
	x1=rx0-delta > 0
	x2=rx0+delta < sz[1]-1
	y1=ry0-delta > 0
	y2=ry0+delta < sz[2]-1
	scale=rscale[0]/60.0
	xdjpg=xdjpg[*,x1:x2,y1:y2]
	sz=size(xdjpg,/dim)
	xrng = ([0+sz[1]/2,0-sz[1]/2]) * scale
	yrng = ([0-sz[2]/2,0+sz[2]/2]) * scale

	loadct,0,/silent

	plotimage,xdjpg,/preserve,color=max(!d.n_colors), $
		imgxrange=xrng, imgyrange=yrng, $
		xtitle='arcmin', ytitle='arcmin', $
		title='urz '+msklab+'Composite'

	ela=[semimajor/60.0,semiminor/60.0,0,0,float(90.-pa,0)]
	tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
		linestyle=0,color=cgcolor('red'), thick=3, noclip=0
	plots,ela[2],ela[3],psym=4,symsi=2.,color=cgcolor('green'),/data, $
		thick=0.5

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.3],[0.5],'No urz Composite Image',charsize=1

endelse


;;;;;;;;;;;;;;;;;;;;;
; text

!p.charsize=0.8

plot,[0,100],[0,100],xs=4,ys=4,/nodata, pos = [0.085, 0.01, 0.98, 0.33]

plots,[0.06,.98,.98,0.06,0.06],[.32,.32,.01,0.01,0.32],/norm, $
	color=cgcolor('Charcoal')

plots,[0.06,0.98],[0.235,0.235],/norm,color=cgcolor('Charcoal')

xyouts,0,85,id,charsi=1.5
xyouts,0,75,'SDSS',charsi=1.25
xyouts,30,85,'R.A. [J2K, deg]:  '+strn(ra_cen,format='(f12.6)')
xyouts,30,75,'DEC  [J2K, deg]:  '+strn(dec_cen,format='(f12.6)')
xyouts,70,85,'Gal. Lon [l, deg]:  '+strn(l,format='(f12.6)')
xyouts,70,75,'Gal. Lat [b, deg]:  '+strn(b,format='(f12.6)')

xyouts,0,60,"SEMIMAJOR:  "+strn(semimajor/60,format='(f6.2)')+"'"
xyouts,0,50,'RATIO (a/b):  '+strn(semimajor/semiminor,format='(f6.2)')
xyouts,0,40,'P.A. [deg]:  '+strn(pa,format='(f6.2)')
xyouts,0,30,'TYP:  '+typ

xyouts,70,60,'GAL E(B-V):  '+strn(ebv,format='(f5.3)')
xyouts,70,50,'GAL A!du!n:  '+strn(a_u,format='(f5.3)')
xyouts,70,40,'GAL A!dg!n:  '+strn(a_g,format='(f5.3)')
xyouts,70,30,'GAL A!dr!n:  '+strn(a_r,format='(f5.3)')
xyouts,70,20,'GAL A!di!n:  '+strn(a_i,format='(f5.3)')
xyouts,70,10,'GAL A!dz!n:  '+strn(a_z,format='(f5.3)')

xyouts,30,60,'APR: (u-g)!ddr!n = ' + $
	strn(((uaf_mag-a_u)-(gaf_mag-a_g)),format='(f5.2)') + $
	' '+greek("plus_minus", /append) + $
	' '+strn(sqrt(uaf_mag_e^2+gaf_mag_e^2)>0.01,format='(f4.2)')

xyouts,30,50,'APR: (g-r)!ddr!n = ' + $
	strn(((gaf_mag-a_g)-(raf_mag-a_r)),format='(f5.2)') + $
	' '+greek("plus_minus", /append) + $
	' '+strn(sqrt(gaf_mag_e^2+raf_mag_e^2)>0.01,format='(f4.2)')

xyouts,30,40,'APR: (r-i)!ddr!n = ' + $
	strn(((raf_mag-a_r)-(iaf_mag-a_i)),format='(f5.2)') + $
	' '+greek("plus_minus", /append) + $
	' '+strn(sqrt(raf_mag_e^2+iaf_mag_e^2)>0.01,format='(f4.2)')

xyouts,30,30,'APR: (i-z)!ddr!n = ' + $
	strn(((iaf_mag-a_i)-(zaf_mag-a_z)),format='(f5.2)') + $
	' '+greek("plus_minus", /append) + $
	' '+strn(sqrt(iaf_mag_e^2+zaf_mag_e^2)>0.01,format='(f4.2)')

xyouts,77,2,systime(),charsize=.7

!p.multi=0
ms_ps_end

;;;;;;;;;;;;;;;;;;;;;;
; convert to pdf & jpg
if keyword_set(verbose) then print,'convert ',form='($,a)'

p=outpath+'/'
name=id+'_sdss_colors'
spawn,'convert +matte -density 196 -resample 72 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'rm '+p+name+'.ps'


;;;;;;;;;;;;
; all done
if keyword_set(verbose) then print,'Done.'

return
end
