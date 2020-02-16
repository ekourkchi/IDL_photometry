pro wise_plotcolors, id, pathtoprofile=pathtoprofile, type=type, $
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
read_radprof,id,bands[b3], w4tot_a, w4tot_mag, w4tot_mag_e, w4ann_mu, w4ann_mu_e, $
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
if now1band and now3band and now3band and now4band then return
; fix missing values
if now1band then begin
	w1tf_mag=!values.f_nan
	w1tf_mag_e=!values.f_nan
	w1af_mag=!values.f_nan
	w1af_mag_e=!values.f_nan
endif
if now2band then begin
	w2tf_mag=!values.f_nan
	w2tf_mag_e=!values.f_nan
	w2af_mag=!values.f_nan
	w2af_mag_e=!values.f_nan
endif
if now3band then begin
	w3tf_mag=!values.f_nan
	w3tf_mag_e=!values.f_nan
	w3af_mag=!values.f_nan
	w3af_mag_e=!values.f_nan
endif
if now4band then begin
	w4tf_mag=!values.f_nan
	w4tf_mag_e=!values.f_nan
	w4af_mag=!values.f_nan
	w4af_mag_e=!values.f_nan
endif
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
	AD2XY, w1ra_cen[0] ,w1dec_cen[0], w1astr, w1x0, w1y0
	sz = [sxpar(w1hdr,'NAXIS1'),sxpar(w1hdr,'NAXIS2')]
	getrot,w1hdr,w1rot,w1cdelt
	as_pix = abs(w1cdelt[0])*3600.
endif else begin
	print,'No image data for '+id+', returning.'
	return
endelse

;;;;;;;;;;;;;;;;;;;;;;
; start profile plots

if keyword_set(verbose) then $
	print,'plot wise colors ... ',form='($,a)'

filename=outpath+'/'+id+'_wise_colors.ps'
ms_ps_start, filename=filename, xsize=7,ysize=10,/inch,$
          /color,/true,bits=8, xoffset=0.75, yoffset = 0.75

!p.multi=[0,2,3]
!p.charsize=1.2

psize=0.5
;
; clean total mags
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
;
; clean annular mags
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
;
; galactic coords to index dust maps
euler, ra_cen, dec_cen, l, b, 1
ebv=dust_getval(l,b,/interp)			; SFD98 E(B-V)
a_w1  = glga_getextin(ebv,'w1',yuan13=yuan13)	; get extinction for w1
a_w2  = glga_getextin(ebv,'w2',yuan13=yuan13)	; get extinction for w2
a_w3  = glga_getextin(ebv,'w3',yuan13=yuan13)	; get extinction for w3
a_w4  = glga_getextin(ebv,'w4',yuan13=yuan13)	; get extinction for w4
;
; color colors
cs = [cgcolor('forest green'), cgcolor('black'), $
      cgcolor('org5'), cgcolor('brown')]

; plot total color
nn=where(finite(w1tot_mag) and finite(w2tot_mag), cnt1finite)
if cnt1finite gt 0 then begin
	op1color = (w1tot_mag-a_w1) - (w2tot_mag-a_w2)
	op1cerr  = sqrt(w1tot_mag_e^2 + w2tot_mag_e^2)
endif
nn=where(finite(w2tot_mag) and finite(w3tot_mag), cnt2finite)
if cnt2finite gt 0 then begin
	op2color = (w2tot_mag-a_w2) - (w3tot_mag-a_w3)
	op2cerr  = sqrt(w2tot_mag_e^2 + w3tot_mag_e^2)
endif
nn=where(finite(w3tot_mag) and finite(w4tot_mag), cnt3finite)
if cnt3finite gt 0 then begin
	op3color = (w3tot_mag-a_w3) - (w4tot_mag-a_w4)
	op3cerr  = sqrt(w3tot_mag_e^2 + w4tot_mag_e^2)
endif
nn=where(finite(w1tot_mag) and finite(w3tot_mag), cnt4finite)
if cnt4finite gt 0 then begin
	op4color = (w1tot_mag-a_w1) - (w3tot_mag-a_w3)
	op4cerr  = sqrt(w1tot_mag_e^2 + w3tot_mag_e^2)
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
xl = !x.crange[1] - (!x.crange[1]-!x.crange[0]) * 0.045
yl = !y.crange[1]
legend,['(w1-w2)!Ddr!N','(w2-w3)!Ddr!N','(w3-w4)!Ddr!N','(w1-w3)!Ddr!N'], $
	box=0, charsi=0.7, textcolors=cs, pos=[xl,yl]
;
; plot annular color
nn=where(finite(w1ann_mu) and finite(w2ann_mu), cnt1finite)
if cnt1finite gt 0 then begin
	op1color = (w1ann_mu-a_w1) - (w2ann_mu-a_w2)
	op1cerr  = sqrt(w1ann_mu_e^2 + w2ann_mu_e^2)
endif
nn=where(finite(w2ann_mu) and finite(w3ann_mu), cnt2finite)
if cnt2finite gt 0 then begin
	op2color = (w2ann_mu-a_w2) - (w3ann_mu-a_w3)
	op2cerr  = sqrt(w2ann_mu_e^2 + w3ann_mu_e^2)
endif
nn=where(finite(w3ann_mu) and finite(w4ann_mu), cnt3finite)
if cnt3finite gt 0 then begin
	op3color = (w3ann_mu-a_w3) - (w4ann_mu-a_w4)
	op3cerr  = sqrt(w3ann_mu_e^2 + w4ann_mu_e^2)
endif
nn=where(finite(w1ann_mu) and finite(w3ann_mu), cnt4finite)
if cnt4finite gt 0 then begin
	op4color = (w1ann_mu-a_w1) - (w3ann_mu-a_w3)
	op4cerr  = sqrt(w1ann_mu_e^2 + w3ann_mu_e^2)
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
; read in w1w2w3 jpg images

xdjpgfile =jpgpath+'/'+id+'_wx.jpg'

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

	minimsz = 30. / w1scale	; 30 arcsec is minimum image size
	delta=(semimajor/w1scale * 1.01) > minimsz 
	sz=size(xdjpg,/dim)
	x1=w1x0-delta > 0
	x2=w1x0+delta < sz[1]-1
	y1=w1y0-delta > 0
	y2=w1y0+delta < sz[2]-1
	scale=w1scale[0]/60.0
	xdjpg=xdjpg[*,x1:x2,y1:y2]
	sz=size(xdjpg,/dim)
	xrng = ([0+sz[1]/2,0-sz[1]/2]) * scale
	yrng = ([0-sz[2]/2,0+sz[2]/2]) * scale

	loadct,0,/silent

	plotimage,xdjpg,/preserve,color=max(!d.n_colors), $
		imgxrange=xrng, imgyrange=yrng, $
		xtitle='arcmin', ytitle='arcmin', $
		title='w1w2w3 '+msklab+'Composite'

	ela=[semimajor/60.0,semiminor/60.0,0,0,float(90.-pa,0)]
	tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
		linestyle=0,color=cgcolor('red'), thick=3, noclip=0
	plots,ela[2],ela[3],psym=4,symsi=2.,color=cgcolor('green'),/data, $
		thick=0.5

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.3],[0.5],'No w1w2w3 Composite Image',charsize=1

endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in w1w2w4 jpg images

xdjpgfile =jpgpath+'/'+id+'_ww.jpg'

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

	minimsz = 30. / w1scale	; 30 arcsec is minimum image size
	delta=(semimajor/w1scale * 1.01) > minimsz 
	sz=size(xdjpg,/dim)
	x1=w1x0-delta > 0
	x2=w1x0+delta < sz[1]-1
	y1=w1y0-delta > 0
	y2=w1y0+delta < sz[2]-1
	scale=w1scale[0]/60.0
	xdjpg=xdjpg[*,x1:x2,y1:y2]
	sz=size(xdjpg,/dim)
	xrng = ([0+sz[1]/2,0-sz[1]/2]) * scale
	yrng = ([0-sz[2]/2,0+sz[2]/2]) * scale

	loadct,0,/silent

	plotimage,xdjpg,/preserve,color=max(!d.n_colors), $
		imgxrange=xrng, imgyrange=yrng, $
		xtitle='arcmin', ytitle='arcmin', $
		title='w1w2w4 '+msklab+'Composite'

	ela=[semimajor/60.0,semiminor/60.0,0,0,float(90.-pa,0)]
	tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
		linestyle=0,color=cgcolor('red'), thick=3, noclip=0
	plots,ela[2],ela[3],psym=4,symsi=2.,color=cgcolor('green'),/data, $
		thick=0.5

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.3],[0.5],'No w1w2w4 Composite Image',charsize=1

endelse


;;;;;;;;;;;;;;;;;;;;;
; text

!p.charsize=0.8

plot,[0,100],[0,100],xs=4,ys=4,/nodata, pos = [0.085, 0.01, 0.98, 0.33]

plots,[0.06,.98,.98,0.06,0.06],[.32,.32,.01,0.01,0.32],/norm, $
	color=cgcolor('Charcoal')

plots,[0.06,0.98],[0.235,0.235],/norm,color=cgcolor('Charcoal')

xyouts,0,85,id,charsi=1.5
xyouts,0,75,'WISE',charsi=1.25
xyouts,30,85,'R.A. [J2K, deg]:  '+strn(ra_cen,format='(f12.6)')
xyouts,30,75,'DEC  [J2K, deg]:  '+strn(dec_cen,format='(f12.6)')
xyouts,70,85,'Gal. Lon [l, deg]:  '+strn(l,format='(f12.6)')
xyouts,70,75,'Gal. Lat [b, deg]:  '+strn(b,format='(f12.6)')

xyouts,0,60,"SEMIMAJOR:  "+strn(semimajor/60,format='(f6.2)')+"'"
xyouts,0,50,'RATIO (a/b):  '+strn(semimajor/semiminor,format='(f6.2)')
xyouts,0,40,'P.A. [deg]:  '+strn(pa,format='(f6.2)')
xyouts,0,30,'TYP:  '+typ

xyouts,70,60,'GAL E(B-V):  '+strn(ebv,format='(f5.3)')
xyouts,70,50,'GAL A!dw1!n:  '+strn(a_w1,format='(f5.3)')
xyouts,70,40,'GAL A!dw2!n:  '+strn(a_w2,format='(f5.3)')
xyouts,70,30,'GAL A!dw3!n:  '+strn(a_w3,format='(f5.3)')
xyouts,70,20,'GAL A!dw4!n:  '+strn(a_w4,format='(f5.3)')

xyouts,30,60,'APR: (w1-w2)!ddr!n = ' + $
	strn(((w1af_mag-a_w1)-(w2af_mag-a_w2)),format='(f5.2)') + $
	' '+greek("plus_minus", /append) + $
	' '+strn(sqrt(w1af_mag_e^2+w2af_mag_e^2)>0.01,format='(f4.2)')

xyouts,30,50,'APR: (w2-w3)!ddr!n = ' + $
	strn(((w2af_mag-a_w2)-(w3af_mag-a_w3)),format='(f5.2)') + $
	' '+greek("plus_minus", /append) + $
	' '+strn(sqrt(w2af_mag_e^2+w3af_mag_e^2)>0.01,format='(f4.2)')

xyouts,30,40,'APR: (w3-w4)!ddr!n = ' + $
	strn(((w3af_mag-a_w3)-(w4af_mag-a_w4)),format='(f5.2)') + $
	' '+greek("plus_minus", /append) + $
	' '+strn(sqrt(w3af_mag_e^2+w4af_mag_e^2)>0.01,format='(f4.2)')

xyouts,30,30,'APR: (w1-w3)!ddr!n = ' + $
	strn(((w1af_mag-a_w1)-(w3af_mag-a_w3)),format='(f5.2)') + $
	' '+greek("plus_minus", /append) + $
	' '+strn(sqrt(w1af_mag_e^2+w3af_mag_e^2)>0.01,format='(f4.2)')

xyouts,77,2,systime(),charsize=.7

!p.multi=0
ms_ps_end

;;;;;;;;;;;;;;;;;;;;;;
; convert to pdf & jpg
if keyword_set(verbose) then print,'convert ',form='($,a)'

p=outpath+'/'
name=id+'_wise_colors'
spawn,'convert +matte -density 196 -resample 72 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'rm '+p+name+'.ps'


;;;;;;;;;;;;
; all done
if keyword_set(verbose) then print,'Done.'

return
end
