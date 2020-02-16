pro twomass_plotcolors, id, pathtoprofile=pathtoprofile, type=type, $
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

bands = ['j','h','k']
nband = n_elements(bands)

if not keyword_set(type) then typ = '-' else typ = strtrim(type,2)
if not keyword_set(pathtoprofile) then pathtoprofile='./'
if not keyword_set(jpgpath) then jpgpath='./'
b0 = 0
b1 = 1
b2 = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in profile data

read_radprof,id,bands[b0], jtot_a, jtot_mag, jtot_mag_e, jann_mu, jann_mu_e, $
	jra_cen, jdec_cen, jsemimajor, jsemiminor, jpa, jscale, $
	jtf_mag, jtf_mag_e, jaf_mag, jaf_mag_e, jmu_bg, jmu_bg_e, $
	jskyradius_in, jskyradius_out, jtf_a, jsymag, jsymag_e, jsyma, $
	pathtoprofile=pathtoprofile
read_radprof,id,bands[b1], htot_a, htot_mag, htot_mag_e, hann_mu, hann_mu_e, $
	hra_cen, hdec_cen, hsemimajor, hsemiminor, hpa, hscale, $
	htf_mag, htf_mag_e, haf_mag, haf_mag_e, hmu_bg, hmu_bg_e, $
	hskyradius_in, hskyradius_out, htf_a, hsymag, hsymag_e, hsyma, $
	pathtoprofile=pathtoprofile
read_radprof,id,bands[b2], ktot_a, ktot_mag, ktot_mag_e, kann_mu, kann_mu_e, $
	kra_cen, kdec_cen, ksemimajor, ksemiminor, kpa, kscale, $
	ktf_mag, ktf_mag_e, kaf_mag, kaf_mag_e, kmu_bg, kmu_bg_e, $
	kskyradius_in, kskyradius_out, ktf_a, ksymag, ksymag_e, ksyma, $
	pathtoprofile=pathtoprofile
; check inputs
if jtot_a[0] eq -1 then nojband = (1 eq 1) else nojband = (1 eq 0)
if htot_a[0] eq -1 then nohband = (1 eq 1) else nohband = (1 eq 0)
if ktot_a[0] eq -1 then nokband = (1 eq 1) else nokband = (1 eq 0)
if nokband then begin
	print,'No primary image photometry for '+id+', returning.'
	return
endif
if nohband and nojband and nokband then return
; fix missing values
if nojband then begin
	jtf_mag=!values.f_nan
	jtf_mag_e=!values.f_nan
	jaf_mag=!values.f_nan
	jaf_mag_e=!values.f_nan
endif
if nohband then begin
	htf_mag=!values.f_nan
	htf_mag_e=!values.f_nan
	haf_mag=!values.f_nan
	haf_mag_e=!values.f_nan
endif
if nokband then begin
	ktf_mag=!values.f_nan
	ktf_mag_e=!values.f_nan
	kaf_mag=!values.f_nan
	kaf_mag_e=!values.f_nan
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set some useful values to
; set up plots regardless of band

 a=ktot_a
 tot_mag=ktot_mag
 ann_mu=kann_mu 
 ra_cen=kra_cen[0]
 dec_cen=kdec_cen[0]
 semimajor=ksemimajor[0]
 semiminor=ksemiminor[0]
 pa=kpa[0]>0.
 r_ap=semimajor
 r_ski=kskyradius_in
 r_sko=kskyradius_out

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in headers

if keyword_set(intfile) and file_exist(intfile) then begin
	khdr=headfits(intfile,ext=0)
	extast,khdr,kastr
	AD2XY, kra_cen[0] ,kdec_cen[0], kastr, kx0, ky0
	sz = [sxpar(khdr,'NAXIS1'),sxpar(khdr,'NAXIS2')]
	getrot,khdr,krot,kcdelt
	as_pix = abs(kcdelt[0])*3600.
endif else begin
	print,'No image data for '+id+', returning.'
	return
endelse

;;;;;;;;;;;;;;;;;;;;;;
; start profile plots

if keyword_set(verbose) then $
	print,'plot 2mass colors ... ',form='($,a)'

filename=outpath+'/'+id+'_2mass_colors.ps'
ms_ps_start, filename=filename, xsize=7,ysize=10,/inch,$
          /color,/true,bits=8, xoffset=0.75, yoffset = 0.75

!p.multi=[0,2,3]
!p.charsize=1.2

psize=0.5
;
; clean total mags
if not nojband then begin
	bad=where(jtot_mag lt 0., nbad)
	if nbad gt 0 then jtot_mag[bad] = !values.f_nan
endif else jtot_mag = !values.f_nan
if not nohband then begin
	bad=where(htot_mag lt 0., nbad)
	if nbad gt 0 then htot_mag[bad] = !values.f_nan
endif else htot_mag = !values.f_nan
if not nokband then begin
	bad=where(ktot_mag lt 0., nbad)
	if nbad gt 0 then ktot_mag[bad] = !values.f_nan
endif else ktot_mag = !values.f_nan
;
; clean annular mags
if not nojband then begin
	bad=where(jann_mu lt 0., nbad)
	if nbad gt 0 then jann_mu[bad] = !values.f_nan
endif else jann_mu = !values.f_nan
if not nohband then begin
	bad=where(hann_mu lt 0., nbad)
	if nbad gt 0 then hann_mu[bad] = !values.f_nan
endif else hann_mu = !values.f_nan
if not nokband then begin
	bad=where(kann_mu lt 0., nbad)
	if nbad gt 0 then kann_mu[bad] = !values.f_nan
endif else kann_mu = !values.f_nan
;
; galactic coords to index dust maps
euler, ra_cen, dec_cen, l, b, 1
ebv=dust_getval(l,b,/interp)			; SFD98 E(B-V)
a_j  = glga_getextin(ebv,'j',yuan13=yuan13)	; get extinction for j
a_h  = glga_getextin(ebv,'h',yuan13=yuan13)	; get extinction for h
a_k  = glga_getextin(ebv,'k',yuan13=yuan13)	; get extinction for k
;
; color colors
cs = [cgcolor('forest green'), cgcolor('black'), cgcolor('org5')]

; plot total color
nn=where(finite(jtot_mag) and finite(htot_mag), cnt1finite)
if cnt1finite gt 0 then begin
	op1color = (jtot_mag-a_j) - (htot_mag-a_h)
	op1cerr  = sqrt(jtot_mag_e^2 + htot_mag_e^2)
endif
nn=where(finite(htot_mag) and finite(ktot_mag), cnt2finite)
if cnt2finite gt 0 then begin
	op2color = (htot_mag-a_h) - (ktot_mag-a_k)
	op2cerr  = sqrt(htot_mag_e^2 + ktot_mag_e^2)
endif
nn=where(finite(jtot_mag) and finite(ktot_mag), cnt3finite)
if cnt3finite gt 0 then begin
	op3color = (jtot_mag-a_j) - (ktot_mag-a_k)
	op3cerr  = sqrt(jtot_mag_e^2 + ktot_mag_e^2)
endif

apidx = where(a le r_ap)
if total([cnt1finite,cnt2finite,cnt3finite]) gt 0 then begin
	min=min([op1color[apidx],op2color[apidx],op3color[apidx]],/nan)
	max=max([op1color[apidx],op2color[apidx],op3color[apidx]],/nan)
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
endif 
oplot, [0,1e4],[0,0], color=cgcolor('black')
oplot,[r_ap/60,r_ap/60],[-100,100], color=cgcolor('red')
xl = !x.crange[1] - (!x.crange[1]-!x.crange[0]) * 0.035
yl = !y.crange[1]
legend,['(j-h)!Ddr!N','(h-k)!Ddr!N','(j-k)!Ddr!N'], $
	box=0, charsi=0.9, textcolors=cs, pos = [xl,yl]
;
; plot annular color
nn=where(finite(jann_mu) and finite(hann_mu), cnt1finite)
if cnt1finite gt 0 then begin
	op1color = (jann_mu-a_j) - (hann_mu-a_h)
	op1cerr  = sqrt(jann_mu_e^2 + hann_mu_e^2)
endif
nn=where(finite(hann_mu) and finite(kann_mu), cnt2finite)
if cnt2finite gt 0 then begin
	op2color = (hann_mu-a_h) - (kann_mu-a_k)
	op2cerr  = sqrt(hann_mu_e^2 + kann_mu_e^2)
endif
nn=where(finite(jann_mu) and finite(kann_mu), cnt3finite)
if cnt3finite gt 0 then begin
	op3color = (jann_mu-a_j) - (kann_mu-a_k)
	op3cerr  = sqrt(jann_mu_e^2 + kann_mu_e^2)
endif

if total([cnt1finite,cnt2finite,cnt3finite]) gt 0 then begin
	min=min([op1color[apidx],op2color[apidx],op3color[apidx]],/nan)
	max=max([op1color[apidx],op2color[apidx],op3color[apidx]],/nan)
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
endif 
oplot, [0,1e4],[0,0], color=cgcolor('black')
oplot,[r_ap/60,r_ap/60],[-100,100], color=cgcolor('red')

;;;;;;;;;;;;;;;;;;;;;;
; start image plots

;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in jhk jpg images

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
	;	if keyword_set(maskimgfile) then begin
	;		if keyword_set(verbose) then $
	;			print,'mask1 ',form='($,a)'
	;		mskimg = glga_getmask(maskimgfile,sz,rastr,as_pix)
	;		maskidx = where(mskimg ge 1, nmaskidx)
	;		if nmaskidx gt 0 then begin
	;			xdjpgr = xdjpg[0, *, *]
	;			xdjpgg = xdjpg[1, *, *]
	;			xdjpgb = xdjpg[2, *, *]
	;			xdjpgr[maskidx] = 244
	;			xdjpgg[maskidx] = 164
	;			xdjpgb[maskidx] = 96
	;			xdjpg[0, *, *] = xdjpgr
	;			xdjpg[1, *, *] = xdjpgg 
	;			xdjpg[2, *, *] = xdjpgb
	;			msklab = 'Masked '
	;		endif
	; 	endif

	minimsz = 30. / kscale	; 30 arcsec is minimum image size
	delta=(semimajor/kscale * 1.01) > minimsz 
	sz=size(xdjpg,/dim)
	x1=kx0-delta > 0
	x2=kx0+delta < sz[1]-1
	y1=ky0-delta > 0
	y2=ky0+delta < sz[2]-1
	scale=kscale[0]/60.0
	xdjpg=xdjpg[*,x1:x2,y1:y2]
	sz=size(xdjpg,/dim)
	xrng = ([0+sz[1]/2,0-sz[1]/2]) * scale
	yrng = ([0-sz[2]/2,0+sz[2]/2]) * scale

	loadct,0,/silent

	plotimage,xdjpg,/preserve,color=max(!d.n_colors), $
		imgxrange=xrng, imgyrange=yrng, $
		xtitle='arcmin', ytitle='arcmin', $
		title='jhk '+msklab+'Composite'

	ela=[semimajor/60.0,semiminor/60.0,0,0,float(90.-pa,0)]
	tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
		linestyle=0,color=cgcolor('red'), thick=3, noclip=0
	plots,ela[2],ela[3],psym=4,symsi=2.,color=cgcolor('green'),/data, $
		thick=0.5

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.3],[0.5],'No jhk Composite Image',charsize=1

endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in jhk jpg images

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

	minimsz = 30. / kscale	; 30 arcsec is minimum image size
	delta=(semimajor/kscale * 1.01) > minimsz 
	sz=size(xdjpg,/dim)
	x1=kx0-delta > 0
	x2=kx0+delta < sz[1]-1
	y1=ky0-delta > 0
	y2=ky0+delta < sz[2]-1
	scale=kscale[0]/60.0
	xdjpg=xdjpg[*,x1:x2,y1:y2]
	sz=size(xdjpg,/dim)
	xrng = ([0+sz[1]/2,0-sz[1]/2]) * scale
	yrng = ([0-sz[2]/2,0+sz[2]/2]) * scale

	loadct,0,/silent

	plotimage,xdjpg,/preserve,color=max(!d.n_colors), $
		imgxrange=xrng, imgyrange=yrng, $
		xtitle='arcmin', ytitle='arcmin', $
		title='jhk '+msklab+'Composite'

	ela=[semimajor/60.0,semiminor/60.0,0,0,float(90.-pa,0)]
	tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
		linestyle=0,color=cgcolor('red'), thick=3, noclip=0
	plots,ela[2],ela[3],psym=4,symsi=2.,color=cgcolor('green'),/data, $
		thick=0.5

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.3],[0.5],'No jhk Composite Image',charsize=1

endelse


;;;;;;;;;;;;;;;;;;;;;
; text

!p.charsize=0.8

plot,[0,100],[0,100],xs=4,ys=4,/nodata, pos = [0.085, 0.01, 0.98, 0.33]

plots,[0.06,.98,.98,0.06,0.06],[.32,.32,.01,0.01,0.32],/norm, $
	color=cgcolor('Charcoal')

plots,[0.06,0.98],[0.235,0.235],/norm,color=cgcolor('Charcoal')

xyouts,0,85,id,charsi=1.5
xyouts,0,75,'2MASS',charsi=1.25
xyouts,30,85,'R.A. [J2K, deg]:  '+strn(ra_cen,format='(f12.6)')
xyouts,30,75,'DEC  [J2K, deg]:  '+strn(dec_cen,format='(f12.6)')
xyouts,70,85,'Gal. Lon [l, deg]:  '+strn(l,format='(f12.6)')
xyouts,70,75,'Gal. Lat [b, deg]:  '+strn(b,format='(f12.6)')

xyouts,0,60,"SEMIMAJOR:  "+strn(semimajor/60,format='(f6.2)')+"'"
xyouts,0,50,'RATIO (a/b):  '+strn(semimajor/semiminor,format='(f6.2)')
xyouts,0,40,'P.A. [deg]:  '+strn(pa,format='(f6.2)')
xyouts,0,30,'TYP:  '+typ

xyouts,70,60,'GAL E(B-V):  '+strn(ebv,format='(f5.3)')
xyouts,70,50,'GAL A!dj!n:  '+strn(a_j,format='(f5.3)')
xyouts,70,40,'GAL A!dh!n:  '+strn(a_h,format='(f5.3)')
xyouts,70,30,'GAL A!dk!n:  '+strn(a_k,format='(f5.3)')

xyouts,30,60,'APR: (j-h)!ddr!n = ' + $
	strn(((jaf_mag-a_j)-(haf_mag-a_h)),format='(f5.2)') + $
	' '+greek("plus_minus", /append) + $
	' '+strn(sqrt(jaf_mag_e^2+haf_mag_e^2)>0.01,format='(f4.2)')

xyouts,30,50,'APR: (h-k)!ddr!n = ' + $
	strn(((haf_mag-a_h)-(kaf_mag-a_k)),format='(f5.2)') + $
	' '+greek("plus_minus", /append) + $
	' '+strn(sqrt(haf_mag_e^2+kaf_mag_e^2)>0.01,format='(f4.2)')

xyouts,30,40,'APR: (j-k)!ddr!n = ' + $
	strn(((jaf_mag-a_j)-(kaf_mag-a_k)),format='(f5.2)') + $
	' '+greek("plus_minus", /append) + $
	' '+strn(sqrt(jaf_mag_e^2+kaf_mag_e^2)>0.01,format='(f4.2)')

xyouts,77,2,systime(),charsize=.7

!p.multi=0
ms_ps_end

;;;;;;;;;;;;;;;;;;;;;;
; convert to pdf & jpg
if keyword_set(verbose) then print,'convert ',form='($,a)'

p=outpath+'/'
name=id+'_2mass_colors'
spawn,'convert +matte -density 196 -resample 72 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'rm '+p+name+'.ps'


;;;;;;;;;;;;
; all done
if keyword_set(verbose) then print,'Done.'

return
end
