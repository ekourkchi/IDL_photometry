pro galex_plotcolors, id, pathtoprofile=pathtoprofile, type=type, $
    fintfile=fintfile, nintfile=nintfile, maskimgfile=maskimgfile, $
    uvjpgpath=uvjpgpath, outpath=outpath, verbose=verbose, yuan13=yuan13
;
; note: set yuan13 keyword to use Yuan et al. 2013 empirical extintion coeff's
;	Yuan coeff's assume SFD98 over-estimates E(B-V) by 14% (see table 2)
;
; Warning this procedure spawns convert and ps2pdf14
; (without any checking) to convert postscript output 
; to pdf and jpeg files. Can just comment out if postscript
; is OK. 
;


if not keyword_set(pathtoprofile) then pathtoprofile='./'
if not keyword_set(type) then typ = '-' else typ = strtrim(type,2)
if not keyword_set(uvjpgpath) then uvjpgpath='./'
if not keyword_set(outpath) then outpath='./'

;;;;;;;;;;;;;;;;;;;;;;;
; read in profile data

;nuv

read_radprof,id,'NUV', ntot_a, ntot_mag, ntot_mag_e, nann_mu, nann_mu_e, $
	nra_cen, ndec_cen, nsemimajor, nsemiminor, npa, nscale, $
	ntf_mag, ntf_mag_e, naf_mag, naf_mag_e, nmu_bg, nmu_bg_e, $
	nskyradius_in, nskyradius_out, ntf_a, nyf_mag, nyf_mag_e, nyf_a, $
	tot_int=ntot_int, pathtoprofile=pathtoprofile
; check input data
if ntot_a[0] eq -1 then return

;fuv

read_radprof,id,'FUV', ftot_a, ftot_mag, ftot_mag_e, fann_mu, fann_mu_e, $
	fra_cen, fdec_cen, fsemimajor, fsemiminor, fpa, fscale, $
	ftf_mag, ftf_mag_e, faf_mag, faf_mag_e, fmu_bg, fmu_bg_e, $
	fskyradius_in, fskyradius_out, ftf_a, fyf_mag, fyf_mag_e, fyf_a, $
	tot_int=ftot_int, pathtoprofile=pathtoprofile
; check input data
if ftot_a[0] eq -1 then return

pre='GALEX_PLOTCOLORS: '

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set some useful values to
; set up plots regardless of band

a=ntot_a
tot_mag=ntot_mag
ann_mu=nann_mu 
ra_cen=nra_cen[0]
dec_cen=ndec_cen[0]
semimajor=nsemimajor[0]
semiminor=nsemiminor[0]
pa=npa[0]>0.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in uv int headers

nintdata=0
if keyword_set(nintfile) and file_exist(nintfile) then begin
	nhdr=headfits(nintfile,ext=0)
	extast,nhdr,nastr
	AD2XY, ra_cen[0] ,dec_cen[0], nastr, nx0, ny0
	nintdata=1
	sz = [sxpar(nhdr,'NAXIS1'),sxpar(nhdr,'NAXIS2')]
	astr = nastr
	getrot,nhdr,nrot,ncdelt
	as_pix = abs(ncdelt[0])*3600.
endif

fintdata=0
if keyword_set(fintfile) and file_exist(fintfile) then begin
	fhdr=headfits(fintfile,ext=0)
	extast,fhdr,fastr
	AD2XY, ra_cen[0] ,dec_cen[0], fastr, fx0, fy0
	fintdata=1
	if not nintdata then begin
		sz = [sxpar(fhdr,'NAXIS1'),sxpar(fhdr,'NAXIS2')]
		astr = fastr
		getrot,fhdr,frot,fcdelt
		as_pix = abs(fcdelt[0])*3600.
	endif
endif

;;;;;;;;;;;;;;;;;;;;;;
; start profile plots

if keyword_set(verbose) then $
	print,'plot galex colors ... ',form='($,a)'

filename=outpath+'/'+id+'_galex_colors.ps'
ms_ps_start, filename=filename, xsize=7,ysize=10,/inch,$
          /color,/true,bits=8,xoffset=0.75, yoffset = 0.75

!p.multi=[0,2,3]
!p.charsize=1.2

psize=0.5
ncol=cgcolor('magenta')
fcol=cgcolor('cyan')
;
; clean total mags
bad=where(ftot_mag lt 0., nbad)
if nbad gt 0 then ftot_mag[bad] = !values.f_nan
bad=where(ntot_mag lt 0., nbad)
if nbad gt 0 then ntot_mag[bad] = !values.f_nan
;
; clean annular mags
bad=where(fann_mu lt 0., nbad)
if nbad gt 0 then fann_mu[bad] = !values.f_nan
bad=where(nann_mu lt 0., nbad)
if nbad gt 0 then nann_mu[bad] = !values.f_nan
;
; galactic coords to index dust maps
euler, ra_cen, dec_cen, l, b, 1
ebv=dust_getval(l,b,/interp)			; SFD98 E(B-V)
A_NUV=glga_getextin(ebv,'NUV',yuan13=yuan13)	; get extinction for NUV
A_FUV=glga_getextin(ebv,'FUV',yuan13=yuan13)	; get extinction for FUV
;
; plot total color
fn=where(finite(ftot_mag),cntffinite)
nn=where(finite(ntot_mag),cntnfinite)

if cntffinite gt 0 and cntnfinite gt 0 then begin
	uvtotcolor = (ftot_mag-A_FUV) - (ntot_mag-A_NUV)
	uvtotcerr = sqrt(ftot_mag_e^2 + ntot_mag_e^2)
	min=min(uvtotcolor,/nan)
	max=max(uvtotcolor,/nan)
	plot, [a/60.], [uvtotcolor], psym=-sym(1, psize = psize),$
		ytitle = 'Total (FUV-NUV)!ddr!n', xtitle = 'R!da!n [arcmin]',$
		xr = [0, max(a/60.)+0.2], /xs, $
		yr = [min-0.05, max+0.05], /ys, /nodata
	oploterror, [a/60.],[uvtotcolor],[uvtotcerr*0],[uvtotcerr],$
		psym=-sym(1, psize = psize), /nohat
	oplot,[fyf_a[0],fyf_a[0]]/60.,[-100,100],line=1,color=cgcolor('cyan')
	oplot,[nyf_a[0],nyf_a[0]]/60.,[-100,100],line=1,color=cgcolor('magenta')
endif else begin
	plot, [a/60.], [a/60.], psym=-sym(1, psize = psize),$
		ytit = 'Total (FUV-NUV)!ddr!n', xtit='R!da!n [arcmin]',$
		xr = [0, max(a/60.)+0.1], /xs, $
		yr = [-1, 1], /ys, /nodata
endelse
oplot, [0,1e4],[0,0], color=cgcolor('black')
oplot, [semimajor,semimajor]/60.,[-100,100], color=cgcolor('red')
xl = !x.crange[1] - (!x.crange[1]-!x.crange[0]) * 0.035
yl = !y.crange[1]
legend,['ASY!DNUV!N','ASY!DFUV!N','APR'],box=0,charsize=0.7, pos=[xl,yl], $
	textcolor=[cgcolor('magenta'),cgcolor('cyan'),cgcolor('red')]
;
; plot annular color
fn=where(finite(fann_mu),cntffinite)
nn=where(finite(nann_mu),cntnfinite)

if cntffinite gt 0 and cntnfinite gt 0 then begin
	uvanncolor = (fann_mu-A_FUV) - (nann_mu-A_NUV)
	uvanncerr = sqrt(fann_mu_e^2 + nann_mu_e^2)
	min=min(uvanncolor,/nan)
	max=max(uvanncolor,/nan)
	plot, [a/60.], [uvanncolor], psym=-sym(1, psize = psize),$
		ytit = 'Annular (FUV-NUV)!ddr!n', xtit='R!da!n [arcmin]',$
		xr = [0, max(a/60.)+0.2], /xs, $
		yr = [min-0.05, max+0.05],/ys,/nodata
	oploterror, [a/60.],[uvanncolor],[uvanncerr*0],[uvanncerr],$
		psym=-sym(1, psize = psize), /nohat
	oplot,[fyf_a[0],fyf_a[0]]/60.,[-100,100],line=1,color=cgcolor('cyan')
	oplot,[nyf_a[0],nyf_a[0]]/60.,[-100,100],line=1,color=cgcolor('magenta')
endif else begin
	plot, [a/60.], [a/60.], psym=-sym(1, psize = psize),$
	    ytit = 'Annular (FUV-NUV)!ddr!n', xtit='R!da!n [arcmin]',$
	    xr = [0, max(a/60.)+0.1], /xs, $
	    yr = [-1, 1],/ys,/nodata
endelse
oplot, [0,1e4],[0,0], color=cgcolor('black')
oplot, [semimajor,semimajor]/60.,[-100,100], color=cgcolor('red')

;;;;;;;;;;;;;;;;;;;;;;
; start image plots

;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in jpg image

xdjpgfile =uvjpgpath+'/'+id+'_FUVNUV.jpg'

xdpix=0
if file_test(xdjpgfile) then begin
	read_jpeg,xdjpgfile,xdjpg,/true
	xdmask = xdjpg
	xdpic=1

endif

if xdpic then begin

	minimsz = 30. / nscale	; 30 arcsec is minimum image size
	delta=(nsemimajor/nscale * 1.01) > minimsz
	sz = size(xdjpg,/dim)
	x1=nx0-delta > 0
	x2=nx0+delta < sz[1]-1
	y1=ny0-delta > 0
	y2=ny0+delta < sz[2]-1
	scale=nscale[0]/60.0
	xdjpg=xdjpg[*,x1:x2,y1:y2]
	sz = size(xdjpg,/dim)
	xrng = ([0+sz[1]/2,0-sz[1]/2]) * scale
	yrng = ([0-sz[2]/2,0+sz[2]/2]) * scale
	loadct,0,/silent
	plotimage,xdjpg,/preserve,color=max(!d.n_colors),$
		imgxrange=xrng, imgyrange=yrng, $
		xtitle='arcmin', ytitle='arcmin', title='FUV,NUV Composite'
	
	ratio=nsemiminor/nsemimajor

	ela=[nsemimajor/60.0,nsemiminor/60.0,0,0,float(90.-npa,0)]
	tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
		linestyle=0,color=cgcolor('red'), thick=3, noclip=0

	eln=[nyf_a[0]/60.0,nyf_a[0]*ratio/60.0,0,0,float(90.-npa,0)]
	tvellipse,eln[0],eln[1],eln[2],eln[3],eln[4],/data,$
		linestyle=1,color=cgcolor('magenta'), thick=3, noclip=0

	elf=[fyf_a[0]/60.0,fyf_a[0]*ratio/60.0,0,0,float(90.-fpa,0)]
	tvellipse,elf[0],elf[1],elf[2],elf[3],elf[4],/data,$
		linestyle=1,color=cgcolor('cyan'), thick=3, noclip=0

	plots,ela[2],ela[3],psym=4,symsi=2.,color=cgcolor('green'),/data, $
		thick=0.5

endif else begin

	plot,[0,1],[0,1],/nodata,xs=4,ys=4
	xyouts,[0.32],[0.5],'No Composite Image',charsize=1

endelse

; COLOR (make use of last derived values)

if xdpic then begin
;
;mask the 2color image
	if keyword_set(maskimgfile) and (nintdata or fintdata) then begin
		if keyword_set(verbose) then $
			print,'mask ',form='($,a)'
		mskimg = glga_getmask(maskimgfile,sz,astr,as_pix)
		maskidx = where(mskimg ge 1,nmaskidx)
		if nmaskidx gt 0 then begin
			xdjpgr = xdmask[0, *, *]
			xdjpgg = xdmask[1, *, *]
			xdjpgb = xdmask[2, *, *]
			xdjpgr[maskidx] = 244
			xdjpgg[maskidx] = 164
			xdjpgb[maskidx] = 96
			xdmask[0, *, *] = xdjpgr
			xdmask[1, *, *] = xdjpgg 
			xdmask[2, *, *] = xdjpgb
		endif
	endif

	minimsz = 30. / nscale	; 30 arcsec is minimum image size
	delta=(nskyradius_out/nscale * 1.01) > minimsz
	sz = size(xdmask,/dim)
	x1=nx0-delta > 0
	x2=nx0+delta < (sz[1]-1)
	y1=ny0-delta > 0
	y2=ny0+delta < (sz[2]-1)
	xdmask=xdmask[*,x1:x2,y1:y2]
	sz = size(xdmask,/dim)
	scale=nscale[0]/60.0
	xrng = ([0+sz[1]/2,0-sz[1]/2]) * scale
	yrng = ([0-sz[2]/2,0+sz[2]/2]) * scale

	loadct,0,/silent

	plotimage,xdmask,/preserve,$
	imgxrange=xrng, imgyrange=yrng, xtitle='arcmin', ytitle='arcmin', $
		title='FUV,NUV Masked Composite'

	tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
		linestyle=0,color=cgcolor('red'), thick=3, noclip=0

	tvellipse,eln[0],eln[1],eln[2],eln[3],eln[4],/data,$
		linestyle=1,color=cgcolor('magenta'), thick=3, noclip=0

	tvellipse,elf[0],elf[1],elf[2],elf[3],elf[4],/data,$
		linestyle=1,color=cgcolor('cyan'), thick=3, noclip=0

	plots,ela[2],ela[3],psym=4,symsi=2.,color=cgcolor('green'),/data, $
		thick=0.5


	eli=[nskyradius_out/60.0,ratio*nskyradius_out/60.0,0,0,float(90.-npa,0)]
	tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
	  linestyle=5,color=cgcolor('white'), thick=1, noclip=0

	elo=[nskyradius_in/60.0,ratio*nskyradius_in/60.0,0,0,float(90.-npa,0)]
	tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
	  linestyle=5,color=cgcolor('white'), thick=1, noclip=0

endif else begin

	plot,[0,1],[0,1],/nodata,xs=4,ys=4
	xyouts,[0.3],[0.5],'No Masked Image',charsize=1

endelse

;;;;;;;;;;;;;;;;;;
; text

!p.charsize=0.8

plot,[0,100],[0,100],xs=4,ys=4,/nodata, pos = [0.085, 0.01, 0.98, 0.33]

plots,[0.06,.98,.98,0.06,0.06],[.32,.32,.01,0.01,0.32],/norm, $
	color=cgcolor('Charcoal')

plots,[0.06,0.98],[0.235,0.235],/norm,color=cgcolor('Charcoal')

xyouts,0,85,id,charsi=1.5
xyouts,0,75,'GALEX',charsi=1.25
xyouts,30,85,'R.A. [J2K, deg]:  '+strn(ra_cen,format='(f12.6)')
xyouts,30,75,'DEC  [J2K, deg]:  '+strn(dec_cen,format='(f12.6)')
xyouts,70,85,'Gal. Lon [l, deg]:  '+strn(l,format='(f12.6)')
xyouts,70,75,'Gal. Lat [b, deg]:  '+strn(b,format='(f12.6)')

if nintdata then ntime=SXPAR(nhdr,'EXPTIME')>0 else ntime=-999
if fintdata then ftime=SXPAR(fhdr,'EXPTIME')>0 else ftime=-999
xyouts,0,60,'FUV EXPT [s]:  '+strn(ftime)
xyouts,0,50,'NUV EXPT [s]:  '+strn(ntime)
xyouts,0,40,"SEMIMAJOR:  "+strn(semimajor/60,format='(f6.2)')+"'"
xyouts,0,30,'RATIO (a/b):  '+strn(semimajor/semiminor,format='(f6.2)')
xyouts,0,20,'P.A. [deg]:  '+strn(pa,format='(f6.2)')
xyouts,0,10,'TYP:  '+typ

xyouts,70,60,'GAL E(B-V):  '+strn(ebv,format='(f5.3)')
xyouts,70,50,'GAL A!dNUV!n:  '+strn(A_NUV,format='(f5.3)')
xyouts,70,40,'GAL A!dFUV!n:  '+strn(A_FUV,format='(f5.3)')

xyouts,30,60,'APR: (FUV-NUV)!ddr!n = ' + $
	strn(((faf_mag-A_FUV)-(naf_mag-A_NUV)),format='(f5.2)') + $
	' '+greek("plus_minus", /append) + $
	strn(sqrt(faf_mag_e^2+naf_mag_e^2)>0.01,format='(f4.2)')

xyouts,30,50,'ASY: (FUV-NUV)!ddr!n = ' + $
	strn(((fyf_mag-A_FUV)-(nyf_mag-A_NUV)),format='(f5.2)') + $
	' '+greek("plus_minus", /append) + $
	strn(sqrt(fyf_mag_e^2+nyf_mag_e^2)>0.01,format='(f4.2)')

xyouts,30,40, 'FUV R!dASY!n = '+string(fyf_a/60.,format='(f5.2)')+"'"
xyouts,30,30, 'NUV R!dASY!n = '+string(nyf_a/60.,format='(f5.2)')+"'"

xyouts,77,2.0,systime(),charsize=.7

!p.multi=0
ms_ps_end

;;;;;;;;;;;;;;;;;;;;;;
; convert to pdf & jpg
if keyword_set(verbose) then print,'convert ',form='($,a)'

p=outpath+'/'
name=id+'_galex_colors'
spawn,'convert +matte -density 196 -resample 72 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'rm '+p+name+'.ps'

;;;;;;;;;;;;
; all done
if keyword_set(verbose) then print,'Done.'

return
end
