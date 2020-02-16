pro galex_plotradprof, id, pathtoprofile=pathtoprofile, type=type, $
    fintfile=fintfile, nintfile=nintfile, maskimgfile=maskimgfile, $
    uvjpgpath=uvjpgpath, dssfile=dssfile ,outpath=outpath, verbose=verbose, $
    yuan13=yuan13
;
; note: set yuan13 keyword to use Yuan et al. 2013 empirical extintion coeff's
;	Yuan coeff's assume SFD98 over-estimates E(B-V) by 14% (see table 2)
;
; Warning this procedure spawns convert and ps2pdf14
; (without any checking) to convert postscript output 
; to pdf and jpeg files. Can just comment out if postscript
; is OK. 
;

nofuv=0
nonuv=0
nintdata=0
fintdata=0
xdpic=0
nfpic=0
fdpic=0
dssdata=0

if not keyword_set(pathtoprofile) then pathtoprofile='./'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; files names from galex_radprof

ntotfile =pathtoprofile+'/'+id+'_NUV_totprofile.dat'
nannfile =pathtoprofile+'/'+id+'_NUV_annprofile.dat'
nbgfile  =pathtoprofile+'/'+id+'_NUV_background.dat'
nellfile =pathtoprofile+'/'+id+'_NUV_ellipsepar.dat'
nasymflxfile =pathtoprofile+'/'+id+'_NUV_asymptotic.dat'
naperflxfile =pathtoprofile+'/'+id+'_NUV_aperture.dat'

ftotfile =pathtoprofile+'/'+id+'_FUV_totprofile.dat'
fannfile =pathtoprofile+'/'+id+'_FUV_annprofile.dat'
fbgfile  =pathtoprofile+'/'+id+'_FUV_background.dat'
fellfile =pathtoprofile+'/'+id+'_FUV_ellipsepar.dat'
fasymflxfile =pathtoprofile+'/'+id+'_FUV_total.dat'
faperflxfile =pathtoprofile+'/'+id+'_FUV_aperture.dat'

;;;;;;;;;;;;;;;;;;;;;;;
; read in profile data

;nuv

if file_exist(ntotfile) and file_exist(nannfile)and $
   file_exist(nbgfile) and file_exist(nellfile) then begin

 readcol, ntotfile, ntot_a, ntot_mag, ntot_mag_e, ntot_int, ntot_int_e, $
          ntot_int_e_cnt, ntot_int_e_bg, ntot_bg, ntot_npixels,/silent

 readcol, nasymflxfile, nyf_a, nyf_mag, nyf_mag_e, nyf_int, nyf_int_e, $
          nyf_int_e_cnt, nyf_int_e_bg, nyf_bg, nyf_npixels,/silent

 readcol, naperflxfile, naf_a, naf_mag, naf_mag_e, naf_int, naf_int_e, $
          naf_int_e_cnt, naf_int_e_bg, naf_bg, naf_npixels,/silent

 readcol, nannfile, nann_a, nann_mu, nann_mu_e, nann_int, nann_int_e, $
          nann_int_e_cnt, nann_int_e_bg, nann_bg, nann_npixels,/silent

 readcol, nellfile, nra_cen, ndec_cen, nsemimajor, nsemiminor, npa,/silent

 readcol, nbgfile, nbg, nbg_e, nmu_bg, nmu_bg_e, nscale, $
          nskyradius_in, nskyradius_out,/silent

endif else nonuv=1

;fuv

if file_exist(ftotfile) and file_exist(fannfile)and $
   file_exist(fbgfile) and file_exist(fellfile) then begin

 readcol, ftotfile, ftot_a, ftot_mag, ftot_mag_e, ftot_int, ftot_int_e, $
          ftot_int_e_cnt, ftot_int_e_bg, ftot_bg, ftot_npixels,/silent

 readcol, fasymflxfile, fyf_a, fyf_mag, fyf_mag_e, fyf_int, fyf_int_e, $
          fyf_int_e_cnt, fyf_int_e_bg, fyf_bg, fyf_npixels,/silent

 readcol, faperflxfile, faf_a, faf_mag, faf_mag_e, faf_int, faf_int_e, $
          faf_int_e_cnt, faf_int_e_bg, faf_bg, faf_npixels,/silent

 readcol, fannfile, fann_a, fann_mu, fann_mu_e, fann_int, fann_int_e, $
          fann_int_e_cnt, fann_int_e_bg, fann_bg, fann_npixels,/silent

 readcol, fellfile, fra_cen, fdec_cen, fsemimajor, fsemiminor, fpa,/silent

 readcol, fbgfile, fbg, fbg_e, fmu_bg, fmu_bg_e, fscale, $
          fskyradius_in, fskyradius_out,/silent

endif else nofuv=1

pre='GALEX_PLOTRADPROF: '
if nofuv and nonuv then begin
 print,pre+'Can not find any profile data.'
 return
endif

if keyword_set(verbose) then $
	print,'galex ',form='($,a)'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set some useful values to
; set up plots regardless of band

if ~nonuv then begin
 a=ntot_a
 tot_mag=ntot_mag
 ann_mu=nann_mu 
 ra_cen=nra_cen[0]
 dec_cen=ndec_cen[0]
 semimajor=nsemimajor[0]
 semiminor=nsemiminor[0]
 pa=npa[0]>0.
endif else begin
 a=ftot_a
 tot_mag=ftot_mag
 ann_mu=fann_mu
 ra_cen=fra_cen[0]
 dec_cen=fdec_cen[0]
 semimajor=fsemimajor[0]
 semiminor=fsemiminor[0]
 pa=fpa[0]>0.
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in uv int headers

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

if keyword_set(dssfile) and file_exist(dssfile) then begin
 dss = mrdfits(dssfile, 0, dhdr,/silent)
 extast, dhdr, dssastr
 AD2XY, ra_cen,dec_cen, dssastr, dssx0, dssy0
 getrot,dhdr,rot_dss,cdelt_dss
 dssscale = abs(cdelt_dss[0])*3600.
 dssdata=1
endif
 
;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in jpg images

xdjpgfile =uvjpgpath+'/'+id+'_FUVNUV.jpg'
ndjpgfile =uvjpgpath+'/'+id+'_NUV.jpg'
fdjpgfile =uvjpgpath+'/'+id+'_FUV.jpg'

if keyword_set(xdjpgfile) then begin
 read_jpeg,xdjpgfile,xdjpg,/true
 xdpic=1

 ;mask the 2color image

 if keyword_set(maskimgfile) and (nintdata or fintdata) then begin
  if keyword_set(verbose) then $
	 print,'masking ... ',form='($,a)'
  mskimg = glga_getmask(maskimgfile,sz,astr,as_pix)
  maskidx = where(mskimg ge 1,nmaskidx)
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

endif

if keyword_set(ndjpgfile) and not nonuv then begin
 read_jpeg,ndjpgfile,ndjpg,/true
 sz=size(ndjpg,/dim)
 ndimg=bytarr(sz[1],sz[2])
 ndimg[*,*]=ndjpg[0,*,*]
 ndpic=1
endif else ndpic=0

if keyword_set(fdjpgfile) and not nofuv then begin
 read_jpeg,fdjpgfile,fdjpg,/true
 sz=size(fdjpg,/dim)
 fdimg=bytarr(sz[1],sz[2])
 fdimg[*,*]=fdjpg[0,*,*]
 fdpic=1
endif else fdpic=0

;;;;;;;;;;;;;;;;;;;;;;
; start profile plots

if keyword_set(verbose) then $
	print,'profiles ',form='($,a)'

filename=outpath+'/'+id+'_galex_profile.ps'
ms_ps_start, filename=filename, xsize=10.5,ysize=8.,/inch,$
          /color,bit=8, xoffset=0.15, yoffset = 10.75,/landscape

!p.multi=[0,2,3]
!p.charsize=2

setplotcolors
psize=0.75
ncol=!magenta
fcol=!cyan

;growth curve

if ~nofuv then begin
 bad=where(ftot_mag lt 0., nbad)
 if nbad gt 0 then ftot_mag[bad] = !values.f_nan
 fminmax=minmax(ftot_mag, /nan)
 min=fminmax[0]
 max=fminmax[1]
endif
if ~nonuv then begin
 bad=where(ntot_mag lt 0., nbad)
 if nbad gt 0 then ntot_mag[bad] = !values.f_nan
 nminmax=minmax(ntot_mag, /nan)
 min=nminmax[0]
 max=nminmax[1]
endif
if ~nofuv and ~nonuv then begin
 min=min([fminmax[0],nminmax[0]])
 max=max([fminmax[1],nminmax[1]])
endif
ydel = max-min
yr0 = max+ydel*0.1
yr1 = min-ydel*0.2

xdel = nskyradius_out[0]/60.
xr0 = min(a/60.) - xdel*0.02
xr1 = nskyradius_out[0]/60. + xdel*0.02

; growth curve linear x-axis
plot,[a/60],[tot_mag],yr=[yr0, yr1],/ys,/nodata,$
    xr = [xr0,xr1], /xs,$
    xtitle='R!da!n [arcmin]',ytitle='Growth Curve [AB mag]'
;
; make sure we are including the sky annulus in the plot
if a[n_elements(a)-1] gt nskyradius_in[0] then begin
	oplot,[nskyradius_in[0],nskyradius_out[0]]/60., $
		[min-ydel*0.1,min-ydel*0.1], color=!blue
	xyouts,nskyradius_in[0]/60.,min-ydel*0.12,'SKY ANN', $
		color=!blue, charsize=0.5
endif

if ~noNUV then begin
 oploterror, [a/60.],[ntot_mag],[ntot_mag_e*0],[ntot_mag_e],$
  color=ncol,errcolor=ncol,psym=-sym(1, psize = psize), /nohat
 oplot,!x.crange,[nyf_mag,nyf_mag],linesty=1,color=ncol
 oplot,[nyf_a,nyf_a]/60.,[yr0,nyf_mag],line=2,color=ncol,thick=4
endif

if ~noFUV then begin
 oploterror, [a/60.],[ftot_mag],[ftot_mag_e*0],[ftot_mag_e],$
  color=fcol,errcolor=fcol,psym=-sym(1, psize = psize), /nohat
 oplot,!x.crange,[fyf_mag,fyf_mag],linesty=1,color=fcol
 oplot,[fyf_a,fyf_a]/60.,[yr0,fyf_mag],line=2,color=fcol,thick=4
endif

 oplot,[semimajor,semimajor]/60.,[yr0, yr1],thick=4,color=!red
 xx = semimajor/60. - abs(!x.crange[1]-!x.crange[0])*0.09
 yy = !y.crange[0] - abs(!y.crange[0]-!y.crange[1])*0.05
 xyouts,xx,yy,'APR',color=!red,charsiz=0.8


; growth curve log x-axis
plot,[a/60],[tot_mag],yr=[yr0, yr1],/ys,/nodata,$
    xr = [0.08, xr1], /xs, /xlog,$
    xtit='R!da!n [arcmin]',ytit='Growth Curve [AB mag]'

xyouts,0.1,yr1+ydel*0.14,'ASY',color=!blue,charsiz=0.8

if ~noNUV then begin
 oploterror, [a/60.],[ntot_mag],[ntot_mag_e*0],[ntot_mag_e],$
  color=ncol,errcolor=ncol,psym=-sym(1, psize = psize), /nohat
 oplot,[0.08,1.e9],[nyf_mag,nyf_mag],linesty=1,color=ncol
endif

if ~nofUV then begin
 oploterror, [a/60.],[ftot_mag],[ftot_mag_e*0],[ftot_mag_e],$
  color=fcol,errcolor=fcol,psym=-sym(1, psize = psize), /nohat
 oplot,[0.08,1.e9],[fyf_mag,fyf_mag],linesty=1,color=fcol
endif

 oplot,[semimajor/60,semimajor/60],[yr0, yr1],thick=4,color=!red

;annular surface brightness

if ~nofuv then begin
 bad=where(fann_mu lt 0., nbad)
 if nbad gt 0 then fann_mu[bad] = !values.f_nan
 fminmax=minmax(fann_mu, /nan)
 min=fminmax[0]
 max=fminmax[1]
endif
if ~nonuv then begin
 bad=where(nann_mu lt 0., nbad)
 if nbad gt 0 then nann_mu[bad] = !values.f_nan
 nminmax=minmax(nann_mu, /nan)
 min=nminmax[0]
 max=nminmax[1]
endif
if ~nofuv and ~nonuv then begin
 min=min([fminmax[0],nminmax[0]])
 max=max([fminmax[1],nminmax[1]])
endif

; annular surface brightness linear x-axis
plot,[a/60],[ann_mu],yr=[max+0.5, min-0.5],/ys,/nodata,$
    xr = [xr0,xr1], /xs,$
    xtit='R!da!n [arcmin]',$
    ytit=greek('mu',/append)+' [AB mag/arcsec!u2!n]'

if ~noNUV then $
 oploterror, [a/60.],[nann_mu],[nann_mu_e*0],[nann_mu_e],$
  color=ncol,errcolor=ncol,psym=-sym(1, psize = psize), /nohat

if ~noFUV then $
 oploterror, [a/60.],[fann_mu],[fann_mu_e*0],[fann_mu_e],$
  color=fcol,errcolor=fcol,psym=-sym(1, psize = psize), /nohat

 oplot,[semimajor/60,semimajor/60],[max+0.5,min-0.5],thick=4,color=!red
 xx = semimajor/60. - abs(!x.crange[1]-!x.crange[0])*0.09
 yy = !y.crange[1] + abs(!y.crange[0]-!y.crange[1])*0.1
 xyouts,xx,yy,'APR',color=!red,charsiz=0.8

legend,['NUV','FUV'],textcolors=[ncol,fcol],/bot,/left,box=0,charsize=0.8

; annular surface brightness log x-axis
plot,[a/60],[ann_mu],yr=[max+0.5, min-0.5],/ys,/nodata,$
    xr = [0.08, xr1], /xs, /xlog, $
    xtit='R!da!n [arcmin]',$
    ytit=greek('mu',/append)+' [AB mag/arcsec!u2!n]'

if ~noNUV then $
 oploterror, [a/60.],[nann_mu],[nann_mu_e*0],[nann_mu_e],$
  color=ncol,errcolor=ncol,psym=-sym(1, psize = psize), /nohat

if ~noFUV then $
 oploterror, [a/60.],[fann_mu],[fann_mu_e*0],[fann_mu_e],$
  color=fcol,errcolor=fcol,psym=-sym(1, psize = psize), /nohat

 oplot,[semimajor/60,semimajor/60],[max+0.5,min-0.5],thick=4,color=!red

legend,['NUV','FUV'],textcolors=[ncol,fcol],/bot,/left,box=0,charsize=0.8


; text

!p.charsize=1.0

plot,[0,100],[0,100],xs=4,ys=4,/nodata, pos = [0.085, 0.01, 0.98, 0.33]

xyouts,0,90,id,charsize=1.5
if nintdata then ntime=SXPAR(nhdr,'EXPTIME')>0 else ntime=-999
if fintdata then ftime=SXPAR(fhdr,'EXPTIME')>0 else ftime=-999
xyouts,0,80,'FUV EXPTIME [s]:  '+strn(ftime)
xyouts,0,70,'NUV EXPTIME [s]:  '+strn(ntime)
xyouts,0,60,'R.A.  [J2K]:  '+strn(ra_cen[0],format='(f12.6)')
xyouts,0,50,'DEC [J2K]:  '+strn(dec_cen[0],format='(f12.6)')
xyouts,0,40,'SEMIMAJOR [arcmin]:  '+strn(semimajor/60,format='(f6.2)')
xyouts,0,30,'RATIO (a/b):  '+strn(semimajor/semiminor,format='(f6.2)')
xyouts,0,20,'P.A. [deg]:  '+strn(pa[0],format='(f6.2)')
if keyword_set(type) then $
	xyouts,0,10,'TYP:  '+strtrim(type,2)

if nofuv then begin
	faf_mag=!values.f_nan 
	faf_mag_e=0
endif

xyouts,28,80,'APR:  m(FUV!do!n) = ' + $
 strn(faf_mag>0.<99.,format='(f6.2)')+ ' '+greek("plus_minus", /append)+$
 ' '+strn(faf_mag_e>0.01,format='(F5.2)')

if nonuv then begin
	naf_mag=!values.f_nan 
	naf_mag_e=0
endif

xyouts,28,70,'APR:  m(NUV!do!n) = ' + $
 strn(naf_mag>0.<99.,format='(f6.2)')+ ' '+greek("plus_minus", /append)+$
 ' '+strn(naf_mag_e>0.01,format='(F5.2)')

uvcolor=(faf_mag-naf_mag)>(-99.)<99.
err=sqrt(faf_mag_e^2 + naf_mag_e^2)>0.01<9.
xyouts,28,60,'APR:  FUV-NUV = '+$
 strn(uvcolor,format='(f6.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(err,format='(F5.2)')


if nofuv then begin
	fyf_mag=!values.f_nan 
	fyf_mag_e=0
	fyf_a=!values.f_nan
endif

xyouts,28,40,'ASY:  m(FUV!do!n) = ' + $
 strn(fyf_mag>0.<99.,format='(f6.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(fyf_mag_e>0.01,format='(F5.2)') $
 + '  R!da!n = '+strn(fyf_a/60.,format='(f8.2)')+"'"

if nonuv then begin
	nyf_mag=!values.f_nan 
	nyf_mag_e=0
	nyf_a=!values.f_nan
endif

xyouts,28,30,'ASY:  m(NUV!do!n) = ' + $
 strn(nyf_mag>0.<99.,format='(f6.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(nyf_mag_e>0.01,format='(F5.2)') $
 + '  R!da!n = '+strn(nyf_a/60.,format='(f8.2)')+"'"

uvcolor=(fyf_mag-nyf_mag)>(-99.)<99.
err=sqrt(fyf_mag_e^2 + nyf_mag_e^2)>0.01<9.
xyouts,28,20,'ASY:  FUV-NUV = '+$
 strn(uvcolor,format='(f6.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(err,format='(F5.2)')

if nofuv then begin
	fmu_bg=!values.f_nan 
	fmu_bg_e=0
endif

xyouts,75,80,greek('mu',/append)+'!DBG!N FUV:  '+$
 strn(fmu_bg>0.<99.,format='(f6.2)') +$
 ' '+greek("plus_minus", /append)+' '+strn(fmu_bg_e>0.<99.,format='(F5.2)')

if nonuv then begin
	nmu_bg=!values.f_nan 
	nmu_bg_e=0
endif

xyouts,75,70,greek('mu',/append)+'!DBG!N NUV:  '+$
 strn(nmu_bg>0.<99.,format='(f6.2)') +$
 ' '+greek("plus_minus", /append)+' '+strn(nmu_bg_e>0.<99.,format='(F5.2)')

; galactic coords to index dust maps
euler, ra_cen, dec_cen, l, b, 1
ebv=dust_getval(l,b,/interp)			; SFD98 E(B-V)
A_NUV=glga_getextin(ebv,'NUV',yuan13=yuan13)	; get extinction for NUV
A_FUV=glga_getextin(ebv,'FUV',yuan13=yuan13)	; get extinction for FUV

xyouts,75,60,'GAL E(B-V):  '+strn(ebv,format='(f5.3)')
xyouts,75,50,'GAL A!dFUV!n:  '+strn(A_FUV,format='(f5.3)')
xyouts,75,40,'GAL A!dNUV!n:  '+strn(A_NUV,format='(f5.3)')

xyouts,0.96,0.077,systime(),charsize=.7,/norm,ori=90.

; read in qa
qafile = strmid(nintfile,0,strpos(nintfile,'_')+1) +'qa.txt'
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
name=id+'_galex_profile'
spawn,'convert +matte -density 196 -resample 72 -rotate -90 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'rm '+p+name+'.ps'

;;;;;;;;;;;;;;;;;;;;;;
; start image plots

if keyword_set(verbose) then $
	print,'images ',form='($,a)'

filename=outpath+'/'+id+'_galex_images.ps'
ms_ps_start, filename=filename, xsize=10.5,ysize=8.,/inch,$
          /color,bit=8, xoffset=0.15, yoffset = 11,/landscape

!p.charsize=0.75
!p.multi=[0,2,2]


; FUV 

if fdpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(fdimg,/dim)
 delta = fskyradius_out[0]/60. * 1.01 > 0.5
 scale = fscale[0]/60.0 ; arcmins
 ximgrng = ([fx0,fx0-sz[0]]) * scale
 yimgrng = ([0-fy0, sz[1]-fy0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, fdimg, /preserve, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='FUV'

 setplotcolors

 ela=[fsemimajor/60.0,fsemiminor/60.0,0,0,float(90.-fpa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3,noclip=0

 ratio=fsemiminor/fsemimajor

 eli=[fskyradius_out/60.0,ratio*fskyradius_out/60.0,0,0,float(90.-fpa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1,noclip=0

 elo=[fskyradius_in/60.0,ratio*fskyradius_in/60.0,0,0,float(90.-fpa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1,noclip=0

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No FUV Image',charsize=1

endelse

; NUV 

if ndpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(ndimg,/dim)
 delta = nskyradius_out[0]/60. * 1.01 > 0.5
 scale = nscale[0]/60.0 ; arcmins
 ximgrng = ([nx0,(nx0-sz[0])]) * scale
 yimgrng = ([0-ny0, sz[1]-ny0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, ndimg, /preserve, color=max(!d.n_colors), $
	 xran=xpltrng,yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='NUV'

 setplotcolors

 ela=[nsemimajor/60.0,nsemiminor/60.0,0,0,float(90.-npa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3,noclip=0

 ratio=nsemiminor/nsemimajor

 eli=[nskyradius_out/60.0,ratio*nskyradius_out/60.0,0,0,float(90.-npa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1,noclip=0

 elo=[nskyradius_in/60.0,ratio*nskyradius_in/60.0,0,0,float(90.-npa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1,noclip=0

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No NUV Image',charsize=1

endelse

; DSS

if dssdata then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz_dss = size(dss,/dim)
 dssximgrng = ([dssx0,(dssx0-sz_dss[0])]) * dssscale[0]/60.
 dssyimgrng = ([0-dssy0, sz_dss[1]-dssy0]) * dssscale[0]/60.
 dss = asinh(dss*1.)
 meanclip, dss, dss_sky, dss_sigma
 min = dss_sky*0.98
 max = max(dss) < (dss_sky+(20*dss_sigma))
 plotimage, dss, /preserve, range=[min,max], color=max(!d.n_colors), $
	xran=xpltrng, yran=ypltrng, imgxrange=dssximgrng, imgyrange=dssyimgrng,$
	xtitle='arcmin', ytitle='arcmin', title='POSS2-R'

 setplotcolors

 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3,noclip=0
 plots,ela[2],ela[3],psym=4,symsi=2.,thick=0.5,color=!green,/data

 if ~nonuv then ratio=nsemiminor/nsemimajor else ratio=fsemiminor/fsemimajor
 if ~nonuv then skyradius_out=nskyradius_out else skyradius_out=fskyradius_out
 if ~nonuv then skyradius_in=nskyradius_in else skyradius_in=fskyradius_in

 eli=[skyradius_out/60.0,ratio*skyradius_out/60.0,0,0,float(90.-pa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1,noclip=0

 elo=[skyradius_in/60.0,ratio*skyradius_in/60.0,0,0,float(90.-pa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1,noclip=0

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.31],[0.5],'No POSS2 Image',charsize=1

endelse

; COLOR (make use of last derived values)

if xdpic then begin

 loadct,0,/silent

 plotimage, xdjpg, /preserve, $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='Masked Composite'

 setplotcolors

 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3,noclip=0
 plots,ela[2],ela[3],psym=4,symsi=2.,thick=0.5,color=!green,/data

 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!gray, thick=1,noclip=0

 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!gray, thick=1,noclip=0

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
name=id+'_galex_images'
spawn,'convert +matte -density 196 -resample 72 -rotate -90 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'rm '+p+name+'.ps'

if keyword_set(verbose) then $
	print,'Done.'

;;;;;;;;;;;;
; all done

return
end
