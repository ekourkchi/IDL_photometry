pro irac_plotradprof, id, pathtoprofile=pathtoprofile, type=type, $
    intfile=intfile, maskimgfile=maskimgfile, yuan13=yuan13, $
    jpgpath=jpgpath, dssfile=dssfile ,outpath=outpath, verbose=verbose
;
; note: set yuan13 keyword to use Yuan et al. 2013 empirical extintion coeff's
;	Yuan coeff's assume SFD98 over-estimates E(B-V) by 14% (see table 2)
;
; Warning this procedure spawns convert and ps2pdf14
; (without any checking) to convert postscript output 
; to pdf and jpeg files. Can just comment out if postscript
; is OK. 
;

bands = ['3p6um','4p5um']
pands = ['3.6um','4.5um']
nband = n_elements(bands)
nonuv=0
nintdata=0
xdpic=0
dssdata=0

if not keyword_set(pathtoprofile) then pathtoprofile='./'

pathtofits=repstr(pathtoprofile,'photometry','irac/fits')

if not keyword_set(jpgpath) then jpgpath='./'
b0 = 0
b1 = 1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in profile data

b0file=pathtoprofile+'/'+id+'_'+bands[0]+'_totprofile.dat'
b1file=pathtoprofile+'/'+id+'_'+bands[1]+'_totprofile.dat'

if not file_exist(b0file) and not file_exist(b1file) then begin
 print,'IRAC_PLOTRADPROF: No profile data, exiting'
 return
endif

read_radprof,id,bands[b0], gtot_a, gtot_mag, gtot_mag_e, gann_mu, gann_mu_e, $
	gra_cen, gdec_cen, gsemimajor, gsemiminor, gpa, gscale, $
	gtf_mag, gtf_mag_e, gaf_mag, gaf_mag_e, gmu_bg, gmu_bg_e, $
	gskyradius_in, gskyradius_out, gtf_a, gsymag, gsymag_e, gsyma, $
	pathtoprofile=pathtoprofile
read_radprof,id,bands[b1], ntot_a, ntot_mag, ntot_mag_e, nann_mu, nann_mu_e, $
	nra_cen, ndec_cen, nsemimajor, nsemiminor, npa, nscale, $
	ntf_mag, ntf_mag_e, naf_mag, naf_mag_e, nmu_bg, nmu_bg_e, $
	nskyradius_in, nskyradius_out, ntf_a, rsymag, rsymag_e, rsyma, $
	pathtoprofile=pathtoprofile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in asymptotic mag

swasymfile=pathtoprofile+'/'+id+'_'+bands[0]+'_asymptotic.dat'
readcol, swasymfile, swasy_a, swasy_mag, swasy_mag_e, swasy_int, swasy_int_e

lwasymfile=pathtoprofile+'/'+id+'_'+bands[1]+'_asymptotic.dat'
readcol, lwasymfile, lwasy_a, lwasy_mag, lwasy_mag_e, lwasy_int, lwasy_int_e

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in 50% mag

sw50file=pathtoprofile+'/'+id+'_'+bands[0]+'_50percent.dat'
readcol, sw50file, sw50_a, sw50_mag, sw50_mag_e, sw50_int, sw50_int_e

lw50file=pathtoprofile+'/'+id+'_'+bands[1]+'_50percent.dat'
readcol, lw50file, lw50_a, lw50_mag, lw50_mag_e, lw50_int, lw50_int_e

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in 80% mag

sw80file=pathtoprofile+'/'+id+'_'+bands[0]+'_80percent.dat'
readcol, sw80file, sw80_a, sw80_mag, sw80_mag_e, sw80_int, sw80_int_e

lw80file=pathtoprofile+'/'+id+'_'+bands[1]+'_80percent.dat'
readcol, lw80file, lw80_a, lw80_mag, lw80_mag_e, lw80_int, lw80_int_e

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in 90% mag

sw90file=pathtoprofile+'/'+id+'_'+bands[0]+'_90percent.dat'
readcol, sw90file, sw90_a, sw90_mag, sw90_mag_e, sw90_int, sw90_int_e

lw90file=pathtoprofile+'/'+id+'_'+bands[1]+'_90percent.dat'
readcol, lw90file, lw90_a, lw90_mag, lw90_mag_e, lw90_int, lw90_int_e

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in background data

swbgfile=pathtoprofile+id+'_'+bands[b0]+'_background.dat'
readcol, swbgfile, swbg, swbg_e, swmu_bg, swmu_bg_e, swscale, $
         swskyradius_in, swskyradius_out,/silent

lwbgfile=pathtoprofile+id+'_'+bands[b1]+'_background.dat'
readcol, lwbgfile, lwbg, lwbg_e, lwmu_bg, lwmu_bg_e, lwscale, $
         lwskyradius_in, lwskyradius_out,/silent

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in exptime

swfile=pathtofits+id+'_'+bands[b0]+'.fits*'
lwfile=pathtofits+id+'_'+bands[b1]+'.fits*'
if file_exist(swfile) then hdr1=headfits(swfile,ext=0)
if file_exist(swfile) then time1=SXPAR(hdr1,'mean_exp')>0 else time1=-999
if file_exist(lwfile) then hdr2=headfits(lwfile,ext=0)
if file_exist(lwfile) then time2=SXPAR(hdr2,'mean_exp')>0 else time2=-999


; check inputs
if gtot_a[0] eq -1 then nogband = (1 eq 1) else nogband = (1 eq 0)
if ntot_a[0] eq -1 then norband = (1 eq 1) else norband = (1 eq 0)
if nogband and norband then return
; fix missing values
if nogband then begin
	gtf_mag=!values.f_nan
	gtf_mag_e=!values.f_nan
	gaf_mag=!values.f_nan
	gaf_mag_e=!values.f_nan
endif
if norband then begin
	ntf_mag=!values.f_nan
	ntf_mag_e=!values.f_nan
	naf_mag=!values.f_nan
	naf_mag_e=!values.f_nan
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set some useful values to
; set up plots regardless of band

 if norband then begin
	a=gtot_a
	tot_mag=gtot_mag
	ann_mu=gann_mu 
	ra_cen=gra_cen[0]
	dec_cen=gdec_cen[0]
	semimajor=gsemimajor[0]
	semiminor=gsemiminor[0]
	pa=gpa[0]>0.
 endif else begin
	a=ntot_a
	tot_mag=ntot_mag
	ann_mu=nann_mu 
	ra_cen=nra_cen[0]
	dec_cen=ndec_cen[0]
	semimajor=nsemimajor[0]
	semiminor=nsemiminor[0]
	pa=npa[0]>0.
 endelse
 b=a*(semiminor/semimajor)
 r=sqrt(a*b)
 r_ap=sqrt(semimajor*semiminor)
 ator = sqrt(semiminor/semimajor)
 sw_r50 = sqrt(sw50_a*sw50_a*(semiminor/semimajor))
 sw_r80 = sqrt(sw80_a*sw80_a*(semiminor/semimajor))
 sw_r90 = sqrt(sw90_a*sw90_a*(semiminor/semimajor))
 lw_r50 = sqrt(lw50_a*lw50_a*(semiminor/semimajor))
 lw_r80 = sqrt(lw80_a*lw80_a*(semiminor/semimajor))
 lw_r90 = sqrt(lw90_a*lw90_a*(semiminor/semimajor))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in headers

if keyword_set(intfile) and file_exist(intfile) then begin
 nhdr=headfits(intfile,ext=0)
 extast,nhdr,nastr
 AD2XY, nra_cen[0] ,ndec_cen[0], nastr, nx0, ny0
 nintdata=1
endif

if keyword_set(dssfile) and file_exist(dssfile) then begin
 dss = mrdfits(dssfile, 0, dhdr)
 extast, dhdr, dssastr
 AD2XY, ra_cen,dec_cen, dssastr, dssx0, dssy0
 getrot,dhdr,rot_dss,cdelt_dss
 dssscale = abs(cdelt_dss[0])*3600.
 dssdata=1
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in jpg images

xdjpgfile =jpgpath+'/'+id+'_3p6um4p5um.jpg'

if file_test(xdjpgfile) then begin
 read_jpeg,xdjpgfile,xdjpg,/true
 xdpic=1

 ;mask the 3color image

 if keyword_set(maskimgfile) then begin
  mskimg = glga_getmask(maskimgfile,sz,astr,as_pix)
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

gdjpgfile =jpgpath+'/'+id+'_3p6um.jpg'

if file_test(gdjpgfile) then begin
 read_jpeg,gdjpgfile,gdjpg,/true
 sz=size(gdjpg,/dim)
 gdimg=bytarr(sz[1],sz[2])
 gdimg[*,*]=gdjpg[0,*,*]
 gdpic=1
endif else gdpic=0

rdjpgfile =jpgpath+'/'+id+'_4p5um.jpg'

if file_test(rdjpgfile) then begin
 read_jpeg,rdjpgfile,rdjpg,/true
 sz=size(rdjpg,/dim)
 rdimg=bytarr(sz[1],sz[2])
 rdimg[*,*]=rdjpg[0,*,*]
 rdpic=1
endif else rdpic=0

;;;;;;;;;;;;;;;;;;;;;;
; start profile plots


filename=outpath+'/'+id+'_irac_profile.ps'
ms_ps_start, filename=filename, xsize=10.5,ysize=8.,/inch,$
          /color,bit=8, xoffset=0.15, yoffset = 10.75,/landscape

!p.multi=[0,2,3]

setplotcolors
psize=0.75
chrsz=2

cs = [!orange, !red]

;growth curve
if not nogband then begin
	bad=where(gtot_mag lt 0., nbad)
	if nbad gt 0 then gtot_mag[bad] = !values.f_nan
endif else gtot_mag = !values.f_nan
if not norband then begin
	bad=where(ntot_mag lt 0., nbad)
	if nbad gt 0 then ntot_mag[bad] = !values.f_nan
endif else rtot_mag = !values.f_nan

nminmax=minmax([gtot_mag,ntot_mag], /nan)
min=nminmax[0]
max=nminmax[1]

plot,[r/60],[tot_mag],yr=[max+0.5, min-0.5],/ys,/nodata,$
    xr = [0, max(r/60.)+0.2], /xs,$
    xtit='R!de!n [arc minutes]',$
    ytit='Growth Curve [AB mag]',charsize=chrsz

if not nogband then $
 oploterror, [r/60.],[gtot_mag],[gtot_mag_e*0],[gtot_mag_e],$
  color=cs[b0],errcolor=cs[b0],psym=-sym(1, psize = psize), /nohat

if not nogband and swasy_mag gt 0 then begin
 oplot,[0,max(r/60.)+0.2],[swasy_mag,swasy_mag],color=cs[b0]
 oplot,[0,sw_r50/60],[sw50_mag,sw50_mag],color=cs[b0],linestyle=3
 oplot,[sw_r50/60,sw_r50/60],[100,sw50_mag],color=cs[b0],linestyle=3
 oplot,[0,sw_r80/60],[sw80_mag,sw80_mag],color=cs[b0],linestyle=4
 oplot,[sw_r80/60,sw_r80/60],[100,sw80_mag],color=cs[b0],linestyle=4
 oplot,[0,sw_r90/60],[sw90_mag,sw90_mag],color=cs[b0],linestyle=2
 oplot,[sw_r90/60,sw_r90/60],[100,sw90_mag],color=cs[b0],linestyle=2
endif

if not norband then $
 oploterror, [r/60.],[ntot_mag],[ntot_mag_e*0],[ntot_mag_e],$
  color=cs[b1],errcolor=cs[b1],psym=-sym(1, psize = psize), /nohat

if not norband and lwasy_mag gt 0 then begin
 oplot,[0,max(r/60.)+0.2],[lwasy_mag,lwasy_mag],color=cs[b1]
 oplot,[0,lw_r50/60],[lw50_mag,lw50_mag],color=cs[b1],linestyle=3
 oplot,[lw_r50/60,lw_r50/60],[100,lw50_mag],color=cs[b1],linestyle=3
 oplot,[0,lw_r80/60],[lw80_mag,lw80_mag],color=cs[b1],linestyle=4
 oplot,[lw_r80/60,lw_r80/60],[100,lw80_mag],color=cs[b1],linestyle=4
 oplot,[0,lw_r90/60],[lw90_mag,lw90_mag],color=cs[b1],linestyle=2
 oplot,[lw_r90/60,lw_r90/60],[100,lw90_mag],color=cs[b1],linestyle=2
endif

 oplot,[r_ap/60,r_ap/60],[max+0.5,min-0.5],line=1,thick=4,color=!red
 xx = r_ap/60. - abs(!x.crange[1]-!x.crange[0])*0.09
 yy = !y.crange[0] - abs(!y.crange[0]-!y.crange[1])*0.05
 xyouts,xx,yy,'APR',color=!red,charsiz=0.8

legend,pands[[b0,b1]],textcolors=cs[[b0,b1]],/bot,/right,box=0, $
	charsize=1.0

;annular surface brightness
if not nogband then begin
	bad=where(gann_mu lt 0., nbad)
	if nbad gt 0 then gann_mu[bad] = !values.f_nan
endif else gann_mu = !values.f_nan
if not norband then begin
	bad=where(nann_mu lt 0., nbad)
	if nbad gt 0 then nann_mu[bad] = !values.f_nan
endif else nann_mu = !values.f_nan

nminmax=minmax([gann_mu,nann_mu], /nan)
min=nminmax[0]
max=nminmax[1]

plot,[r/60],[ann_mu],yr=[max+0.5, min-0.5],/ys,/nodata,$
    xr = [0, max(r/60.)+0.2], /xs,$
    xtit='R!de!n [arc minutes]',$
    ytit=greek('mu',/append)+' [AB mag/arcsec!u2!n]',charsize=chrsz

if not nogband then $
 oploterror, [r/60.],[gann_mu],[gann_mu_e*0],[gann_mu_e],$
  color=cs[b0],errcolor=cs[b0],psym=-sym(1, psize = psize), /nohat

if not norband then $
 oploterror, [r/60.],[nann_mu],[nann_mu_e*0],[nann_mu_e],$
  color=cs[b1],errcolor=cs[b1],psym=-sym(1, psize = psize), /nohat

 oplot,[r_ap/60,r_ap/60],[max+0.5,min-0.5],line=1,thick=4,color=!red
 xx = r_ap/60. + abs(!x.crange[1]-!x.crange[0])*0.03
 yy = !y.crange[1] + abs(!y.crange[0]-!y.crange[1])*0.1
 xyouts,xx,yy,'APR',color=!red,charsiz=0.8

legend,pands[[b0,b1]],textcolors=cs[[b0,b1]],/bot,/left,box=0, $
	charsize=1.0

; total color

;fix errors (one should read n values instead)!!!!!!!!!!
zeropoint = 3631000.0
gtot_flux = zeropoint*10^(-0.4*gtot_mag)
gtot_fluxerr = gtot_mag_e*gtot_flux*alog(10)/2.5
ntot_flux = zeropoint*10^(-0.4*ntot_mag)
ntot_fluxerr = gtot_mag_e*ntot_flux*alog(10)/2.5
flxerr=sqrt(ntot_fluxerr^2 + gtot_fluxerr^2)
magerr=2.5/alog(10)*(flxerr/(ntot_flux+gtot_flux))


nn=where(finite(ntot_mag) and finite(gtot_mag), cntnfinite)

if cntnfinite gt 0 then begin
 min=min(gtot_mag - ntot_mag,/nan)
 max=max(gtot_mag - ntot_mag,/nan)
 yttl = 'Total '+pands[b0]+'-'+pands[b1]
 plot, [r/60.], [gtot_mag - ntot_mag], psym=-sym(1, psize = psize),$
     charsize=chrsz, ytit = yttl, xtit='R!de!n [arc minutes]',$
     xr = [0, max(r/60.)+0.2], /xs, $
     yr = [min-0.05, max+0.05],/ys,/nodata
 ;err=sqrt(gtot_mag_e^2 + ntot_mag_e^2)
 err=magerr

 oploterror, [r/60.],[gtot_mag - ntot_mag],[err*0],[err],$
  psym=-sym(1, psize = psize), /nohat
endif 
oplot, [0,1e4],[0,0], color=!black

; annular color

nn=where(finite(gann_mu) and finite(nann_mu),cntnfinite)

if cntnfinite gt 0 then begin
 min=min(gann_mu - nann_mu,/nan)
 max=max(gann_mu - nann_mu,/nan)
 yttl = 'Annular '+pands[b0]+'-'+pands[b1]
 plot, [r/60.], [gann_mu - nann_mu], psym=-sym(1, psize = psize),$
     charsize=chrsz, ytit = yttl, xtit='R!de!n [arc minutes]',$
     xr = [0, max(r/60.)+0.2], /xs, $
     yr = [min-0.05, max+0.05],/ys,/nodata
 err=sqrt(gann_mu_e^2 + nann_mu_e^2)
 oploterror, [r/60.],[gann_mu - nann_mu],[err*0],[err],$
  psym=-sym(1, psize = psize), /nohat
endif else begin
 plot, [r/60.], [r/60.], psym=-sym(1, psize = psize),$
     charsize=chrsz, ytit = yttl, xtit='R!de!n [arc minutes]',$
     xr = [0, max(r/60.)+0.1], /xs, $
     yr = [-1, 1],/ys,/nodata
endelse
oplot, [0,1e4],[0,0], color=!black

; text

chrsz=0.8
plot,[0,100],[0,100],xs=4,ys=4,/nodata, pos = [0.085, 0.01, 0.98, 0.33]

xyouts,0,90,'ID:  '+id,charsize=chrsz
if keyword_set(type) then $
 xyouts,0,80,'TYP:  '+strtrim(type,2),charsiz=chrsz
xyouts,0,70,pands[b0]+' EXPTIME [s]:  '+strn(time1,format='(f5.1)'),charsize=chrsz
xyouts,0,60,pands[b1]+' EXPTIME [s]:  '+strn(time2,format='(f5.1)'),charsize=chrsz

xyouts,0,50,'R.A.  [J2K]:  '+strn(ra_cen[0],format='(f12.6)'),charsize=chrsz
xyouts,0,40,'DEC [J2K]:  '+strn(dec_cen[0],format='(f12.6)'),charsize=chrsz
xyouts,0,30,'SEMIMAJOR [arcmin]:  '+strn(semimajor/60,format='(f6.2)')$
 ,charsize=chrsz
xyouts,0,20,'RATIO (a/b):  '+strn(semimajor/semiminor,format='(f6.2)'),$
 charsize=chrsz  
xyouts,0,10,'P.A. [deg]:  '+strn(pa[0],format='(f6.2)'),charsize=chrsz

euler, ra_cen, dec_cen, l, b, 1
ebv=dust_getval(l,b,/interp)
a_3p6um  = glga_getextin(ebv,'3p6um',yuan13=yuan13)
a_4p5um  = glga_getextin(ebv,'4p5um',yuan13=yuan13)

a_0 = a_3p6um
a_1 = a_4p5um

xyouts,26,80,'ASY:  m('+pands[b0]+'!do!n)='$
 +strn(swasy_mag,format='(f5.2)')+'   m('+pands[b0]+'!ddr!n)='$
 +strn(swasy_mag-a_0,format='(f5.2)')$
 +'  ('+greek("plus_minus", /append)+$
 strn(swasy_mag_e>0.01,format='(F4.2)')+')',charsize=chrsz

xyouts,26,70,'ASY:  m('+pands[b1]+'!do!n)='$
 +strn(lwasy_mag,format='(f5.2)')+'   m('+pands[b1]+'!ddr!n)='$
 +strn(lwasy_mag-a_1,format='(f5.2)')$
 +'  ('+greek("plus_minus", /append)+$
 strn(lwasy_mag_e>0.01,format='(F4.2)')+')',charsize=chrsz

xyouts,26,60,'80%:  m('+pands[b0]+'!do!n)='$
 +strn(sw80_mag,format='(f5.2)')+'   m('+pands[b0]+'!ddr!n)='$
 +strn(sw80_mag-a_0,format='(f5.2)')$
 +'  ('+greek("plus_minus", /append)+$
 strn(sw80_mag_e>0.01,format='(F4.2)')+')',charsize=chrsz

xyouts,26,50,'80%:  m('+pands[b0]+'!do!n)='$
 +strn(lw80_mag,format='(f5.2)')+'   m('+pands[b0]+'!ddr!n)='$
 +strn(lw80_mag-a_0,format='(f5.2)')$
 +'  ('+greek("plus_minus", /append)+$
 strn(lw80_mag_e>0.01,format='(F4.2)')+')',charsize=chrsz


xyouts,26,40,'APR:  m('+pands[b0]+'!do!n)='$
 +strn(gaf_mag,format='(f5.2)')+'   m('+pands[b0]+'!ddr!n)='$
 +strn(gaf_mag-a_0,format='(f5.2)')$
 +'  ('+greek("plus_minus", /append)+$
 strn(gaf_mag_e>0.01,format='(F4.2)')+')',charsize=chrsz

xyouts,26,30,'APR:  m('+pands[b1]+'!do!n)='$
 +strn(naf_mag,format='(f5.2)')+'   m('+pands[b1]+'!ddr!n)='$
 +strn(naf_mag-a_1,format='(f5.2)')$
 +'  ('+greek("plus_minus", /append)+$
 strn(naf_mag_e>0.01,format='(F4.2)')+')',charsize=chrsz

uvcolor=swasy_mag-lwasy_mag
err=sqrt(swasy_mag_e^2 + lwasy_mag_e^2)>0.01
xyouts,26,20,'ASY:  '+pands[b0]+'-'+pands[b1]+'='+$
 strn(uvcolor,format='(f5.2)')$
 +'  ('+greek("plus_minus", /append)+strn(err,$
 format='(F4.2)')+')',charsize=chrsz

uvcolor=gaf_mag-naf_mag
err=sqrt(gaf_mag_e^2 + naf_mag_e^2)>0.01
xyouts,26,10,'APR:  '+pands[b0]+'-'+pands[b1]+'='+$
 strn(uvcolor,format='(f5.2)')$
 +'  ('+greek("plus_minus", /append)+strn(err,$
 format='(F4.2)')+')',charsize=chrsz


xyouts,73,80,'GAL E(B-V):  '+strn(ebv,format='(f5.3)'),charsize=chrsz
xyouts,73,70,'GAL A!d'+pands[b0]+'!n:  '+strn(A_3p6um,format='(f5.3)'), $
	charsize=chrsz
xyouts,73,60,'GAL A!d'+pands[b1]+'!n:  '+strn(A_4p5um,format='(f5.3)'), $
	charsize=chrsz


xyouts,73,40,greek('mu',/append)+' BG '+pands[b0]+':  ' $
 +strn(swmu_bg,format='(f5.2)')$
 +'  ('+greek("plus_minus", /append)+strn(swmu_bg_e,format='(F4.2)')+')'$
 ,charsize=chrsz
xyouts,73,30,greek('mu',/append)+' BG '+pands[b1]+':  ' $
 +strn(lwmu_bg,format='(f5.2)')$
 +'  ('+greek("plus_minus", /append)+strn(lwmu_bg_e,format='(F4.2)')+')'$
 ,charsize=chrsz


xyouts,85,3,systime(),charsize=.7

ms_ps_end

;;;;;;;;;;;;;;;;;;;;;;
; convert to pdf & jpg

p=outpath+'/'
name=id+'_irac_profile'
;spawn,'convert -density 150 -rotate -90 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'convert -density 150 '+p+name+'.pdf '+p+name+'.jpg'
spawn,'convert -resize 75% '+p+name+'.jpg '+p+name+'.jpg'
spawn,'rm '+p+name+'.ps'

;;;;;;;;;;;;;;;;;;;;;;
; start image plots

filename=outpath+'/'+id+'_irac_images.ps'
ms_ps_start, filename=filename, xsize=10.5,ysize=8.,/inch,$
          /color,bit=8, xoffset=0.15, yoffset = 11,/landscape

!p.charsize=1.0
chrsz=0.75

!p.multi=[0,2,2]

minimsz = 30. / nscale	; 30 arcsec is minimum image size

; 3p6um

if gdpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(gdimg,/dim)
 delta = nskyradius_out[0]/60. * 1.01 > 0.5
 scale = nscale[0]/60.0 ; arcmins
 ximgrng = ([nx0,(nx0-sz[0])]) * scale
 yimgrng = ([0-ny0, sz[1]-ny0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, gdimg, /preserve, color=max(!d.n_colors), $
	 xran=xpltrng,yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='3.6um'

 setplotcolors

 ela=[nsemimajor/60.0,nsemiminor/60.0,0,0,float(90.+npa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=nsemiminor/nsemimajor

 eli=[nskyradius_out/60.0,ratio*nskyradius_out/60.0,0,0,float(90.+npa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[nskyradius_in/60.0,ratio*nskyradius_in/60.0,0,0,float(90.+npa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1

 el50=[sw50_a/60.0,ratio*sw50_a/60.0,0,0,float(90.+npa,0)]
 tvellipse,el50[0],el50[1],el50[2],el50[3],el50[4],/data,$
   linestyle=2,color=!green, thick=1

 el80=[sw80_a/60.0,ratio*sw80_a/60.0,0,0,float(90.+npa,0)]
 tvellipse,el80[0],el80[1],el80[2],el80[3],el80[4],/data,$
   linestyle=2,color=!blue, thick=1

 el90=[sw90_a/60.0,ratio*sw90_a/60.0,0,0,float(90.+npa,0)]
 tvellipse,el90[0],el90[1],el90[2],el90[3],el90[4],/data,$
   linestyle=2,color=!dmagenta, thick=1

 elasy=[swasy_a/60.0,ratio*swasy_a/60.0,0,0,float(90.+npa,0)]
 tvellipse,elasy[0],elasy[1],elasy[2],elasy[3],elasy[4],/data,$
   linestyle=2,color=!cyan, thick=1

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No 3p6um Image',charsize=1

endelse


; 4p5um

if rdpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(rdimg,/dim)
 delta = nskyradius_out[0]/60. * 1.01 > 0.5
 scale = nscale[0]/60.0 ; arcmins
 ximgrng = ([nx0,(nx0-sz[0])]) * scale
 yimgrng = ([0-ny0, sz[1]-ny0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, rdimg, /preserve, color=max(!d.n_colors), $
	 xran=xpltrng,yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='3.6um'

 setplotcolors

 ela=[nsemimajor/60.0,nsemiminor/60.0,0,0,float(90.+npa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=nsemiminor/nsemimajor

 eli=[nskyradius_out/60.0,ratio*nskyradius_out/60.0,0,0,float(90.+npa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[nskyradius_in/60.0,ratio*nskyradius_in/60.0,0,0,float(90.+npa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1

 el50=[lw50_a/60.0,ratio*lw50_a/60.0,0,0,float(90.+npa,0)]
 tvellipse,el50[0],el50[1],el50[2],el50[3],el50[4],/data,$
   linestyle=2,color=!green, thick=1

 el80=[lw80_a/60.0,ratio*lw80_a/60.0,0,0,float(90.+npa,0)]
 tvellipse,el80[0],el80[1],el80[2],el80[3],el50[4],/data,$
   linestyle=2,color=!blue, thick=1

 el90=[lw90_a/60.0,ratio*lw90_a/60.0,0,0,float(90.+npa,0)]
 tvellipse,el90[0],el90[1],el90[2],el90[3],el90[4],/data,$
   linestyle=2,color=!dmagenta, thick=1

 elasy=[lwasy_a/60.0,ratio*lwasy_a/60.0,0,0,float(90.+npa,0)]
 tvellipse,elasy[0],elasy[1],elasy[2],elasy[3],elasy[4],/data,$
   linestyle=2,color=!cyan, thick=1

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No 4p5um Image',charsize=1

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
   linestyle=1,color=!red, thick=3
 plots,ela[2],ela[3],psym=4,symsi=0.5,color=!red,/data

 if ~norband then ratio=nsemiminor/nsemimajor else ratio=gsemiminor/gsemimajor
 if ~norband then skyradius_out=nskyradius_out else skyradius_out=gskyradius_out
 if ~norband then skyradius_in=nskyradius_in else skyradius_in=gskyradius_in

 eli=[skyradius_out/60.0,ratio*skyradius_out/60.0,0,0,float(90.+pa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 eli=[skyradius_out/60.0,ratio*skyradius_out/60.0,0,0,float(90.+pa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[skyradius_in/60.0,ratio*skyradius_in/60.0,0,0,float(90.+pa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1

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
   linestyle=1,color=!red, thick=3
 plots,ela[2],ela[3],psym=4,symsi=0.5,color=!red,/data

 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!gray, thick=1

 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!gray, thick=1

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
name=id+'_irac_images'
;spawn,'convert -density 150 -rotate -90 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'convert -density 150 '+p+name+'.pdf '+p+name+'.jpg'
spawn,'convert -resize 75% '+p+name+'.jpg '+p+name+'.jpg'
spawn,'rm '+p+name+'.ps'

;;;;;;;;;;;;
; all done

return
end
