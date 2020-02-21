pro acs_plotradprof, id, pathtoprofile=pathtoprofile, type=type, $
    intfile=intfile, maskimgfile=maskimgfile, fpath=fpath, $
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

bands = ['F606W','F814W']
nband = n_elements(bands)

if not keyword_set(type) then typ = '-' else typ = strtrim(type,2)
if not keyword_set(pathtoprofile) then pathtoprofile='./'
if not keyword_set(jpgpath) then jpgpath='./'
b0 = 0
b1 = 1


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
; if norband then begin
; 	rtf_mag=!values.f_nan
; 	rtf_mag_e=!values.f_nan
; 	raf_mag=!values.f_nan
; 	raf_mag_e=!values.f_nan
; 	rsyma=!values.f_nan
; 	rsymag=!values.f_nan
; 	rsymag_e=!values.f_nan
; endif
; if noiband then begin
; 	itf_mag=!values.f_nan
; 	itf_mag_e=!values.f_nan
; 	iaf_mag=!values.f_nan
; 	iaf_mag_e=!values.f_nan
; 	isyma=!values.f_nan
; 	isymag=!values.f_nan
; 	isymag_e=!values.f_nan
; endif
; if nozband then begin
; 	ztf_mag=!values.f_nan
; 	ztf_mag_e=!values.f_nan
; 	zaf_mag=!values.f_nan
; 	zaf_mag_e=!values.f_nan
; 	zsyma=!values.f_nan
; 	zsymag=!values.f_nan
; 	zsymag_e=!values.f_nan
; endif
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
; read in headers

if keyword_set(intfile) and file_exist(intfile) then begin
 ghdr=headfits(intfile,ext=0)
 extast,ghdr,gastr
 AD2XY, gra_cen[0] ,gdec_cen[0], gastr, gx0, gy0
 sz = [sxpar(ghdr,'NAXIS1'),sxpar(ghdr,'NAXIS2')]
 getrot,ghdr,grot,gcdelt
 as_pix = abs(gcdelt[0])*3600.
endif else begin
	print,'No image data for '+id+', returning.'
	print, 'File not found: '+intfile
	return
endelse

if keyword_set(verbose) then $
	print,'acs ',form='($,a)'
;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in jpg images

xdjpgfile =jpgpath+'/'+id+'_color.jpg'

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

band_texts = ['F606W', 'F814W']
udjpgfile =jpgpath+'/'+id+'_F606W.jpg'

if file_test(udjpgfile) then begin
 
 read_jpeg,udjpgfile,udjpg,/true
 udpic=1
 
 udjpgr = udjpg[0, *, *]
 udjpgg = udjpg[1, *, *]
 udjpgb = udjpg[2, *, *]
 
 if keyword_set(fpath) then begin
    bandflag=fpath+id+'_F606W_flag.fits'
    if file_test(bandflag) then begin
      bmask = mrdfits(bandflag,0,header,/fscale,/silent)
      flagidx = where(bmask ge 1, nflagidx)
      if nflagidx gt 0 then begin
         udjpgr[flagidx] = 0
         udjpgg[flagidx] = 0
         udjpgb[flagidx] = 0
         udjpg[0, *, *] = 255-udjpgr
         udjpg[1, *, *] = 255-udjpgg
         udjpg[2, *, *] = 255-udjpgb
      endif
    endif  
 endif
 
;   sz=size(udjpg,/dim)
;   udimg=bytarr(sz[1],sz[2])
;   udimg[*,*]=udjpg[0,*,*]
  udimg = udjpg

endif else udpic=0



gdjpgfile =jpgpath+'/'+id+'_F606W.jpg'

if file_test(udjpgfile) then begin
 
 read_jpeg,gdjpgfile,gdjpg,/true
 gdpic=1
 
 gdjpgr = gdjpg[0, *, *]
 gdjpgg = gdjpg[1, *, *]
 gdjpgb = gdjpg[2, *, *]
 
 if keyword_set(fpath) then begin
    bandflag=fpath+id+'_F814W_flag.fits'
    if file_test(bandflag) then begin
      bmask = mrdfits(bandflag,0,header,/fscale,/silent)
      flagidx = where(bmask ge 1, nflagidx)
      if nflagidx gt 0 then begin
         gdjpgr[flagidx] = 0
         gdjpgg[flagidx] = 0
         gdjpgb[flagidx] = 0
         gdjpg[0, *, *] = 255-gdjpgr
         gdjpg[1, *, *] = 255-gdjpgg
         gdjpg[2, *, *] = 255-gdjpgb
      endif
    endif  
 endif
 

  gdimg = gdjpg

endif else udpic=0




;;;;;;;;;;;;;;;;;;;;;;
; start profile plots

if keyword_set(verbose) then $
	print,'profiles ',form='($,a)'

filename=outpath+'/'+id+'_acs_profile.ps'
ms_ps_start, filename=filename, xsize=10.5,ysize=8.,/inch,$
          /color,bit=8, xoffset=0.15, yoffset = 10.75,/landscape

!p.multi=[0,2,3]
!p.charsize=2

setplotcolors
psize=0.75

cs = [!blue, !red]

;growth curve
if not nouband then begin
	bad=where(utot_mag lt 0., nbad)
	if nbad gt 0 then utot_mag[bad] = !values.f_nan
endif else utot_mag = !values.f_nan
if not nogband then begin
	bad=where(gtot_mag lt 0., nbad)
	if nbad gt 0 then gtot_mag[bad] = !values.f_nan
endif else gtot_mag = !values.f_nan


; nminmax=minmax([utot_mag,gtot_mag,rtot_mag,itot_mag,ztot_mag], /nan)

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


 oplot,[r_ap,r_ap]/60.,[yr0, yr1],thick=4,color=!red

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
plot,[a/60],[ann_mu],yr=[max+0.5, min-0.5],/ys,/nodata,$
    xr = [xr0, xr1], /xs,$
    xtit='R!da!n [arcmin]',$
    ytit=greek('mu',/append)+' [AB mag/arcsec!u2!n]'
; 
if not nouband then $
 oploterror, [a/60.],[uann_mu],[uann_mu_e*0],[uann_mu_e],$
  color=cs[b0],errcolor=cs[b0],psym=-sym(1, psize = psize), /nohat

if not nogband then $
 oploterror, [a/60.],[gann_mu],[gann_mu_e*0],[gann_mu_e],$
  color=cs[b1],errcolor=cs[b1],psym=-sym(1, psize = psize), /nohat


 oplot,[r_ap,r_ap]/60.,[max+0.5,min-0.5],thick=4,color=!red
 xx = r_ap/60. - abs(!x.crange[1]-!x.crange[0])*0.09
 yy = !y.crange[1] + abs(!y.crange[0]-!y.crange[1])*0.1
 xyouts,xx,yy,'APR',color=!red,charsiz=0.8

 

legend,band_texts[[b0,b1]],textcolors=cs[[b0,b1]],/bot,/left, $
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


 oplot,[r_ap/60,r_ap/60],[max+0.5,min-0.5],thick=4,color=!red

legend,band_texts[[b0,b1]],textcolors=cs[[b0,b1]],/bot,/left, $
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




xyouts,25,80,'APR:  m('+band_texts[b0]+') = ' + $
 strn(uaf_mag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(uaf_mag_e>0.01,format='(F4.2)')

xyouts,25,70,'APR:  m('+band_texts[b1]+') = ' + $
 strn(gaf_mag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(gaf_mag_e>0.01,format='(F4.2)')


xyouts,55,80,'ASY:  m('+band_texts[b0]+') = ' + $
 strn(usymag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(usymag_e>0.01,format='(F4.2)') + $
 '    R!da!n = '+strn(usyma/60.,format='(f8.2)')+"'"

xyouts,55,70,'ASY:  m('+band_texts[b1]+') = ' + $
 strn(gsymag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(gsymag_e>0.01,format='(F4.2)') + $
 '    R!da!n = '+strn(gsyma/60.,format='(f8.2)')+"'"

opcolor=((uaf_mag-a_u)-(gaf_mag-a_g))>(-99.)<99.
err=sqrt(uaf_mag_e^2 + gaf_mag_e^2)>0.01<9.
xyouts,25,30,'APR:  ('+band_texts[b0]+'-'+band_texts[b1]+')!d0!n = '+$
 strn(opcolor,format='(f5.2)') + $
 ' '+greek("plus_minus", /append)+' '+strn(err,$
 format='(F4.2)')

; opcolor=((gaf_mag-a_g)-(raf_mag-a_r))>(-99.)<99.
; err=sqrt(gaf_mag_e^2 + raf_mag_e^2)>0.01<9.
; xyouts,25,20,'APR:  ('+bands[b1]+'-'+bands[b2]+')!d0!n = '+$
;  strn(opcolor,format='(f5.2)') + $
;  ' '+greek("plus_minus", /append)+' '+strn(err,$
;  format='(F4.2)')
; 
; 
; 
; opcolor=((usymag-a_u)-(gsymag-a_g))>(-99.)<99.
; err=sqrt(usymag_e^2 + gsymag_e^2)>0.01<9.
; xyouts,50,30,'ASY:  ('+bands[b0]+'-'+bands[b1]+')!d0!n = '+$
;  strn(opcolor,format='(f5.2)') + $
;  ' '+greek("plus_minus", /append)+' '+strn(err,$
;  format='(F4.2)')
; 
; opcolor=((gsymag-a_g)-(rsymag-a_r))>(-99.)<99.
; err=sqrt(gsymag_e^2 + rsymag_e^2)>0.01<9.
; xyouts,50,20,'ASY:  ('+bands[b1]+'-'+bands[b2]+')!d0!n = '+$
;  strn(opcolor,format='(f5.2)') + $
;  ' '+greek("plus_minus", /append)+' '+strn(err,$
;  format='(F4.2)')
; 


; xyouts,83,80,'GAL E(B-V):  '+strn(ebv,format='(f5.3)')
; xyouts,83,70,'GAL A!d'+bands[b0]+'!n:  '+strn(A_u,format='(f5.3)')
; xyouts,83,60,'GAL A!d'+bands[b1]+'!n:  '+strn(A_g,format='(f5.3)')


; ; xyouts,83,20,greek('mu',/append)+'!DBG!N '+bands[b1]+':  ' + $
; ;  strn(gmu_bg>0.<99.,format='(f8.2)') + $
; ;  ' '+greek("plus_minus", /append)+' '+strn(gmu_bg_e>0.<99.,format='(F8.2)')
; ;  

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
name=id+'_acs_profile'
spawn,'convert +matte -density 196 -resample 72 -rotate -90 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'rm '+p+name+'.ps'

;;;;;;;;;;;;;;;;;;;;;;
; start image plots

if keyword_set(verbose) then $
	print,'images ',form='($,a)'

filename=outpath+'/'+id+'_acs_images.ps'
ms_ps_start, filename=filename, xsize=10.5,ysize=8.,/inch,$
          /color,bit=8, xoffset=0.15, yoffset = 11,/landscape

!p.charsize=1.0
chrsz=0.75

!p.multi=[0,2,2]

minimsz = 30. / gscale	; 30 arcsec is minimum image size

; I

if gdpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(gdimg,/dim)
 sz = sz[1:2]
 
 delta = gskyradius_out[0]/60. * 1.01 > 0.5
 scale=gscale[0]/60.0
 ximgrng = ([gx0,gx0-sz[0]]) * scale
 yimgrng = ([0-gy0,sz[1]-gy0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, gdimg, /preserve, charsize=chrsz, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='F814W'

 setplotcolors

 ela=[gsemimajor/60.0,gsemiminor/60.0,0,0,float(90.-gpa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 ratio=gsemiminor/gsemimajor

 eli=[gskyradius_out/60.0,ratio*gskyradius_out/60.0,0,0,float(90.-gpa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[gskyradius_in/60.0,ratio*gskyradius_in/60.0,0,0,float(90.-gpa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No F814W Image',charsize=1

endelse

; V

if udpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(udimg,/dim)
 sz = sz[1:2]
 
 
 delta = uskyradius_out[0]/60. * 1.01 > 0.5
 scale=gscale[0]/60.0
 ximgrng = ([gx0,gx0-sz[0]]) * scale
 yimgrng = ([0-gy0,sz[1]-gy0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, udimg, /preserve, charsize=chrsz, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='F606W'

 setplotcolors

 ela=[gsemimajor/60.0,gsemiminor/60.0,0,0,float(90.-gpa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3

 print, "V:", ela
 
 
 ratio=gsemiminor/gsemimajor

 eli=[gskyradius_out/60.0,ratio*gskyradius_out/60.0,0,0,float(90.-gpa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1

 elo=[gskyradius_in/60.0,ratio*gskyradius_in/60.0,0,0,float(90.-gpa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No F606W Image',charsize=1

endelse



; COLOR (make use of last derived values)

if xdpic then begin

 loadct,0,/silent

 plotimage, xdjpg, /preserve, charsize=chrsz, $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='Masked Composite'

 setplotcolors
 ratio=gsemiminor/gsemimajor
 ela=[gsemimajor/60.0,gsemiminor/60.0,0,0,float(90.-gpa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3
 plots,ela[2],ela[3],psym=4,symsi=2.,thick=0.5,color=!green,/data
 eli=[gskyradius_out/60.0,ratio*gskyradius_out/60.0,0,0,float(90.-gpa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!gray, thick=1
 elo=[gskyradius_in/60.0,ratio*gskyradius_in/60.0,0,0,float(90.-gpa,0)]
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
name=id+'_acs_images'
; spawn,'convert +matte -density 196 -resample 72 -rotate -90 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
; spawn,'rm '+p+name+'.ps'

if keyword_set(verbose) then $
	print,'Done.'

;;;;;;;;;;;;
; all done

return
end
