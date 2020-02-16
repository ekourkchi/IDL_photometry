pro twomass_plotradprof, id, pathtoprofile=pathtoprofile, type=type, $
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
	pathtoprofile=pathtoprofile, annuli_size=annuli_size
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
if nojband then begin
	print,'No primary image photometry for '+id+', returning.'
	return
endif
if nojband and nohband and nokband then return
; fix missing values
if nojband then begin
	jtf_mag=!values.f_nan
	jtf_mag_e=!values.f_nan
	jaf_mag=!values.f_nan
	jaf_mag_e=!values.f_nan
	jsyma=!values.f_nan
	jsymag=!values.f_nan
	jsymag_e=!values.f_nan
endif
if nohband then begin
	htf_mag=!values.f_nan
	htf_mag_e=!values.f_nan
	haf_mag=!values.f_nan
	haf_mag_e=!values.f_nan
	hsyma=!values.f_nan
	hsymag=!values.f_nan
	hsymag_e=!values.f_nan
endif
if nokband then begin
	ktf_mag=!values.f_nan
	ktf_mag_e=!values.f_nan
	kaf_mag=!values.f_nan
	kaf_mag_e=!values.f_nan
	ksyma=!values.f_nan
	ksymag=!values.f_nan
	ksymag_e=!values.f_nan
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set some useful values to
; set up plots regardless of band

a=jtot_a
tot_mag=jtot_mag
ann_mu=jann_mu 
ra_cen=jra_cen[0]
dec_cen=jdec_cen[0]
semimajor=jsemimajor[0]
semiminor=jsemiminor[0]
pa=jpa[0]>0.
r_ap=semimajor
r_ski=jskyradius_in
r_sko=jskyradius_out

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in headers

if keyword_set(intfile) and file_exist(intfile) then begin
	jhdr=headfits(intfile,ext=0)
	extast,jhdr,jastr
	AD2XY, jra_cen[0] ,jdec_cen[0], jastr, jx0, jy0
	sz = [sxpar(jhdr,'NAXIS1'),sxpar(jhdr,'NAXIS2')]
	getrot,jhdr,jrot,jcdelt
	as_pix = abs(jcdelt[0])*3600.
endif else begin
	print,'No image data for '+id+', returning.'
	return
endelse

if keyword_set(verbose) then $
	print,'2mass ',form='($,a)'
;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in jpg images

xdjpgfile =jpgpath+'/'+id+'_jhk.jpg'

if file_test(xdjpgfile) then begin
 read_jpeg,xdjpgfile,xdjpg,/true
 xdpic=1

 ;mask the 3color image

 if keyword_set(maskimgfile) then begin
  if keyword_set(verbose) then $
	 print,'masking ... ',form='($,a)'
  mskimg = glga_getmask(maskimgfile,sz,jastr,as_pix)
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

jdjpgfile =jpgpath+'/'+id+'_j.jpg'

if file_test(jdjpgfile) then begin
 read_jpeg,jdjpgfile,jdjpg,/true
 sz=size(jdjpg,/dim)
 jdimg=bytarr(sz[1],sz[2])
 jdimg[*,*]=jdjpg[0,*,*]
 jdpic=1
endif else jdpic=0

hdjpgfile =jpgpath+'/'+id+'_h.jpg'

if file_test(hdjpgfile) then begin
 read_jpeg,hdjpgfile,hdjpg,/true
 sz=size(hdjpg,/dim)
 hdimg=bytarr(sz[1],sz[2])
 hdimg[*,*]=hdjpg[0,*,*]
 hdpic=1
endif else hdpic=0

kdjpgfile =jpgpath+'/'+id+'_k.jpg'

if file_test(kdjpgfile) then begin
 read_jpeg,kdjpgfile,kdjpg,/true
 sz=size(kdjpg,/dim)
 kdimg=bytarr(sz[1],sz[2])
 kdimg[*,*]=kdjpg[0,*,*]
 kdpic=1
endif else kdpic=0

;;;;;;;;;;;;;;;;;;;;;;
; start profile plots

if keyword_set(verbose) then $
	print,'profiles ',form='($,a)'

filename=outpath+'/'+id+'_2mass_profile.ps'
ms_ps_start, filename=filename, xsize=10.5,ysize=8.,/inch,$
          /color,bit=8, xoffset=0.15, yoffset = 10.75,/landscape

!p.multi=[0,2,3]
!p.charsize=2

setplotcolors
psize=0.75

cs = [!green, !orange, !red]

;growth curve
if not nojband then begin
	bad=where(jtot_mag lt 0., nbad)
	if nbad gt 0 then jtot_mag[bad] = !values.f_nan
endif else jtot_mag = !values.f_nan
if not nohband then begin
	bad=where(htot_mag lt 0., nbad)
	if nbad gt 0 then htot_mag[bad] = !values.f_nan
endif else rtot_mag = !values.f_nan
if not nokband then begin
	bad=where(ktot_mag lt 0., nbad)
	if nbad gt 0 then ktot_mag[bad] = !values.f_nan
endif else ktot_mag = !values.f_nan

nminmax=minmax([jtot_mag,htot_mag,ktot_mag], /nan)
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
    xtitle='R!da!n [arcmin]',ytitle='Growth Curve [AB mag]'
;
; make sure we are including the sky annulus in the plot
if a[n_elements(a)-1] gt r_ski then begin
	oplot,[r_ski,r_sko]/60.,[min-ydel*0.1,min-ydel*0.1], color=!blue
	xyouts,r_ski/60.,min-ydel*0.12,'SKY ANN', $
		color=!blue, charsize=0.5
endif

if not nojband then begin
 oploterror, [a/60.],[jtot_mag],[jtot_mag_e*0],[jtot_mag_e],$
  color=cs[b0],errcolor=cs[b0],psym=-sym(1, psize = psize), /nohat
 oplot,!x.crange,[jsymag,jsymag],linesty=1,color=cs[b0]
 oplot,[jsyma,jsyma]/60.,[yr0,jsymag],linesty=2,color=cs[b0],thick=4
endif

if not nohband then begin
 oploterror, [a/60.],[htot_mag],[htot_mag_e*0],[htot_mag_e],$
  color=cs[b1],errcolor=cs[b1],psym=-sym(1, psize = psize), /nohat
 oplot,!x.crange,[hsymag,hsymag],linesty=1,color=cs[b1]
 oplot,[hsyma,hsyma]/60.,[yr0,hsymag],linesty=2,color=cs[b1],thick=4
endif

if not nokband then begin
 oploterror, [a/60.],[ktot_mag],[ktot_mag_e*0],[ktot_mag_e],$
  color=cs[b2],errcolor=cs[b2],psym=-sym(1, psize = psize), /nohat
 oplot,!x.crange,[ksymag,ksymag],linesty=1,color=cs[b2]
 oplot,[ksyma,ksyma]/60.,[yr0,ksymag],linesty=2,color=cs[b2],thick=4
endif

oplot,[r_ap,r_ap]/60.,[yr0, yr1],thick=4,color=!red
xx = r_ap/60. - abs(!x.crange[1]-!x.crange[0])*0.09
yy = !y.crange[0] - abs(!y.crange[0]-!y.crange[1])*0.05
xyouts,xx,yy,'APR',color=!red,charsiz=0.8

; growth curve log x-axis

plot,[a/60],[tot_mag],yr=[yr0, yr1],/ys,/nodata,$
    xr = [xlr0, xr1], /xs, /xlog, $
    xtitle='R!da!n [arcmin]',ytitle='Growth Curve [AB mag]'

xyouts,0.1,yr1+ydel*0.14,'ASY',color=!blue,charsiz=0.8

if not nojband then begin
 oploterror, [a/60.],[jtot_mag],[jtot_mag_e*0],[jtot_mag_e],$
  color=cs[b0],errcolor=cs[b0],psym=-sym(1, psize = psize), /nohat
 oplot,[xlr0,1.e9],[jsymag,jsymag],linesty=1,color=cs[b0]
endif

if not nohband then begin
 oploterror, [a/60.],[htot_mag],[htot_mag_e*0],[htot_mag_e],$
  color=cs[b1],errcolor=cs[b1],psym=-sym(1, psize = psize), /nohat
 oplot,[xlr0,1.e9],[hsymag,hsymag],linesty=1,color=cs[b1]
endif

if not nokband then begin
 oploterror, [a/60.],[ktot_mag],[ktot_mag_e*0],[ktot_mag_e],$
  color=cs[b2],errcolor=cs[b2],psym=-sym(1, psize = psize), /nohat
 oplot,[xlr0,1.e9],[ksymag,ksymag],linesty=1,color=cs[b2]
endif

 oplot,[r_ap,r_ap]/60.,[yr0, yr1],thick=4,color=!red

;annular surface brightness
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

nminmax=minmax([jann_mu,hann_mu,kann_mu], /nan)
min=nminmax[0]
max=nminmax[1]
; annular surface brightness linear x-axis
plot,[a/60],[ann_mu],yr=[max+0.5, min-0.5],/ys,/nodata,$
    xr = [xr0, xr1], /xs,$
    xtitle='R!da!n [arcmin]',$
    ytitle=greek('mu',/append)+' [AB mag/arcsec!u2!n]'

if not nojband then $
 oploterror, [a/60.],[jann_mu],[jann_mu_e*0],[jann_mu_e],$
  color=cs[b0],errcolor=cs[b0],psym=-sym(1, psize = psize), /nohat

if not nohband then $
 oploterror, [a/60.],[hann_mu],[hann_mu_e*0],[hann_mu_e],$
  color=cs[b1],errcolor=cs[b1],psym=-sym(1, psize = psize), /nohat

if not nokband then $
 oploterror, [a/60.],[kann_mu],[kann_mu_e*0],[kann_mu_e],$
  color=cs[b2],errcolor=cs[b2],psym=-sym(1, psize = psize), /nohat

 oplot,[r_ap,r_ap]/60.,[max+0.5,min-0.5],thick=4,color=!red
 xx = r_ap/60. - abs(!x.crange[1]-!x.crange[0])*0.09
 yy = !y.crange[1] + abs(!y.crange[0]-!y.crange[1])*0.1
 xyouts,xx,yy,'APR',color=!red,charsiz=0.8

legend,bands[[b0,b1,b2]],textcolors=cs[[b0,b1,b2]],/bot,/left,box=0, $
	charsize=1.0

; annular surface brightness log x-axis
plot,[a/60],[ann_mu],yr=[max+0.5, min-0.5],/ys,/nodata,$
    xr = [xlr0, xr1], /xs, /xlog, $
    xtitle='R!da!n [arcmin]',$
    ytitle=greek('mu',/append)+' [AB mag/arcsec!u2!n]'

if not nojband then $
 oploterror, [a/60.],[jann_mu],[jann_mu_e*0],[jann_mu_e],$
  color=cs[b0],errcolor=cs[b0],psym=-sym(1, psize = psize), /nohat

if not nohband then $
 oploterror, [a/60.],[hann_mu],[hann_mu_e*0],[hann_mu_e],$
  color=cs[b1],errcolor=cs[b1],psym=-sym(1, psize = psize), /nohat

if not nokband then $
 oploterror, [a/60.],[kann_mu],[kann_mu_e*0],[kann_mu_e],$
  color=cs[b2],errcolor=cs[b2],psym=-sym(1, psize = psize), /nohat

 oplot,[r_ap,r_ap]/60.,[max+0.5,min-0.5],thick=4,color=!red

legend,bands[[b0,b1,b2]],textcolors=cs[[b0,b1,b2]],/bot,/left,box=0, $
	charsize=1.0

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
a_j  = glga_getextin(ebv,'j',yuan13=yuan13)
a_h  = glga_getextin(ebv,'h',yuan13=yuan13)
a_k  = glga_getextin(ebv,'k',yuan13=yuan13)

xyouts,28,80,'APR:  m('+bands[b0]+') = ' + $
 strn(jaf_mag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(jaf_mag_e>0.01,format='(F4.2)')

xyouts,28,70,'APR:  m('+bands[b1]+') = ' + $
 strn(haf_mag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(haf_mag_e>0.01,format='(F4.2)')

xyouts,28,60,'APR:  m('+bands[b2]+') = ' + $
 strn(kaf_mag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(kaf_mag_e>0.01,format='(F4.2)')

nircolor=((jaf_mag-a_j)-(kaf_mag-a_k))>(-99.)<99.
err=sqrt(jaf_mag_e^2 + kaf_mag_e^2)>0.01<9.
xyouts,28,50,'APR:  ('+bands[b0]+'-'+bands[b2]+')!d0!n = '+$
 strn(nircolor,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(err,format='(F4.2)')

xyouts,28,40,'ASY:  m('+bands[b0]+') = ' + $
 strn(jsymag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(jsymag_e>0.01,format='(F4.2)') + $
 ' R!da!n = '+strn(jsyma/60.,format='(f8.2)')+"'"

xyouts,28,30,'ASY:  m('+bands[b1]+') = ' + $
 strn(hsymag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(hsymag_e>0.01,format='(F4.2)') + $
 ' R!da!n = '+strn(hsyma/60.,format='(f8.2)')+"'"

xyouts,28,20,'ASY:  m('+bands[b2]+') = ' + $
 strn(ksymag>0.<99.,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(ksymag_e>0.01,format='(F4.2)') + $
 ' R!da!n = '+strn(ksyma/60.,format='(f8.2)')+"'"

nircolor=((jsymag-a_j)-(ksymag-a_k))>(-99.)<99.
err=sqrt(jsymag_e^2 + ksymag_e^2)>0.01<9.
xyouts,28,10,'ASY:  ('+bands[b0]+'-'+bands[b2]+')!d0!n = '+$
 strn(nircolor,format='(f5.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(err,format='(F4.2)')

xyouts,75,80,'GAL E(B-V):  '+strn(ebv,format='(f5.3)')
xyouts,75,70,'GAL A!d'+bands[b0]+'!n:  '+strn(a_j,format='(f5.3)')
xyouts,75,60,'GAL A!d'+bands[b1]+'!n:  '+strn(a_h,format='(f5.3)')
xyouts,75,50,'GAL A!d'+bands[b2]+'!n:  '+strn(a_k,format='(f5.3)')

xyouts,75,40,greek('mu',/append)+' BG '+bands[b0]+':  ' + $
 strn(jmu_bg>0.<99.,format='(f8.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(jmu_bg_e>0.01,format='(F8.2)')

xyouts,75,30,greek('mu',/append)+' BG '+bands[b1]+':  ' + $
 strn(hmu_bg>0.<99.,format='(f8.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(hmu_bg_e>0.01,format='(F8.2)')

xyouts,75,20,greek('mu',/append)+' BG '+bands[b2]+':  ' + $
 strn(kmu_bg>0.<99.,format='(f8.2)')+' '+greek("plus_minus", /append)+$
 ' '+strn(kmu_bg_e>0.01,format='(F8.2)')


xyouts,0.96,0.077,systime(),charsize=.7,/norm,ori=90.

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
name=id+'_2mass_profile'
spawn,'convert +matte -density 196 -resample 72 -rotate -90 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'rm '+p+name+'.ps'

;;;;;;;;;;;;;;;;;;;;;;
; start image plots

if keyword_set(verbose) then $
	print,'images ',form='($,a)'

filename=outpath+'/'+id+'_2mass_images.ps'
ms_ps_start, filename=filename, xsize=10.5,ysize=8.,/inch,$
          /color,bit=8, xoffset=0.15, yoffset = 11,/landscape

!p.charsize=1.0
chrsz=0.75

!p.multi=[0,2,2]

minimsz = 30. / jscale	; 30 arcsec is minimum image size

; j

if jdpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(jdimg,/dim)
 delta = jskyradius_out[0]/60. * 1.01 > 0.5
 scale = jscale[0]/60.0 ; arcmins
 ximgrng = ([jx0,jx0-sz[0]]) * scale
 yimgrng = ([0-jy0, sz[1]-jy0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, jdimg, /preserve, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='j'

 setplotcolors

 ela=[jsemimajor/60.0,jsemiminor/60.0,0,0,float(90.-jpa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3, noclip=0

 ratio=jsemiminor/jsemimajor

 eli=[jskyradius_out/60.0,ratio*jskyradius_out/60.0,0,0,float(90.-jpa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1, noclip=0

 elo=[jskyradius_in/60.0,ratio*jskyradius_in/60.0,0,0,float(90.-jpa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1, noclip=0

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No j Image',charsize=1

endelse


; h

if hdpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(hdimg,/dim)
 delta = hskyradius_out[0]/60. * 1.01 > 0.5
 scale = hscale[0]/60.0 ; arcmins
 ximgrng = ([jx0,jx0-sz[0]]) * scale
 yimgrng = ([0-jy0, sz[1]-jy0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, hdimg, /preserve, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='h'

 setplotcolors

 ela=[hsemimajor/60.0,hsemiminor/60.0,0,0,float(90.-hpa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3, noclip=0

 ratio=hsemiminor/hsemimajor

 eli=[hskyradius_out/60.0,ratio*hskyradius_out/60.0,0,0,float(90.-hpa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1, noclip=0

 elo=[hskyradius_in/60.0,ratio*hskyradius_in/60.0,0,0,float(90.-hpa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1, noclip=0

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No h Image',charsize=1

endelse


; k

if kdpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(kdimg,/dim)
 delta = kskyradius_out[0]/60. * 1.01 > 0.5
 scale = kscale[0]/60.0 ; arcmins
 ximgrng = ([jx0,jx0-sz[0]]) * scale
 yimgrng = ([0-jy0, sz[1]-jy0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, kdimg, /preserve, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='k'

 setplotcolors

 ela=[ksemimajor/60.0,ksemiminor/60.0,0,0,float(90.-kpa,0)]
 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3, noclip=0

 ratio=ksemiminor/ksemimajor

 eli=[kskyradius_out/60.0,ratio*kskyradius_out/60.0,0,0,float(90.-kpa,0)]
 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!black, thick=1, noclip=0

 elo=[kskyradius_in/60.0,ratio*kskyradius_in/60.0,0,0,float(90.-kpa,0)]
 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!black, thick=1, noclip=0

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No k Image',charsize=1

endelse


; COLOR (make use of last derived values)

if xdpic then begin

 loadct,0,/silent

 plotimage, xdjpg, /preserve, $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='Masked Composite'

 setplotcolors

 tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
   linestyle=1,color=!red, thick=3, noclip=0
 plots,ela[2],ela[3],psym=4,symsi=2.,thick=0.5,color=!green,/data

 tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
   linestyle=2,color=!gray, thick=1, noclip=0

 tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
   linestyle=2,color=!gray, thick=1, noclip=0

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
name=id+'_2mass_images'
spawn,'convert +matte -density 196 -resample 72 -rotate -90 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'
spawn,'rm '+p+name+'.ps'

if keyword_set(verbose) then $
	print,'Done.'

;;;;;;;;;;;;
; all done

return
end
