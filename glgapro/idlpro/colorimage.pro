

Pro colorimage




pathtoprofile = '/home/ehsan/db_esn/Luca/298D/photometry'
id = 'dwarf'
fpath = '/home/ehsan/db_esn/Luca/298D/acs/fits/'
intfile=fpath+id+'_I.fits'
outpath = '/home/ehsan/db_esn/Luca/298D/plots/'
filename=outpath+id+'_acs_image.esn.ps'
jpgpath='/home/ehsan/db_esn/Luca/298D/acs/jpg'

bands = ['V', 'I']
nband = n_elements(bands)


if not keyword_set(type) then typ = '-' else typ = strtrim(type,2)
if not keyword_set(pathtoprofile) then pathtoprofile='./'
if not keyword_set(jpgpath) then jpgpath='./'
b0 = 0
b1 = 1


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


 ghdr=headfits(intfile,ext=0)
 extast,ghdr,gastr
 AD2XY, gra_cen[0] ,gdec_cen[0], gastr, gx0, gy0
 sz = [sxpar(ghdr,'NAXIS1'),sxpar(ghdr,'NAXIS2')]
 getrot,ghdr,grot,gcdelt
 as_pix = abs(gcdelt[0])*3600.
          
          
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in jpg images

xdjpgfile =jpgpath+'/'+id+'_VI.jpg'

if file_test(xdjpgfile) then begin
 read_jpeg,xdjpgfile,xdjpg,/true
 xdpic=1
endif else xdpic=0

udjpgfile =jpgpath+'/'+id+'_V.jpg'

if file_test(udjpgfile) then begin
 
 read_jpeg,udjpgfile,udjpg,/true
 udpic=1
 udimg = udjpg
endif else udpic=0

gdjpgfile =jpgpath+'/'+id+'_I.jpg'

if file_test(udjpgfile) then begin
 
 read_jpeg,gdjpgfile,gdjpg,/true
 gdpic=1
 gdimg = gdjpg
endif else gdpic=0



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
 delta = 0.4
;  delta = gskyradius_out[0]/60. * 1.01 > 0.5
 scale=gscale[0]/60.0
 ximgrng = ([gx0,gx0-sz[0]]) * scale
 yimgrng = ([0-gy0,sz[1]-gy0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, gdimg, /preserve, charsize=chrsz, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='I'

 setplotcolors

 ela=[gsemimajor/60.0,gsemiminor/60.0,0,0,float(90.-gpa,0)]
;  tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
;    linestyle=1,color=!red, thick=3

 ratio=gsemiminor/gsemimajor

 eli=[gskyradius_out/60.0,ratio*gskyradius_out/60.0,0,0,float(90.-gpa,0)]
;  tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
;    linestyle=2,color=!black, thick=1

 elo=[gskyradius_in/60.0,ratio*gskyradius_in/60.0,0,0,float(90.-gpa,0)]
;  tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
;    linestyle=2,color=!black, thick=1

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No I Image',charsize=1

endelse

; V

if udpic then begin

 loadct,0,/silent & gamma_ct, 2 & rct

 sz = size(udimg,/dim)
 sz = sz[1:2]
 
 delta = 0.4
 
;  delta = uskyradius_out[0]/60. * 1.01 > 0.5
 scale=gscale[0]/60.0
 ximgrng = ([gx0,gx0-sz[0]]) * scale
 yimgrng = ([0-gy0,sz[1]-gy0]) * scale
 xpltrng = [delta,0-delta]
 ypltrng = [0-delta,delta]
 plotimage, udimg, /preserve, charsize=chrsz, color=max(!d.n_colors), $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='V'

 setplotcolors

 ela=[gsemimajor/60.0,gsemiminor/60.0,0,0,float(90.-gpa,0)]
;  tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
;    linestyle=1,color=!red, thick=3

 print, "V:", ela
 
 
 ratio=gsemiminor/gsemimajor

 eli=[gskyradius_out/60.0,ratio*gskyradius_out/60.0,0,0,float(90.-gpa,0)]
;  tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
;    linestyle=2,color=!black, thick=1

 elo=[gskyradius_in/60.0,ratio*gskyradius_in/60.0,0,0,float(90.-gpa,0)]
;  tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
;    linestyle=2,color=!black, thick=1

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.32],[0.5],'No V Image',charsize=1

endelse


; COLOR (make use of last derived values)

if xdpic then begin

 loadct,0,/silent

 plotimage, xdjpg, /preserve, charsize=chrsz, $
	 xran=xpltrng, yran=ypltrng, imgxrange=ximgrng, imgyrange=yimgrng, $
	 xtitle='arcmin', ytitle='arcmin', title='Composite'

 setplotcolors
 ratio=gsemiminor/gsemimajor
 ela=[gsemimajor/60.0,gsemiminor/60.0,0,0,float(90.-gpa,0)]
;  tvellipse,ela[0],ela[1],ela[2],ela[3],ela[4],/data,$
;    linestyle=1,color=!red, thick=3
;  plots,ela[2],ela[3],psym=4,symsi=2.,thick=0.5,color=!green,/data
 eli=[gskyradius_out/60.0,ratio*gskyradius_out/60.0,0,0,float(90.-gpa,0)]
;  tvellipse,eli[0],eli[1],eli[2],eli[3],eli[4],/data,$
;    linestyle=2,color=!gray, thick=1
 elo=[gskyradius_in/60.0,ratio*gskyradius_in/60.0,0,0,float(90.-gpa,0)]
;  tvellipse,elo[0],elo[1],elo[2],elo[3],elo[4],/data,$
;    linestyle=2,color=!gray, thick=1

endif else begin

 plot,[0,1],[0,1],/nodata,xs=4,ys=4
 xyouts,[0.3],[0.5],'No Composite Image',charsize=1

endelse



ms_ps_end
;;;;;;;;;;;;;;;;;;;;;;
; convert to pdf & jpg

p=outpath+'/'
name=id+'_acs_image.esn'
spawn,'convert +matte -density 196 -resample 72 -rotate -90 '+p+name+'.ps '+p+name+'.jpg'
spawn,'ps2pdf14 '+p+name+'.ps '+p+name+'.pdf'


;;;;;;;;;;;;;;;;;;;;;;
END
