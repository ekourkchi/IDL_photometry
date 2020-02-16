;+
; NAME:
;   
;   irac_radprof
;
; PURPOSE:
;
;   generate radial profile of Spitzer IRAC data
;
; CATEGORY:
;
;   radial photometry
;
; CALLING SEQUENCE:
;
;   irac_radprof, id, ra, dec, $
;                  majordiam, minordiam, pa, $
;                  intfile, $
;                  zeropoint=zeropoint, diam_units=diam_units,$
;                  ellipsefile=ellipsefile, $ 
;                  maskimgfile=maskimgfile,annuli_size=annuli_size,$
;                  extendtoskyannulus=extendtoskyannulus,$
;                  outpath=outpath, /plotsky, /verbose, status=status
;
; INPUTS:
;
;    id: scalar, string  
;    ra: scalar, decimal degree (J2K) 
;    dec:scalar, decimal degree (J2K) 
;    majordiam: scalar, specifiy units with diam_units keyword (default='arcsec')
;    minordiam: scalar, specifiy units with diam_units keyword (default='arcsec') 
;    pa: scalar, position angle in degrees East of North 
;    intfile: scalar string, IRAC intensity map in mJy
;
; OPTIONAL INPUTS/KEYWORDS:
;
;    zeropoint=zeropoint: scalar (defaults: 3631,000.0 mJy)
;    diam_units=diam_units :'pixel'or 'arcmin' or 'arcsec' (default='arcsec')
;    ellipsefile=ellipsefile : QA file of modified ellipse parameters
;    maskimgfile=maskimgfile: QA file of idl indicies of masked pixels
;    annuli_size=annuli_size: annuli increment in arcseconds
;    extendtoskyannulus=extendtoskyannulus: take profile to skyannulus-1pix
;    outpath=outpath: (default=current directory)
;    /plotsky: plots distribution of sky counts (via poissonsky) 
;    /verbose: progress messages
;
; OUTPUTS:
;
;    total profile    = ID_BAND_totprofile.dat  BAND = u, g, r, i, or z
;    annular profile  = ID_BAND_annprofile.dat  "
;    background       = ID_BAND_background.dat  "
;    ellipse params   = ID_BAND_ellipsepar.dat  "
;    best total flux  = ID_BAND_total.dat       "
;    aperture flux    = ID_BAND_aperture.dat    "
;    limiting mags    = ID_BAND_limitingmags.dat    "
;    status=status, returns 0 if all OK
;
; OPTIONAL OUTPUTS:
; COMMON BLOCKS:
; SIDE EFFECTS:
;
;    writes out 6 files
;
; RESTRICTIONS:
; PROCEDURE:
;
;   calls: poissonsky, posmeanclip, ellint, and many astrolib routines
;
; EXAMPLE:
;
;   band='3p6um'
;   filebase='~/ngc4656/data/irac/fits/NGC4656B'
;   intfile=filebase+'_'+band+'.fit*'
;   ra=191.060      
;   dec=32.2784
;   id='NGC4656B'
;   majordiam= 2*2.84809
;   minordiam= 2*1.41348
;   pa=38
;
;   resolve_routine,'irac_radprof'
;
;   irac_radprof,id,ra,dec,$
;                 majordiam,minordiam,pa,$
;                 intfile, $
;                 diam_units='arcmin',/verbose
;
; MODIFICATION HISTORY:
;
;         v1.0 M. Seibert 9/10/2008
;	  v1.0.1 D. Neill 09apr09
;-

pro irac_radprof, $
    id, ra, dec, $
    majordiam, minordiam, pa, $
    intfile, $
    zeropoint=zeropoint,verbose=verbose,diam_units=diam_units,$
    ellipsefile=ellipsefile, $
    maskimgfile=maskimgfile, annuli_size=annuli_size, $
    extendtoskyannulus=extendtoskyannulus,$
    outpath=outpath, plotsky=plotsky, status=status, yuan13=yuan13
;
; note: set yuan13 keyword to use Yuan et al. 2013 empirical extintion coeff's
;	Yuan coeff's assume SFD98 over-estimates E(B-V) by 14% (see table 2)
;
ver='1.0.1'
pre='IRAC_RADPROF: '
status = -1

if not file_exist(intfile) then begin 
   print,pre+'Can not find the input file.'
   return
endif

;
; steradians per arcsec^2
stpa = 1.D0 / ( (180.D0/!DPI)^2 * 3600.D0^2 )
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; read input files

 int=mrdfits(intfile,0,inthdr,/silent) * stpa * 1.e9	; mJy / arcsec^2
 rute=strmid(intfile,0,strpos(intfile,'fit',/reverse_search)-1)
 elist=file_search(rute+'_sigma.fit*',count=nfe)
 if nfe eq 1 then $
 	err=mrdfits(elist[0],0,errhdr,/silent) * stpa * 1.e9 $
 else	begin
 	err=-1.
	if keyword_set(verbose) then print,pre,'No err image found'
 endelse

 ; basic tests
 t=where(finite(int),nt)
 if nt le 0 then begin
	 print,pre+'Bad image.'
	 return
 endif
 if max(int(where(finite(int)))) le 0. then begin
	 print,pre+'Bad image.'
	 return
 endif
 sz=size(int,/dim)

 band=strmid(rute,strlen(rute)-5,5)

 if keyword_set(verbose) then print, pre+'processing '+band+' data'

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; set zeropoint mags to defualt values

 if not keyword_Set(zeropoint) then zeropoint = 3631000.0
 if keyword_set(verbose) then print, pre+'using zeropint flux (mJy): ',zeropoint

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; get pixel scale and center

 extast,inthdr,astr
 AD2XY, ra, dec, astr, x0, y0
 getrot,inthdr,rot,cdelt
 as_pix = abs(cdelt[0])*3600.
 if keyword_set(verbose) then print,pre+'scale (arcsec/pixel) =',as_pix
 ; convert images to unit per pixel
 int = int * as_pix^2
 err = err * as_pix^2


 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; get Galactic dust extinction mags

 euler, ra, dec, l, b, 1
 ebv=dust_getval(l,b,/interp)
 A_Galactic = glga_getextin(ebv,band,yuan13=yuan13)
 if keyword_set(verbose) then print,pre+'Gal. extinction mag = ',A_Galactic

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; convert diameters to pixel units
 ; if diam_units keyword not set then
 ;  assumes diameters are arcseconds

 if not keyword_Set(diam_units) then diam_units = 'arcsec'
 case  diam_units of
  'arcmin': scale = 60./as_pix
  'pixel':  scale = 1.0
  'arcsec': scale = 1.0/as_pix
  else: scale=999
 endcase
 
 if scale eq 999 then begin
    print,pre+'diam_units unknown. try arcmin, arcsec, or pixel.'
    return
 endif

 mjdiam=majordiam*scale	; convert to pixels
 mndiam=minordiam*scale ; convert to pixels
 diam=mjdiam
 axratio=mjdiam/mndiam
 semimajor=mjdiam*0.5
 semiminor=mndiam*0.5

 if keyword_set(verbose) and diam_units ne 'pixel' then $
  print,pre+'converted diam_units from '+diam_units+' to pixel'

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;use corrected ellipse info from QA (pixel units)

 if file_exist(ellipsefile) then begin
  readcol,ellipsefile,majdiam_as,mindiam_as,el_ra,el_dec,pa_, $
	  majdiam_px, mindiam_px, x0_,y0_, astrom_bool, el_as_pix, $
	  format='f,f,d,d,f,f,f,f,i,f',/silent
  if astrom_bool ne 1 then begin
	if keyword_set(verbose) then begin
		print,pre+'No astrometry, using default ellipse parameters'
		print,pre+'Major diameter (arcmin): ',majordiam
	endif
  endif else begin
	ad2xy,el_ra,el_dec,astr,x0_,y0_
	ra=el_ra
	dec=el_dec
  	x0=x0_[0]
  	y0=y0_[0]
  	pa=pa_[0]
  	diam=majdiam_as[0]/as_pix
  	mjdiam=majdiam_as[0]/as_pix
	mndiam=mindiam_as[0]/as_pix
	axratio=mjdiam[0]/mndiam[0]
	semimajor=mjdiam*0.5
	semiminor=mndiam*0.5
	if keyword_set(verbose) then begin
		print,pre+'Using ellipse info: ',ellipsefile
		print,pre+'Major diameter (arcmin): ',majdiam_as[0]/60.
	endif
  endelse
 endif else if keyword_set(verbose) then $
	 print,pre+'Major diameter (arcmin): ',majordiam

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; establish sky annuli sizes

 minskywid = 27.0 / as_pix	; 27 arcsec is minimum sky width
 factor=1.5
 skyinner=factor*semimajor
 skyouter=semimajor*sqrt(1+factor^2) > (skyinner+minskywid)  
 skysemimajors=[skyinner,skyouter]

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; read in mask file currently read
 ; list of idl image indicies

 maskidx=-1

 if keyword_set(maskimgfile) then begin 
  if keyword_set(verbose) then print,pre+'Using maskimgfile: ',maskimgfile
  mskimg = glga_getmask(maskimgfile,sz,astr,as_pix,/update)
  maskidx = where(mskimg ge 1)
  delvarx,mskimg,/free
 endif

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; make mask image for photometry where good pixels=1

 dist_ellipse, im, sz, x0, y0, axratio, pa
 mask=bytarr(sz)+1b
 nanpixs=where(finite(int,/nan),countnan)
 if countnan gt 0 then mask[nanpixs]=0b
 ;mask[where(finite(int,/nan) eq 1)]=0b
 if maskidx[0] ge 0 then $
 	mask[maskidx]=0b
 mask[where(im gt skyouter)]=0b
 mask[where(im lt skyinner)]=0b
 delvarx, im, /free
 in=where(mask eq 1b, sky_npix)
  

 ;;;;;;;;;;;;;;;;;;;;;;;;
 ; measure sky in annulus


 ;take first stab at sky
 ;MMM,int[in],skyint,skysig,nsky=nsky
 ;skyint=avg(int[in])
 ;skysig=stdev(int[in])
 meanclip, int[in], skyint, skysig, clipsig=10  
 skyidx=in

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; determine large scale sky variations

 xy=get_indices(skyidx,sz)
 dy=xy[1,*]-y0 & dx=xy[0,*]-x0 
 theta = asin(dx/(sqrt( dx^2 + dy^2 ))) ; -pi/2 to pi/2

 ;break into 8 sections

 s0=where(dy gt 0 and theta ge !pi/4, counts1)
 s1=where(dy gt 0 and theta ge 0 and theta lt !pi/4, counts2)   
 s2=where(dy le 0 and theta ge !pi/4, counts3)
 s3=where(dy le 0 and theta ge 0 and theta lt !pi/4, counts4) 
 s4=where(dy gt 0 and theta le -!pi/4, counts5)
 s5=where(dy gt 0 and theta le 0 and theta gt -!pi/4, counts6)   
 s6=where(dy le 0 and theta le -!pi/4, counts7)
 s7=where(dy le 0 and theta le 0 and theta gt -!pi/4, counts8) 

 s={sky:!values.f_nan,skyerr:!values.f_nan}
 s=replicate(s,8)

 for j=0,7 do begin    
  result=execute('sub = s'+strn(j))  
  if n_elements(sub) gt 19 then begin
   MMM,int[skyidx[sub]],sjsky,sjskyerr
   s[j].sky = sjsky
   s[j].skyerr = sjskyerr
  endif
 endfor

 if keyword_set(verbose) then begin
	print, pre+'Subsection sky/px and stdev:'
	forprint, strarr(n_elements(s))+pre,indgen(n_elements(s))+1, $
		s.sky,s.skyerr, format='(a,i8,2f13.8)'
 endif

 g=where(finite(s.sky) eq 1, ng)
 if ng gt 1 then begin
 	mom=moment(s.sky,/nan)
 	meannskyint = mom[0]
 	skyintspatialerr = sqrt(mom[1])
 endif else if ng eq 1 then begin
	meannskyint = s[g[0]].sky
	skyintspatialerr = s[g[0]].skyerr
 endif else if ng lt 1 then begin
	meannskyint = -999.99
	skyintspatialerr = skysig       ; set it to something
 endif

 if keyword_set(verbose) then begin
   print,pre+'Sub <sky> and sigma_<sky>[DN]: ', $
	meannskyint, skyintspatialerr, format='(a,2f13.8)'
   print,pre+'Img sky, sigma_sky   [mJy/px]: ', $
	skyint, skysig, format='(a,2f13.8)'
 endif

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;if we assume that the counting error and spatial error 
 ; are unrelated we can add them in quadrature
 ; but the spatial may have a component of the counting error
 ; will take the max error (spatial always) 
 ; will account for sky counting error later

 ;skyinterr =sqrt(skysig^2 + skyintspatialerr^2)
 ;skyinterr =max([skysig,skyintspatialerr])

 ;always make sky error large scale spatial error for IRAC galaxies
 skyinterr = skyintspatialerr 

 ;musky = -2.5*alog10( (skyinterr*5./(as_pix^2)) / zeropoint )	; 5 sigma lim
 musky = -2.5*alog10( (skyint/(as_pix^2)) / zeropoint )	
 if finite(musky) then begin
	muskyerr = 2.5/alog(10)*(skyinterr/skyint)
 endif else begin
	musky = -99.
	muskyerr = -9.
 endelse
 if keyword_set(verbose) then begin
	print, pre+'final sky/px, skyerr/pix [DN]: ', $
		skyint, skyinterr, format='(a,2f13.8)'
	print, pre+'final musky, muskyerr [mag/as^2]: ', $
		musky, muskyerr, format='(a,2f10.3)'
 endif

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; start radial profile

 ;radial profile increments (default = 6 arcsec)
 if not keyword_set(annuli_size) then annuli_size = 6.0
 inc = 1.0*annuli_size / as_pix
 max_annuli = semimajor 
 if keyword_set(extendtoskyannulus) then max_annuli = skyouter ;skyinner-1
 n_annuli=ceil(max_annuli/inc)>1
 a=(findgen(n_annuli)+1)*inc
 a = a < max_annuli > inc
 dups=where(a ge max_annuli,countdups)
 if countdups gt 1 then a=a[0:dups[0]]
 
 ;run just to get pixel counts w/out masking
 ;accounts for fractional pixels
 ellint, int, x0,y0,axratio,a,pa,$
   tot, npix,annulus, npixannulus

 ;pixel mask for photometry
 mask=intarr(sz)+1
 if maskidx[0] ge 0 then $
	 mask[maskidx]=0

 ;fix for nan's in image that are not part of mask
 nanpixs=where(finite(int,/nan),countnan)
 if countnan gt 0 then mask[nanpixs]=0

 ;measure intensity
 ellint, int, x0,y0,axratio,a,pa,$
   tot_int, npix_int, annulus_int, npixannulus_int, mask=mask

 ;pixel correction is only appropriate for annulii
 ;otherwise avg surf brightness is used each time.
 annulus_int = annulus_int * npixannulus/npixannulus_int
 tot_int = total(annulus_int,/cumulative,/nan)

 ;measure sky intensity
 ellint, fltarr(sz)*0+skyint, x0,y0,axratio,a,pa,$
   tot_skyint, npix_skyint, annulus_skyint, npixannulus_skyint, mask=mask

 ;pixel correction is only appropriate for annulii
 ;otherwise avg surf brightness is used each time.
 annulus_skyint = annulus_skyint * npixannulus/npixannulus_skyint
 tot_skyint = total(annulus_skyint,/cumulative,/nan)

 ;correct intensities for masked regions
 ;tot_int = tot_int * npix/npix_int 
 ;tot_skyint = tot_skyint * npix/npix_skyint 
 ;annulus_int = annulus_int * npixannulus/npixannulus_int
 ;annulus_skyint = annulus_skyint * npixannulus/npixannulus_skyint

 ;skysubtracted intensity values
 ss_int = tot_int - tot_skyint
 ss_annulus_int = annulus_int - annulus_skyint
 
 ;;err includes sky and other sources of noise (not true for IRAC)
 ;ADD pixel errors in quadrature
 if n_elements(err) gt 1 then begin

	good = where(err gt 0.)
	bad  = where(err le 0., nbad)
	if nbad gt 0 then err[bad] = 0.

 	;ellint, err, x0, y0, axratio, a, pa, $
   	;	tot_int_cnterr, npix_int_cnterr, $
   	;	annulus_int_cnterr, npixannulus_int_cnterr, mask=mask

  	ellint, err^2, x0, y0, axratio, a, pa, $
   		tot_int_cnterr, npix_int_cnterr, $
   		annulus_int_cnterr, npixannulus_int_cnterr, mask=mask
 
        tot_int_cnterr = sqrt(tot_int_cnterr)
        annulus_int_cnterr = sqrt(annulus_int_cnterr)

 endif else begin

	errt = (int - int) + skyinterr	; sigma/pix

 	ellint, errt, x0, y0, axratio, a, pa, $
   		tot_int_cnterr, npix_int_cnterr, $
   		annulus_int_cnterr, npixannulus_int_cnterr, mask=mask

 endelse

 ;account for root npix reduction in error
 ;tot_int_cnterr = tot_int_cnterr/sqrt(npix_int_cnterr)
 ;annulus_int_cnterr = annulus_int_cnterr/sqrt(npixannulus_int_cnterr)

 ;sky subtracted errors already included (not for IRAC)
 ;ss_int_err = tot_int_cnterr
 ;ss_annulus_int_err = annulus_int_cnterr
 ss_int_err = sqrt(tot_int_cnterr^2 + (skyinterr*npix)^2)
 ss_annulus_int_err = sqrt(annulus_int_cnterr^2 + (skyinterr*npixannulus)^2)

 ;magnitudes and errors               
 mag = -2.5*alog10(ss_int/zeropoint)
 mag_err = 2.5/alog(10)*(ss_int_err/ss_int)

 ;limiting mags
 lmag3 = -2.5*alog10(ss_int_err*3./zeropoint)
 lmag5 = -2.5*alog10(ss_int_err*5./zeropoint)
 
 ;annular surface brightness and errors   
 area=npixannulus*(as_pix)^2
 mu_annulus = -2.5*alog10( (ss_annulus_int/area) / zeropoint )
 mu_annulus_err = 2.5/alog(10)*(ss_annulus_int_err/ss_annulus_int)

 ;asymptotic mag
 dm=shift(mag,1) - mag
 dr=shift(a*as_pix,1) - a*as_pix
 n=n_elements(dm)
 good=indgen(n)
 good=where(good ge floor(0.25*n) and $ ;ss_annulus_int gt 0 and $
	 finite(dm) and mag_err gt 0,ngood)
 if ngood gt 2 then begin
	 x=dm[good]/dr[good]
	 y=mag[good]
	 erry=mag_err[good]
	 fit=linfit(x,y,measure_errors=erry,sigma=sigma)
	 asymag=fit[0]
	 asymag_err=sigma[0]
         iasy=value_to_index(mag,asymag)
         asymap=a[iasy]
	 ;dif = abs(mag-asymag)
	 ;w = where(dif eq min(dif[where(dif ge 0.)]),nw)
	 ;if nw ge 1 then $
	 ;	 asymap = a[w[0]] $
	 ;else	 asymap = a[n_elements(a)-1]
 endif else begin
	 asymag=-99.999
	 asymag_err=-9.999
	 asymap=-9.999
 endelse

 ;print apr mag and error
 good = where(a le (semimajor+1.) and ss_annulus_int gt 0, cntgood)
 if cntgood gt 0 then begin
	best=where(a eq max(a[good]))
	print,pre+"a['],mag,merr[Vega]:", $
		a[best]*as_pix/60.,mag[best],mag_err[best],form='(a,f8.3,2f9.4)'
 endif else print,'No good mags.'

if asymag gt 0 then begin
 asy_int = 1e3 * 1d23 * 10.^(-0.4*(asymag+48.6))
 i20=value_to_index(ss_int,asy_int*0.2)
 ap20=a[i20]
 i50=value_to_index(ss_int,asy_int*0.5)
 ap50=a[i50]
 i80=value_to_index(ss_int,asy_int*0.8)
 ap80=a[i80]
 i90=value_to_index(ss_int,asy_int*0.9)
 ap90=a[i90]
endif else begin
 ap20=-99
 ap50=-99
 ap80=-99
 ap90=-99
endelse


 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; write out tables

 ; write out total profile

 v=['a','mag','mag_e','int','int_e','int_e_cnt','int_e_bg','bg','pixels']
 u=['[arcsec]','[AB]','[AB]','[mJy]','[mJy]','[mJy]','[mJy]',$
    '[mJy]','[N]']


 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_totprofile.dat'
 openw,lun,outfile,/get_lun
 printf,lun,'# IRAC_RADPROF v'+ver+': total profile'
 printf,lun,'# Date: '+systime()
 printf,lun,'# ID  : '+id
 printf,lun,'# Band: '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')
 printf,lun,'# Intensity map: ',intfile,$
  format='(a17,a'+strcompress(strlen(intfile))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img     : ',maskimgfile,$
  format='(a17,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'# Additional info in ellipse parameter file and background file.'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 array_clean,-99.999,0.,50.,mag,mag_err,/nan
 for i=0,n_elements(a)-1 do $
 printf,lun,a[i]*as_pix,mag[i], mag_err[i], ss_int[i], ss_int_err[i], $
  tot_int_cnterr[i], sqrt(npix[i])*skyinterr,tot_skyint[i], npix[i],$
  format=$
  '(f12.3,x,f7.3,x,f7.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'total profile written to ',outfile
 
 ; write out annular profile

 v=['a','mu','mu_e','int','int_e','int_e_cnt','int_e_bg','bg','pixels']
 u=['[arcsec]','[mag/as^2]','[mag/as^2]','[mJy]','[mJy]','[mJy]',$
    '[mJy]','[mJy]','[N]']

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_annprofile.dat'
 openw,lun,outfile,/get_lun
 printf,lun,'# IRAC_RADPROF v'+ver+': annular profile'
 printf,lun,'# Date: '+systime()
 printf,lun,'# ID  : '+id
 printf,lun,'# Band: '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')
 printf,lun,'# Intensity map: ',intfile,$
  format='(a17,a'+strcompress(strlen(intfile))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img     : ',maskimgfile,$
  format='(a17,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'# Additional info in ellipse parameter file and background file.'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a11,x,a11,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a11,x,a11,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 array_clean,-99.999,0.,50.,mu_annulus,mu_annulus_err,/nan
 for i=0,n_elements(a)-1 do $
 printf,lun,a[i]*as_pix, mu_annulus[i], mu_annulus_err[i], $
            ss_annulus_int[i], ss_annulus_int_err[i], $
            annulus_int_cnterr[i], sqrt(npixannulus[i])*skyinterr,$
            annulus_skyint[i], npixannulus[i],$
  format=$
  '(f12.3,x,f11.3,x,f11.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'annular profile written to ',outfile

 ; write out a total flux file

 v=['a','mag','mag_e','int','int_e','int_e_cnt','int_e_bg','bg','pixels']
 u=['[arcsec]','[AB]','[AB]','[mJy]','[mJy]','[mJy]','[mJy]',$
    '[mJy]','[N]']

 sn=ss_annulus_int/ss_annulus_int_err
 good = where(sn ge 2,cntgood)
 if cntgood gt 0 then $
	 best=where(a eq max(a[good])) $
 else	 best=[n_elements(a)-1L]
 ff=-99.999

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_total.dat' 
 openw,lun,outfile,/get_lun
 printf,lun,'# IRAC_RADPROF v'+ver+': profile parameters'
 printf,lun,'# Date: '+systime()
 printf,lun,'# ID  : '+id
 printf,lun,'# Band: '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 if cntgood gt 0 then $
 printf,lun,a[best]*as_pix, mag[best], mag_err[best], $
  ss_int[best], ss_int_err[best], $
  tot_int_cnterr[best],sqrt(npix[best])*skyinterr,tot_skyint[best], npix[best],$
  format=$
  '(f12.3,x,f7.3,x,f7.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)' $
 else $
 printf,lun,ff,ff,ff,ff,ff,ff,ff,ff,ff,  format=$
  '(f12.3,x,f7.3,x,f7.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'total flux written to ',outfile


 ; write out asymtotic magnitudes

 v=['a','mag','mag_e','int','int_e','int_e_cnt','int_e_bg','bg','pixels']
 u=['[arcsec]','[AB]','[AB]','[mJy]','[mJy]','[mJy]','[mJy]',$
    '[mJy]','[N]']

 best=where(a eq asymap)
 if best[0] lt 0 then best=[n_elements(a)-1L]
 ff=-99.999

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_asymptotic.dat' 
 openw,lun,outfile,/get_lun
 printf,lun,'# IRAC_RADPROF v'+ver+': profile parameters'
 printf,lun,'# Date: '+systime()
 printf,lun,'# ID  : '+id
 printf,lun,'# Band: '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 array_clean,-99.999,0.,50.,asymag,asymag_err,/nan
 ;printf,lun,a[best]*as_pix, asymag, asymag_err, $
 printf,lun,a[best]*as_pix, mag[best], mag_err[best], $
  ss_int[best], ss_int_err[best], $
  tot_int_cnterr[best],sqrt(npix[best])*skyinterr,tot_skyint[best], npix[best],$
  format=$
  '(f12.3,x,f7.3,x,f7.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'asymptotic flux written to ',outfile

 ; write out 50% light magnitudes

 v=['a','mag','mag_e','int','int_e','int_e_cnt','int_e_bg','bg','pixels']
 u=['[arcsec]','[AB]','[AB]','[mJy]','[mJy]','[mJy]','[mJy]',$
    '[mJy]','[N]']

 best=where(a eq ap50)
 if best[0] lt 0 then best=[n_elements(a)-1L]
 ff=-99.999

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_50percent.dat' 
 openw,lun,outfile,/get_lun
 printf,lun,'# IRAC_RADPROF v'+ver+': profile parameters'
 printf,lun,'# Date: '+systime()
 printf,lun,'# ID  : '+id
 printf,lun,'# Band: '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,a[best]*as_pix, mag[best], mag_err[best], $
  ss_int[best], ss_int_err[best], $
  tot_int_cnterr[best],sqrt(npix[best])*skyinterr,tot_skyint[best], npix[best],$
  format=$
  '(f12.3,x,f7.3,x,f7.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'50% flux written to ',outfile

 ; write out 80% light magnitudes

 v=['a','mag','mag_e','int','int_e','int_e_cnt','int_e_bg','bg','pixels']
 u=['[arcsec]','[AB]','[AB]','[mJy]','[mJy]','[mJy]','[mJy]',$
    '[mJy]','[N]']

 best=where(a eq ap80)
 if best[0] lt 0 then best=[n_elements(a)-1L]
 ff=-99.999

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_80percent.dat' 
 openw,lun,outfile,/get_lun
 printf,lun,'# IRAC_RADPROF v'+ver+': profile parameters'
 printf,lun,'# Date: '+systime()
 printf,lun,'# ID  : '+id
 printf,lun,'# Band: '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,a[best]*as_pix, mag[best], mag_err[best], $
  ss_int[best], ss_int_err[best], $
  tot_int_cnterr[best],sqrt(npix[best])*skyinterr,tot_skyint[best], npix[best],$
  format=$
  '(f12.3,x,f7.3,x,f7.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'80% flux written to ',outfile

 ; write out 80% light magnitudes

 v=['a','mag','mag_e','int','int_e','int_e_cnt','int_e_bg','bg','pixels']
 u=['[arcsec]','[AB]','[AB]','[mJy]','[mJy]','[mJy]','[mJy]',$
    '[mJy]','[N]']

 best=where(a eq ap90)
 if best[0] lt 0 then best=[n_elements(a)-1L]
 ff=-99.999

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_90percent.dat' 
 openw,lun,outfile,/get_lun
 printf,lun,'# IRAC_RADPROF v'+ver+': profile parameters'
 printf,lun,'# Date: '+systime()
 printf,lun,'# ID  : '+id
 printf,lun,'# Band: '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,a[best]*as_pix, mag[best], mag_err[best], $
  ss_int[best], ss_int_err[best], $
  tot_int_cnterr[best],sqrt(npix[best])*skyinterr,tot_skyint[best], npix[best],$
  format=$
  '(f12.3,x,f7.3,x,f7.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'90% flux written to ',outfile


 ; write out 20% light magnitudes

 v=['a','mag','mag_e','int','int_e','int_e_cnt','int_e_bg','bg','pixels']
 u=['[arcsec]','[AB]','[AB]','[mJy]','[mJy]','[mJy]','[mJy]',$
    '[mJy]','[N]']

 best=where(a eq ap20)
 if best[0] lt 0 then best=[n_elements(a)-1L]
 ff=-99.999

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_20percent.dat' 
 openw,lun,outfile,/get_lun
 printf,lun,'# IRAC_RADPROF v'+ver+': profile parameters'
 printf,lun,'# Date: '+systime()
 printf,lun,'# ID  : '+id
 printf,lun,'# Band: '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,a[best]*as_pix, mag[best], mag_err[best], $
  ss_int[best], ss_int_err[best], $
  tot_int_cnterr[best],sqrt(npix[best])*skyinterr,tot_skyint[best], npix[best],$
  format=$
  '(f12.3,x,f7.3,x,f7.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'20% flux written to ',outfile



 ; write out a aperture flux file

 v=['a','mag','mag_e','int','int_e','int_e_cnt','int_e_bg','bg','pixels']
 u=['[arcsec]','[AB]','[AB]','[mJy]','[mJy]','[mJy]','[mJy]',$
    '[mJy]','[N]']

 ; get best within a pixel
 good = where(a le (semimajor+1.) and ss_annulus_int gt 0, cntgood)
 if cntgood gt 0 then best=where(a eq max(a[good]))
 ff=-99.999

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_aperture.dat' 
 openw,lun,outfile,/get_lun
 printf,lun,'# IRAC_RADPROF v'+ver+': profile parameters'
 printf,lun,'# Date: '+systime()
 printf,lun,'# ID  : '+id
 printf,lun,'# Band: '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 if cntgood gt 0 then $
 printf,lun,a[best]*as_pix, mag[best], mag_err[best], $
  ss_int[best], ss_int_err[best], $
  tot_int_cnterr[best],sqrt(npix[best])*skyinterr,tot_skyint[best], npix[best],$
  format=$
  '(f12.3,x,f7.3,x,f7.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)' $
 else $
 printf,lun,ff,ff,ff,ff,ff,ff,ff,ff,ff,  format=$
  '(f12.3,x,f7.3,x,f7.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'aperture flux written to ',outfile
 

 ; write out ellipse parameters

 v=['ra_cen','dec_cen','semimajor','semiminor','p.a.']
 u=['[j2k]','[j2k]','[arcsec]','[arcsec]','[deg]']

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_ellipsepar.dat' 
 openw,lun,outfile,/get_lun
 printf,lun,'# IRAC_RADPROF v'+ver+': profile parameters'
 printf,lun,'# Date: '+systime()
 printf,lun,'# ID  : '+id
 printf,lun,'# Band: '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')
 printf,lun,'# Intensity map: ',intfile,$
  format='(a17,a'+strcompress(strlen(intfile))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img     : ',maskimgfile,$
  format='(a17,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],$
   format='(a1,a11,x,a12,x,a12,x,a12,x,a8)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],$
   format='(a1,a11,x,a12,x,a12,x,a12,x,a8)'
 printf,lun, ra, dec, semimajor*as_pix, semiminor*as_pix, pa, $
   format='(f12.6,x,f12.6,x,f12.3,x,f12.3,x,f8.3)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'ellipse param info written to ', $
  outfile


 ; write out background  parameters

 v=['bg','bg_e','mu_bg', 'mu_bg_e', 'scale','skyradius_in','skyradius_out']
 u=['[mJy]','[mJy]','[mag/as^2]','[mag/as^2]','[as/pixel]','[as]','[as]']

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_background.dat' 
 openw,lun,outfile,/get_lun
 printf,lun,'# IRAC_RADPROF v'+ver+': background'
 printf,lun,'# Date: '+systime()
 printf,lun,'# ID  : '+id
 printf,lun,'# Band: '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')
 printf,lun,'# Intensity map: ',intfile,$
  format='(a17,a'+strcompress(strlen(intfile))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img     : ',maskimgfile,$
  format='(a17,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],$
   format='(a1,a11,x,a12,x,a12,x,a12,x,a12,x,a14,x,a14)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],$
   format='(a1,a11,x,a12,x,a12,x,a12,x,a12,x,a14,x,a14)'
 array_clean,-99.999,0.,50.,musky,muskyerr,/nan
 printf,lun, skyint, skyinterr,musky, muskyerr,as_pix,$
              skyinner*as_pix,skyouter*as_pix,$
 format='(e12.5,x,e12.5,x,f12.3,x,f12.3,x,f12.3,x,f14.3,x,f14.3)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'background info written to ',outfile

 ; write out limiting mags

 v=['a','lmag_sn3','lmag_sn5']
 u=['[arcsec]','[AB]','[AB]']

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_limitingmags.dat'
 openw,lun,outfile,/get_lun
 printf,lun,'# IRAC_RADPROF v'+ver+': limiting mags'
 printf,lun,'# Date: '+systime()
 printf,lun,'# ID  : '+id
 printf,lun,'# Band: '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')
 printf,lun,'# Intensity map: ',intfile,$
  format='(a17,a'+strcompress(strlen(intfile))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img     : ',maskimgfile,$
  format='(a17,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'# Additional info in ellipse parameter file and background file.'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2], format='(a1,a11,x,a11,x,a11)'
 printf,lun,'#',u[0],u[1],u[2], format='(a1,a11,x,a11,x,a11)'
 array_clean,-99.999,0.,50.,lmag3,lmag5,/nan
 for i=0,n_elements(a)-1 do $
 printf,lun,a[i]*as_pix,lmag3[i], lmag5[i], format='(f12.3,x,f11.3,x,f11.3)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'limiting mags written to ',outfile
 

 ;all done
 if keyword_set(verbose) then print,pre+'finis.' 

 status = 0

 return
end
