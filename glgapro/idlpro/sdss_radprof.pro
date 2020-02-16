; $Id: sdss_radprof.pro,v 1.48 2014/03/03 19:23:35 neill Exp $
;+
; NAME:
;   
;   sdss_radprof
;
; PURPOSE:
;
;   generate radial profile of SDSS data
;
; CATEGORY:
;
;   radial photometry
;
; CALLING SEQUENCE:
;
;   sdss_radprof, id, ra, dec, $
;                  majordiam, minordiam, pa, $
;                  intfile, $
;                  zeropoint=zeropoint, diam_units=diam_units,$
;                  ellipsefile=ellipsefile, $ 
;                  maskimgfile=maskimgfile,annuli_size=annuli_size,$
;                  extendtoskyannulus=extendtoskyannulus,$
;                  outpath=outpath, /verbose, status=status
;
; INPUTS:
;
;    id: scalar, string  
;    ra: scalar, decimal degree (J2K) 
;    dec:scalar, decimal degree (J2K) 
;    majordiam: scalar, specifiy units with diam_units keyword (default='arcsec')
;    minordiam: scalar, specifiy units with diam_units keyword (default='arcsec') 
;    pa: scalar, position angle in degrees East of North 
;    intfile: scalar string, SDSS intensity map in mJy
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
;    /verbose: progress messages
;
; OUTPUTS:
;
;    total profile    = ID_BAND_totprofile.dat   BAND = u, g, r, i, or z
;    annular profile  = ID_BAND_annprofile.dat   "
;    background       = ID_BAND_background.dat   "
;    ellipse params   = ID_BAND_ellipsepar.dat   "
;    best total flux  = ID_BAND_total.dat        "
;    aperture flux    = ID_BAND_aperture.dat     "
;    asymptotic flux  = ID_BAND_asymptotic.dat   "
;    limiting mags    = ID_BAND_limitingmags.dat "
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
;   calls: ellint, and many astrolib routines
;
; EXAMPLE:
;
;   band='g'
;   filebase='~/ngc4656/data/sdss/fits/NGC4656B'
;   intfile=filebase+'_'+band+'.fit*'
;   ra=191.060      
;   dec=32.2784
;   id='NGC4656B'
;   majordiam= 2*2.84809
;   minordiam= 2*1.41348
;   pa=38
;
;   resolve_routine,'sdss_radprof'
;
;   sdss_radprof,id,ra,dec,$
;                 majordiam,minordiam,pa,$
;                 intfile, $
;                 diam_units='arcmin',/verbose
;
; MODIFICATION HISTORY:
;
;	2008-SEP-09	M. Seibert
;	2009-APR-09	D. Neill
;	2012-MAY-31	D. Neill - fixed noise model
;-

pro sdss_radprof, $
    id, ra, dec, $
    majordiam, minordiam, pa, $
    intfile, $
    zeropoint=zeropoint,verbose=verbose,diam_units=diam_units,$
    ellipsefile=ellipsefile, $
    maskimgfile=maskimgfile, bandmask=bandmask, annuli_size=annuli_size, $
    extendtoskyannulus=extendtoskyannulus, auxbase=auxbase, $
    outpath=outpath, status=status, yuan13=yuan13
;
; note: set yuan13 keyword to use Yuan et al. 2013 empirical extintion coeff's
;	Yuan coeff's assume SFD98 over-estimates E(B-V) by 14% (see table 2)
;
ver=repstr('$Revision: 1.48 $ $Date: 2014/03/03 19:23:35 $','$','')
pre='SDSS_RADPROF: '
status = -1

int_fi = file_info(intfile)
if not int_fi.exists then begin 
   print,pre+'Cannot find the input file.'
   return
endif

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; read input files

 int=mrdfits(int_fi.name,0,inthdr,/fscale,/silent)
 rute=strmid(intfile,0,strpos(intfile,'fit',/reverse_search)-1)
 elist=file_search(rute+'_ivar.fit*',count=nfe)
 
 bbb = strpos(intfile,'.fit')
 aaa = strpos(intfile,'_', /REVERSE_SEARCH)
 band_f = strmid(intfile,aaa+1,bbb-aaa-1) 
 
 if nfe eq 1 then $
	var=1./( mrdfits(elist[0],0,errhdr,/fscale,/silent) ) $
 else	begin
	var=-1.
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

 band=strmid(rute,strlen(rute)-1,1)

 if keyword_set(verbose) then print, pre+'processing '+band+' data'

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; set zeropoint to default values

 if not keyword_set(zeropoint) then zeropoint = 3631000.0
 if keyword_set(verbose) then print, pre+'using zeropoint flux (mJy): ',zeropoint

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; get pixel scale and center

 extast,inthdr,astr
 AD2XY, ra, dec, astr, x0, y0
 getrot,inthdr,rot,cdelt
 as_pix = abs(cdelt[0])*3600.
 if keyword_set(verbose) then print,pre+'scale (arcsec/pixel) =',as_pix


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

 if not keyword_set(diam_units) then diam_units = 'arcsec'
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

 minskywid = 27. / as_pix	; 27 arcsec is minimum sky width
 factor= 1.5     ;  1.5
 skyinner=factor*semimajor
 skyouter=semimajor*sqrt(1+factor^2) > (skyinner+minskywid)  
 skysemimajors=[skyinner,skyouter]

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; read in mask file currently read
 ; list of idl image indicies

 maskidx=-1
 maskidx_band=-1

 if keyword_set(maskimgfile) then begin 
  if keyword_set(verbose) then print,pre+'Using maskimgfile: '+maskimgfile
  mskimg = glga_getmask(maskimgfile,sz,astr,as_pix,/update)
  maskidx = where(mskimg ge 1)
  delvarx,mskimg,/free
 endif

 if keyword_set(bandmask) and file_test(bandmask) then begin 
  if keyword_set(verbose) then print,pre+'Using maskimgfile: '+bandmask
  mskimg2 = glga_getmask(bandmask,sz,astr,as_pix,/update)
  maskidx_band = where(mskimg2 ge 1)
  delvarx,mskimg2,/free
 endif  
 
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; make mask image for photometry where good pixels=1

 dist_ellipse, im, sz, x0, y0, axratio, pa
 mask=bytarr(sz)+1b
 if maskidx[0] ge 0 then $
 	mask[maskidx]=0b
 
 ; taking care of masking of each individual mask
 if maskidx_band[0] ge 0 then $   
 	mask[maskidx_band]=0b 	
 	
 loco=where(im gt skyouter,nloco)
 if nloco gt 0 then mask[loco]=0b
 loci=where(im lt skyinner,nloci)
 if nloci gt 0 then mask[loci]=0b
 delvarx, im, /free
 in=where(mask eq 1b, sky_npix)
 if sky_npix le 0 then begin
	 sky_npix = n_elements(int)
	 in = lindgen(sky_npix)	; use whole image
 endif
  

 ; @ehsan: Taking care of NaN values
 indx_nan = where(finite(int, /NAN), no_nan)
 if no_nan gt 0 then begin 
    
    mask[indx_nan]=0b
 
    bmap = bytarr(sz)
    bmap[indx_nan]=1b
    badmap = auxbase+'_'+band_f+'_badmap.fits.gz'
    writefits,badmap,bmap,/compress
    
 endif
  
  
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; measure sky and sky uncertainty in annulus


 ;take first stab at sky
 ;MMM,int[in],skyint,skysig,nsky=nsky
 meanclip,int[in],skyint,skysig,clipsig=10
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
	skyintspatialerr = skysig	; set it to something
 endif

 if keyword_set(verbose) then begin
	print,pre+'Sub <sky> and sigma_<sky>   [mJy/px]: ', $
		meannskyint, skyintspatialerr, format='(a,2f13.8)'
	print,pre+'Img sky, sigma_sky          [mJy/px]: ', $
		skyint, skysig, format='(a,2f13.8)'
 endif
 
;  skyint = 1.01*skyint     ; ehsan-edit
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
 ;if we assume that the counting error and spatial error 
 ; are unrelated we can add them in quadrature
 ; but the spatial may have a component of the counting error
 ; will take the max error (spatial always) 
 ; will account for sky counting error later

 ;skyinterr =sqrt(skysig^2 + skyintspatialerr^2)
 skyinterr =max([skysig,skyintspatialerr])
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
 vertices = auxbase+'_stv_background_vertices.dat'
 background = auxbase+'_stv_background_output.dat'   
  
 backindx = -1
 if file_test(background) then begin
    readcol,background,rndx,form='l',/silent  
    backindx = rndx
 
 mask=intarr(sz)+1
 if maskidx[0] ge 0 then $
 	mask[maskidx]=0
 
 mask[backindx] = mask[backindx]+1
 in=where(mask eq 2, sky_npix)
 my_sky = int[in]
 meanclip,my_sky,esn_sky,esn_skysig,clipsig=10
  
 if keyword_set(verbose) then begin
        print, ''
        print, '---------------' 
	print,pre+'Region final sky/px, skyerr/px     [mJy/px]: ', $
		esn_sky, esn_skysig, format='(a,2f13.8)'
	print, '---------------' 
 endif  
 
   skyint = esn_sky
   skyinterr = esn_skysig
   
 endif   ; end-if sky background region has been set manually
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


 musky = -2.5*alog10( (skyint/(as_pix^2)) / zeropoint )
 if finite(musky) then begin
	muskyerr = 2.5/alog(10)*(skyinterr/skyint)
 endif else begin
	musky = -99.
	muskyerr = -9.
 endelse
 if keyword_set(verbose) then begin
	print,pre+'final sky/px, skyerr/px     [mJy/px]: ', $
		skyint, skyinterr, format='(a,2f13.8)'
	print,pre+'final musky, muskyerr  [AB mag/as^2]: ', $
		musky, muskyerr, format='(a,2f10.3)'
 endif
  
  
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; start radial profile

 ;radial profile increments (default = 3 arcsec)
 if not keyword_set(annuli_size) then annuli_size = 3.0
 inc = 1.0*annuli_size / as_pix
 max_annuli = semimajor
 if keyword_set(extendtoskyannulus) then max_annuli = skyouter
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

 if maskidx_band[0] ge 0 then $   
 	mask[maskidx_band]=0
 if no_nan gt 0 then mask[indx_nan]=0
 
 
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

 ;skysubtracted intensity values
 ss_int = tot_int - tot_skyint
 ss_annulus_int = annulus_int - annulus_skyint
 
 ; var is variance image output from swarp
 ; it includes sky and other sources of noise
 if n_elements(var) gt 1 then begin

	good = where(var gt 0.)
	bad  = where(var lt 0., nbad)
	if nbad gt 0 then var[bad] = 0.

 	ellint, var, x0, y0, axratio, a, pa, $
   		tot_int_var, npix_int_var, $
   		annulus_int_var, npixannulus_int_var, mask=mask
 
 endif else begin

	vart = (int - int) + skyinterr^2	; variance/pix

 	ellint, vart, x0, y0, axratio, a, pa, $
   		tot_int_var, npix_int_var, $
   		annulus_int_var, npixannulus_int_var, mask=mask

 endelse

 ;account for all sources of noise

 ;sky subtracted errors
 ss_int_err = sqrt( tot_int_var + (npix_int^2/npix_skyint) * skyinterr^2 )
 ss_annulus_int_err = sqrt( annulus_int_var + $
	 (npixannulus_int^2/npixannulus_skyint) * skyinterr^2 )

 ;magnitudes and errors               
 totarea = npix*(as_pix)^2
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
 good=where(good ge floor(0.25*n) and ss_annulus_int gt 0 and $
	 finite(dm) and mag_err gt 0.,ngood)
 if ngood gt 2 then begin
	 x=dm[good]/dr[good]
	 y=mag[good]
	 erry=mag_err[good]
	 fit=linfit(x,y,measure_errors=erry,sigma=sigma)
	 asymag=fit[0]
	 asymag_err=sigma[0]
	 dif = abs(mag-asymag)
	 w = where(dif eq min(dif[where(dif ge 0.)]),nw)
	 if nw ge 1 then $
		 asymap = a[w[0]] $
	 else	 asymap = a[n_elements(a)-1]
 endif else begin
	 asymag=-99.999
	 asymag_err=-9.999
	 asymap=-9.999
 endelse

 ;print apr mag and error
 good = where(a le (semimajor+1.) and ss_annulus_int gt 0, cntgood)
 if cntgood gt 0 then begin
	best=where(a eq max(a[good]))
	print,pre+"a['],mag,merr[AB]:", $
		a[best]*as_pix/60.,mag[best],mag_err[best],form='(a,f8.3,2f9.4)'
 endif else print,'No good mags.'

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; write out tables

 ; write out total profile

 v=['a','mag','mag_e','int','int_e','int_e_src','int_e_bg','bg','area']
 u=['[arcsec]','[AB]','[AB]','[mJy]','[mJy]','[mJy]','[mJy]',$
    '[mJy]','[as^2]']


 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_totprofile.dat'
 openw,lun,outfile,/get_lun
 printf,lun,'# SDSS_RADPROF '+ver+': total profile'
 printf,lun,'# Date : '+systime()
 printf,lun,'# ID   : '+id
 printf,lun,'# Band : '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')+' mag'
 printf,lun,'# PxScl:'+string(as_pix,form='(f7.3)')+' asec/px'
 printf,lun,'# AnnSz:'+string(annuli_size,form='(f7.3)')+' asec'
 printf,lun,'# Intensity map (mJy): ',int_fi.name,$
  format='(a,a'+strcompress(strlen(int_fi.name))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img           : ',maskimgfile,$
  format='(a,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'# Additional info in ellipse parameter file and background file.'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 array_clean,-99.999,0.,50.,mag,mag_err,/nan
 for i=0,n_elements(a)-1 do $
 printf,lun,a[i]*as_pix,mag[i], mag_err[i], ss_int[i], ss_int_err[i], $
  sqrt(tot_int_var[i]), sqrt(npix[i])*skyinterr,tot_skyint[i], totarea[i],$
  format=$
  '(f12.3,x,f7.3,x,f7.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'total profile written to     : '+outfile
 
 ; write out annular profile

 v=['a','mu','mu_e','int','int_e','int_e_src','int_e_bg','bg','area']
 u=['[arcsec]','[mag/as^2]','[mag/as^2]','[mJy]','[mJy]','[mJy]',$
    '[mJy]','[mJy]','[as^2]']

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_annprofile.dat'
 openw,lun,outfile,/get_lun
 printf,lun,'# SDSS_RADPROF '+ver+': annular profile'
 printf,lun,'# Date : '+systime()
 printf,lun,'# ID   : '+id
 printf,lun,'# Band : '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')+' mag'
 printf,lun,'# PxScl:'+string(as_pix,form='(f7.3)')+' asec/px'
 printf,lun,'# AnnSz:'+string(annuli_size,form='(f7.3)')+' asec'
 printf,lun,'# Intensity map (mJy): ',int_fi.name,$
  format='(a,a'+strcompress(strlen(int_fi.name))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img           : ',maskimgfile,$
  format='(a,a'+strcompress(strlen(maskimgfile))+')'
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
            sqrt(annulus_int_var[i]), sqrt(npixannulus[i])*skyinterr,$
            annulus_skyint[i], area[i],$
  format=$
  '(f12.3,x,f11.3,x,f11.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'annular profile written to   : '+outfile

 ; write out a total flux file

 v=['a','mag','mag_e','int','int_e','int_e_src','int_e_bg','bg','area']
 u=['[arcsec]','[AB]','[AB]','[mJy]','[mJy]','[mJy]','[mJy]',$
    '[mJy]','[as^2]']

 sn=ss_annulus_int/ss_annulus_int_err
 good = where(sn ge 2,cntgood)
 if cntgood gt 0 then $
	 best=where(a eq max(a[good])) $
 else	 best=[n_elements(a)-1L]
 ff=-99.999

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_total.dat' 
 openw,lun,outfile,/get_lun
 printf,lun,'# SDSS_RADPROF '+ver+': total magnitude'
 printf,lun,'# Date : '+systime()
 printf,lun,'# ID   : '+id
 printf,lun,'# Band : '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')+' mag'
 printf,lun,'# PxScl:'+string(as_pix,form='(f7.3)')+' asec/px'
 printf,lun,'# AnnSz:'+string(annuli_size,form='(f7.3)')+' asec'
 printf,lun,'# Intensity map (mJy): ',int_fi.name,$
  format='(a,a'+strcompress(strlen(int_fi.name))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img           : ',maskimgfile,$
  format='(a,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'# Additional info in ellipse parameter file and background file.'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 if cntgood gt 0 then $
 printf,lun,a[best]*as_pix, mag[best], mag_err[best], $
  ss_int[best], ss_int_err[best], $
  sqrt(tot_int_var[best]),sqrt(npix[best])*skyinterr,tot_skyint[best],totarea[best],$
  format=$
  '(f12.3,x,f7.3,x,f7.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)' $
 else $
 printf,lun,ff,ff,ff,ff,ff,ff,ff,ff,ff,  format=$
  '(f12.3,x,f7.3,x,f7.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'total flux written to        : '+outfile


 ; write out asymtotic magnitudes

 v=['a','mag','mag_e','int','int_e','int_e_src','int_e_bg','bg','area']
 u=['[arcsec]','[AB]','[AB]','[mJy]','[mJy]','[mJy]','[mJy]',$
    '[mJy]','[as^2]']

 best=where(a eq asymap)
 if best[0] lt 0 then best=[n_elements(a)-1L]
 ff=-99.999

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_asymptotic.dat' 
 openw,lun,outfile,/get_lun
 printf,lun,'# SDSS_RADPROF '+ver+': asymptotic magnitude'
 printf,lun,'# Date : '+systime()
 printf,lun,'# ID   : '+id
 printf,lun,'# Band : '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')+' mag'
 printf,lun,'# PxScl:'+string(as_pix,form='(f7.3)')+' asec/px'
 printf,lun,'# AnnSz:'+string(annuli_size,form='(f7.3)')+' asec'
 printf,lun,'# Intensity map (mJy): ',int_fi.name,$
  format='(a,a'+strcompress(strlen(int_fi.name))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img           : ',maskimgfile,$
  format='(a,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'# Additional info in ellipse parameter file and background file.'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 array_clean,-99.999,0.,50.,asymag,asymag_err,/nan
 printf,lun,a[best]*as_pix, asymag, asymag_err, $
  ss_int[best], ss_int_err[best], $
  sqrt(tot_int_var[best]),sqrt(npix[best])*skyinterr,tot_skyint[best],totarea[best],$
  format=$
  '(f12.3,x,f7.3,x,f7.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'asymptotic flux written to   : '+outfile


 ; write out a aperture flux file

 v=['a','mag','mag_e','int','int_e','int_e_src','int_e_bg','bg','area']
 u=['[arcsec]','[AB]','[AB]','[mJy]','[mJy]','[mJy]','[mJy]',$
    '[mJy]','[as^2]']

 ; get best within a pixel
 good = where(a le (semimajor+1.) and ss_annulus_int gt 0, cntgood)
 if cntgood gt 0 then best=where(a eq max(a[good]))
 ff=-99.999

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_aperture.dat' 
 openw,lun,outfile,/get_lun
 printf,lun,'# SDSS_RADPROF '+ver+': aperture magnitude'
 printf,lun,'# Date : '+systime()
 printf,lun,'# ID   : '+id
 printf,lun,'# Band : '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')+' mag'
 printf,lun,'# PxScl:'+string(as_pix,form='(f7.3)')+' asec/px'
 printf,lun,'# AnnSz:'+string(annuli_size,form='(f7.3)')+' asec'
 printf,lun,'# Intensity map (mJy): ',int_fi.name,$
  format='(a,a'+strcompress(strlen(int_fi.name))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img           : ',maskimgfile,$
  format='(a,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'# Additional info in ellipse parameter file and background file.'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 if cntgood gt 0 then $
 printf,lun,a[best]*as_pix, mag[best], mag_err[best], $
  ss_int[best], ss_int_err[best], $
  sqrt(tot_int_var[best]),sqrt(npix[best])*skyinterr,tot_skyint[best],totarea[best],$
  format=$
  '(f12.3,x,f7.3,x,f7.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)' $
 else $
 printf,lun,ff,ff,ff,ff,ff,ff,ff,ff,ff,  format=$
  '(f12.3,x,f7.3,x,f7.3,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5,x,e12.5)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'aperture flux written to     : '+outfile
 

 ; write out ellipse parameters

 v=['ra_cen','dec_cen','semimajor','semiminor','p.a.']
 u=['[j2k]','[j2k]','[arcsec]','[arcsec]','[deg]']

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_ellipsepar.dat' 
 openw,lun,outfile,/get_lun
 printf,lun,'# SDSS_RADPROF '+ver+': ellipse parameters'
 printf,lun,'# Date : '+systime()
 printf,lun,'# ID   : '+id
 printf,lun,'# Band : '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')+' mag'
 printf,lun,'# PxScl:'+string(as_pix,form='(f7.3)')+' asec/px'
 printf,lun,'# AnnSz:'+string(annuli_size,form='(f7.3)')+' asec'
 printf,lun,'# Intensity map (mJy): ',int_fi.name,$
  format='(a,a'+strcompress(strlen(int_fi.name))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img           : ',maskimgfile,$
  format='(a,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'# Additional info in ellipse parameter file and background file.'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],$
   format='(a1,a11,x,a12,x,a12,x,a12,x,a8)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],$
   format='(a1,a11,x,a12,x,a12,x,a12,x,a8)'
 printf,lun, ra, dec, semimajor*as_pix, semiminor*as_pix, pa, $
   format='(f12.6,x,f12.6,x,f12.3,x,f12.3,x,f8.3)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'ellipse param info written to: '+outfile


 ; write out background  parameters

 v=['bg','bg_e','mu_bg', 'mu_bg_e', 'scale','skyradius_in','skyradius_out']
 u=['[mJy/as^2]','[mJy/as^2]','[mag/as^2]','[mag/as^2]','[as/pixel]','[as]','[as]']

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_background.dat' 
 openw,lun,outfile,/get_lun
 printf,lun,'# SDSS_RADPROF '+ver+': background'
 printf,lun,'# Date : '+systime()
 printf,lun,'# ID   : '+id
 printf,lun,'# Band : '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')+' mag'
 printf,lun,'# PxScl:'+string(as_pix,form='(f7.3)')+' asec/px'
 printf,lun,'# AnnSz:'+string(annuli_size,form='(f7.3)')+' asec'
 printf,lun,'# Intensity map (mJy): ',int_fi.name,$
  format='(a,a'+strcompress(strlen(int_fi.name))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img           : ',maskimgfile,$
  format='(a,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'# Additional info in ellipse parameter file and background file.'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],$
   format='(a1,a11,x,a12,x,a12,x,a12,x,a12,x,a14,x,a14)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],$
   format='(a1,a11,x,a12,x,a12,x,a12,x,a12,x,a14,x,a14)'
 array_clean,-99.999,0.,50.,musky,muskyerr,/nan
 printf,lun, skyint/as_pix^2, skyinterr/as_pix^2, musky, muskyerr, as_pix,$
              skyinner*as_pix, skyouter*as_pix,$
 format='(e12.5,x,e12.5,x,f12.3,x,f12.3,x,f12.3,x,f14.3,x,f14.3)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'background info written to   : '+outfile

 ; write out limiting mags

 v=['a','lmag_sn3','lmag_sn5']
 u=['[arcsec]','[AB]','[AB]']

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_limitingmags.dat'
 openw,lun,outfile,/get_lun
 printf,lun,'# SDSS_RADPROF '+ver+': limiting mags'
 printf,lun,'# Date : '+systime()
 printf,lun,'# ID   : '+id
 printf,lun,'# Band : '+band
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')+' mag'
 printf,lun,'# PxScl:'+string(as_pix,form='(f7.3)')+' asec/px'
 printf,lun,'# AnnSz:'+string(annuli_size,form='(f7.3)')+' asec'
 printf,lun,'# Intensity map (mJy): ',int_fi.name,$
  format='(a,a'+strcompress(strlen(int_fi.name))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img           : ',maskimgfile,$
  format='(a,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'# Additional info in ellipse parameter file and background file.'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2], format='(a1,a11,x,a11,x,a11)'
 printf,lun,'#',u[0],u[1],u[2], format='(a1,a11,x,a11,x,a11)'
 array_clean,-99.999,0.,50.,lmag3,lmag5,/nan
 for i=0,n_elements(a)-1 do $
 printf,lun,a[i]*as_pix,lmag3[i], lmag5[i], format='(f12.3,x,f11.3,x,f11.3)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'limiting mags written to     : '+outfile

 ;all done
 if keyword_set(verbose) then print,pre+'finis.' 

 status = 0

 return
end
