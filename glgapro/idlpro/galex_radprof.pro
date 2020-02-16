; $Id: galex_radprof.pro,v 1.39 2014/03/03 19:25:42 neill Exp $
;+
; NAME:
;   
;   galex_radprof
;
; PURPOSE:
;
;   generate radial profile of GALEX data
;
; CATEGORY:
;
;   radial photometry
;
; CALLING SEQUENCE:
;
;   galex_radprof, id, ra, dec, $
;                  majordiam, minordiam, pa, $
;                  intfile, rrhrfile, cntfile, $
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
;    intfile: scalar string, GALEX intensity map  
;    rrhrfile: scalar string, GALEX hi-res relative response map 
;    cntfile: scalalr string, GALEX count map 
;
; OPTIONAL INPUTS/KEYWORDS:
;
;    zeropoint=zeropoint: scalar (defaults: NUV=20.08, FUV=18.82)
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
;    total profile    = ID_BAND_totprofile.dat   BAND = NUV or FUV
;    annular profile  = ID_BAND_annprofile.dat   "
;    background       = ID_BAND_background.dat   "
;    ellipse params   = ID_BAND_ellipsepar.dat   "
;    best total flux  = ID_BAND_total.dat        "
;    aperture flux    = ID_BAND_aperture.dat     "
;    asymptotic flux  = ID_BAND_asymptotic.dat   "
;    limiting mags    = ID_BAND_limitingmags.dat "
;    status=status: returns 0 if all OK
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
;   calls: poissonsky, ellint, and many astrolib routines
;
; EXAMPLE:
;
;   band='NUV'
;   filebase='~/ngc4656/data/galex/fits/NGC4656B'
;   intfile=filebase+'-'+band+'.fit*'
;   rrhrfile=filebase+'-'+band+'_rr.fit*'
;   cntfile=filebase+'-'+band+'_cnt.fit*'
;   ra=191.060      
;   dec=32.2784
;   id='NGC4656B'
;   majordiam= 2*2.84809
;   minordiam= 2*1.41348
;   pa=38
;
;   resolve_routine,'galex_radprof'
;
;   galex_radprof,id,ra,dec,$
;                 majordiam,minordiam,pa,$
;                 intfile,rrhrfile,cntfile,$
;                 diam_units='arcmin',/verbose
;
; MODIFICATION HISTORY:
;
;	2008-SEP-09	M. Seibert
;	2014-FEB-28	D. Neill - update for standard outputs
;
;-

pro galex_radprof, $
    id, ra, dec, $
    majordiam, minordiam, pa, $
    intfile, rrhrfile, cntfile, $
    zeropoint=zeropoint,verbose=verbose,diam_units=diam_units,$
    ellipsefile=ellipsefile, $
    maskimgfile=maskimgfile, annuli_size=annuli_size, $
    extendtoskyannulus=extendtoskyannulus,$
    outpath=outpath, plotsky=plotsky, status=status, yuan13=yuan13
;
; note: set yuan13 keyword to use Yuan et al. 2013 empirical extintion coeff's
;	Yuan coeff's assume SFD98 over-estimates E(B-V) by 14% (see table 2)
;
ver=repstr('$Revision: 1.39 $ $Date: 2014/03/03 19:25:42 $','$','')
pre='GALEX_RADPROF: '
status=-1

int_fi = file_info(intfile)
rrh_fi = file_info(rrhrfile)
cnt_fi = file_info(cntfile)
if not int_fi.exists or $
   not rrh_fi.exists or $
   not cnt_fi.exists then begin 
   print,pre+'Can not find one of the input files.'
   return
endif

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; read input files

 int=mrdfits(int_fi.name,0,inthdr,/silent)
 cnt=mrdfits(cnt_fi.name,0,cnthdr,/silent)
 rrhr=mrdfits(rrh_fi.name,0,rrhrhdr,/silent)
 sz=size(rrhr,/dim)

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; check input files are for same band

 intband=sxpar(inthdr,'BAND')
 cntband=sxpar(cnthdr,'BAND')
 rrhrband=sxpar(rrhrhdr,'BAND')

 if cntband ne intband or rrhrband ne intband then begin
  print, pre+'ERROR: input files are not for the same pass-band.'
  return
 endif

 case intband of
  1: band='NUV'
  2: band='FUV'
  else: band='???'
 endcase

 if band eq '???' then begin
  print,pre+'unknown pass-band.'
  return
 endif

 if keyword_set(verbose) then print, pre+'processing '+band+' data'

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; set zeropoint mags to defualt values

 if not keyword_set(zeropoint) then begin
  if band eq 'NUV' then zeropoint = '20.08'
  if band eq 'FUV' then zeropoint = '18.82'
 endif
 if keyword_set(verbose) then print, pre+'using zeropint mag: ',zeropoint

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; convert intensities from cnt/s to mJy
 mjc = 10.^(0.4 * (16.4 - zeropoint))

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
 ebv=dust_getval(l,b,/interp)				; SFD98 E(B-V)
 A_Galactic = glga_getextin(ebv,band,yuan13=yuan13)	; MW extinction
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

 mjdiam=majordiam*scale ; convert to pixels
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

 minskywid = 27.0 / as_pix	; 27 arcsec is minimum sky width
 factor=1.5
 skyinner=factor*semimajor
 skyouter=semimajor*sqrt(1+factor^2) > (skyinner+minskywid)  
 skysemimajors=[skyinner,skyouter]

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; read in mask file currently read
 ; list of idl image indicies
 ; then sets relative response map values to zero

 if keyword_set(verbose) then print,pre+'Using maskimgfile: '+maskimgfile
 mskimg = glga_getmask(maskimgfile,sz,astr,as_pix,/update)
 maskidx = where(mskimg ge 1,count)
 if count gt 0 then rrhr[maskidx]=0
 delvarx,mskimg,maskidx,/free

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; make mask image for photometry where good pixels=1

 dist_ellipse, im, sz, x0, y0, axratio, pa
 mask=bytarr(sz)*0b
 loco=where(rrhr gt 0,nloco) ; only good repsonse regions
 if nloco gt 0 then mask[loco]=1b
 loco=where(im gt skyouter,nloco)
 if nloco gt 0 then mask[loco]=0b
 loci=where(im lt skyinner,nloci)
 if nloci gt 0 then mask[loci]=0b
 delvarx, im, /free
 in=where(mask eq 1b, sky_npix)
 if sky_npix le 0 then begin
	 if keyword_set(verbose) then print,pre+'No good sky pixels, returning'
	 return
 endif


 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; measure sky and sky uncertainty in annulus


 ;take first stab at sky
 poissonsky,cnt[in],skymean=skycntmean,skysigma=skycntsigma,$
           skysdom=skycntsdom, used=used,clipstart=5,plot=plotsky
 if used[0] eq -1 then begin
	 print,pre+'No good sky pixels, returning'
	 return
 endif
 skyidx=in[used] 

 ;mask 7.5x7.5 arcsec (5x5 pixel) area around any 
 ; 3 sigma pixel and measure sky again 
 ; max of 3 iterations
 iter=0
 if band eq 'FUV' then nsig=3.0
 if band eq 'NUV' then nsig=2.0
 hi=where(cnt[in] ge (skycntmean + nsig*skycntsigma),counthi )
 if keyword_set(verbose) then $
    print,pre+'3sig pix region mask: ',iter,skycntmean,skycntsigma,counthi, $
    	format='(a,i4,2f13.8,i6)'
 while counthi gt 0 do begin  
  iter=iter+1
  msize=3 ; arcsec   
  msize_pix = round(msize/as_pix)
  xy=get_indices(in[hi],sz)
  x1=xy[0,*]-msize_pix > 0
  x2=xy[0,*]+msize_pix < (sz[0]-1)
  y1=xy[1,*]-msize_pix > 0
  y2=xy[1,*]+msize_pix < (sz[1]-1)
  for k=0l,n_elements(hi)-1 do mask[x1[k]:x2[k],y1[k]:y2[k]]=0b    
  in=where(mask eq 1b, sky_npix)
  if sky_npix le 30 then goto, skipout
  ;take second stab 
  poissonsky,cnt[in],skymean=skycntmean,skysigma=skycntsigma,$
             skysdom=skycntsdom, used=used,clipstart=5,plot=plotsky;,/verb
  if skycntmean eq 0 or skycntsigma eq 0 then goto,skipout 
  if used[0] eq -1 then begin
	 print,pre+'No good sky pixels, returning'
	 return
  endif
  skyidx=in[used]
  hi=where(cnt[in] ge (skycntmean + nsig*skycntsigma),counthi)
  if keyword_set(verbose) then $
     print,pre+'3sig pix region mask: ',iter,skycntmean,skycntsigma,counthi, $
     	format='(a,i4,2f13.8,i6)'
  if iter ge 3 then counthi=0
 endwhile
 skipout:

 ;compute sky intensity per pixel
 ;  <rate>=total(counts)/total(time) is better than total(rate)/Npixels
 ;  when exptime (rel. resp.) varies from pixel to pixel and will be the same
 ;  if exptime is uniform  
 ;truth be told this is typically not too different from the mean and 
 ; sdom of intensity 

 skyint = total(int[skyidx]*rrhr[skyidx])/total(rrhr[skyidx])
 
 ;compute sky error from counting statistics
   
 skyintcounterr = sqrt(total(int[skyidx]*rrhr[skyidx])) / total(rrhr[skyidx])

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
  if n_elements(sub) gt 1 then begin
   s[j].sky = total(int[skyidx[sub]]*rrhr[skyidx[sub]])/ $
              total(rrhr[skyidx[sub]])
   s[j].skyerr = sqrt(total(int[skyidx[sub]]*rrhr[skyidx[sub]])) / $
                 total(rrhr[skyidx[sub]])
  endif
 endfor

 if keyword_set(verbose) then begin
	print, pre+'Subsection sky and counting errors:'
	forprint, strarr(n_elements(s))+pre,indgen(n_elements(s))+1, $
		mjc*s.sky,mjc*s.skyerr, format='(a,i8,2f13.8)'
 endif

 g=where(finite(s.sky) eq 1, ng)
 if ng gt 1 then $
 	mom=moment(s.sky,/nan) $
 else if ng eq 1 then $
 	mom = [s[g].sky,0.] $
 else if ng lt 1 then $
 	mom = [!values.f_nan,!values.f_nan]

 meannskyint = mom[0]
 skyintspatialerr = sqrt(mom[1])

 if keyword_set(verbose) then begin
	print,pre+'Sub <sky> and sigma_<sky>   [mJy/px]: ', $
		mjc*meannskyint, mjc*skyintspatialerr, format='(a,2f13.8)'
	print,pre+'Img sky, sigma_sky          [mJy/px]: ', $
		mjc*skyint, mjc*skyintcounterr, format='(a,2f13.8)'
 endif

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;if we assume that the counting error and spatial error 
 ; are unrelated we can add them in quadrature
 ; but the spatial may have a component of the counting error
 ; will take the max error (spatial always) 
 ; will account for sky counting error later

 ;skyerr =sqrt(skyintcounterr^2 + skyintspatialerr^2)
 skyinterr =max([skyintcounterr,skyintspatialerr])

 musky = zeropoint - 2.5*alog10(skyint/(as_pix^2))
 if finite(musky) then begin
	muskyerr = 2.5/alog(10)*(skyinterr/skyint)
 endif else begin
	musky = -99.
	muskyerr = -9.
 endelse
 if keyword_set(verbose) then begin
	print,pre+'final sky/px, skyerr/px     [mJy/px]: ', $
		mjc*skyint, mjc*skyinterr, format='(a,2f13.8)'
   	print,pre+'final musky, muskyerr  [AB mag/as^2]: ', $
		musky, muskyerr, format='(a,2f10.3)'
 endif
 
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; new sky method

 mask=intarr(sz)*0
 mask[where(rrhr gt 0)]=1 ; only good repsonse regions   

 ellint, int, x0,y0,axratio,[skyinner,skyouter],pa,$
   tot_skyint, npix_skyint, annulus_skyint, npixannulus_skyint, mask=mask
 skyint=annulus_skyint[1]/npixannulus_skyint[1]

 musky = zeropoint - 2.5*alog10(skyint/(as_pix^2))
 if finite(musky) then begin
	muskyerr = 2.5/alog(10)*skyinterr/skyint
 endif else begin
	 musky = -99.
	 muskyerr = -9.
 endelse
 if keyword_set(verbose) then begin
	print,pre+'FINAL sky/px, skyerr/px     [mJy/px]: ', $
		mjc*skyint, mjc*skyinterr, format='(a,2f13.8)'
   	print,pre+'FINAL musky, muskyerr     [mag/as^2]: ', $
		musky, muskyerr, format='(a,2f10.3)'
 endif
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; start radial profile

 ;radial profile increments (default = 6 arcsec)
 if not keyword_set(annuli_size) then annuli_size = 6.0
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
 ;note rel. resp. map has been set to zero 
 ;where we want to mask 
 mask=intarr(sz)*0
 mask[where(rrhr gt 0)]=1 ; only good repsonse regions   

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

 ;measure counterrors in the form of counts/time^2
 ; thus rate count errors are sqrt of integrated values
 ellint, int/(rrhr > 1e-20), x0,y0,axratio,a,pa,$
   tot_int_cnterr, npix_int_cnterr, $
   annulus_int_cnterr, npixannulus_int_cnterr, mask=mask
 
 tot_int_cnterr = sqrt(tot_int_cnterr)
 annulus_int_cnterr = sqrt(annulus_int_cnterr)

 ;skysubtracted intensity values
 ss_int = tot_int - tot_skyint
 ss_annulus_int = annulus_int - annulus_skyint
 
 ;intensity errors = counting errors + sky errors
 ss_int_err = sqrt(tot_int_cnterr^2 + (sqrt(npix)*skyinterr)^2)
 ss_annulus_int_err = sqrt(annulus_int_cnterr^2+(sqrt(npixannulus)*skyinterr)^2)

 ;magnitudes and errors               
 totarea = npix*(as_pix)^2
 mag = zeropoint - 2.5*alog10(ss_int)
 mag_err = 2.5/alog(10)*ss_int_err/ss_int

 ;limiting mags
 lmag3 = zeropoint - 2.5*alog10(ss_int_err*3.)
 lmag5 = zeropoint - 2.5*alog10(ss_int_err*5.)
 
 ;annular surface brightness and errors   
 area=npixannulus*(as_pix)^2
 mu_annulus = zeropoint - 2.5*alog10(ss_annulus_int/area)
 mu_annulus_err = 2.5/alog(10)*ss_annulus_int_err/ss_annulus_int

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
 printf,lun,'# GALEX_RADPROF '+ver+': total profile'
 printf,lun,'# Date : '+systime()
 printf,lun,'# ID   : '+id
 printf,lun,'# Band : '+band
 printf,lun,'# MagZp:'+string(zeropoint,form='(f7.3)')+' ABmag'
 printf,lun,'# INmJy:'+string(mjc,form='(g9.3)')+' mJy/cnt/s'
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')+' mag'
 printf,lun,'# PxScl:'+string(as_pix,form='(f7.3)')+' asec/px'
 printf,lun,'# AnnSz:'+string(annuli_size,form='(f7.3)')+' asec'
 printf,lun,'# Intensity map (cnt/s): ',int_fi.name,$
  format='(a,a'+strcompress(strlen(int_fi.name))+')'
 printf,lun,'# Rel. resp. map (s)   : ',rrh_fi.name,$
  format='(a,a'+strcompress(strlen(rrh_fi.name))+')'
 printf,lun,'# Count map (cnt)      : ',cnt_fi.name,$
  format='(a,a'+strcompress(strlen(cnt_fi.name))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img             : ',maskimgfile,$
  format='(a,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'# Additional info in ellipse parameter file and background file.'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 array_clean,-99.999,0.,50.,mag,mag_err,/nan
 for i=0,n_elements(a)-1 do $
 printf,lun,a[i]*as_pix,mag[i], mag_err[i], mjc*ss_int[i], mjc*ss_int_err[i], $
  mjc*tot_int_cnterr[i], mjc*sqrt(npix[i])*skyinterr,tot_skyint[i], totarea[i],$
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
 printf,lun,'# GALEX_RADPROF '+ver+': annular profile'
 printf,lun,'# Date : '+systime()
 printf,lun,'# ID   : '+id
 printf,lun,'# Band : '+band
 printf,lun,'# MagZp:'+string(zeropoint,form='(f7.3)')+' ABmag'
 printf,lun,'# INmJy:'+string(mjc,form='(g9.3)')+' mJy/cnt/s'
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')+' mag'
 printf,lun,'# PxScl:'+string(as_pix,form='(f7.3)')+' asec/px'
 printf,lun,'# AnnSz:'+string(annuli_size,form='(f7.3)')+' asec'
 printf,lun,'# Intensity map (cnt/s): ',int_fi.name,$
  format='(a,a'+strcompress(strlen(int_fi.name))+')'
 printf,lun,'# Rel. resp. map (s)   : ',rrh_fi.name,$
  format='(a,a'+strcompress(strlen(rrh_fi.name))+')'
 printf,lun,'# Count map (cnt)      : ',cnt_fi.name,$
  format='(a,a'+strcompress(strlen(cnt_fi.name))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img             : ',maskimgfile,$
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
            mjc*ss_annulus_int[i], mjc*ss_annulus_int_err[i], $
            mjc*annulus_int_cnterr[i], mjc*sqrt(npixannulus[i])*skyinterr,$
            mjc*annulus_skyint[i], area[i],$
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
 printf,lun,'# GALEX_RADPROF '+ver+': total magnitude'
 printf,lun,'# Date : '+systime()
 printf,lun,'# ID   : '+id
 printf,lun,'# Band : '+band
 printf,lun,'# MagZp:'+string(zeropoint,form='(f7.3)')+' ABmag'
 printf,lun,'# INmJy:'+string(mjc,form='(g9.3)')+' mJy/cnt/s'
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')+' mag'
 printf,lun,'# PxScl:'+string(as_pix,form='(f7.3)')+' asec/px'
 printf,lun,'# AnnSz:'+string(annuli_size,form='(f7.3)')+' asec'
 printf,lun,'# Intensity map (cnt/s): ',int_fi.name,$
  format='(a,a'+strcompress(strlen(int_fi.name))+')'
 printf,lun,'# Rel. resp. map (s)   : ',rrh_fi.name,$
  format='(a,a'+strcompress(strlen(rrh_fi.name))+')'
 printf,lun,'# Count map (cnt)      : ',cnt_fi.name,$
  format='(a,a'+strcompress(strlen(cnt_fi.name))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img             : ',maskimgfile,$
  format='(a,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'# Additional info in ellipse parameter file and background file.'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 if cntgood gt 0 then $
 printf,lun,a[best]*as_pix, mag[best], mag_err[best], $
  mjc*ss_int[best], mjc*ss_int_err[best], $
  mjc*tot_int_cnterr[best], mjc*sqrt(npix[best])*skyinterr, $
  mjc*tot_skyint[best], totarea[best],$
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
 printf,lun,'# GALEX_RADPROF '+ver+': asymptotic magnitude'
 printf,lun,'# Date : '+systime()
 printf,lun,'# ID   : '+id
 printf,lun,'# Band : '+band
 printf,lun,'# MagZp:'+string(zeropoint,form='(f7.3)')+' ABmag'
 printf,lun,'# INmJy:'+string(mjc,form='(g9.3)')+' mJy/cnt/s'
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')+' mag'
 printf,lun,'# PxScl:'+string(as_pix,form='(f7.3)')+' asec/px'
 printf,lun,'# AnnSz:'+string(annuli_size,form='(f7.3)')+' asec'
 printf,lun,'# Intensity map (cnt/s): ',int_fi.name,$
  format='(a,a'+strcompress(strlen(int_fi.name))+')'
 printf,lun,'# Rel. resp. map (s)   : ',rrh_fi.name,$
  format='(a,a'+strcompress(strlen(rrh_fi.name))+')'
 printf,lun,'# Count map (cnt)      : ',cnt_fi.name,$
  format='(a,a'+strcompress(strlen(cnt_fi.name))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img             : ',maskimgfile,$
  format='(a,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'# Additional info in ellipse parameter file and background file.'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 array_clean,-99.999,0.,50.,asymag,asymag_err,/nan
 printf,lun,a[best]*as_pix, asymag, asymag_err, $
  mjc*ss_int[best], mjc*ss_int_err[best], $
  mjc*tot_int_cnterr[best], mjc*sqrt(npix[best])*skyinterr, $
  mjc*tot_skyint[best], totarea[best],$
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
 printf,lun,'# GALEX_RADPROF '+ver+': aperture magnitude'
 printf,lun,'# Date : '+systime()
 printf,lun,'# ID   : '+id
 printf,lun,'# Band : '+band
 printf,lun,'# MagZp:'+string(zeropoint,form='(f7.3)')+' ABmag'
 printf,lun,'# INmJy:'+string(mjc,form='(g9.3)')+' mJy/cnt/s'
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')+' mag'
 printf,lun,'# PxScl:'+string(as_pix,form='(f7.3)')+' asec/px'
 printf,lun,'# AnnSz:'+string(annuli_size,form='(f7.3)')+' asec'
 printf,lun,'# Intensity map (cnt/s): ',int_fi.name,$
  format='(a,a'+strcompress(strlen(int_fi.name))+')'
 printf,lun,'# Rel. resp. map (s)   : ',rrh_fi.name,$
  format='(a,a'+strcompress(strlen(rrh_fi.name))+')'
 printf,lun,'# Count map (cnt)      : ',cnt_fi.name,$
  format='(a,a'+strcompress(strlen(cnt_fi.name))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img             : ',maskimgfile,$
  format='(a,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'# Additional info in ellipse parameter file and background file.'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],u[7],u[8],$
   format='(a1,a11,x,a7,x,a7,x,a12,x,a12,x,a12,x,a12,x,a12,x,a12)'
 if cntgood gt 0 then $
 printf,lun,a[best]*as_pix, mag[best], mag_err[best], $
  mjc*ss_int[best], mjc*ss_int_err[best], $
  mjc*tot_int_cnterr[best],mjc*sqrt(npix[best])*skyinterr, $
  mjc*tot_skyint[best], totarea[best],$
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
 printf,lun,'# GALEX_RADPROF '+ver+': ellipse parameters'
 printf,lun,'# Date : '+systime()
 printf,lun,'# ID   : '+id
 printf,lun,'# Band : '+band
 printf,lun,'# MagZp:'+string(zeropoint,form='(f7.3)')+' ABmag'
 printf,lun,'# INmJy:'+string(mjc,form='(g9.3)')+' mJy/cnt/s'
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')+' mag'
 printf,lun,'# PxScl:'+string(as_pix,form='(f7.3)')+' asec/px'
 printf,lun,'# AnnSz:'+string(annuli_size,form='(f7.3)')+' asec'
 printf,lun,'# Intensity map (cnt/s): ',int_fi.name,$
  format='(a,a'+strcompress(strlen(int_fi.name))+')'
 printf,lun,'# Rel. resp. map (s)   : ',rrh_fi.name,$
  format='(a,a'+strcompress(strlen(rrh_fi.name))+')'
 printf,lun,'# Count map (cnt)      : ',cnt_fi.name,$
  format='(a,a'+strcompress(strlen(cnt_fi.name))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img             : ',maskimgfile,$
  format='(a,a'+strcompress(strlen(maskimgfile))+')'
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
 printf,lun,'# GALEX_RADPROF '+ver+': background'
 printf,lun,'# Date : '+systime()
 printf,lun,'# ID   : '+id
 printf,lun,'# Band : '+band
 printf,lun,'# MagZp:'+string(zeropoint,form='(f7.3)')+' ABmag'
 printf,lun,'# INmJy:'+string(mjc,form='(g9.3)')+' mJy/cnt/s'
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')+' mag'
 printf,lun,'# PxScl:'+string(as_pix,form='(f7.3)')+' asec/px'
 printf,lun,'# AnnSz:'+string(annuli_size,form='(f7.3)')+' asec'
 printf,lun,'# Intensity map (cnt/s): ',int_fi.name,$
  format='(a,a'+strcompress(strlen(int_fi.name))+')'
 printf,lun,'# Rel. resp. map (s)   : ',rrh_fi.name,$
  format='(a,a'+strcompress(strlen(rrh_fi.name))+')'
 printf,lun,'# Count map (cnt)      : ',cnt_fi.name,$
  format='(a,a'+strcompress(strlen(cnt_fi.name))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img             : ',maskimgfile,$
  format='(a,a'+strcompress(strlen(maskimgfile))+')'
 printf,lun,'#'
 printf,lun,'#',v[0],v[1],v[2],v[3],v[4],v[5],v[6],$
   format='(a1,a11,x,a12,x,a12,x,a12,x,a12,x,a14,x,a14)'
 printf,lun,'#',u[0],u[1],u[2],u[3],u[4],u[5],u[6],$
   format='(a1,a11,x,a12,x,a12,x,a12,x,a12,x,a14,x,a14)'
 array_clean,-99.999,0.,50.,musky,muskyerr,/nan
 printf,lun, mjc*skyint/as_pix^2, mjc*skyinterr/as_pix^2, musky, muskyerr, $
	 as_pix, skyinner*as_pix, skyouter*as_pix,$
 format='(e12.5,x,e12.5,x,f12.3,x,f12.3,x,f12.3,x,f14.3,x,f14.3)'
 free_lun,lun
 if keyword_set(verbose) then print,pre+'background info written to   : '+outfile

 ; write out limiting mags

 v=['a','lmag_sn3','lmag_sn5']
 u=['[arcsec]','[AB]','[AB]']

 if not keyword_set(outpath) then outpath="./"
 outfile=outpath+'/'+id+'_'+band+'_limitingmags.dat'
 openw,lun,outfile,/get_lun
 printf,lun,'# GALEX_RADPROF '+ver+': limiting mags'
 printf,lun,'# Date : '+systime()
 printf,lun,'# ID   : '+id
 printf,lun,'# Band : '+band
 printf,lun,'# MagZp:'+string(zeropoint,form='(f7.3)')+' ABmag'
 printf,lun,'# INmJy:'+string(mjc,form='(g9.3)')+' mJy/cnt/s'
 printf,lun,'# A_Gal:'+string(A_Galactic,form='(f7.3)')+' mag'
 printf,lun,'# PxScl:'+string(as_pix,form='(f7.3)')+' asec/px'
 printf,lun,'# AnnSz:'+string(annuli_size,form='(f7.3)')+' asec'
 printf,lun,'# Intensity map (cnt/s): ',int_fi.name,$
  format='(a,a'+strcompress(strlen(int_fi.name))+')'
 printf,lun,'# Rel. resp. map (s)   : ',rrh_fi.name,$
  format='(a,a'+strcompress(strlen(rrh_fi.name))+')'
 printf,lun,'# Count map (cnt)      : ',cnt_fi.name,$
  format='(a,a'+strcompress(strlen(cnt_fi.name))+')'
 if keyword_set(maskimgfile) then  $
  printf,lun,'# Mask img             : ',maskimgfile,$
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
