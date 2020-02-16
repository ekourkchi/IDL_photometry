PRO ellint, image, x0, y0, ratio ,a ,pa , tot, npix, totdif, npixdif, $
    mask=mask, verbose=verbose, apmean=apmean, apsigma=apsigma
;+
; NAME:
;    ELLINT
; PURPOSE:
;    Integrate flux in concentric elliptical apertures.  User
;    responsible for 
;    subtracting sky first.  See notes (1) and (2) below.
; CALLING SEQUENCE:
;    ellint, image, x0, y0, ratio,a,pa, tot, npix, totdif, npixdif, $
;           mask=mask, /verbose
; INPUT PARAMETERS:
;    image      2-D image array
;    x0, y0     center of circular apertures (need not be integral)
;    ratio      semimajor/semiminor ratio (scalar)
;    a          vector of semimajor axis'
;    pa         position angle scalar
; OUTPUT PARAMETERS (see note 2):
;    tot        total raw flux in each aperture
;    npix       number of pixels in each aperture
;    totdif     total raw flux for each annulus (difference of successive
;               apertures)
;    npixdif    number of pixels in each annulus
; OPTIONAL INPUT KEYWORDS:
;    mask       Must be same size as image.  If a pixel in
;               the mask is 1, the corresponding pixel in the image is
;               counted.  If 0, the corresponding pixel in the image is
;               ignored (in all results).  See note (3).
;    verbose    Use of /verbose will make the routine tell you what stage
;               it's at.  See note (1).
;    apmean     mean value in aperture
;    apsigma    std deviation in aperture
;
; COMMON BLOCKS:
;    none
; NOTES:
;   (1) If you aren't sure you've set up right, use the /verbose keyword,
;       because the routine is fairly slow.
;   (2) For a surface brightness profile, plot totdif/npixdif vs. radius.
;       For an aperture growth curve, plot tot vs. radius.
;   (3) Mask is intended to mask out stars or garbage.  Depending on your
;       application, you might be better off modifying the input image.
; SIDE EFFECTS:  none
; PROCEDURE:   Similar to IDL/UIT/DAOPHOT aperture routine.
; MODIFICATION HISTORY:
;  Integrate in circular aperture.
;  RSH - STX - 20 Aug 1990
;  Modified to ignore some pixels according to mask image.
;  RSH - STX - 19 Sept 90
;  Totals done for annuli rather than discs.
;  RSH - STX - 20 Sept 90
;  Use of mask corrected.  RSH - STX - 27 Sept 90
;  Small change to conserve array space.  RSH - STX - 5 Nov 90
;  Fractional-pixel approximation adopted from Wayne Landsman's version
;     of DAOPHOT APER.  Annuli computed from disks.  RSH - STX - 17 July 91
;  Several bugs fixed.  RSH - STX - 3 Oct 91
;  Spiffed up for UIT IDL library.  RSH - Hughes STX - 11 June 92
;  Speeded up initializaion.  RSH - HSTX - 5 August 1993
;  Converted from circint to ellint - MHS 23 AUG 1999
;  Added aperture mean and sigma options - MHS AUG 2008
;-
on_error,1
IF n_params(0) LT 1 THEN BEGIN
   message,'Calling sequence:  ellint,image,x0,y0,' $
          +'ratio,a,pa,tot,npix,totdif,npixdif,' $
          +'mask=mask,/verbose'
ENDIF
IF NOT keyword_set(verbose) THEN verbose=0
sz       = size(image)
bell     = string(07b)
IF keyword_set(mask) THEN BEGIN
   szm = size(mask)
   IF (szm(1) NE sz(1)) OR (szm(2) NE sz(2)) THEN BEGIN
      message,'Mask image must be same size as flux image.'
   ENDIF
ENDIF
IF verbose NE 0 THEN $ 
   message,/inf,'Beginning routine to integrate in elliptical apertures ...'
nsm       = n_elements(a)
order_test=a 
if nsm ne 1 then order_test = a(1:nsm-1) - a(0:nsm-2)
ww         = where((order_test LE 0),count)
IF count GT 0 THEN $                             
   message,'Semimajors must be specified in ascending order.'
ww  = where((a LT 0.5), count)
IF count GT 0 THEN $
   message,'All semimajors must be 0.5 pixels or greater.'

;Create an elliptical image mask
DIST_ELLIPSE, ell, [sz(1),sz(2)], x0, y0, ratio, pa

tot      = fltarr(nsm)
totdif   = tot
npix     = fltarr(nsm)
npixdif  = npix
apmean   = fltarr(nsm)
apsigma  = fltarr(nsm)

IF verbose NE 0 THEN message,/inf,'Initialization complete.  Now integrating.'
FOR ir=0, nsm-1 DO BEGIN
   IF verbose NE 0 THEN $
      message,/inf,'Now doing aperture of ' $
                   +strtrim(a(ir),2)+' pixel semimajor'
   IF keyword_set(mask) THEN BEGIN
      within = (ell LE a(ir)) AND mask
   ENDIF ELSE within = ell LE a(ir)
   w = where(within, nw)
   if nw gt 2 then begin
   	apmoment=moment(image(w))
   	apmean(ir)=apmoment[0]               
   	apsigma(ir)=sqrt(apmoment[1])
   	thisapr = ell(w)
	;Fraction of pixels to count
	;fractn needs to be double or we run into errors
   	fractn = (a(ir)-thisapr < 1.0 >0.0 )*1.d0 
   	tot(ir) = total(image(w)*fractn)
   	npix(ir)= total(fractn)
   	thisapr=0 & fractn=0 & within=0
   endif else begin
	apmean(ir) = 0.
	apsigma(ir) = 0.
	tot(ir) = 0.
	npix(ir) = 0L
   endelse
ENDFOR
totdif (0)  = tot(0)
npixdif(0) = npix(0)

FOR ir=1,nsm-1 DO BEGIN
   npixdif(ir) = npix(ir) - npix(ir-1)
   totdif(ir)  = tot(ir) - tot(ir-1)
ENDFOR

IF verbose NE 0 THEN message,/inf,bell+'Done.'
RETURN
END

