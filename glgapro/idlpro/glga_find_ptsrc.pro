pro glga_find_ptsrc,im,el,astr,as_pix,psrc,nobj, $
	sdss=sdss,galex=galex,panstarrs=panstarrs, twomass=twomass,wise=wise,irac=irac, acs=acs, $
	sigmult=sigmult,maglim=maglim,galaxy=galaxy,outside=outside, $
	verbose=verbose, zpmag=zpmag
;+
; glga_find_ptsrc - automatically generate a list of point sources to mask
;
; im - intensity image
; el - ellipse info: [a,b,x,y,pa]
; astr - astrometry for image
; as_pix - arcseconds per pixel
; psrc - point source array with [a,d,x,y,rpx,psrc_astr?,1] for each star found
; nobj - number of objects in psrc
; 
; sigmult - sigma multiplier for threshhold (default = 2.5)
; maglim - magnitude faint limit (default = 22)
; outside - exclude galaxy
; galaxy - measure only galaxy
; keywords specify data type
; verbose - extra output
;-
; data type
srvy='galex'
fwhm=4.0
zeropoint = 1.0
zeromag = 20.08
int=im
if keyword_set(sdss) then begin
	srvy='sdss '
	fwhm=4.0
        zeropoint = 3631000.0
	zeromag = 0.0
        int=im
endif
if keyword_set(panstarrs) then begin
	srvy='panstarrs '
	fwhm=4.0
        zeropoint = 3631000.0
	zeromag = 0.0
        int=im
endif
if keyword_set(acs) then begin
	srvy='acs '
	fwhm=4.0
        zeropoint = 3631000.0
	zeromag = 0.0
        int=im
endif
if keyword_set(twomass) then begin
	srvy='2mass'
	fwhm=4.0
	zeropoint = 666700.0
	zeromag = 0.0
        int=im
endif
if keyword_set(wise) then begin
	srvy='wise'
	fwhm=4.0
	zeropoint = 1.0
	zeromag = 22.5
	if keyword_set(zpmag) then $
		zeromag = zpmag
	int=im
endif
if keyword_set(irac) then begin
	srvy='irac'
	fwhm=4.0
        ;convert from MJy/str to mJy to make mag cut
        zeropoint = 3631000.0
	zeromag = 0.0
        stpa = 1.D0 / ( (180.D0/!DPI)^2 * 3600.D0^2 ); steradians per arcsec^2
        int=im*stpa* 1.e9 ; mJy / arcsec^2
        int = int * as_pix^2
endif
;
if keyword_set(sigmult) then $
	sgm = sigmult $
else	sgm = 2.5
if keyword_set(maglim) then $
	maxmag = maglim $
else	maxmag = 22
;
; imrad
minskywid = 27. / as_pix	; minimum sky width 5. for ACS
factor = 1.5			; expand by 50%
skyinner=factor*el[0]
skyouter=el[0]*sqrt(1.+factor^2) < (skyinner+minskywid) ; (skyinner+minskywid) ; 
imrad = skyouter


imrad = skyouter*3.
;
; eldist
sz  = size(im,/dim)
rat = el[0]/el[1]
pa  = el[4]
x   = el[2]
y   = el[3]
dist_ellipse,eldist,sz,x,y,rat,pa
dist_circle,circdist,sz,x,y
;
mask=bytarr(sz)*0b
loc=where(eldist gt skyinner and eldist lt skyouter and finite(im), nloc)
if nloc gt 0 then mask[loc]=1b
in=where(mask eq 1b, sky_npix)
if sky_npix le 0 then in=lindgen(n_elements(im))
meanclip, int[in], sky, skysg, clipsig = 5
if not keyword_set(galex) then begin
	mmm,int[in], mmm_sky, mmm_sigma
	if mmm_sigma gt 0. and finite(mmm_sky) eq 1 then begin
		sky = mmm_sky
		skysg = mmm_sigma
	endif
endif
;
; set limits
roundlim=[-1.5,1.5]
sharplim=[0.1,1.1]
;hmin=sky+sgm*skysg	; intensity thresh
hmin=sgm*skysg		; intensity thresh
;
find,int,x,y,flux,sharp,roundness,hmin,fwhm,roundlim,sharplim,/silent
;
;magnitude cut
if n_elements(x) gt 0 then begin
 mag = -2.5*alog10(flux/zeropoint) + zeromag
 brite=where(mag le maxmag,countbrite)
 print,minmax(mag), 'N < ',maxmag,' = ',countbrite
 if countbrite gt 0 then begin
   x=x[brite]
   y=y[brite]
   flux=flux[brite]
   sharp=sharp[brite]
   roundness=roundness[brite]
 endif else delvarx, x,y,flux,sharp,roundness
endif
;
nobj=n_elements(x)
if nobj le 0 then begin
	if keyword_set(verbose) then $
		print,'No point sources found.'
	return
endif
;
psrc = fltarr(7,nobj)
;
; loop over nobj
for i=0L,nobj-1L do begin
    xp = fix(x[i]+0.5)
    yp = fix(y[i]+0.5)
;
; avoid nucleus
    if circdist[xp,yp] gt 2.*fwhm then begin
	if keyword_set(outside) then begin
		if circdist[xp,yp] le imrad*2.0 and $ ;+4.*fwhm and $
		   eldist[xp,yp] ge el[0]-4.*fwhm then begin
			xy2ad,x[i],y[i],astr,a,d
			psrc_ast = 1
;
; use subim to speed things up
			x0=fix(x[i]-63.)>0
			x1=x0+127<(sz[0]-1)
			y0=fix(y[i]-63.)>0
			y1=y0+127<(sz[1]-1)
			radius = get_mask_radius(im[x0:x1,y0:y1], $
				x[i]-float(x0),y[i]-float(y0),srvy)
			psrc[*,i] = [a,d,x[i],y[i],radius,psrc_ast,3.]
		endif	; eldist[xp,yp] le imrad
	endif else if keyword_set(galaxy) then begin
		if eldist[xp,yp] le el[0]+4.*fwhm then begin
			xy2ad,x[i],y[i],astr,a,d
			psrc_ast = 1
;
; use subim to speed things up
			x0=fix(x[i]-63.)>0
			x1=x0+127<(sz[0]-1)
			y0=fix(y[i]-63.)>0
			y1=y0+127<(sz[1]-1)
			radius = get_mask_radius(im[x0:x1,y0:y1], $
				x[i]-float(x0),y[i]-float(y0),srvy)
			psrc[*,i] = [a,d,x[i],y[i],radius,psrc_ast,4.]
		endif	; eldist[xp,yp] le imrad
	endif else begin
		if circdist[xp,yp] le imrad*2.0 then begin
			xy2ad,x[i],y[i],astr,a,d
			psrc_ast = 1
;
; use subim to speed things up
			x0=fix(x[i]-63.)>0
			x1=x0+127<(sz[0]-1)
			y0=fix(y[i]-63.)>0
			y1=y0+127<(sz[1]-1)

			print, "prsrc size:", size(im)
			print, x0, x1, y0, y1

			radius = get_mask_radius(im[x0:x1,y0:y1], $
				x[i]-float(x0),y[i]-float(y0),srvy)
			psrc[*,i] = [a,d,x[i],y[i],radius,psrc_ast,2.]
		endif	; eldist[xp,yp] le imrad
	endelse
    endif else begin	; skip nucleus
	    if keyword_set(verbose) then $
		    print,'Skipped nuclear pt. source: ',xp,yp
    endelse
endfor	; loop over nobj
;
w=where(psrc[6,*] ge 2., nobj)
if nobj gt 0 then $
	psrc=psrc[*,w] $
else	psrc=-1.
;
return
end
