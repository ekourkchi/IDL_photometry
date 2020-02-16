pro rgb_irac_asinh, fspec=fspec,scales=scales,brite=brite, $
	nonlinearity=nonlinearity, origin=origin, verbose=verbose, $
	quality=quality, outdir=outdir, update=update, greenmix=greenmix, $
	swonly=swonly, lwonly=lwonly
;+
; NAME:
;  rgb_irac_asinh
;
;
; PURPOSE:
;  Generate RGB composite jpeg images from Spitzer IRAC intensity maps.
;  RED=NUV, BLUE=FUV, GREEN=combination of FUV & NUV.
;
;  Images are lambda*flamda and asinh scaled/fit using Wherry,
;  Blanton, Hogg IDL routines of the the Lupton, et al. algorithm.
;
;  Defaults to full resolution images.
;  Optionally creates 1/2, 1/4 and 1/8 resolution images.
;  Will process single files or lists of files as arrays.
;
; CATEGORY:
;  image processing
;
;
; CALLING SEQUENCE:
;  rgb_irac_asinh, fspec=fspec,scales=scales,brite=brite, $
;                     nonlinearity=nonlinearity, origin=origin, $
;                     verbose=verbose, quality=quality
;
; INPUTS:
;          
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;  brite = scale factor applied to all three RGB planes. Default=1
;  verbose = set keyword to complain if an error occurs
;          and give status of jpeg file creation
;  quality = jpeg compression quality for full and half res images.
;            Default = 75. Other resolutions are always quality=100.
;
;  The following keywords are associated with the Wherry et al routines.
;  ----------------------------------------------------------------------
;  scales = scaling applied to red green and blue. Default=[.085,.09,.085]
;  nonlinearity = nonlinearity factor used in asinh scaling. Default=2.5
;  origin = Limits the pixel values of the image to a 'box', so that
;           the colors do not saturate to white but to a specific color. 
;           Default=[0,0,0]
;
;
; OUTPUTS:
;  full scale jpeg -> target.jpg
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;  uses file_exist.pro and astrolib routines
;  assumes Unix platform
;  
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;  M. Seibert 12/24/2004
;  M. Seibert 12/06/2005
;  D. Neill, kludged for SN survey 13/07/2007
;  D. Neill, further kludge for 2mass LGA data 16/08/2007
;
;
;-
; get file list
;if not keyword_set(fspec) then fspec='*_3p6um.fits'
;swfile = file_search(fspec)
;
; get file list
if not keyword_set(fspec) then begin
 fspec='*_3p6um.fits'
 fspec = file_search(fspec)
endif

swfile = fspec
;-----------------------------------------------------
; some quick checks

if keyword_set(verbose) then print, '| Starting asinh image creation' 

if keyword_set(scales) then $
	scls=scales $
else	scls=[400.,300.,300.]
if n_elements(scls) ne 3 then begin $
 if keyword_set(verbose) then $
  print,'*** ERROR: scales must be a 3 element array: [red,green,blue].'
  print,'*** ERROR: default value is [1.,1.,1.]'
 return
endif

if not keyword_set(quality) then quality=75

if not keyword_set(brite) then brite=1.

if not keyword_set(nonlinearity) then nonlinearity=3.5

if not keyword_set(origin) then origin=[0,0,0]

if keyword_set(outdir) then $
	odir = outdir $
else	odir = strarr(n_elements(swfile))

if not keyword_set(greenmix) then greenmix=[.5, .5]
;-----------------------------------------------------
stpa = 1.D0 / ( (180.D0/!DPI)^2 * 3600.D0^2 )

for i=0,n_elements(swfile)-1 do begin 

 lwfile=repstr(swfile[i],'3p6um','4p5um')
 xfile=(reverse(strsplit(swfile[i],'/',/extract)))[0]
 xfile=repstr(xfile,'_3p6um','_3p6um4p5um')
 
 if keyword_set(swonly) then begin
	 re=(1 eq 0)
	 gg=(1 eq 0)
	 be=file_exist(swfile[i])
	 xfile=repstr(xfile,'_3p6um4p5um','_3p6um')
 endif else if keyword_set(lwonly) then begin
	 re=file_exist(lwfile)
	 gg=(1 eq 0)
	 be=(1 eq 0)
	 xfile=repstr(xfile,'_3p6um4p5um','_4p5um')
 endif else begin
	 re=file_exist(lwfile)
	 be=file_exist(swfile[i])
 endelse
 ;
 ; check update
 jname=odir[i]+'/'+gettok(xfile,'.fit')+'.jpg'
 if file_exist(jname) eq 1 and not keyword_set(update) then goto,skipall

 ;if be or re or gg then begin	; if fits files exist
 if be or re then begin	; if fits files exist

  if be then begin
    blue = mrdfits(swfile[i],0,bh,/fscale,/silent) * stpa * 1.e9; mJy / arcsec^2
    sky,blue[where(finite(blue))],bsky,bsig,/silent,/nan
    blue = blue - bsky
    hdr  = bh
  endif
  if re then begin
    red  = mrdfits(lwfile,0,rh,/fscale,/silent) * stpa * 1.e9	; mJy / arcsec^2
    sky,red[where(finite(red))],rsky,rsig,/silent,/nan
    red  = red - rsky
    hdr  = rh
  endif

  if not be and re then begin
	  blue = red
  endif else if not re and be then begin
	  red = blue
  endif

  if re and be then begin ; don't let NaNs cause a problem
   rednan=where(finite(red,/nan) and finite(blue),countrednan)
   bluenan=where(finite(blue,/nan) and finite(red),countbluenan)
   if countrednan gt 0 then red[rednan]=0
   if countbluenan gt 0 then blue[bluenan]=0
  endif

  green = greenmix[0]*blue + greenmix[1]*red

  nx=sxpar(hdr,'NAXIS1')
  ny=sxpar(hdr,'NAXIS2')
  getrot,hdr,rot,cdelt
  as_pix = abs(cdelt[0])*3600.

  ;----------------------------------
  ; build RGB composite

  sz=size(red)

  image=fltarr(sz[1],sz[2],3)
  image[*,*,0]=red * as_pix^2		; convert to mJy/px
  image[*,*,1]=green * as_pix^2
  image[*,*,2]=blue * as_pix^2

  image = nw_scale_rgb(image,scales=scls*brite)
  image = nw_arcsinh_fit(image,nonlinearity=nonlinearity)
  image = nw_fit_to_box(image,origin=origin)
  image = nw_float_to_byte(image)
  
  ;------------------------------------
  ; write file
  write_jpeg,jname,image,true=3,quality=quality
  if keyword_set(verbose) then print, '| Wrote '+jname

  set_plot,'x'

endif else if keyword_set(verbose) then $
	print,'Cannot find any file: ',swfile,' ',lwfile

skipall:

endfor 

end
