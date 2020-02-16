pro rgb_2mass_asinh, fspec=fspec,scales=scales,brite=brite,noskysub=noskysub, $
                     nonlinearity=nonlinearity, origin=origin, $
                     verbose=verbose, quality=quality, outdir=outdir, $
		     update=update, jonly=jonly, honly=honly, konly=konly
;+
; NAME:
;  rgb_2mass_asinh
;
;
; PURPOSE:
;  Generate RGB composite jpeg images from GALEX FUV and NUV intensity maps.
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
;  rgb_2mass_asinh, fspec=fspec,scales=scales,brite=brite, $
;                     nonlinearity=nonlinearity, origin=origin, $
;                     verbose=verbose, quality=quality
;
; INPUTS:
;          
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;  noskysub - set to disable sky subtraction of each image
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
if not keyword_set(fspec) then fspec='*_k.fits.gz'
pfile = file_search(fspec)
;-----------------------------------------------------
; some quick checks

if keyword_set(verbose) then print, '| Starting asinh image creation' 

if keyword_set(scales) then $
	scls=scales $
else	scls=[12.5, 9.5, 10.5]	; k, h, j --> r, g, b
if n_elements(scls) ne 3 then begin $
 if keyword_set(verbose) then $
  print,'*** ERROR: scales must be a 3 element array: [red,green,blue].'
  print,'*** ERROR: default value is [1.,1.,1.]'
 return
endif

if not keyword_set(quality) then quality=75

if not keyword_set(brite) then brite=1.

if not keyword_set(nonlinearity) then nonlinearity=2.5

if not keyword_set(origin) then origin=[0,0,0]

if keyword_set(outdir) then $
	odir = outdir $
else	odir = ''

;-----------------------------------------------------

for i=0,n_elements(pfile)-1 do begin 

 jfile=repstr(pfile[i],'_k','_j')
 hfile=repstr(pfile[i],'_k','_h')
 kfile=pfile[i]
 xfile=odir+'/'+(reverse(strsplit(pfile[i],'/',/extract)))[0]
 xfile=repstr(xfile,'_k','_jhk')

 if keyword_set(konly) then begin
 	re=file_exist(kfile)
	gg=(1 eq 0)
 	be=(1 eq 0)
	xfile=repstr(xfile,'_jhk','_k')
	scls=[scls[0],scls[0],scls[0]]
 endif else if keyword_set(honly) then begin
 	re=(1 eq 0)
 	gg=file_exist(hfile)
 	be=(1 eq 0)
	xfile=repstr(xfile,'_jhk','_h')
	scls=[scls[1],scls[1],scls[1]]
 endif else if keyword_set(jonly) then begin
 	re=(1 eq 0)
	gg=(1 eq 0)
 	be=file_exist(jfile)
	xfile=repstr(xfile,'_jhk','_j')
	scls=[scls[2],scls[2],scls[2]]
 endif else begin
	re=file_exist(kfile)
	gg=file_exist(hfile)
	be=file_exist(jfile)
 endelse
 ;
 ; check update
 pdot = strpos(xfile,'.fit')
 if pdot ge 0 then $
	jname=strmid(xfile,0,pdot)+'.jpg' $
 else	jname=xfile+'.jpg'
 if file_exist(jname) and not keyword_set(update) then goto,skipall

 if be or re or gg then begin  ; if fits files exists

  if be     then  begin
    blue = mrdfits(jfile,0,bh,/fscale,/silent)
    shf = sxpar(bh,'SHFACT')    ; get shrink factor if present
    if shf gt 1.0 then blue = blue / shf^2
    if not keyword_set(noskysub) then begin
    	sky,blue,bsky,/silent
    	blue=blue-bsky
    endif
    hdr = bh
  endif
  if gg     then  begin
    green  = mrdfits(hfile,0,gh,/fscale,/silent)
    shf = sxpar(gh,'SHFACT')    ; get shrink factor if present
    if shf gt 1.0 then green = green / shf^2
    if not keyword_set(noskysub) then begin
    	sky,green,gsky,/silent
    	green=green-gsky
    endif
    hdr = gh
  endif
  if re     then  begin
    red  = mrdfits(kfile,0,rh,/fscale,/silent)
    shf = sxpar(rh,'SHFACT')    ; get shrink factor if present
    if shf gt 1.0 then red = red / shf^2
    if not keyword_set(noskysub) then begin
    	sky,red,rsky,/silent
    	red=red-rsky
    endif
    hdr = rh
  endif
  if not be and re then begin 
    blue = red
  endif else if not be and gg then begin
    blue = green
  endif
  if not re and gg then begin 
    red = green
  endif else if not re and be then begin
    red = blue
  endif
  if not gg and re then begin 
    green = red
  endif else if not gg and be then begin
    green = blue
  endif

  nx=sxpar(hdr,'NAXIS1')
  ny=sxpar(hdr,'NAXIS2')


  ;----------------------------------
  ; build RGB composite

  image=fltarr(nx,ny,3)
  image[*,*,0]=red
  image[*,*,1]=green
  image[*,*,2]=blue

  image = nw_scale_rgb(image,scales=scls*brite)
  image = nw_arcsinh_fit(image,nonlinearity=nonlinearity)
  image = nw_fit_to_box(image,origin=origin)
  image = nw_float_to_byte(image)
  
  ;------------------------------------
  ; write files


  write_jpeg,jname,image,true=3,quality=quality
  if keyword_set(verbose) then print, '| Wrote '+jname

  set_plot,'x'

endif else if keyword_set(verbose) then $
  print,'Cannot find any file:',hfile,jfile,kfile

skipall:

endfor 

end
