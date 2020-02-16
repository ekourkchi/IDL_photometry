pro rgb_galex_asinh, fuvfile, nuvfile, target, $
                     bgsub=bgsub,bgvals=bgvals,$
                     greenmix=greenmix, scales=scales,brite=brite, $
                     nonlinearity=nonlinearity, origin=origin, $
                     thumb=thumb, verbose=verbose, $
                     quality=quality,tif=tif,fuvshift=fuvshift, $
                     nofuvsmooth = nofuvsmooth, nuvsmooth = nuvsmooth, $
                     bgouter = bgouter

;+
; NAME:
;  rgb_galex_asinh
;
;
; PURPOSE:
;  Generate RGB composite jpeg images from GALEX FUV and NUV intensity maps.
;  RED=NUV, BLUE=FUV, GREEN=combination of FUV & NUV.
;
;  Images are lambda*flamda and asinh scaled/fit using Wherry,
;  Blanton, Hogg IDL routines of the the Lupton, et al. algorithm.
;
; CATEGORY:
;  image processing
;
;
; CALLING SEQUENCE:
;  rgb_galex_fields, fuv, nuv, target,
;                   greenmix=greenmix, scales=scales,
;                   brite=brite, nonlinearity=nonlinearity, origin=origin
;
; INPUTS:
;  fuvfile=FUV intensity fits file(s) (..._FUV.fits)
;  nuvfile=NUV intensity fits file(s) (..._NUV.fits)
;  target =name(s) of jpeg files to create and title of image
;          Field name is recommended.
;          
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;  bgsub = subtract a 'global' background level
;  greenmix = blue and red combination level for green where red and
;             blue are the lambda*flamda scaled
;             values. Default=[.2,.8]
;  brite = scale factor applied to all three RGB planes. Default=1
;          if bgsub then brite=(1+(sxpar(hdr,'EXPTIME')/1e4))<3,
;          setting brite will override this.
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
; MODIFICATION HISTORY:
;  M. Seibert 12/24/2004
;  M. Seibert 12/06/2005
;  D. Neill   02/10/2013 - removed unused code
;-

;-----------------------------------------------------
; some quick checks

if keyword_set(verbose) then print, '| Starting asinh image creation' 

if n_params() lt 3 then begin
 print,'*** ERROR: Requires at least 3 parameters.'
 return
endif

if n_elements(nuvfile) ne n_elements(fuvfile) or $
   n_elements(target) ne n_elements(fuvfile) then begin $
  if keyword_set(verbose) then $
    print,'*** ERROR: Image and target lists must have same number of elements.'
  return
endif

if not keyword_set(greenmix) then greenmix=[.2,.8];[.5,.5]
if n_elements(greenmix) ne 2 then begin $
 if keyword_set(verbose) then $
  print,'*** ERROR: greenmix must be a 2 element array: [fraction_FUV,fraction_NUV].'
 return
endif


if not keyword_set(scales) then scales=[.085,.09,.085]
if n_elements(scales) ne 3 then begin $
 if keyword_set(verbose) then $
  print,'*** ERROR: scales must be a 3 element array: [red,green,blue].'
  print,'*** ERROR: default value is [.085,.09,.085].'
 return
endif

if not keyword_set(quality) then quality=75


;-----------------------------------------------------

for i=0,n_elements(fuvfile)-1 do begin 

 be=file_exist(fuvfile[i])
 re=file_exist(nuvfile[i])

 if be or re then begin  ; if fits files exists

  if be     then  begin
    blue = mrdfits(fuvfile[i],0,bh)
    blue= blue*1.4*1525  ; scale as lambda*flambda
    if not keyword_set(nofuvsmooth) then $
     blue = smooth(blue,2,/nan)>0 ; lightly smooth the FUV image
    if keyword_set(fuvshift) then blue=shift(blue,fuvshift) 
  endif
  if re     then  begin
    red  = mrdfits(nuvfile[i],0,rh)
    red=  red*0.206*2297
    if keyword_set(nuvsmooth) then red = smooth(red,2,/nan)>0
  endif
  if re and be then hdr = rh
  if not be then begin 
    blue = red
    hdr  = rh
    greenmix=[0.5,0.5]
    scales=[.085,.085,.085]    
  endif
  if not re then  begin
    red  = blue
    hdr  = bh
    greenmix=[0.5,0.5]
    scales=[.085,.085,.085]
  endif


  if keyword_set(bgsub) then begin

   sz = size(red)
   mindim = sz[1] < sz[2]
   skycirc = fix(0.73*(mindim/2.))

   if not keyword_set(bgvals) and not keyword_set(bgouter) then begin
    if keyword_set(verbose) then print, '| Background subtracting' 
    sky,red,   rsky, rskysigma,circ=skycirc, /silent; quick sky & sigma
    sky,blue,  bsky, bskysigma,circ=skycirc, /silent
   
    if rsky le 0 or bsky le 0 then begin
     if keyword_set(verbose) then print, '| Computing alternative background' 
     dist_circle, circle, mindim, mindim/2, mindim/2
     area=where(circle le skycirc)
     delvarx,circle
     meanclip,red[area],rsky,clip,clipsig=3
     meanclip,blue[area],bsky,clip,clipsig=3
     if keyword_set(verbose) then print, '| ',bsky, rsky
    endif
   endif

   if keyword_set(bgouter) then begin
     dist_circle, circle, mindim, mindim/2, mindim/2
     area=where(circle gt skycirc)
     delvarx,circle
     meanclip,red[area],rsky,clip,clipsig=3
     meanclip,blue[area],bsky,clip,clipsig=3
     if keyword_set(verbose) then print, '| ',bsky, rsky
   endif

   if keyword_set(bgvals) then begin
    bsky=bgvals[0]
    rsky=bgvals[1]
    if keyword_set(verbose) then print, '| ',bsky, rsky
   endif

   rsky = rsky>0
   bsky = bsky>0

   if keyword_set(verbose) then print, '| ',bsky, rsky

   red=(temporary(red)-rsky)>0
   blue=(temporary(blue)-bsky)>0

  endif

  ;----------------------------------
  ; define green color

  green=(greenmix[0]*blue + greenmix[1]*red)


  ;----------------------------------
  ; build RGB composite

  sz=size(red)


  image=fltarr(sz[1],sz[2],3)
  image[*,*,0]=red
  image[*,*,1]=green
  image[*,*,2]=blue



  if not keyword_set(brite) then begin
    brite=1
    if keyword_set(bgsub) then brite=(1+(sxpar(hdr,'EXPTIME')/1e4))<3
  endif 
  if not keyword_set(nonlinearity) then nonlinearity=2.5
  if not keyword_Set(origin) then origin=[0,0,0]

  image = nw_scale_rgb(image,scales=scales*brite)
  image = nw_arcsinh_fit(image,nonlinearity=nonlinearity)
  image = nw_fit_to_box(image,origin=origin)
  image = nw_float_to_byte(image)
  
  set_plot,'x'

  ;------------------------------------
  ; write files

  jname=target[i]+'.jpg'
  write_jpeg,jname,image,true=3,quality=quality
  if keyword_set(verbose) then print, '| Wrote '+jname

endif else if keyword_set(verbose) then $
  print,'Cannot find either file:',fuvfile[i],nuvfile[i] 

endfor 

end
