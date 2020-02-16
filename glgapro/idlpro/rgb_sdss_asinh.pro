pro rgb_sdss_asinh, fspec=fspec,scales=scales,brite=brite, $
	nonlinearity=nonlinearity, origin=origin, verbose=verbose, $
	quality=quality, outdir=outdir, update=update, wide=wide, gri=gri,$
	uonly=uonly, gonly=gonly,ronly=ronly,ionly=ionly,zonly=zonly
;+
; NAME:
;  rgb_sdss_asinh
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
;  rgb_sdss_asinh, fspec=fspec,scales=scales,brite=brite, $
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
;  E. Kourkchi, The linear scaling of the final JPG image is changed to have a better contrast
;  - No constant value (e.g. 0.0002) is added to the background subtracted image
;  - This way more faint objects would be visulaised for a better masking Jun/22/2016
;
;
;-
; get file list
if not keyword_set(fspec) then fspec='*_i.fit*'
pfile = file_search(fspec)
;-----------------------------------------------------
; some quick checks

if keyword_set(verbose) then print, '| Starting asinh image creation' 

if keyword_set(scales) then $
	scls=scales $
else	if keyword_set(wide) then $
	scls=[300.,390.,1200.] $
else	scls=[300.,390.,600.]
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

 ufile=repstr(pfile[i],'_i','_u')
 gfile=repstr(pfile[i],'_i','_g')
 rfile=repstr(pfile[i],'_i','_r')
 ifile=pfile[i]
 zfile=repstr(pfile[i],'_i','_z')
 xfile=odir+'/'+(reverse(strsplit(pfile[i],'/',/extract)))[0]
 if keyword_set(wide) then $
	 xfile=repstr(xfile,'_i','_urz') $
 else    xfile=repstr(xfile,'_i','_gri')

 if keyword_set(ionly) then begin
 	re=file_exist(ifile)
	gg=(1 eq 0)
 	be=(1 eq 0)
	xfile=repstr(xfile,'_gri','_i')
	scls=[scls[0],scls[0],scls[0]]
 endif else if keyword_set(ronly) then begin
 	re=(1 eq 0)
 	gg=file_exist(rfile)
 	be=(1 eq 0)
	xfile=repstr(xfile,'_gri','_r')
	scls=[scls[1],scls[1],scls[1]]
 endif else if keyword_set(gonly) then begin
 	re=(1 eq 0)
	gg=(1 eq 0)
 	be=file_exist(gfile)
	xfile=repstr(xfile,'_gri','_g')
	scls=[scls[2],scls[2],scls[2]]
 endif else if keyword_set(uonly) then begin
 	re=(1 eq 0)
	gg=(1 eq 0)
 	be=file_exist(ufile)
	xfile=repstr(xfile,'_gri','_u')
	scls=[scls[2],scls[2],scls[2]]
 endif else if keyword_set(zonly) then begin
 	re=file_exist(zfile)
	gg=(1 eq 0)
 	be=(1 eq 0)
	xfile=repstr(xfile,'_gri','_z')
	scls=[scls[0],scls[0],scls[0]]
 endif else begin
 	if keyword_set(wide) then $
		re=file_exist(zfile) $
	else	re=file_exist(ifile)
 	gg=file_exist(rfile)
 	if keyword_set(wide) then $
		be=file_exist(ufile) $
	else	be=file_exist(gfile)
 endelse
 ;
 ; check update
 pdot = strpos(xfile,'.fit')
 if pdot ge 0 then $
	jname=strmid(xfile,0,pdot)+'.jpg' $
 else	jname=xfile+'.jpg'
 if file_exist(jname) eq 1 and not keyword_set(update) then goto,skipall

 if be or re or gg then begin  ; if fits files exist

  if be     then  begin
    if keyword_set(wide) or keyword_set(uonly) then $
	    blue = mrdfits(ufile,0,bh,/fscale,/silent) $
    else    blue = mrdfits(gfile,0,bh,/fscale,/silent)
    shf = sxpar(bh,'SHFACT')	; get shrink factor if present
    if shf gt 1.0 then blue = blue / shf^2
    mmm,blue[where(blue gt 0)],sky,sig
    if sig le 0. then $
    	sky,blue[where(blue gt 0)],sky,sig
    blue = (blue - sky) + 0.0002
    hdr = bh
  endif
  if gg     then  begin
    green  = mrdfits(rfile,0,gh,/fscale,/silent)
    shf = sxpar(gh,'SHFACT')	; get shrink factor if present
    if shf gt 1.0 then green = green / shf^2
    mmm,green[where(green gt 0)],sky,sig
    if sig le 0. then $
    	sky,green[where(green gt 0)],sky,sig
    green = (green - sky) + 0.0002
    hdr = gh
  endif
  if re     then  begin
    if keyword_set(wide) or keyword_set(zonly) then $
    	    red  = mrdfits(zfile,0,rh,/fscale,/silent) $
    else    red  = mrdfits(ifile,0,rh,/fscale,/silent)
    shf = sxpar(rh,'SHFACT')	; get shrink factor if present
    if shf gt 1.0 then red = red / shf^2
    mmm,red[where(red gt 0)],sky,sig
    if sig le 0. then $
    	sky,red[where(red gt 0)],sky,sig
    red = (red - sky) + 0.0002
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

  sz=size(red)

  image=fltarr(sz[1],sz[2],3)
  
;;;;;;;;;;;;;;; Added section by E. Kourkchi 
; if keyword_set(gri) then begin
;   a_min = 0.00009 ; 0.00001
;   a_max = 0.2 ; 0.3 default
;   alfa = (a_max/a_min)^(1./255)
;   
;   red   -= 0.0002
;   green -= 0.0002
;   blue  -= 0.0002
;   for p = 0, sz[1]-1 do begin
;      for q = 0, sz[2]-1 do begin
;          
;          if red[p,q] le a_min then red[p,q] = 0.
;          if red[p,q] ge a_max then red[p,q] = 1.
; ;          if red[p,q] gt 0. and red[p,q] le a_max then red[p,q] = alog(red[p,q]/a_min)/alog(alfa)/255.
;          if red[p,q] gt 0. and red[p,q] le a_max then red[p,q] = (red[p,q]-a_min)/(a_max-a_min)
;   
;          if green[p,q] le a_min then green[p,q] = 0.
;          if green[p,q] ge a_max then green[p,q] = 1.
; ;          if green[p,q] gt 0. and green[p,q] le a_max then green[p,q] = alog(green[p,q]/a_min)/alog(alfa)/255.
;          if green[p,q] gt 0. and green[p,q] le a_max then green[p,q] = (green[p,q]-a_min)/(a_max-a_min)
;   
;   
;          if blue[p,q] le a_min then blue[p,q] = 0.
;          if blue[p,q] ge a_max then blue[p,q] = 1.
; ;          if blue[p,q] gt 0. and blue[p,q] le a_max then blue[p,q] = alog(blue[p,q]/a_min)/alog(alfa)/255.  
;          if blue[p,q] gt 0. and blue[p,q] le a_max then blue[p,q] = (blue[p,q]-a_min)/(a_max-a_min)
;   
;   endfor
;       endfor 
; endif
;;;;;;;;;;;;;;; End - added section  

  image[*,*,0]=red
  image[*,*,1]=green
  image[*,*,2]=blue

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
  print,'Cannot find any file:',gfile,' ',rfile,' ',ifile

skipall:

endfor 

end
