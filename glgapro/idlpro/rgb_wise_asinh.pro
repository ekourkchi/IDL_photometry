pro rgb_wise_asinh, fspec=fspec,scales=scales,brite=brite, $
	nonlinearity=nonlinearity, origin=origin, verbose=verbose, $
	quality=quality, outdir=outdir, update=update, wide=wide, $
	w1only=w1only, w2only=w2only, w3only=w3only, w4only=w4only
;+
; NAME:
;  rgb_wise_asinh
;
;
; PURPOSE:
;  Generate RGB composite jpeg images from WISE intensity maps.
;  RED=w3, BLUE=w1, GREEN=w2.
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
;  rgb_wise_asinh, fspec=fspec,scales=scales,brite=brite, $
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
;  D. Neill, yet more of a kludge for WISE data 14/04/2011
;
;
;-
; setup
pre = 'RGB_WISE_ASINH'
;
; get file list
if not keyword_set(fspec) then fspec='*_w1.fit*'
pfile = file_search(fspec)
;-----------------------------------------------------
; some quick checks

if keyword_set(verbose) then print, '| Starting asinh image creation' 

if keyword_set(scales) then $
	scl=scales $
else	scl=[0.07,0.18,0.07]
if n_elements(scl) ne 3 then begin $
 if keyword_set(verbose) then $
  print,'*** ERROR: scales must be a 3 element array: [red,green,blue].'
  print,'*** ERROR: default value is [1.,1.,1.]'
 return
endif

if not keyword_set(quality) then quality=100

if not keyword_set(brite) then brite=1.5

if not keyword_set(nonlinearity) then nonlinearity=2.5

if not keyword_set(origin) then origin=[0,0,0]

if keyword_set(outdir) then $
	odir = outdir $
else	odir = ''

;-----------------------------------------------------

for i=0,n_elements(pfile)-1 do begin 

 w1file=pfile[i]
 w2file=repstr(pfile[i],'_w1','_w2')
 w3file=repstr(pfile[i],'_w1','_w3')
 w4file=repstr(pfile[i],'_w1','_w4')
 xfile=odir+'/'+(reverse(strsplit(pfile[i],'/',/extract)))[0]
 if keyword_set(wide) then $
	 xfile=repstr(xfile,'_w1','_w124') $
 else    xfile=repstr(xfile,'_w1','_w123')
 scls=scl
 if keyword_set(w3only) then begin
 	be=(1 eq 0)
	gg=(1 eq 0)
 	re=file_exist(w3file)
	xfile=repstr(xfile,'_w123','_w3')
	scls=[scl[0],scl[0],scl[0]]
 endif else if keyword_set(w2only) then begin
 	re=(1 eq 0)
 	gg=file_exist(w2file)
 	be=(1 eq 0)
	xfile=repstr(xfile,'_w123','_w2')
	scls=[scl[1],scl[1],scl[1]]
 endif else if keyword_set(w1only) then begin
 	re=(1 eq 0)
	gg=(1 eq 0)
 	be=file_exist(w1file)
	xfile=repstr(xfile,'_w123','_w1')
	scls=[scl[2],scl[2],scl[2]]
 endif else if keyword_set(w4only) then begin
 	re=file_exist(w4file)
	gg=(1 eq 0)
 	be=(1 eq 0)
	xfile=repstr(xfile,'_w123','_w4')
	scls=[scl[0],scl[0],scl[0]]*12.0
 endif else begin
	be=file_exist(w1file)
 	gg=file_exist(w2file)
 	if keyword_set(wide) then begin
		re=file_exist(w4file)
		scls=[scl[0]*12.0,scl[1],scl[2]]
	endif else	re=file_exist(w3file)
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
    blue = mrdfits(w1file,0,bh,/fscale,/silent)
    good = where(blue gt 0,ngood)
    if ngood gt 0 then begin
    	mmm,blue[good],sky,sig
    	if sig le 0. then $
    		sky,blue[good],sky,sig
    	blue = blue - sky
    	hdr = bh
    endif else begin
	print,pre,' ERROR: no good pixels in file: ',w1file, $
		format='(a,a,a)'
	return
    endelse
  endif
  if gg     then  begin
    green  = mrdfits(w2file,0,gh,/fscale,/silent)
    good = where(green gt 0,ngood)
    if ngood gt 0 then begin
    	mmm,green[good],sky,sig
    	if sig le 0. then $
    		sky,green[good],sky,sig
    	green = green - sky
    	hdr = gh
    endif else begin
	print,pre,' ERROR: no good pixels in file: ',w2file, $
		format='(a,a,a)'
	return
    endelse
  endif
  if re     then  begin
    if keyword_set(wide) or keyword_set(w4only) then $
    	  red  = mrdfits(w4file,0,rh,/fscale,/silent) $
    else  red  = mrdfits(w3file,0,rh,/fscale,/silent)
    good = where(red gt 0,ngood)
    if ngood gt 0 then begin
    	mmm,red[good],sky,sig
    	if sig le 0. then $
    		sky,red[good],sky,sig
    	red = red - sky
    	hdr = rh
    endif else begin
	if keyword_set(wide) or keyword_set(w4only) then begin
		print,pre,' ERROR: no good pixels in file: ',w4file, $
			format='(a,a,a)'
	endif else begin
		print,pre,' ERROR: no good pixels in file: ',w3file, $
			format='(a,a,a)'
	endelse
	return
    endelse
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
  print,'Cannot find any file:',w2file,' ',w3file,' ',w4file

skipall:

endfor 

end
