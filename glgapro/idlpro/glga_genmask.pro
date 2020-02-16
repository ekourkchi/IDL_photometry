pro glga_genmask,sz,astr,as_pix,maskfile,verbose=verbose,filter=filter
;+
; glga_genmask - generate mask image files removing foreground point sources
;	and neighboring galaxies and artifacts
;
; sz		- image size [nx,ny]
; astr		- astrometry structure from image header (extast.pro)
; as_pix	- image scale in arcsec per pixel
; ptsrcfile	- point source file list
; roifile	- region of interest index list
; maskfile	- output mask image file
;-
; mask image
gim = glga_getmask(maskfile,sz,astr,as_pix,verbose=verbose,filter=filter)
;
return
end
