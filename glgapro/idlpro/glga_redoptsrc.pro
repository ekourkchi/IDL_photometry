pro glga_redoptsrc,lfile,sdss=sdss,galex=galex,twomass=twomass,wise=wise, $
	irac=irac,verbose=verbose,logfile=logfile
;+
; glga_redoptsrc - redo point source radii (with new psf, algorthm, etc.)
;
; lfile - list of objects with one line per object with these columns:
;	id	- string
;	ra,dec	- degrees
;	majdiam,mindiam	- arcmin
;	pa	- degrees
;	type	- string
;
; keywords:
;	sdss,galex,twomass,wise,irac - data to analyze
;	verbose - print out more details
;	logfile - file to log new point source files
;
;-
; set defaults and check keywords
srvy='galex'
ptail='_NUV.fit*'	; primary file tail
atail='_FUV.fit*'	; alternate file tail
if keyword_set(sdss) then begin
	srvy='sdss'
	ptail='_r.fit*'
	atail='_g.fit*'
endif
if keyword_set(twomass) then begin
	srvy='2mass'
	ptail='_k.fit*'
	atail='_j.fit*'
endif
if keyword_set(irac) then begin
	srvy='irac'
	ptail='_4p5um.fit*'
	atail='_3p6um.fit*'
endif
if keyword_set(wise) then begin
	srvy='wise'
	ptail='_w1.fit*'
	atail='_w2.fit*'
endif
;
; open log file
if keyword_set(logfile) then begin
	openw,ll,logfile,/get_lun
	printf,ll,'# GLGA_REDOPTSRC: ' + systime(0)
	printf,ll,'# '+srvy
	printf,ll,'# NEW POINT SOURCE RADII'
endif
;
; read in sample data
readcol,lfile, id, ra, dec, majdiam, mindiam, pa, type, $
	format='a,d,d,f,f,f,a', /silent, comment='#'
;
; define top level directory
deg = string(floor(ra), format='(i3.3)')+'D'
;
; used to generate file names
id = strcompress(id,/rem)

filebase=!GLGA_ROOT+'data/'+deg+'/'+srvy+'/fits/'+id
;
nloop = n_elements(id)
nmsks = 0L
;
; loop over object list
for i=0L, nloop-1 do begin
;
; print host
  print,i+1,'/',nloop,id[i],deg[i], $
	  format = '(i6,a1,i6,2x,a-25,a5)'
;
; directories
  outdir_aux = !GLGA_ROOT+'data/'+deg[i]+'/aux/'  
;
; generate filenames
  pfile=filebase[i]+ptail	; primary file
  afile=filebase[i]+atail	; alternate file
  psrcfile = outdir_aux+id[i]+'_'+srvy+'_pointsrc.dat'
  roifile = outdir_aux+id[i]+'_'+srvy+'_roi.dat'
  mimgfile = outdir_aux+id[i]+'_'+srvy+'_mask.fits.gz'
;
; check primary and alternate files
  nopfile = 0 &  noafile = 0
  if not file_exist(pfile) then $
	  nopfile = 1
  if not file_exist(afile) then $
	  noafile = 1
  if nopfile then $
	usefile = afile $
  else	usefile = pfile
;
; skip if neither is available
  if nopfile and noafile then begin
	  print,'No images found'
	  goto,skiperr
  endif
;
; check mask, and point source files
  nopsrcfile = 0
  if not file_exist(psrcfile) then $
	  nopsrcfile = 1
;
; skip if no point source file
  if nopsrcfile then $
	  goto, skiperr
;
; read in image data to use
  im = mrdfits(usefile,0,hdr,/fscale,/silent)
  sz = size(im,/dim)
;
; get astrometry
  extast, hdr, astr
  ad2xy, ra[i], dec[i], astr, xn, yn
  getrot, hdr, rot, cdelt
  as_pix = abs(cdelt[0])*3600.
;
; read point sources
  readcol,psrcfile, a, d, xin, yin, rad, psc_astr, $
	format='d,d,f,f,f,i',/silent
  filestamp,psrcfile
  openw,ol,psrcfile,/get_lun
  printf,ol,'# STV centroid output (a,d,x,y,r_asec,astr?): '+systime(0)
;
; convert ra, dec to x,y px
  ad2xy,a,d,astr,xcen,ycen ; convert ra,dec to x,y px
;
; loop over point sources
  for j=0l,n_elements(a)-1 do begin
;
; are we inside the image?
    if xcen[j] ge 0 and xcen[j] lt sz[0] and $
       ycen[j] ge 0 and ycen[j] lt sz[1] then begin
;
; use subim to speed things up
	x0=fix(xcen[j]-63.)>0
	x1=x0+127<(sz[0]-1)
	y0=fix(ycen[j]-63.)>0
	y1=y0+127<(sz[1]-1)
	; get mask radius in pixels
	radius = get_mask_radius(im[x0:x1,y0:y1], $
		xcen[j]-float(x0),ycen[j]-float(y0),srvy)
	; convert to radius in arcsec
	printf,ol,a[j],d[j],xcen[j],ycen[j],radius*as_pix,psc_astr[j], $
			format='(2f13.8,3f9.3,i5)'
     endif else print,'Outside image: ',xcen[j],ycen[j],sz[0],sz[1]
  endfor
  free_lun,ol
  nmsks = nmsks + 1L
;
; re-generate mask list file
  glga_genmask,sz,astr,as_pix,mimgfile,verbose=verbose
;
; log new masks
  if keyword_set(logfile) then $
	printf,ll,id[i],ra[i],dec[i],majdiam[i],mindiam[i],pa[i], $
		type[i],form='(a-25,2f13.8,3f9.3,2x,a)'
;
; go here if we can't generate mask
  skiperr:

endfor	; loop over object list
if keyword_set(logfile) then free_lun,ll
;
print,' '
print,'Wrote ',nmsks,' new point source files'
return
end
