pro glga_batchmask,lfile,sdss=sdss,galex=galex,twomass=twomass,$
	wise=wise,irac=irac,verbose=verbose
;+
; glga_batchmask - generate mask image files for objects listed in lfile
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
if keyword_set(wise) then begin
	srvy='wise'
	ptail='_w1.fit*'
	atail='_w2.fit*'
endif
if keyword_set(irac) then begin
	sryv='irac'
	ptail='_3p6um.fit*'
	atail='_4p5um.fit*'
endif
;
; read in sample data
readcol,lfile, id, ra, dec, format='a,d,d', /silent
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
  usefile = pfile
  if nopfile then usefile = afile
;
; skip if neither is available
  if nopfile and noafile then $
	  goto,skiperr
;
; check pointsource and ROI files
  nopsrcfile = 0 & noroifile = 0
  if not file_exist(psrcfile) then $
	  nopsrcfile = 1
  if not file_exist(roifile) then $
	  noroifile = 1
;
; skip if neither is available
  if nopsrcfile and noroifile then $
	  goto, skiperr
;
; log what happens
  qalogfile = filebase[i] + '_qa.txt'
  qa=glga_read_qa_stat(qalogfile)
;
; set qa status
  if nopsrcfile then $
	  qa.psrc = 0 $
  else	  qa.psrc = 1
  if noroifile then $
	  qa.roi = 0 $
  else	  qa.roi = 1
;
; read in image data to use
  hdr = headfits(usefile,/silent)
  sz = [sxpar(hdr,'NAXIS1'),sxpar(hdr,'NAXIS2')]
;
; get astrometry
  extast, hdr, astr
  ad2xy, ra[i], dec[i], astr, x, y
  getrot, hdr, rot, cdelt
  as_pix = abs(cdelt[0])*3600.
;
; check ROI file if needed
  if not noroifile then begin
	  ; get image dimensions
	  spawn,'grep IMAGE '+roifile,sout
	  ; do we have the requisite record?
	  if strlen(sout[0]) gt 0 then begin
		  rx = 0L & ry = 0L
		  ; us last (most recent) image size
		  for ii=0l,n_elements(sout)-1 do begin
			  cmnt = sout[ii]
			  for j=0,2 do jnk = gettok(cmnt,' ')
			  rx = long(gettok(cmnt,' '))
			  ry = long(gettok(cmnt,' '))
		  endfor
		  ; read in roi data
		  readcol,roifile,indx,form='l',/silent
		  ; do we need to adjust roi?
		  ; we expanded
		  if sz[0] gt rx then begin
			  x0 = (sz[0] - rx) / 2
			  y0 = (sz[1] - ry) / 2
			  x1 = x0 + rx - 1
			  y1 = y0 + ry - 1
			  large = bytarr(sz)
			  small = bytarr(rx,ry)
			  small[indx] = 1b
			  large[x0:x1,y0:y1] = small[*,*]
			  indx = where(large eq 1,nindx)
		  endif
		  ; we shrank
		  if sz[1] lt rx then begin
			  x0 = (rx - sz[0]) / 2
			  y0 = (ry - sz[1]) / 2
			  x1 = x0 + sz[0] - 1
			  y1 = y0 + sz[1] - 1
			  large = bytarr(rx,ry)
			  large[indx] = 1b
			  small = large[x0:x1,y0:y1]
			  indx = where(small eq 1)
		  endif
		  ; re-write file with new img dims
		  print,'Re-writing ROI file: '+roifile
		  openw,lun,roifile,/get_lun
		  printf,lun,'# STV ROI output (image index): '+systime(0)
		  printf,lun,'# IMAGE NX,NY: ',sz[0],sz[1]
		  for ii=0L,n_elements(indx)-1L do $
			printf,lun,indx[ii]
		  free_lun,lun
	  endif
  endif

;
; re-generate mask image file
  glga_genmask,sz,astr,as_pix,mimgfile,verbose=verbose
  qa.mask = 1
  nmsks = nmsks + 1L
;
; update qa log, but don't increment qa count
  glga_write_qa_stat,qa,/batch
;
; go here if we can't generate mask
  skiperr:

endfor	; loop over object list
;
print,' '
print,'Wrote ',nmsks,' mask files'
return
end
