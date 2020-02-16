pro glga_batchptsrc,lfile,sdss=sdss,galex=galex,twomass=twomass,wise=wise, $
	irac=irac,verbose=verbose,logfile=logfile
;+
; glga_batchptsrc - generate mask files for objects listed in lfile
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
zp=0.
srvy='galex'
ptail='_NUV.fit*'	; primary file tail
atail='_FUV.fit*'	; alternate file tail
sigmult = 2.5
maglim = 25.0
if keyword_set(sdss) then begin
	srvy='sdss'
	ptail='_r.fit*'
	atail='_g.fit*'
	sigmult = 2.5
	maglim = 24.0
endif
if keyword_set(twomass) then begin
	srvy='2mass'
	ptail='_k.fit*'
	atail='_j.fit*'
	sigmult = 4.6
	maglim = 18.0
endif
if keyword_set(wise) then begin
	srvy='wise'
	ptail='_w1.fit*'
	atail='_w2.fit*'
	sigmult = 4.6
	maglim = 21.0
endif
if keyword_set(irac) then begin
	srvy='irac'
	ptail='_4p5um.fit*'
	atail='_3p6um.fit*'
	sigmult = 2.5
	maglim = 22.0
endif
;
; open log file
if keyword_set(logfile) then begin
	openw,ll,logfile,/get_lun
	printf,ll,'# GLGA_BATCHPTSRC: ' + systime(0)
	printf,ll,'# '+srvy+ $
		', SIGMULT = ' + string(sigmult,form='(f3.1)') + $
		', MAGLIM = ' + string(maglim, form='(f4.1)')
	printf,ll,'# NEW POINT SOURCE MASKS'
endif
;
; read in sample data
readcol,lfile, id, ra, dec, majdiam, mindiam, pa, type, $
	format='a,d,d,f,f,f,a', /silent, comment='#'
;
; deal with negative sizes
a = where(majdiam le 0.0, count)
if count gt 0 then majdiam[a] = 0.5
if count gt 0 then mindiam[a] = 0.5
a = where(mindiam le 0.0,  count)
if count gt 0 then mindiam[a] = majdiam[a]
;
rat = majdiam/mindiam
;
;convert diam into arcsec and expand by 50%
d=majdiam*60.0*1.5
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
  ellipsefile = outdir_aux+id[i]+'_ellipse.dat'
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
  nomask = 0 & nopsrcfile = 0 & noroifile = 0
  if not file_exist(mimgfile) then $
	  nomask = 1
  if not file_exist(psrcfile) then $
	  nopsrcfile = 1
;
; skip if already masked
  if not nomask then $
	  goto, skiperr
;
; log what happens
  qalogfile = filebase[i] + '_qa.txt'
  qa=glga_read_qa_stat(qalogfile)
;
; set qa status
  qa.psrc = 1
;
; read in image data to use
  im = mrdfits(usefile,0,hdr,/fscale,/silent)
  sz = [sxpar(hdr,'NAXIS1'),sxpar(hdr,'NAXIS2')]
;
; check image: allow up to 2/3 blanks
  blank = where( im eq min(im), nblank)
  if double(nblank) / double(n_elements(im)) gt 0.67 then begin
	  print,'Too many blanks in image'
	  goto, skiperr
  endif
;
; get zero point for WISE data
  if keyword_set(wise) then $
	  zp = sxpar(hdr,'MAGZP')
;
; get astrometry
  extast, hdr, astr
  ad2xy, ra[i], dec[i], astr, xn, yn
  getrot, hdr, rot, cdelt
  as_pix = abs(cdelt[0])*3600.
;
; get ellipse defining photometric apertures
;
; is there a file modifying the default ellipse?
  if file_exist(ellipsefile) then begin
	qa.ellipse_update = 1
	readcol,ellipsefile,majdiam_as,mindiam_as,el_ra,el_dec,pa_, $
		majdiam_px, mindiam_px, x0_,y0_, astrom_bool, el_as_pix, $
		format='f,f,d,d,f,f,f,f,i,f',/silent
	if astrom_bool ne 1 and keyword_set(verbose) then $
		print,'No astrometry, using default ellipse parameters' $
	else begin
		ad2xy,el_ra,el_dec,astr,x0_,y0_
		x=x0_[0]
		y=y0_[0]
		pa[i]=pa_[0]
		d[i]=majdiam_as[0]
		rat[i]=majdiam_as[0]/mindiam_as[0]
		if keyword_set(verbose) then $
			print,'Using ellipse info: ',ellipsefile
	endelse
  endif else begin
	  qa.ellipse_update = 0
	  x=xn
	  y=yn
	  if keyword_set(verbose) then $
	  	print,'No ellipse file, using default ellipse parameters'
  endelse
;
; construct ellipse vector
  el = [0.5*d[i]/as_pix, 0.5*(d[i]/rat[i]/as_pix), x, y, pa[i],xn,yn]
;
; find point sources
  glga_find_ptsrc,im,el,astr,as_pix,asrc,nobj,sigmult=sigmult,maglim=maglim, $
	galex=galex,sdss=sdss,twomass=twomass,wise=wise,irac=irac, $
	zpmag=zp, /verbose
;
; post-processing
  r = sqrt( (x-asrc[2,*])^2 + (y-asrc[3,*])^2 )
  g = where( r gt el[0]/5., gobj)
  if gobj gt 0 then begin
	openw,lun,psrcfile,/get_lun
	for j=0L,gobj-1 do begin
		p = g[j]
	 	if asrc[6,p] gt 0. then $
			printf,lun,asrc[0,p],asrc[1,p],asrc[2,p],asrc[3,p], $
				asrc[4,p]*as_pix,fix(asrc[5,p]), $
					format='(2f13.8,3f9.3,i5)'
	endfor
	free_lun,lun
;
; re-generate mask list file
	glga_genmask,sz,astr,as_pix,mimgfile,verbose=verbose
	qa.mask = 1
	nmsks = nmsks + 1L
;
; log new masks
	if keyword_set(logfile) then $
		printf,ll,id[i],ra[i],dec[i],majdiam[i],mindiam[i],pa[i], $
			type[i],form='(a-25,2f13.8,3f9.3,2x,a)'

  endif
;
; update qa log, but don't increment qa count
  glga_write_qa_stat,qa,/batch
;
; go here if we can't generate mask
  skiperr:

endfor	; loop over object list
if keyword_set(logfile) then free_lun,ll
;
print,' '
print,'Wrote ',nmsks,' mask files'
return
end
