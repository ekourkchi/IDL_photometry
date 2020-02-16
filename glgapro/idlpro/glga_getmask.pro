function glga_getmask,maskfile,sz,astr,as_pix,verbose=verbose,update=update, filter=filter
;+
; glga_getmask - return a byte image mask with foreground point sources
;	and neighboring galaxies and artifacts masked with a value of 1
;
; maskfile	- mask image filename
; sz		- image size [nx,ny]
; astr		- astrometry structure from image header (extast.pro)
; as_pix	- image scale in arcsec per pixel
; ptsrcfile	- point source file list
; roifile	- region of interest index list
; update	- set to force update of mask
;-
;
; get point source and region files
rute = strmid(maskfile,0,strpos(maskfile,'.fit'))
ptsrcfile = repstr(rute,'mask','pointsrc')+'.dat'
roifile = repstr(rute,'mask','roi')+'.dat'
;
; is the mask image current?
img_current = (1 eq 1)	; assume yes
;
; get information on mask, pointsrc, and roi files
mfli = file_info(maskfile)
pfli = file_info(ptsrcfile)
rfli = file_info(roifile)
;
; does the mask image exist?
if mfli.exists then begin
;
; read in header
	ghdr = headfits(maskfile)
;
; have image dimensions changed?
	if sxpar(ghdr,'NAXIS1') ne sz[0] or sxpar(ghdr,'NAXIS2') ne sz[1] then $
		img_current = (1 eq 0)
;
; if it exists, have there been any updates to the pointsrc or roi file?
	if pfli.mtime gt mfli.mtime or rfli.mtime gt mfli.mtime then $
		img_current = (1 eq 0)
endif else	img_current = (1 eq 0)	; mask doesn't exist, not current
;
; if image is current, just read it in
if img_current and not keyword_set(update) then begin
	gim = mrdfits(maskfile,0,/silent)
endif else begin
;
; image is not current or doesn't exist, so regenerate
	gim = bytarr(sz)	; mask image
	did_mask = (1 eq 0)	; no masking yet
	if keyword_set(verbose) then print,' '
;
; mask point sources
	if pfli.exists then begin
;
; read list
		readcol, ptsrcfile, a, d, xin, yin, r_as, psc_astr, $
			format='d,d,f,f,f,i',/silent
;
; clean out baddies
		good = where(finite(a) eq 1 and finite(d) eq 1 and $
			     finite(xin) eq 1 and finite(yin) eq 1, ngood)
;
; any good ones left?
		if ngood gt 0 then begin
			a = a[good]
			d = d[good]
			xin = xin[good]
			yin = yin[good]
			r_as= r_as[good]
			psc_astr = psc_astr[good]
;
; convert radii to pixels
			r_px = r_as / as_pix
;
; max r
			rmax = fix(max(r_px)) + 10
;
; loop over point source list
			if keyword_set(verbose) then $
			     print,'Masking ',ngood,' point sources from: ', $
				ptsrcfile, format='(a,i5,a,a)'

			for i = 0L, ngood-1L do begin
;
; if astrometry good, convert to pix coords
				if psc_astr[i] eq 1 then begin
					AD2XY,a[i],d[i],astr,x,y
;
; otherwise use input x,y
				endif else begin
					x = xin[i]
					y = yin[i]
				endelse
;
; now check position in image
				if x gt 0. and x lt sz[0]-1. and $
				   y gt 0. and y lt sz[1]-1. then begin
;
; get subim limits
					sim = bytarr(2*rmax,2*rmax)
					xs0 = 0
					xs1 = 2*rmax - 1
					ys0 = 0
					ys1 = 2*rmax - 1
					x0 = (fix(x) - rmax) > 0
					y0 = (fix(y) - rmax) > 0
					x1 = (x0 + 2*rmax) - 1
					if x1 ge sz[0] then begin
						delx = x1 - sz[0]
						xs1 = xs1 - delx - 1
						x1 = sz[0] - 1
					endif
					y1 = (y0 + 2*rmax) - 1
					if y1 ge sz[1] then begin
						dely = y1 - sz[1]
						ys1 = ys1 - dely - 1
						y1 = sz[1] - 1
					endif
					sim[xs0:xs1,ys0:ys1] = gim[x0:x1,y0:y1]
;
; get index of masks
					dist_circle,im,2*rmax, $
						x-float(x0),y-float(y0)
					w = where(im le r_px[i])
					sim[w] = 1B
;
; set mask with using 'or' ensuring we don't unflag previously flagged pixels
					gim[x0:x1,y0:y1] = gim[x0:x1,y0:y1] or $
						   	sim[xs0:xs1,ys0:ys1]
;
; end check position in image
				endif
;
			endfor	; loop over point source list
			did_mask = (1 eq 1)	; masking was done
		endif	; any good point sources to mask?
	endif	; pfli.exists

; regions
	if rfli.exists then begin
		if keyword_set(verbose) then $
			print,'Masking regions from: ',roifile
		readcol,roifile,w,format='l',/silent
		gim[w] = 1B
		did_mask = (1 eq 1)
	endif
;
;write out mask list file
;
; are there any to mask?
	if did_mask then begin
		writefits,maskfile,gim,/compress
		if keyword_set(verbose) then $
			print,'Wrote mask img file: ',maskfile
	endif else $
		if keyword_set(verbose) then $
			print,'Nothing to mask'

endelse



if keyword_set(filter) then begin
    
   for i=0, n_elements(filter)-1 do begin
   
      roifile = repstr(rute,'mask','roi')+'_'+filter[i]+'.dat'
      maskfile = rute+'_'+filter[i]+'.fits.gz'

      rfli = file_info(roifile)
      mfli = file_info(maskfile)
      
      gim = bytarr(sz)
      if rfli.exists  then begin
        readcol,roifile,w,format='l',/silent
        gim[w] = 1B
        writefits,maskfile,gim,/compress
      endif else begin
        if mfli.exists then begin
	  spawn, 'rm ' + maskfile
        endif
      endelse

   endfor
   
   
endif




return,gim
end
