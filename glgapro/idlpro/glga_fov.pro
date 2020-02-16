pro glga_fov, lfile, galex=galex,sdss=sdss,twomass=twomass,$
	wise=wise,irac=irac
;+
; glga_fov - test image size against outer sky annulus to see if expansion of
;		the image is required.
;
; list of objects with these columns:
;       id	- string
;       ra,dec  - degrees
;       majdiam,mindiam - arcmin
;       pa      - degrees
;	type	- string
;
;-
; set defaults and check keywords
band = 'NUV'
srvy='galex'		; data type
if keyword_set(sdss) then begin
	band='r'
	srvy='sdss'
endif
if keyword_set(twomass) then begin
	band = 'k'
	srvy='2mass'
endif
if keyword_set(wise) then begin
	band = 'w1'
	srvy='wise'
endif
if keyword_set(irac) then begin
	band = '3p6um'
	srvy='irac'
endif
;
; read in sample data
readcol, lfile, id, ira, idec, mjx, mnx, ipa, ty, format='a,d,d,f,f,f,a'
;
; output file
ofile = strmid(lfile,0,strpos(lfile,'.')) + '_expnd_'+srvy+'.glga'
filestamp,ofile
openw,ol,ofile,/get_lun
printf,ol,'# GLGA_FOV - '+systime(0)
printf,ol,'# input list: '+lfile
;
; define top level directory
deg = string(floor(ira), format='(i3.3)')+'D'
;
; loop over object list
for i=0,n_elements(id)-1 do begin
;
; input files
	profpath=!GLGA_ROOT+'/data/'+deg[i]+'/photometry/'
	ellfile = profpath+id[i]+'_'+band+'_ellipsepar.dat'
	skyfile = profpath+id[i]+'_'+band+'_background.dat'
	fpath=!GLGA_ROOT+'data/'+deg[i]+'/'+srvy+'/fits/'
	imgfile = fpath+id[i]+'_'+band+'.fit*'
	;
	; get input parameters
	if file_test(ellfile,/read) and file_test(skyfile,/read) then begin
		readcol,ellfile,ra,dec,semimajor,semiminor,pa,/silent
		readcol,skyfile,sflx,sflxe,smu,smue,sscale,sradi,srado,/silent
		;
		; output values
		mjxo = semimajor[0]/30.	; convert to diameter in arcmin
		mnxo = semiminor[0]/30.
		;
		; use outer sky annulus as test
		ratio = semiminor/semimajor
		semimajor = srado
		semiminor = srado * ratio
		;
		; get fits file header
		flist = file_search(imgfile,count=nf)
		if nf eq 1 then begin
			hdr = headfits(flist[0])
			extast,hdr,astr
			ad2xy,ra,dec,astr,x0,y0
			sz = [sxpar(hdr,'naxis1'),sxpar(hdr,'naxis2')]
			as_pix = abs(astr.cdelt[0])*3600.d0
			pos_ang = (90. - pa[0])/!RADEG	; to radians
			semimjpx = semimajor[0]/as_pix	; to pixels
			semimnpx = semiminor[0]/as_pix
			npoints = 1000
			phi = 2.*!pi*(findgen(npoints)/(npoints-1))
			cosang = cos(pos_ang)
			sinang = sin(pos_ang)

			x = semimjpx*cos(phi)
			y = semimnpx*sin(phi)

			xp = x0[0] + x*cosang - y*sinang
			yp = y0[0] + x*sinang + y*cosang
			print,i+1,'/',n_elements(id),': ',id[i], $
				sz,minmax(xp),minmax(yp), $
				format='(i6,a1,i6,a,a-25,2x,2i6,4f9.1)'
			;
			; now test dimensions
			if min(xp) lt 0 or min(yp) lt 0 or $
			   max(xp) ge sz[0] or max(yp) ge sz[1] then begin
			   	printf,ol,id[i],ira[i],idec[i], $
					mjxo,mnxo,ipa[i],ty[i], $
					format='(a-25,2f13.8,3f9.3,2x,a)'
			   	print,id[i],ira[i],idec[i], $
					mjxo,mnxo,ipa[i],ty[i], $
					format='(a-25,2f13.8,3f9.3,2x,a)'
			endif
		endif else begin
			print,i+1,'/',n_elements(id),': ',id[i]
			print,'input fits file not found: ',imgfile
		endelse
	endif else begin
		print,i+1,'/',n_elements(id),': ',id[i]
		if not file_test(ellfile,/read) then $
			print,'ellipse file not found: ',ellfile
		if not file_test(skyfile,/read) then $
			print,'background file not found: ',skyfile
	endelse
endfor

free_lun,ol

return
end
