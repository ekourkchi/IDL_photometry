pro glga_batchqaplots, lfile, galex=galex,sdss=sdss,twomass=twomass,$
	wise=wise,irac=irac,verbose=verbose
;+
; glga_batchqaplots - plot QA version of radial profile photometry and images
;	of objects in lfile
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
srvy='galex'		; data type
if keyword_set(sdss) then srvy='sdss'
if keyword_set(twomass) then srvy='2mass'
if keyword_set(wise) then srvy='wise'
if keyword_set(irac) then srvy='irac'
;
; read in sample data
readcol, lfile, id, ra, dec, majdiam, mindiam, pa, type, format='a,d,d,f,f,f,a'
;
; define top level directory
deg = string(floor(ra), format='(i3.3)')+'D'
;
; loop over object list
for i=0,n_elements(id)-1 do begin
	print,i+1,'/',n_elements(id),': ',id[i]
;
; directories
	auxpath=!GLGA_ROOT+'/data/'+deg[i]+'/aux/'
	dsspath=!GLGA_ROOT+'data/'+deg[i]+'/dss/fits/'
	photpath=!GLGA_ROOT+'data/'+deg[i]+'/photometry/'
	plotpath=!GLGA_ROOT+'data/'+deg[i]+'/plots/'
	fpath=!GLGA_ROOT+'data/'+deg[i]+'/'+srvy+'/fits/'
	jpath=!GLGA_ROOT+'data/'+deg[i]+'/'+srvy+'/jpg/'
;
; get mask, dss image file names
	maskimgfile=auxpath+id[i]+'_'+srvy+'_mask.fits.gz'
	dssfile=dsspath+id[i]+'_dss2_red.fits*'
	;
	; create plots
	if strpos(srvy,'galex') ge 0 then begin
		galex_plotradprof, id[i], type=type[i], $
			pathtoprofile=photpath, $
    			fintfile=fpath+id[i]+'_FUV.fit*', $
			nintfile=fpath+id[i]+'_NUV.fit*', $
    			maskimgfile=maskimgfile,$
    			uvjpgpath=jpath, dssfile=dssfile, $
			outpath=plotpath, verbose=verbose, /yuan13
	endif else if strpos(srvy,'sdss') ge 0 then begin
		sdss_plotradprof, id[i], type=type[i], $
			pathtoprofile=photpath, $
			intfile=fpath+id[i]+'_r.fit*', $
			maskimgfile=maskimgfile, $
			jpgpath=jpath, outpath=plotpath, $
			verbose=verbose, /yuan13
	endif else if strpos(srvy,'2mass') ge 0 then begin
		twomass_plotradprof, id[i], type=type[i], $
			pathtoprofile=photpath, $
			intfile=fpath+id[i]+'_h.fit*', $
			maskimgfile=maskimgfile, $
			jpgpath=jpath, outpath=plotpath, $
			verbose=verbose, /yuan13
	endif else if strpos(srvy,'wise') ge 0 then begin
		wise_plotradprof, id[i], type=type[i], $
			pathtoprofile=photpath, $
			intfile=fpath+id[i]+'_w1.fit*', $
			maskimgfile=maskimgfile, $
			jpgpath=jpath, outpath=plotpath, $
			verbose=verbose, /yuan13
	endif else if strpos(srvy,'irac') ge 0 then begin
		irac_plotradprof, id[i], type=type[i], $
			pathtoprofile=photpath, $
			intfile=fpath+id[i]+'_4p5um.fit*', $
			maskimgfile=maskimgfile, dssfile=dssfile, $
			jpgpath=jpath, outpath=plotpath, $
			verbose=verbose, /yuan13
	endif

endfor

return
end
