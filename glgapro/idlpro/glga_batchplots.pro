pro glga_batchplots, lfile, ra=ra, galex=galex,sdss=sdss,twomass=twomass,$
	wise=wise,irac=irac, $
	outpath=outpath, verbose=verbose, update=update, rotatee=rotatee, fromindex=fromindex
;+
; glga_batchplots - plot per band radial profile photometry and images
;	of objects in lfile
;
; list of objects with these columns:
;       id	- string
;       ra,dec  - degrees
;       majdiam,mindiam - arcmin
;       pa      - degrees
;	type	- string
;       
;- Set the keyword 'ra' together with the name of the galaxy for a single galaxy
; Example: IDL> glga_batchplots, 'UGC01329', ra=27., /wise, /verbose, /update
; 
; set defaults and check keywords
bands = ['FUV', 'NUV']
srvy='galex'		; data type
if keyword_set(sdss) then begin
	bands=['u','g','r','i','z']
	srvy='sdss'
endif
if keyword_set(twomass) then begin
	bands = ['j','h','k']
	srvy='2mass'
endif
if keyword_set(wise) then begin
	bands = ['w1','w2','w3','w4']
	srvy='wise'
endif
if keyword_set(irac) then begin
	bands = ['3p6um','4p5um']
	srvy='irac'
endif
nbands = n_elements(bands)


; check input file
if not file_test(lfile) then begin
	print,'GLGA_QA: Error - file not found: ',lfile
	
	if keyword_set(ra) then begin
	   
	   id = [lfile]
	   ra = [ra]
	   dec = [0]
	   majdiam = [1.]
	   mindiam = [1.]
	   pa = [45.]
	   type = ['???']
	
	endif else return
endif


;
; read in sample data
if not keyword_set(ra) then begin
readcol,lfile, id, ra, dec, majdiam, mindiam, pa, type, $
	format='a,d,d,f,f,f,a', /silent
endif
;
; define top level directory
deg = string(floor(ra), format='(i3.3)')+'D'
;
; loop over object list

if not keyword_set(fromindex) then fromindex = 0


for i=fromindex,n_elements(id)-1 do begin
	print,i+1,'/',n_elements(id),': ',id[i]
;
; directories
	auxpath=!GLGA_ROOT+'/data/'+deg[i]+'/aux/'
	dsspath=!GLGA_ROOT+'data/'+deg[i]+'/dss/fits/'
	photpath=!GLGA_ROOT+'data/'+deg[i]+'/photometry/'
	if keyword_set(outpath) then $
		plotpath=outpath $
	else	plotpath=!GLGA_ROOT+'data/'+deg[i]+'/plots/'
	fpath=!GLGA_ROOT+'data/'+deg[i]+'/'+srvy+'/fits/'
	jpath=!GLGA_ROOT+'data/'+deg[i]+'/'+srvy+'/jpg/'
;
; get mask image file name
	maskimgfile=auxpath+id[i]+'_'+srvy+'_mask.fits.gz'
	;
	; loop over bands
	for j=0,nbands-1 do begin
	    if file_test(fpath+id[i]+'_'+bands[j]+'.fit*') eq 1 then begin
	    
	        bandmaskimgfile=auxpath+id[i]+'_'+srvy+'_mask_'+bands[j]+'.fits.gz'
	        
		glga_plotradprof, id[i], bands[j], type=type[i], $
			survey=strupcase(srvy), pathtoprofile=photpath, $
			intfile=fpath+id[i]+'_'+bands[j]+'.fit*', $
			maskimgfile=maskimgfile, bandmaskimgfile=bandmaskimgfile, outpath=plotpath, $
			verbose=verbose, /yuan13, update=update, rotatee=rotatee
	    endif
	endfor

endfor

return
end
