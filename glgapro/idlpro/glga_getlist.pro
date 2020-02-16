pro glga_getlist, photometry=photometry, images=images, jpgs=jpgs, qa=qa, $
	dss=dss, missing=missing, outdated=outdated, unfinished=unfinished, $
	galex=galex,sdss=sdss,twomass=twomass,irac=irac,wise=wise
;+
; glga_getlist - scan the data directory and create a list based on selection
;
; KEYWORDS:
;
;-
; set defaults and check keywords
bands = ['fuv','nuv']	; image filters
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
;
; log file
logfile='list.txt'
filestamp,logfile,/arch
openw,ll,logfile,/get_lun
printf,ll,'# GLGA_GETLIST: '+systime(0)
printf,ll,'# '+strupcase(srvy)+', '+strn(nbands)+' BANDS'
if keyword_set(photometry) then printf,ll,'# PHOT ',form='($,a)'
if keyword_set(images) then 	printf,ll,'# IMGS ',form='($,a)'
if keyword_set(jpgs) then 	printf,ll,'# JPGS ',form='($,a)'
if keyword_set(qa) then 	printf,ll,'# QA   ',form='($,a)'
if keyword_set(dss) then 	printf,ll,'# DSS  ',form='($,a)'
if keyword_set(missing) then 	printf,ll,'missing'
if keyword_set(outdated) then 	printf,ll,'outdated'
if keyword_set(unfinished) then printf,ll,'unfinished'
;
; loop over data directory
for d = 0,359 do begin
;
; define top level directory
	deg = string(d, format='(i3.3)')+'D/'
;
; directories
	auxpath=!GLGA_ROOT+'data/'+deg+'aux/'
	dsspath=!GLGA_ROOT+'data/'+deg+'dss/fits/'
	photpath=!GLGA_ROOT+'data/'+deg+'photometry/'
	plotpath=!GLGA_ROOT+'data/'+deg+'plots/'
	fpath=!GLGA_ROOT+'data/'+deg+srvy+'/fits/'
	jpath=!GLGA_ROOT+'data/'+deg+srvy+'/jpg/'
;
; photometry missing
	if keyword_set(photometry) and keyword_set(missing) then begin
		flist = file_search(fpath+'*.fit*', count=nf)
		if nf gt 0 then begin
			id = extract_ids(flist,nf)
			if nf gt 0 then begin
			    for i=0,nf-1 do begin
				plist = file_search(photpath+id[i] + $
						'_'+bands[0]+'_ellipsepar.dat',$
						count=np)
				if np le 0 then begin
					print,'ID: ',id[i]
					printf,ll,id[i]
				endif
			    endfor
		        endif
		endif
	endif
;
; jpegs missing
	if keyword_set(jpgs) and keyword_set(missing) then begin
		flist = file_search(fpath+'*.fit*', count=nf)
		if nf gt 0 then begin
			id = extract_ids(flist,nf)
			if nf gt 0 then begin
			    for i=0,nf-1 do begin
				jlist = file_search(jpath+id[i]+'_*.jpg', $
					count=nj)
				if nj le 0 then begin
					print,'ID: ',id[i]
					printf,ll,id[i]
				endif
			    endfor
		        endif
		endif
	endif
;
; get other files
	;maskimgfile=auxpath+id[i]+'_'+srvy+'_mask.fits.gz'
	;ellipsefile=auxpath+id[i]+'_ellipse.dat'
	;dssfile=dsspath+id[i]+'_dss2_red.fits*'
	;for j=0,nbands-1 do begin
	;endfor	; loop over bands

endfor

return
end
