pro glga_write_qa_stat,qa,batch=batch, hasNote=hasNote
;+
; glga_write_qa_stat - write out the qa_stat structure in the file qa.file
;
; INPUTS:
;	qa    - qa status struct see qa_stat__define.pro
;
; KEYWORDS:
;	batch - do not increment QA count, ask for note, or change login info
;
;-
cmnts = [ $
	'# last QA timestamp (unix seconds systime(1))', $
	'# primary image time stamp (0 = missing)', $
	'# secondary image time stamp (0 = missing)', $
	'# jpeg image time stamp (0 = missing)', $
	'# phot plots time stamp (0 = missing)', $
	'# image plots time stamp (0 = missing)', $
	'# alt wave image time stamp (0 = missing)', $
	'# alt wave image type (GALEX, SDSS, 2MASS, etc.)', $
	'# username of last person to perform QA', $
	'# machine name on which QA was last done', $
	'', $	; QA note
	'# QA file full path', $
	'# count of QA runs', $
	'# QA completed?', $
	'# did QA require the jpg image?', $
	'# error preventing QA?', $
	'# unresolved QA mask question?', $
	'# galaxy involved with bright star?', $
	'# multiple galaxy?', $
	'# FOV needs expanding?', $
	'# updated ellipse info?', $
	'# ROI used?', $
	'# point sources masked?', $
	'# mask file generated?', $
	'# galaxy involved in band1 (FUV, g, j) image edge?', $
	'# galaxy involved in band1 (FUV, g, j) S/N gradient?', $
	'# galaxy involved in band1 (FUV, g, j) artifact?', $
	'# band1 (FUV, g, j) other error?', $
	'# band1 (FUV, g, j) data missing?', $
	'# galaxy involved in band2 (NUV, r, h) image edge?', $
	'# galaxy involved in band2 (NUV, r, h) S/N gradient?', $
	'# galaxy involved in band2 (NUV, r, h) artifact?', $
	'# band3 (NUV, r, h) other error?', $
	'# band2 (NUV, r, h) data missing?', $
	'# galaxy involved in band3 (i, k) image edge?', $
	'# galaxy involved in band3 (i, k) S/N gradient?', $
	'# galaxy involved in band3 (i, k) artifact?', $
	'# band3 (i, k) other error?', $
	'# band3 (i, k) data missing?', $
	'# Quality of photometry (0-5)', $
	'# Disturbed and/or distorted (confused)?', $
	'# Has tail or trail?', $
	'# NOT a Spiral Galaxy?', $
	'# Face-on Galaxy?', $
	'# Too Faint?', $
	'# Crowded Field?', $
	'# Needed a lot of Masking?', $
	'# Where to modify background level (semimajor axis in arcmin), (0 = no modifications)' $
	]
;
; in batch mode do not change these items
if not keyword_set(batch) then begin
	qa.nqa = qa.nqa + 1
;
; add login info
	linfo = get_login_info()
	qa.user_name = linfo.user_name
	qa.machine_name = linfo.machine_name
;
; get comment
    if  not keyword_set(hasNote) then begin
	cm=''
	read,'Enter QA note: ',cm
	if strlen(cm) gt 0 then $
		qa.note = cm
    endif
    
endif


; get qa struct tags
tags = tag_names(qa)
maxl = max(strlen(tags))
ntags= n_elements(tags)
;
; open qa file
openw,ol,qa.file,/get_lun
;
; loop over tags and write out
for i=0,ntags-1 do begin
	ty=size(qa.(i),/type)
	if strcmp(tags[i],'NOTE') eq 1 then begin
		ofmt = '(a-'+strn(maxl)+',a1,a,a)'
	endif else if strcmp(tags[i],'FILE') eq 1 or $
		      strcmp(tags[i],'MACHINE_NAME') eq 1 then begin
		ofmt = '(a-'+strn(maxl)+',a1,a,2x,a)'
	endif else if ty eq 7 then begin
		ofmt = '(a-'+strn(maxl)+',a1,a19,2x,a)'
	
	endif else if strcmp(tags[i],'BACKRADIUS') eq 1 then begin
                ofmt = '(a-'+strn(maxl)+',a1,F19.2,2x,a)'
	
	endif else ofmt = '(a-'+strn(maxl)+',a1,i19,2x,a)'
	printf,ol,tags[i],'=',qa.(i),cmnts[i],format=ofmt
endfor
;
free_lun,ol
;
return
end
