function glga_read_qa_stat,qfil,stat=stat,verbose=verbose
;+
; glga_read_qa_stat - return current qa status from qfil
;
; INPUTS:
;	qfil - qa file with name: <host>_qa.txt
;
; RETURNS:
;	qa structure with the tags defined in qa_stat__define.pro
;
; KEYWORDS:
;	stat -	0 if no file exists
;		1 if file exists, format current, and qa finished
;		2 if file exists, format current, and qa unfinished
;		3 if file exists, format old
;		4 if file exists, format unreadable
;-
; default qa structure
qa={qa_stat}
qa.file=qfil
stat=0
;
; get file time stamp
finfo=file_info(qfil)
if not finfo.exists then begin
        qa.quality = 5
	return,qa	; blank structure
endif else begin
    qa.ts=finfo.mtime
;
; count keyword value pairs
    rec=''
    neq=0
    openr,il,qfil,/get_lun
    while not eof(il) do begin
	readf,il,rec
	if strpos(rec,'=') ge 0 then neq = neq + 1
    endwhile
    free_lun,il
;
    if neq gt 0 then begin
;
; new format for QA files
	tags = tag_names(qa)
	openr,il,qfil,/get_lun
	while not eof(il) do begin
		readf,il,rec
		if strpos(rec,'=') ge 0 then begin
			key=strupcase(strtrim(gettok(rec,'='),2))
			val=strtrim(gettok(rec,'#'),2)
			it =where(strcmp(tags,key) eq 1,nit)
			if nit le 0 then begin	; no keyword match
				if keyword_set(verbose) then $
					print,'Unknown keyword: ',key
			endif else if nit gt 1 then begin ; multiple match
				if keyword_set(verbose) then $
					print,'Ambiguous keyword: ',key
			endif else begin	; a single match
			    if strcmp(key,'TS') ne 1 and  $
			       strcmp(key,'FILE') ne 1 then begin
				ty=size(qa.(it),/type)
				if ty gt 0 then $
					qa.(it)=fix(val,type=ty) $
				else	if keyword_set(verbose) then $
				print,'Type mismatch for keyword,value: ', $
					key,' ',val
			    endif
			endelse
		endif	; rec contains '='
	endwhile	; not eof(il)
	free_lun,il
	if qa.complete then $
		stat = 1 $
	else	stat = 2
;
; oldest format with 29 cols
    endif else begin
;
	readcol,qfil,jnk,cols,form='(a,i)',/silent
	if max(cols) lt 32 then begin
		readcols,qfil,nqa,comp,mq,bs,mul,fov,err,pim,sim,jim,ppm,ipm, $
			ell,roi,psrc,mask,edge1,sg1,art1,dm1,edge2,sg2,art2, $
			dm2,edge3,sg3,art3,dm3,file, $
	form='i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,a',/silent
;
; check if it woiked
		if n_elements(nqa) eq 1 then begin
			qa.nqa=nqa
			qa.complete=comp
			qa.uncertain_mask=mq
			qa.bright_star=bs
			qa.multiple=mul
			qa.fov_expand=fov
			qa.error=err
			;qa.no_primary_img=pim
			;qa.no_secondary_img=sim
			;qa.no_jpg_img=jim
			;qa.no_phot_plots=ppm
			;qa.no_image_plots=ipm
			qa.ellipse_update=ell
			qa.roi=roi
			qa.psrc=psrc
			qa.mask=mask
			qa.band1_edge=edge1
			qa.band1_sn_grad=sg1
			qa.band1_artifact=art1
			qa.band1_missing=dm1
			qa.band2_edge=edge2
			qa.band2_sn_grad=sg2
			qa.band2_artifact=art2
			qa.band2_missing=dm2
			qa.band3_edge=edge3
			qa.band3_sn_grad=sg3
			qa.band3_artifact=art3
			qa.band3_missing=dm3
			stat=3
		endif else begin
			print,'GLGA_READ_QA_STAT - Unreadable file: ',qfil
			stat=4
		endelse
;
; older format with 32 columns
	endif else begin
		readcols,qfil,nqa,comp,mq,bs,mul,fov,err,pim,sim,jim,ppm,ipm, $
			ell,roi,psrc,mask, $
			edge1,sg1,art1,dm1,ot1, $
			edge2,sg2,art2,dm2,ot2, $
			edge3,sg3,art3,dm3,ot3,file, $
  form='i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,a',/silent
		if n_elements(nqa) eq 1 then begin
			qa.nqa=nqa
			qa.complete=comp
			qa.uncertain_mask=mq
			qa.bright_star=bs
			qa.multiple=mul
			qa.fov_expand=fov
			qa.error=err
			;qa.no_primary_img=pim
			;qa.no_secondary_img=sim
			;qa.no_jpg_img=jim
			;qa.no_phot_plots=ppm
			;qa.no_image_plots=ipm
			qa.ellipse_update=ell
			qa.roi=roi
			qa.psrc=psrc
			qa.mask=mask
			qa.band1_edge=edge1
			qa.band1_sn_grad=sg1
			qa.band1_artifact=art1
			qa.band1_missing=dm1
			qa.band1_other=ot1
			qa.band2_edge=edge2
			qa.band2_sn_grad=sg2
			qa.band2_artifact=art2
			qa.band2_missing=dm2
			qa.band2_other=ot2
			qa.band3_edge=edge3
			qa.band3_sn_grad=sg3
			qa.band3_artifact=art3
			qa.band3_missing=dm3
			qa.band3_other=ot3
			stat=3
		endif else begin
			print,'GLGA_READ_QA_STAT - Unreadable file: ',qfil
			stat=4
		endelse
	endelse
    endelse
endelse
;
return,qa
end
