function glga_get_qa_flags, qa
;+
; glga_get_qa_flags - get flags from user pertaining to QA
;
; INPUT:
;	qa	- qa_stat structure read in with glga_read_qa_stat.pro
;
; RETURNS:
;	a string summarizing the qa status
;-
;
status=''
q=''
read,'Enter <cr> if QA complete with no issues, else set flags: ',q
if q eq '' then begin
	qa.complete = 1
	qa.uncertain_mask = 0
	qa.fov_expand = 0
	qa.band1_edge = 0
	qa.band1_sn_grad = 0
	qa.band1_artifact = 0
	qa.band1_missing = 0
	qa.band2_edge = 0
	qa.band2_sn_grad = 0
	qa.band2_artifact = 0
	qa.band2_missing = 0
	qa.band3_edge = 0
	qa.band3_sn_grad = 0
	qa.band3_artifact = 0
	qa.band3_missing = 0
	return,'QA complete'
endif else begin
	print,'Enter all flags to set, missing ones are cleared'
	print,'    C - complete in spite of flags'
	print,'    M - multiple galaxy'
	print,'    B - bright star involved'
	print,'    Q - open question about masking'
	print,'    V - expand FOV'
	print,'    D - data issue: edge, S/N gradient, missing, artifact'
	read,'Enter upper/lower case in any order (e.g. MBV): ',q
	q=strupcase(q)
	if strpos(q,'C') ge 0 then begin
		qa.complete = 1
		status = 'QA complete'
	endif else qa.complete = 0
	if strpos(q,'M') ge 0 then begin
		qa.multiple = 1
		status = status + ' multiple'
	endif else qa.multiple = 0
	if strpos(q,'B') ge 0 then begin
		qa.bright_star = 1
		status = status + ' bright star'
	endif else qa.bright_star = 0
	if strpos(q,'Q') ge 0 then begin
		qa.uncertain_mask = 1
		status = status + ' mask question'
	endif else qa.uncertain_mask = 0
	if strpos(q,'V') ge 0 then begin
		qa.fov_expand = 1
		status = status + ' needs fov expansion'
	endif else qa.fov_expand = 0
	if strpos(q,'D') ge 0 then begin
		print,'Enter all flags to set, missing ones cleared'
		print,'    1 - BAND1 (FUV, g, j) edge'
		print,'    2 - BAND1 (FUV, g, j) S/N gradient'
		print,'    3 - BAND1 (FUV, g, j) artifact'
		print,'    4 - BAND1 (FUV, g, j) data missing'
		print,'    5 - BAND2 (NUV, r, h) edge'
		print,'    6 - BAND2 (NUV, r, h) S/N gradient'
		print,'    7 - BAND2 (NUV, r, h) artifact'
		print,'    8 - BAND2 (NUV, r, h) data missing'
		print,'    9 - BAND3 (i, k) edge'
		print,'    A - BAND3 (i, k) S/N gradient'
		print,'    B - BAND3 (i, k) artifact'
		print,'    C - BAND3 (i, k) data missing'
		read,'Enter in any order (e.g. 1524C): ',q
		if strpos(q,'1') ge 0 then begin
			qa.band1_edge = 1
			status = status + ' band1 edge'
		endif else qa.band1_edge = 0
		if strpos(q,'2') ge 0 then begin
			qa.band1_sn_grad = 1
			status = status + ' band1 S/N grad'
		endif else qa.band1_sn_grad = 0
		if strpos(q,'3') ge 0 then begin
			qa.band1_artifact = 1
			status = status + ' band1 artifact'
		endif else qa.band1_artifact = 0
		if strpos(q,'4') ge 0 then begin
			qa.band1_missing = 1
			status = status + ' band1 missing'
		endif else qa.band1_missing = 0
		if strpos(q,'5') ge 0 then begin
			qa.band2_edge = 1
			status = status + ' band2 edge'
		endif else qa.band2_edge = 0
		if strpos(q,'6') ge 0 then begin
			qa.band2_sn_grad = 1
			status = status + ' band2 S/N grad'
		endif else qa.band2_sn_grad = 0
		if strpos(q,'7') ge 0 then begin
			qa.band2_artifact = 1
			status = status + ' band2 artifact'
		endif else qa.band2_artifact = 0
		if strpos(q,'8') ge 0 then begin
			qa.band2_missing = 1
			status = status + ' band2 missing'
		endif else qa.band2_missing = 0
		if strpos(q,'9') ge 0 then begin
			qa.band3_edge = 1
			status = status + ' band3 edge'
		endif else qa.band3_edge = 0
		if strpos(q,'A') ge 0 then begin
			qa.band3_sn_grad = 1
			status = status + ' band3 S/N grad'
		endif else qa.band3_sn_grad = 0
		if strpos(q,'B') ge 0 then begin
			qa.band3_artifact = 1
			status = status + ' band3 artifact'
		endif else qa.band3_artifact = 0
		if strpos(q,'C') ge 0 then begin
			qa.band3_missing = 1
			status = status + ' band3 missing'
		endif else qa.band3_missing = 0
	endif else begin
		qa.band1_edge = 0
		qa.band1_sn_grad = 0
		qa.band1_artifact = 0
		qa.band1_missing = 0
		qa.band2_edge = 0
		qa.band2_sn_grad = 0
		qa.band2_artifact = 0
		qa.band2_missing = 0
		qa.band3_edge = 0
		qa.band3_sn_grad = 0
		qa.band3_artifact = 0
		qa.band3_missing = 0
	endelse
endelse
;
return,status
;
end
