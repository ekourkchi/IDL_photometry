pro qa_stat__define
;+
; qa_stat__define - define qa_stat structure
;
; TAGS:
;	TS 			- QA timestamp like systime(1)
;	TS_PRIMARY_IMG		- primary image timestamp (0.0 = missing)
;	TS_SECONDARY_IMG	- secondary image timestamp (0.0 = missing)
;	TS_JPG_IMG		- jpeg image timestamp (0.0 = missing)
;	TS_PHOT_PLOTS		- photometry plots timestamp (0.0 = missing)
;	TS_IMAGE_PLOTS		- image plots timestamp (0.0 = missing)
;	TS_ALTIM_PLOTS		- alt wave image plots timestamp (0.0 = missing)
;	ALTIM_TYPE		- (str) alt wave image type (GALEX, SDSS, etc.)
;	UNAME 			- (str) username of last person to do QA
;	NOTE			- (str) QA note
;	FILE 			- (str) full path and filename of qa file
;	NQA   			- (int) # of times QA has been run on this host
;	--FLAGS (integer) below: 1 = yes, 0 = no
;	COMPLETE 		- QA complete?
;	REQUIRE_JPG		- did QA require jpg image?
;	ERROR 			- error preventing QA?
;	UNCERTAIN_MASK 		- uncertain about masking?
;	BRIGHT_STAR 		- host involved with bright star?
;	MULTIPLE 		- host multiple?
;	FOV_EXPAND 		- FOV needs expansion?
;	ELLIPSE_UPDATE 		- updated ellipse info?
;	ROI			- ROI used?
;	PSRC			- point sources marked?
;	MASK			- mask file generated?
;	BAND1_EDGE		- galaxy on edge of band1 (FUV, g, j)?
;	BAND1_SN_GRAD		- galaxy on S/N gradient in band1 (FUV, g, j)?
;	BAND1_ARTIFACT		- galaxy on artifact in band1 (FUV, g, j)
;	BAND1_OTHER		- band1 (FUV, g, j) other error?
;	BAND1_MISSING		- band1 (FUV, g, j) data missing?
;	BAND2_EDGE		- galaxy on edge of band2 (NUV, r, h)?
;	BAND2_SN_GRAD		- galaxy on S/N gradient in band2 (NUV, r, h)?
;	BAND2_ARTIFACT		- galaxy on artifact in band2 (NUV, r, h)
;	BAND2_OTHER		- band2 (NUV, r, h) other error?
;	BAND2_MISSING		- band2 (NUV, r, h) data missing?
;	BAND3_EDGE		- galaxy on edge of band3 (i, k)?
;	BAND3_SN_GRAD		- galaxy on S/N gradient in band3 (i, k)?
;	BAND3_ARTIFACT		- galaxy on artifact in band3 (i, k)
;	BAND3_OTHER		- band3 (i, k) other error?
;	BAND3_MISSING		- band3 (i, k) data missing?
;
;-
;	
tmp = { qa_stat, $
	ts:0LL, $		; last QA timestamp (unix seconds systime(1))
	ts_primary_img:0LL, $	; primary image time stamp (0 = missing)
	ts_secondary_img:0LL, $	; secondary image time stamp (0 = missing)
	ts_jpg_img:0LL, $	; jpeg image time stamp (0 = missing)
	ts_phot_plots:0LL, $	; photometry plots time stamp (0 = missing)
	ts_image_plots:0LL, $	; image plots time stamp (0 = missing)
	ts_altim_plots:0LL, $	; alt wave image time stamp (0 = missing)
	altim_type:'', $	; alt wave image type (GALEX, SDSS, 2MASS, etc.)
	user_name:'', $		; username of last person to perform QA
	machine_name:'', $	; machine name on which QA was last done
	note:'', $		; QA note
	file:'', $		; filename of QA file: <host>_qa.txt
	nqa:0, $		; count of QA runs
	complete:0, $		; QA complete?
	require_jpg:0, $	; did QA require the jpg image?
	error:0, $		; error preventing QA?
	uncertain_mask:0, $	; unresolved mask question?
	bright_star:0, $	; galaxy involved with bright star?
	multiple:0, $		; galaxy mutiple?
	fov_expand:0, $		; FOV needs expanding?
	ellipse_update:0, $	; updated ellipse info?
	roi:0, $		; ROI used?
	psrc:0, $		; point sources masked?
	mask:0, $		; mask file generated?
	band1_edge:0, $		; galaxy involved in band1 (FUV, g, j) edge?
	band1_sn_grad:0, $	; galaxy involved in band1 (FUV, g, j) S/N grad?
	band1_artifact:0, $	; galaxy involved in band1 (FUV, g, j) artifact?
	band1_other:0, $	; band1 (FUV, g, j) other error?
	band1_missing:0, $	; band1 (FUV, g, j) data missing?
	band2_edge:0, $		; galaxy involved in band2 (NUV, r, h) edge?
	band2_sn_grad:0, $	; galaxy involved in band2 (NUV, r, h) S/N grad?
	band2_artifact:0, $	; galaxy involved in band2 (NUV, r, h) artifact?
	band2_other:0, $	; band2 (NUV, r, h) other error?
	band2_missing:0, $	; band2 (NUV, r, h) data missing?
	band3_edge:0, $		; galaxy involved in band3 (i, k) edge?
	band3_sn_grad:0, $	; galaxy involved in band3 (i, k) S/N grad?
	band3_artifact:0, $	; galaxy involved in band3 (i, k) artifact?
	band3_other:0, $	; band3 (i, k) other error?
	band3_missing:0, $	; band3 (i, k) data missing?
	quality:5, $            ; Quality of photometry (0-5)
	disturbed:0, $          ; Disturbed and/or distorted (confused)           
	trail:0, $              ; Has tail or trail 
	not_spiral:0, $         ; NOT a Spiral Galaxy
	face_on:0, $            ; Face-on Galaxy
	faint:0, $              ; Too Faint
	crowded:0, $            ; Crowded Field
	over_masked:0, $        ; Needed a lot of Masking
	backradius:0.0D $       ; Where to modify background level (semimajor axis in arcmin), (0 = no modifications)
	
	}
end
