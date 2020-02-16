function extract_surveys, inlist, nsurveys
;+
; extract_surveys - extract the galaxy surveys from filenames
;-
nsurveys = n_elements(inlist)
srvs = ''
if nsurveys eq 1 then begin
	sta = strsplit(inlist,'/',/extract)
	if n_elements(sta) ge 6 then $
		srvs = sta[5] $
	else	nsurveys = 0
endif else if nsurveys gt 1 then begin
	srvs = strarr(nsurveys)
	for i=0l,nsurveys-1l do begin
		sta = strsplit(inlist[i],'/',/extract)
		if n_elements(sta) ge 6 then $
			srvs[i] = sta[5]
	endfor
	;
	; clean up
	good = where(strlen(srvs) gt 0, ngood)
	if ngood gt 0 then begin
		srvs = srvs[good]
		srvs = srvs[sort(srvs)]
		srvs = srvs[uniq(srvs)]
		nsurveys = n_elements(srvs)
	endif else begin
		srvs = ''
		nsurveys = 0
	endelse
endif
;
return,srvs
end
