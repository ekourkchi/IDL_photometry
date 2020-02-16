function extract_ids, inlist, nids, degdir=degdir
;+
; extract_ids - extract the galaxy ids from glga full filenames
;-
nids = n_elements(inlist)
id = ''
degdir = ''
if nids eq 1 then begin
	fdecomp,inlist,disk,dir,rute
	id = gettok(rute,'_')
	degdir = strmid(dir,0,strpos(dir,'D/')+2)
endif else if nids gt 1 then begin
	id = strarr(nids)
	degdir = strarr(nids)
	for i=0l,nids-1l do begin
		fdecomp,inlist[i],disk,dir,rute
		idi = gettok(rute,'_')
		id[i] = idi
		degdir[i] = strmid(dir,0,strpos(dir,'D/')+2)
	endfor
	;
	; return a unique list of hosts
	s = sort(id)
	id = id[s]
	degdir = degdir[s]
	u = uniq(id)
	id = id[u]
	degdir = degdir[u]
	nids = n_elements(id)
endif
;
return,id
end
