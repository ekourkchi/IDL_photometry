pro glga_plotgrade,lfile,sdss=sdss,galex=galex,twomass=twomass,wise=wise, $
	start=start
;+
; glga_plotgrade - grade the plots: 1 - OK, 2 - Questionable, 3 - junk!
;
; lfile - list of objects with one line per object with these columns:
;	id
;	ra,dec	- degrees
;	majdiam,mindiam	- arcmin
;	pa	- degrees
;
; keywords:
;	sdss,galex,twomass,wise, - plots to see
;
;-
; set defaults and check keywords
;
filts = ['FUV','NUV','u','g','r','i','z','j','h','k','w1','w2','w3','w4']
if keyword_set(galex) then begin
	filts = ['FUV','NUV']
endif
if keyword_set(sdss) then begin
	filts = ['u','g','r','i','z']
endif
if keyword_set(twomass) then begin
	filts = ['j','h','k']
endif
if keyword_set(wise) then begin
	filts = ['w1','w2','w3','w4']
endif
nfilts = n_elements(filts)
;
; check input file
if not file_test(lfile) then begin
	print,'GLGA_PLOTGRADE: Error - file not found: ',lfile
	return
endif
;
; read in sample data
readcol,lfile, id, ra, dec, majdiam, mindiam, pa, type, $
	format='a,d,d,f,f,f,a', /silent, comment='#'
;
; define top level directory
deg = string(floor(ra), format='(i3.3)')+'D'
;
; used to generate file names
id = strcompress(id,/rem)

filebase=!GLGA_ROOT+'data/'+deg+'/plots/'+id
;
; query variable
q=''
;
; starting index
i0=0L
if keyword_set(start) then begin
	if size(start,/type) eq 7 and strtrim(start,2) ne '' then begin
		w=where(strcmp(id,strtrim(start,2)) eq 1, nw)
		if nw gt 1 then begin
			print,'GLGA_PLOTGRADE: Error - ambiguous id: ',start
			return
		endif else if nw lt 1 then begin
			print,'GLGA_PLOTGRADE: Error - not found: ',start
			return
		endif else i0 = w[0]
	endif else if (size(start,/type) ge 12 and size(start,/type) le 15) or $
		(size(start,/type) ge 1 and size(start,/type) le 3) then $
		i0 = start - 1L
endif
nloop = n_elements(id)
;
; get display size
device,get_screen_size=screen_size
;
; loop over object list
for i=i0, nloop-1 do begin
;
; print host
	print,i+1,'/',nloop,id[i],deg[i],' start browse...', $
		format = '(i6,a1,i6,2x,a-25,a5,a)'
;
; reset window index
	wi = 0
;
; reset plotted boolean
	plotted = ( 1 eq 0 )
;
; loop over filters
	for f = 0, nfilts-1 do begin
		photplot = filebase[i] + '_' + filts[f] + '.jpg'
;
; radial plots
		finfo=file_info(photplot)
		if finfo.exists then begin
			read_jpeg, photplot, jplot, /true
			sz = size(jplot)
			dsize = min([sz[2]*nfilts, screen_size[0]]) + 2
			xp = wi * ((dsize-sz[2])/(nfilts-1))
			window, wi, xsize = sz[2], ysize = sz[3], $
				xpos = xp, ypos = screen_size[1], $
		    		title=id[i] + ' '+filts[f]+ $
				' Photometry ('+strn(i+1)+'/'+strn(nloop)+')'
			tvscl, jplot, /true
			plotted = ( 1 eq 1 )
		endif
		wi += 1
	endfor	; loop over filts
;
; read input
	if not plotted then begin
		print,'no plots, skipping...'
;
; we got some plots to look at
	endif else begin
		read,'what? (<cr> - next, g - grade, q - quit): ',q
;
; shall we grade?
		if strupcase(strmid(q,0,1)) eq 'G' then begin
;
			print,'enter filter=<grade>, <cr> to quit'
			q = ''
			read,': ',q
			while strlen(q) gt 0 do begin
				fpick = gettok(q,'=')	; filter
				gr = fix(q)		; grade
				pick = where(strpos(filts,fpick) ge 0, np)
				if np eq 1 and gr ge 1 and gr le 3 then begin
					glist = file_search(filebase[i] + '_' $
						    + fpick + '_grade?', $
						    count=ngf)
				; delete old grade
					if ngf gt 0 then $
						for j=0,ngf-1 do $
							file_delete,glist[j]
				; create new grade (default grade is 1)
					if gr gt 1 then begin
					  tcmd = 'touch '+filebase[i] + '_' + $
						fpick + '_grade' + strn(gr)
					  spawn,tcmd
					endif
				endif else $
					print,'bad filter/grade: ',fpick,'/',gr
				read,': ',q
			endwhile	; still grading by filter
;
		endif	; grading
;
; clean up windows
		while !d.window ge 0 do wdelete
;
; check for stop GRAD-ing
		if strupcase(strmid(q,0,1)) eq 'Q' then goto, done
	endelse	; we got some plots to look at
;
; print status
	print,' '
	q = ''

endfor	; loop over object list
;
; go here if done GRAD-ing
done:
;
; clean up
while !d.window ge 0 do wdelete
print, 'Finis PLOTGRADE'
end
