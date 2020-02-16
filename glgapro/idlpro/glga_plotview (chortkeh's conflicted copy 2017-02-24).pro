pro glga_plotview,lfile,sdss=sdss,galex=galex,twomass=twomass,wise=wise, $
	start=start, pick_file=pick_file,skip=skip
;+
; glga_plotview - look at all the plots!
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
	print,'GLGA_PLOTVIEW: Error - file not found: ',lfile
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
			print,'GLGA_PLOTVIEW: Error - ambiguous id: ',start
			return
		endif else if nw lt 1 then begin
			print,'GLGA_PLOTVIEW: Error - not found: ',start
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
; set up pick file
if keyword_set(pick_file) then begin
	openw,pl,pick_file,/get_lun,/append
	printf,pl,'# GLGA_PLOTVIEW Selected Objects - '+systime(0)
	printf,pl,'# Selected from: '+lfile
	free_lun,pl
endif
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
print, photplot
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
	if keyword_set(skip) and not plotted then $
		print,'no plots, skipping...' $
	else	read,'what? (<cr> - next, s - select, q - quit): ',q
;
; clean up windows
	while !d.window ge 0 do wdelete
;
; check for stop VIEW-ing
	if strupcase(strmid(q,0,1)) eq 'Q' then goto, done
;
; shall we select?
	if strupcase(strmid(q,0,1)) eq 'S' then $
		pick_obj = (1 eq 1) $
	else	pick_obj = (1 eq 0)
;
; picked?
	if pick_obj and keyword_set(pick_file) then begin
		openw,pl,pick_file,/get_lun,/append
		printf,pl,id[i],ra[i],dec[i],majdiam[i],mindiam[i], $
			pa[i],type[i], format='(a-25,2f13.8,3f9.3,2x,a)'
		free_lun,pl
		print,'Wrote ',id[i],' to ',pick_file
	endif
;
; print status
	print,' '
	q = ''

endfor	; loop over object list
;
; go here if done VIEW-ing
done:
;
; clean up
while !d.window ge 0 do wdelete
print, 'Finis PLOTVIEW'
end
