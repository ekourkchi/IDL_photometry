pro glga_deliver,lfile,photometry=photometry,plots=plots,aux=aux,jpg=jpg,qa=qa,$
	filter_plots=filter_plots, $
	images=images,galex=galex,sdss=sdss,twomass=twomass,wise=wise,dss=dss
;+
; glga_deliver - gather GLGA products in one place
;-
pdirs = ['photometry/','plots/','aux/','plots/']
ddirs = ['2mass/','sdss/','galex/','wise/']
ptail=[['_[jhk]_*.dat','_2mass_*.*','_{ellipse,2mass_*}.dat','_[jhk].*'], $
       ['_[ugriz]_*.dat','_sdss_*.*','_{ellipse,sdss_*}.dat','_[ugriz].*'], $
       ['_?UV_*.dat','_galex_*.*','_{ellipse,galex_*}.dat','_?UV.*'],$
       ['_w?_*.dat','_wise_*.*','_{ellipse,wise_*}.dat','_w?.*']]
dtail=['_[jhk].fit*','_[ugriz].fit*','_?UV*.fit*','_w?*.fit*']
;
pdo = [0,0,0,0]
ddo = [0,0,0,0]
;
; check keywords
if keyword_set(photometry) then pdo[0] = 1
if keyword_set(plots) then pdo[1] = 1
if keyword_set(aux) then pdo[2] = 1
if keyword_set(filter_plots) then pdo[3] = 1
if keyword_set(twomass) then ddo[0] = 1
if keyword_set(sdss) then ddo[1] = 1
if keyword_set(galex) then ddo[2] = 1
if keyword_set(wise) then ddo[3] = 1
;
; read list
readcol,lfile,id,ra,form='a,d',comment='#',/silent
ngal = n_elements(id)
;
; get output
fdecomp,lfile,disk,dir,rute
ofil = rute+'_cp.csh'
openw,ol,ofil,/get_lun
;
; loop over gals
for k = 0,ngal-1 do begin
	bdir = !GLGA_ROOT + 'data/' + glga_degdir(ra[k]) + '/'
	;
	; loop over data type
	for i = 0,3 do begin
    	; is this data type selected?
	    if ddo[i] ne 0 then begin
		;
		; loop over product
		for j = 0,3 do begin
	    	; if this product selected?
	    	    if pdo[j] ne 0 then begin
			fs = bdir+pdirs[j]+id[k]+ptail[j,i]
			fl = file_search(fs, count=nf)
			print,'Found: ',nf,fs,form='(a,i5,2x,a)'
			if nf gt 0 then for m=0,nf-1 do $
				printf,ol,'cp '+fl[m]+' .'
			if nf le 0 then $
				print,'Nothing found for: ',id[k]
		    endif
		endfor	; loop over product
		;
		; check image data
		if keyword_set(images) then begin
			fs = bdir+ddirs[i]+'fits/'+id[k]+dtail[i]
			fl = file_search(fs, count=nf)
			print,'Found: ',nf,fs,form='(a,i5,2x,a)'
			if nf gt 0 then for m=0,nf-1 do $
				printf,ol,'cp '+fl[m]+' .'
			if nf le 0 then $
				print,'No images found for: ',id[k]
		endif
		;
		; check qa files
		if keyword_set(qa) then begin
			fs = bdir+ddirs[i]+'fits/'+id[k]+'_qa.txt'
			fl = file_search(fs, count=nf)
			print,'Found: ',nf,fs,form='(a,i5,2x,a)'
			if nf gt 0 then for m=0,nf-1 do $
				printf,ol,'cp '+fl[m]+' .'
			if nf le 0 then $
				print,'No qa files found for: ',id[k]
		endif
		;
		; check jpg data
		if keyword_set(jpg) then begin
			fs = bdir+ddirs[i]+'jpg/'+id[k]+'_*.jpg'
			fl = file_search(fs, count=nf)
			print,'Found: ',nf,fs,form='(a,i5,2x,a)'
			if nf gt 0 then for m=0,nf-1 do $
				printf,ol,'cp '+fl[m]+' .'
			if nf le 0 then $
				print,'No jpgs found for: ',id[k]
		endif
		;
		; check dss data
		if keyword_set(dss) then begin
			fs = bdir+'dss/fits/'+id[k]+'_*.fit*'
			fl = file_search(fs, count=nf)
			print,'Found: ',nf,fs,form='(a,i5,2x,a)'
			if nf gt 0 then for m=0,nf-1 do $
				printf,ol,'cp '+fl[m]+' .'
			if nf le 0 then $
				print,'No dss images found for: ',id[k]
		endif
	    endif	; is this data type selected?
	endfor	; loop over data type
endfor	; loop over gals
;
free_lun,ol
return
end
