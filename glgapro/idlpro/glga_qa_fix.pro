pro glga_qa_fix,fspec,sdss=sdss,galex=galex,twomass=twomass,wise=wise
;+
; glga_qa_fix - read in old style qa files and write out new ones
;
; fspec - full file spec, e.g. !GLGA_ROOT+'data/*D/*/fits/*_qa.txt'
;
; keywords:
;	sdss,galex,twomass,csp - data to analyze
;
;-
;
flist=file_search(fspec,count=nf)
print,'Examining ',nf,' files.',form='(a,i5,a)'
;
for i=0,nf-1 do begin
	print,string(13B),i+1,'/',nf,flist[i],form='($,a1,i5,a1,i5,2x,a-64)'
	qa = glga_read_qa_stat(flist[i],stat=stat)
	if stat gt 2 then begin	; not current format
;
; get head of QA file spec
	    head=strmid(flist[i],0,strpos(flist[i],'qa.txt')-1)
;
; get type
	    if strpos(head,'/galex/') ge 0 then begin
		pfil=head+'_NUV.fit*'
		afil=head+'_FUV.fit*'
		jfil=repstr(head,'/fits/','/jpg/')+'_FUVNUV.jpg'
		phed=repstr(head,'galex/fits/','plots/')
		phfil=phed+'_galex_profile.jpg'
		imfil=phed+'_galex_images.jpg'
		ai1fil=phed+'_sdss_images.jpg'
		ai2fil=phed+'_wise_images.jpg'
		reqjfile=1
	    endif else if strpos(head,'/sdss/') ge 0 then begin
		pfil=head+'_r.fit*'
		afil=head+'_g.fits'
		jfil=repstr(head,'/fits/','/jpg/')+'_gri.jpg'
		phed=repstr(head,'sdss/fits/','plots/')
		phfil=phed+'_sdss_profile.jpg'
		imfil=phed+'_sdss_images.jpg'
		ai1fil=phed+'_galex_images.jpg'
		ai2fil=phed+'_wise_images.jpg'
		reqjfile=1
	    endif else if strpos(head,'/2mass/') ge 0 then begin
		pfil=head+'_k.fit*'
		afil=head+'_j.fit*'
		jfil=repstr(head,'/fits/','/jpg/')+'_jhk.jpg'
		phed=repstr(head,'2mass/fits/','plots/')
		phfil=phed+'_2mass_profile.jpg'
		imfil=phed+'_2mass_images.jpg'
		ai1fil=phed+'_galex_images.jpg'
		ai2fil=phed+'_wise_images.jpg'
		reqjfile=1
	    endif else if strpos(head,'/wise/') ge 0 then begin
		pfil=head+'_w1.fit*'
		afil=head+'_w2.fit*'
		jfil=repstr(head,'/fits/','/jpg/')+'_w123.jpg'
		phed=repstr(head,'wise/fits/','plots/')
		phfil=phed+'_wise_profile.jpg'
		imfil=phed+'_wise_images.jpg'
		ai1fil=phed+'_galex_images.jpg'
		ai2fil=phed+'_sdss_images.jpg'
		reqjfile=1
	    endif else begin
		pfil=''
		print,'non-sensical QA file: ',flist[i], '  ', fspec
	    endelse
;
; make sure we got good type
	    if strlen(pfil) gt 0 then begin
;
; QA note
		qa.note='GLGA_QA_FIX run on '+systime(0)
;
; require jpg?
		qa.require_jpg = reqjfile
;
; get time stamps for input files
;
; primary and secondary images
		imlis=file_search(pfil, count=nim)
		if nim gt 0 then begin
			finfo=file_info(imlis[0])
			qa.ts_primary_img = finfo.mtime
		endif
		imlis=file_search(afil, count=nim)
		if nim gt 0 then begin
			finfo=file_info(imlis[0])
			qa.ts_secondary_img = finfo.mtime
		endif
;
; jpg image
		finfo=file_info(jfil)
		qa.ts_jpg_img = finfo.mtime
;
; phot and image plots
		finfo=file_info(phfil)
		qa.ts_phot_plots = finfo.mtime
		finfo=file_info(imfil)
		qa.ts_image_plots = finfo.mtime
;
; alt wave image
		finfo=file_info(ai1fil)
		if not finfo.exists then $
			finfo=file_info(ai2fil)
		qa.ts_altim_plots = finfo.mtime
		if finfo.exists then begin
			if strpos(finfo.name,'/galex/') ge 0 then $
				qa.altim_type = 'GALEX' $
			else	if strpos(finfo.name,'/sdss/') then $
				qa.altim_type = 'SDSS' $
			else	if strpos(finfo.name,'/2mass/') then $
				qa.altim_type = '2MASS'
		endif

;
; write it out and give it the right time stamp
		glga_write_qa_stat,qa,/batch
		dar = bin_date(systime(0,qa.ts))
		ts = string(dar[0],form='(i4)')
		for j=1,4 do ts = ts + string(dar[j],form='(i02)')
		ts = ts + '.' + string(dar[5],form='(i02)') + ' '
		cmd = 'touch -am -t '+ts+qa.file
		spawn,cmd
	    endif	; got good type
	endif	; stat gt 2
endfor
print,' '
print,'Done.'
;
return
end
