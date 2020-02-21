pro glga_measure, id, ra, dec, majdiam, mindiam, pa, type, $
        galex=galex,sdss=sdss,panstarrs = panstarrs, twomass=twomass,irac=irac,wise=wise, acs=acs, $
	verbose=verbose, annuli_size=annuli_size, logfile=logfile, $
	qa_run=qa_run, background_radius=background_radius
;+
; glga_measure - run radial profile photometry on objects in list
;
;       id	- string
;       ra,dec  - degrees
;       majdiam,mindiam - arcmin
;       pa      - degrees
;	type	- string
;
;-
;
; log file
if keyword_set(logfile) then begin
	filestamp,logfile,/arch
	openw,ll,logfile,/get_lun
	printf,ll,'# GLGA_BATCHPROFILE: '+systime(0)
	printf,ll,'# List: '+lfile
	printf,ll,'# NGAL: '+strn(n_elements(ra))
endif
;
; set defaults and check keywords
bands = ['FUV','NUV']	; image filters
srvy='galex'		; data type
if keyword_set(sdss) then begin
	bands=['u','g','r','i','z']
	srvy='sdss'
endif
if keyword_set(panstarrs) then begin
	bands=['g','r','i','z']
	srvy='panstarrs'
	print, bands
	print, srvy
endif
if keyword_set(acs) then begin
	bands=['F606W','F814W']
	srvy='acs'
	print, bands
	print, srvy
endif
if keyword_set(twomass) then begin
	bands = ['j','h','k']
	srvy='2mass'
endif
if keyword_set(wise) then begin
	bands = ['w1','w2','w3','w4']
	srvy='wise'
endif
if keyword_set(irac) then begin
	bands = ['3p6um','4p5um']
	srvy='irac'
endif
nbands = n_elements(bands)
if keyword_set(logfile) then $
	printf,ll,'# '+strupcase(srvy)+' '+strn(nbands)
;
; deal with negative sizes
a = where(majdiam le 0.0, count)
if count gt 0 then majdiam[a] = 0.5
if count gt 0 then mindiam[a] = 0.5
a = where(mindiam le 0.0,  count)
if count gt 0 then mindiam[a] = majdiam[a]
;
; expand by 50%
majordiam=majdiam*1.5
minordiam=mindiam*1.5
;
; deal with negative PA
w=where(pa lt 0. and pa gt -360., nw)
if nw gt 0 then pa[w] = pa[w] + 180.
w=where(pa lt 0. and pa gt -180., nw)
if nw gt 0 then pa[w] = pa[w] + 180.
;
; define top level directory
deg = string(floor(ra), format='(i3.3)')+'D'
;
; loop over object list
for i=0,n_elements(id)-1 do begin
	if keyword_set(qa_run) then $
		i0 = qa_run $
	else	i0 = i+1
	print,i0,'/',n_elements(id),': ',id[i]
	stat=0
;
; directories
	auxpath=!GLGA_ROOT+'data/'+deg[i]+'/aux/'
	dsspath=!GLGA_ROOT+'data/'+deg[i]+'/dss/fits/'
	photpath=!GLGA_ROOT+'data/'+deg[i]+'/photometry/'
	plotpath=!GLGA_ROOT+'data/'+deg[i]+'/plots/'
	fpath=!GLGA_ROOT+'data/'+deg[i]+'/'+srvy+'/fits/'
	jpath=!GLGA_ROOT+'data/'+deg[i]+'/'+srvy+'/jpg/'
        auxbase = auxpath+id[i]+'_'+srvy
        
; get object data
	r=ra[i]
	d=dec[i]
	p=pa[i]>0.
	mjdiam=majordiam[i]
	mndiam=minordiam[i]
;
; get object file names
	maskimgfile = auxpath+id[i]+'_'+srvy+'_mask.fits.gz'
	
	ellipsefile=auxpath+id[i]+'_ellipse.dat'
	dssfile=dsspath+id[i]+'_dss2_red.fits*'
	for j=0,nbands-1 do begin
		
		bandmask    = auxpath+id[i]+'_'+srvy+'_mask_'+bands[j]+'.fits.gz'
		
		; GALEX data
		if strpos(srvy,'galex') ge 0 then begin
			galex_radprof,id[i],r,d,$
 				mjdiam,mndiam,p,$
				fpath+id[i]+'_'+bands[j]+'.fits*',$
				fpath+id[i]+'_'+bands[j]+'_rr.fits*',$
				fpath+id[i]+'_'+bands[j]+'_cnt.fits*',$
				annuli_size=annuli_size, $
 				ellipsefile=ellipsefile,$
 				diam_units='arcmin', maskimgfile=maskimgfile,$
 				outpath=photpath, status=pstat, $
 				/verbose, /extend, /yuan13
		;
		; SDSS data
		endif else if strpos(srvy,'sdss') ge 0 then begin
			sdss_radprof,id[i],r,d, $
				mjdiam,mndiam,p, $
				fpath+id[i]+'_'+bands[j]+'.fits*', $
				annuli_size=annuli_size, $
				ellipsefile=ellipsefile, $
				diam_units='arcmin', maskimgfile=maskimgfile, bandmask=bandmask,$
				outpath=photpath,status=pstat, auxbase=auxbase,$
				/verbose, /extend, /yuan13
				
		; PanSTARRS
		endif else if strpos(srvy,'panstarrs') ge 0 then begin
		        
; 		        print, "Ehsan band:", bands[j]
			panstarrs_radprof,id[i],r,d, $
				mjdiam,mndiam,p, $
				fpath+id[i]+'_'+bands[j]+'.fits*', $
				annuli_size=annuli_size, $
				ellipsefile=ellipsefile, $
				diam_units='arcmin', maskimgfile=maskimgfile, bandmask=bandmask,$
				outpath=photpath,status=pstat, auxbase=auxbase, $
				/verbose, /extend, /yuan13
		; ACS
		endif else if strpos(srvy,'acs') ge 0 then begin
		        
		        bandflag=fpath+id[i]+'_'+bands[j]+'_flag.fits*'
; 		        print, "Ehsan band:", bands[j]
			acs_radprof,id[i],r,d, $
				mjdiam,mndiam,p, $
				fpath+id[i]+'_'+bands[j]+'.fits*', $
				annuli_size=annuli_size, $
				ellipsefile=ellipsefile, filter=bands[j],  $
				diam_units='arcmin', maskimgfile=maskimgfile,$
				bandflag=bandflag, outpath=photpath,status=pstat, $
				/verbose, /extend, /yuan13		
		
		; 2MASS data
		endif else if strpos(srvy,'2mass') ge 0 then begin
			twomass_radprof,id[i],r,d, $
				mjdiam,mndiam,p, $
				fpath+id[i]+'_'+bands[j]+'.fits*', $
				annuli_size=annuli_size, $
				ellipsefile=ellipsefile, $
				diam_units='arcmin', maskimgfile=maskimgfile,$
				outpath=photpath,status=pstat, $
				/verbose, /extend, /yuan13
		;
		; WISE data
		endif else if strpos(srvy,'wise') ge 0 then begin
			wise_radprof,id[i],r,d, $
				mjdiam,mndiam,p, $
				fpath+id[i]+'_'+bands[j]+'.fits*', $
				annuli_size=annuli_size, $
				ellipsefile=ellipsefile, $
				diam_units='arcmin', maskimgfile=maskimgfile, bandmask=bandmask,$
				outpath=photpath,status=pstat, auxbase=auxbase, background_radius=background_radius, $
				/verbose, /extend, /yuan13
		;
		; Spitzer IRAC data
		endif else if strpos(srvy,'irac') ge 0 then begin
			irac_radprof,id[i],r,d, $
				mjdiam,mndiam,p, $
				fpath+id[i]+'_'+bands[j]+'.fits*', $
				annuli_size=annuli_size, $
				ellipsefile=ellipsefile, $
				diam_units='arcmin', maskimgfile=maskimgfile,$
				outpath=photpath,status=pstat, $
				/verbose, /extend, /yuan13
		endif
		;
		; accumulate status
		if pstat lt 0 then $
			stat = stat + 1
	endfor	; loop over bands
	print,i0,'/',n_elements(id),': ',id[i],'  DONE'
	;
	; create plots
	if stat lt nbands then begin
	    if strpos(srvy,'galex') ge 0 then begin
		galex_plotradprof, id[i], type=type[i], $
			pathtoprofile=photpath, $
    			fintfile=fpath+id[i]+'_FUV.fits*', $
			nintfile=fpath+id[i]+'_NUV.fits*', $
    			maskimgfile=maskimgfile,$
    			uvjpgpath=jpath, dssfile=dssfile, $
			outpath=plotpath, /yuan13
	    endif else if strpos(srvy,'sdss') ge 0 then begin
		sdss_plotradprof, id[i], type=type[i], $
			pathtoprofile=photpath, $
			intfile=fpath+id[i]+'_r.fits*', $
			maskimgfile=maskimgfile, auxbase=auxbase,$
			jpgpath=jpath, outpath=plotpath, /yuan13
	    endif else if strpos(srvy,'panstarrs') ge 0 then begin
		panstarrs_plotradprof, id[i], type=type[i], $
			pathtoprofile=photpath, $
			intfile=fpath+id[i]+'_r.fits*', $
			maskimgfile=maskimgfile, auxbase=auxbase, $
			jpgpath=jpath, outpath=plotpath, /yuan13			
	    endif else if strpos(srvy,'acs') ge 0 then begin
	        
		acs_plotradprof, id[i], type=type[i], $
			pathtoprofile=photpath, $
			intfile=fpath+id[i]+'_F814W.fits*', $
			maskimgfile=maskimgfile, fpath=fpath, $
			jpgpath=jpath, outpath=plotpath, /yuan13			    
	    endif else if strpos(srvy,'2mass') ge 0 then begin
		twomass_plotradprof, id[i], type=type[i], $
			pathtoprofile=photpath, $
			intfile=fpath+id[i]+'_j.fits*', $
			maskimgfile=maskimgfile, $
			jpgpath=jpath, outpath=plotpath, /yuan13
	    endif else if strpos(srvy,'wise') ge 0 then begin
		wise_plotradprof, id[i], type=type[i], $
			pathtoprofile=photpath, $
			intfile=fpath+id[i]+'_w1.fits*', $
			maskimgfile=maskimgfile, auxbase=auxbase, $
			jpgpath=jpath, outpath=plotpath, /yuan13
	    endif else if strpos(srvy,'irac') ge 0 then begin
		irac_plotradprof, id[i], type=type[i], $
			pathtoprofile=photpath, $
			intfile=fpath+id[i]+'_4p5um.fits*', $
			maskimgfile=maskimgfile, dssfile=dssfile, $
			jpgpath=jpath, outpath=plotpath, /yuan13
	    endif
	endif else begin
		print,'Photometry not done for: '+id[i]
		if keyword_set(logfile) then $
		    printf,ll,id[i],'No photometry performed',form='(a-25,2x,a)'
	endelse

endfor

if keyword_set(logfile) then free_lun,ll

return
end
