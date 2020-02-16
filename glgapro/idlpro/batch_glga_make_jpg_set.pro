pro batch_glga_make_jpg_set, lfile, verbose=verbose,sdss=sdss,twomass=twomass,$
        galex=galex, wise=wise, panstarrs=panstarrs, update=update, start=start
;+
; glga make_jpg_set - make the set of jpg file needed for the GLGA
;
; set UPDATE keyword to overwrite existing files
;-

readcol, lfile, id, ra, dec, majdiam, mindiam, pa, type, format='a,d,d,f,f,f,a', /silent

deg = string(floor(ra), format='(i3.3)')+'D'
d0 = 0L
if keyword_set(start) then d0 = start-1L

for d=d0,n_elements(id)-1 do begin

print, "List row: ", d, '  ', id[d]

;
; 2MASS
if keyword_set(twomass) then begin
  fspec=!GLGA_ROOT+'data/'+deg[d]+'/2mass/fits/'+id[d]+'_k.fit*'
  outdir=!GLGA_ROOT+'data/'+deg[d]+'/2mass/jpg/'
  rgb_2mass_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose
  rgb_2mass_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/jonly
  rgb_2mass_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/honly
  rgb_2mass_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/konly
;
; SDSS
endif else if keyword_set(sdss) then begin
  fspec=!GLGA_ROOT+'data/'+deg[d]+'/sdss/fits/'+id[d]+'_i.fit*'
  outdir=!GLGA_ROOT+'data/'+deg[d]+'/sdss/jpg/'
  rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose
  rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/wide
  rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/uonly
  rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/gonly
  rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/ronly
  rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/ionly
  rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/zonly
;
; WISE
endif else if keyword_set(wise) then begin
  fspec=!GLGA_ROOT+'data/'+deg[d]+'/wise/fits/'+id[d]+'_w1.fit*'
  outdir=!GLGA_ROOT+'data/'+deg[d]+'/wise/jpg/'
  rgb_wise_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose
  rgb_wise_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/wide
  rgb_wise_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/w1only
  rgb_wise_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/w2only
  rgb_wise_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/w3only
  rgb_wise_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/w4only
;
; Pan_STARRS
endif else if keyword_set(panstarrs) then begin
  fspec=!GLGA_ROOT+'data/'+deg[d]+'/panstarrs/fits/'+id[d]+'_z.fit*'
  outdir=!GLGA_ROOT+'data/'+deg[d]+'/panstarrs/jpg/'
  rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose
  rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/gonly
  rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/ronly
  rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/ionly
  rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/zonly
; GALEX
endif else if keyword_set(galex) then begin
  indir=!GLGA_ROOT+'data/'+deg[d]+'/galex/fits/'
  outdir=!GLGA_ROOT+'data/'+deg[d]+'/galex/jpg/'
  fuvfile= indir+id[d]+'_FUV.fits.gz'
  nuvfile= indir+id[d]+'_NUV.fits.gz'

; NUV 
  nuvsmooth=0
  if file_exist(nuvfile) then begin
    hdr=headfits(nuvfile)
    ntime=sxpar(hdr,'EXPTIME')
    if ntime le 300 then $
	    nuvsmooth = 1 $
    else    nuvsmooth = 0
  endif

  target = outdir+id[d]+'_NUV'
  if file_test(target+'.jpg') eq 0 or keyword_set(update) then $
      rgb_galex_asinh, fuvfile+'XX', nuvfile, target, $
	/bgsub,/bgouter, /verbose, brite = 2,  nuvsmooth = nuvsmooth $
  else print,'Skipping NUV: ',target

;
; FUV
  target = outdir+id[d]+'_FUV'
  if file_test(target+'.jpg') eq 0 or keyword_set(update) then $
      rgb_galex_asinh, fuvfile, nuvfile+'XX', target, $
	/bgsub,/bgouter, /verbose, brite = 2,  nuvsmooth = 0 $
  else print,'Skipping FUV: ',target

;
; COMPOSITE
  target = outdir+id[d]+'_FUVNUV'
  if file_test(target+'.jpg') eq 0 or keyword_set(update) then $
      rgb_galex_asinh, fuvfile, nuvfile, target, $
	/bgsub,/bgouter, /verbose, brite = 2,  nuvsmooth = nuvsmooth $
  else print,'Skipping Composite: ',target


endif
;
if keyword_set(verbose) then $
	print,d+1L,'/',n_elements(id),' ',id[d],' Done.'
endfor 
;
return
end
