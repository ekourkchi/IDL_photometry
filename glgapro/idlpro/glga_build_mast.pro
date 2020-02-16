pro glga_build_mast, lfile, xd=xd, $
build_dss=build_dss, $
uvoutpath=uvoutpath,$
dssoutpath=dssoutpath,$
summarylog=summarylog,$
errorlog=errorlog,$
sample_id=sample_id,$
forceupdate=forceupdate,$
exclude_all_post_csp_data=exclude_all_post_csp_data,$
use_mains=use_mains,$
stagepath=stagepath

; exclude_post_csp_data = exclude all data take after 
; the coarse sun pointing event on 5/4/2010 04:30:00
; this requires that visit level data be used for 
; building products. Note the worst uncorrected data 
; is always excluded (123 images) 5/4-6/24/2010.
;
; mains = use coadds (mains) to build instead of visits
; which will be faster and require few files to stage.
; default is to use visit level data.


readcol,lfile, id, ra, dec, majordiam , format='a,d,d,f', comment='#',/silent

;;;;;;;;;
;defaults

;subdirectory designation (000D to 359D)
deg = string(floor(ra), format='(i3.3)')+'D'

if not keyword_set(uvoutpath) then $
 uvoutpath = !GLGA_ROOT+'data/'+deg+'/galex/fits/'

if not keyword_set(dssoutpath) then $
 dssoutpath = !GLGA_ROOT+'data/'+deg+'/dss/fits/'

if not keyword_set(sample_id) then get_date,sample_id

if not keyword_set(summarylog) then $
 summarylog = !GLGA_ROOT+'work/logs/'+sample_id+'_summary.txt'

if not keyword_set(errorlog) then $
 errorlog = !GLGA_ROOT+'work/logs/'+sample_id+'_errors.txt'

if not keyword_set(xd) then $
 xd = 4 ; extract data over the area xd*majordiam
;
; 10 arcmin is minimum size
msize=(majordiam*xd)>10.
bg = where(majordiam ge 50., nbg)
if nbg gt 0 then msize[bg] = majordiam[bg] * 1.5
bg = where(majordiam lt 50. and majordiam ge 30., nbg)
if nbg gt 0 then msize[bg] = majordiam[bg] * 2.0
bg = where(majordiam lt 30 and majordiam ge 20., nbg)
if nbg gt 0 then msize[bg] = majordiam[bg] * 2.5
bg = where(majordiam lt 20 and majordiam ge 10., nbg)
if nbg gt 0 then msize[bg] = majordiam[bg] * 3.0


;default to visits
mastpathfile=!GLGA_ROOT+'work/GR67_MAST_visits.fits.gz'

if keyword_set(use_mains) then $
 mastpathfile=!GLGA_ROOT+'work/GR67_MAST_mains.fits.gz'


if not keyword_set(forceupdate) then $
 forceupdate = 0

;;;;;;;;;;;;;;;;;;;;;;;;
; open log files

if not file_exist(summarylog) then begin
 openw, slog, summarylog, /get_lun, /append
 printf, slog, $
  "# date, id, FUVvisits, NUVvisits, medrr_ftime, medrr_ntime, ra, dec, majordiam, path"
 free_lun, slog
endif

if not file_exist(errorlog) then begin
 openw, elog, errorlog, /get_lun, /append
 printf, elog, "# i, id, issue"
 free_lun, elog
endif


;;;;;;;;;;;;;;;;;;;;;;;;;
;read mast path file

g=mrdfits(mastpathfile,1, /SILENT)
G.NUV_EXPTIME = G.NUV_EXPTIME > 0
G.FUV_EXPTIME = G.FUV_EXPTIME > 0
g.path = strcompress(g.path,/REM)
g.tile = strcompress(g.tile,/REM)
g.base = strcompress(g.base,/REM)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;FILTER ON USABLE IMAGE DATA

A = WHERE( g.nuv_exptime gt 0 or g.fuv_exptime gt 0, count)
g = g[a]

;allways exclude worst coarse sunpoint uncorrected data
; if not using mains (only 123 visits)
if not keyword_set(use_mains) then begin  
 datetime=g.obs_date+'t'+g.obs_time
 a=where(datetime lt '2010-05-04t04:30:00' or datetime ge '2010-06-24t00:00:00')
 g = g[a]
 delvarx,datetime,/free
endif

if keyword_set(exclude_all_post_csp_data) and $
   not keyword_set(use_mains) then begin  
 datetime=g.obs_date+'t'+g.obs_time
 a=where(datetime lt '2010-05-04t04:30:00')
 g = g[a]
 delvarx,datetime,/free
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;loop through the batch from input catalog

for i = 0, n_elements(ra)-1 do begin

 print, '** '+id[i]+': '+strcompress(string(i+1),/rem)+' of '+$
        strcompress(string(n_elements(ra)),/rem)


 gcirc, 1, ra[i]/15., dec[i], g.ra/15., g.dec,  dist

 ; check for good overlap first 
 a = where(dist/60. le 36 and $
           (g.fuv_exptime gt 10 or g.nuv_exptime gt 10),  count)

 ; now get all fuv and nuv in extended region

 f = where(dist/60. le 35 + msize[i] and $
            g.fuv_exptime gt 10,  countfuv)
 n = where(dist/60. le 35 + msize[i] and $
            g.nuv_exptime gt 10,  countnuv)

 ; if there is nothing within 36 arcmin then covereage is too poor
 if count eq 0 then begin
  countfuv=0
  countnuv=0 
 endif


 uv:

 outpath = uvoutpath[i]

 ;;;;;;;;;;
 ;build NUV

 nuvstatus=0

 if countnuv gt 0 then begin

  detailfile = outpath+strcompress(id[i],/rem)+'_detail_nuv.txt'
  nuvintfile=outpath+strcompress(id[i],/rem)+'_NUV.fits.gz'

  ; see if updates are needed 
  if file_exist(nuvintfile) and not forceupdate then begin
    print,'** NUV Coadd already built. Checking data...'
    curhdr=headfits(nuvintfile)
    curtime=sxpar(curhdr,'exptime')
    if total(g[n].nuv_exptime) le curtime then begin
      print,'** No need to update NUV.'
      nuvstatus=1
      goto, skip_nuv_update
    endif
    print,'** NUV Update needed.'
  endif 

  get_galex_from_mast, g[n].base, g[n].path, $ 
    stagepath=stagepath,$ 
    stagesubdir=strcompress(id[i],/rem), $
    /nuv, /cnt, /rr

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'Extracting/building NUV: '+strcompress(countnuv)+' images'

  pathlist=strarr(n_elements(n))+stagepath+'/'+strcompress(id[i],/rem)+'/'

  galex_como, strcompress(id[i],/rem), ra[i], dec[i], msize[i], $
    strcompress(g[n].base, /rem), strcompress(pathlist, /rem), $  
    outdir = outpath, /nuv, status=nuvstatus


  print,'NUV status = '+strcompress(nuvstatus)


  ;rename files to be consistem with glga format 
  if nuvstatus then begin
    file1=outpath+strcompress(id[i],/rem)+'-nd-int.fits.gz'
    file2=outpath+strcompress(id[i],/rem)+'_NUV.fits.gz'
    if file_exist(file1) then $
     spawn,'mv '+file1+' '+file2
    file1=outpath+strcompress(id[i],/rem)+'-nd-cnt.fits.gz'
    file2=outpath+strcompress(id[i],/rem)+'_NUV_cnt.fits.gz'
    if file_exist(file1) then $
     spawn,'mv '+file1+' '+file2
    file1=outpath+strcompress(id[i],/rem)+'-nd-rrhr.fits.gz'
    file2=outpath+strcompress(id[i],/rem)+'_NUV_rr.fits.gz'
    if file_exist(file1) then $
     spawn,'mv '+file1+' '+file2
  endif

  ;update detail file 
  if nuvstatus then begin
   openw, lun, detailfile, /get_lun
   !textunit = lun
   printf, lun, "# NUV visit level data suggested for coadd"
   printf, lun,$ 
    "# tile, base, fov_dist('), ntime,  path"
   forprint, strcompress(g[n].tile, /rem), strcompress(g[n].base, /rem),  $
        dist[n]/60, g[n].nuv_exptime, $
        strcompress(g[n].path, /rem), $
        format = 'a-35,1x,a-45,1x,f6.2,1x,f10.2,1x,a', $
        textout = 5, /silent, /nocom
   free_lun, lun
  endif else begin
   openw, elog, errorlog, /get_lun, /append
   printf, elog, i, id[i]+' NUV coadd fail',format='(i5,2x,a)'
   free_lun, elog
  endelse



  skip_nuv_update:

 endif 

 ;end NUV build
 ;;;;;;;;;;;;;;

 ;;;;;;;;;;;
 ;;;;;;;;;;;
 ;build FUV

 fuvstatus=0

 if countfuv gt 0 then begin

  detailfile = outpath+strcompress(id[i],/rem)+'_detail_fuv.txt'
  fuvintfile=outpath+strcompress(id[i],/rem)+'_FUV.fits.gz'

  ; see if updates are needed 
  if file_exist(fuvintfile) and not forceupdate then begin
    print,'** FUV Coadd already built. Checking data...'
    curhdr=headfits(fuvintfile)
    curtime=sxpar(curhdr,'exptime')
    if total(g[n].fuv_exptime) le curtime then begin
      print,'** No need to update FUV.'
      fuvstatus=1
      goto, skip_fuv_update
    endif
    print,'** FUV Update needed.'
  endif 

  get_galex_from_mast, g[f].base, g[f].path, $ 
    stagepath=stagepath,$ 
    stagesubdir=strcompress(id[i],/rem), $
    /fuv, /cnt, /rr

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'Extracting/building FUV: '+strcompress(countfuv)+' images'

  pathlist=strarr(n_elements(f))+stagepath+'/'+strcompress(id[i],/rem)+'/'

  galex_como, strcompress(id[i],/rem), ra[i], dec[i], msize[i], $
    strcompress(g[f].base, /rem), strcompress(pathlist, /rem), $  
    outdir = outpath, /fuv, status=fuvstatus



  ;rename files to be consistem with glga format 
  if nuvstatus then begin
    file1=outpath+strcompress(id[i],/rem)+'-fd-int.fits.gz'
    file2=outpath+strcompress(id[i],/rem)+'_FUV.fits.gz'
    if file_exist(file1) then $
     spawn,'mv '+file1+' '+file2
    file1=outpath+strcompress(id[i],/rem)+'-fd-cnt.fits.gz'
    file2=outpath+strcompress(id[i],/rem)+'_FUV_cnt.fits.gz'
    if file_exist(file1) then $
     spawn,'mv '+file1+' '+file2
    file1=outpath+strcompress(id[i],/rem)+'-fd-rrhr.fits.gz'
    file2=outpath+strcompress(id[i],/rem)+'_FUV_rr.fits.gz'
    if file_exist(file1) then $
     spawn,'mv '+file1+' '+file2
  endif

  ;update detail file 

  print,'FUV status = '+strcompress(fuvstatus)

  if fuvstatus then begin
   openw, lun, detailfile, /get_lun
   !textunit = lun
   printf, lun, "# FUV visit level data suggested for coadd"
   printf, lun,$ 
    "# tile, base, fov_dist('), ftime,  path"
   forprint, strcompress(g[f].tile, /rem), strcompress(g[f].base, /rem),  $
        dist[f]/60, g[f].fuv_exptime, $
        strcompress(g[f].path, /rem), $
        format = 'a-35,1x,a-45,1x,f6.2,1x,f10.2,1x,a', $
        textout = 5, /silent, /nocom
   free_lun, lun

  endif else begin
   openw, elog, errorlog, /get_lun, /append
   printf, elog, i, id[i]+' FUV coadd fail',format='(i5,2x,a)'
   free_lun, elog
  endelse

  skip_fuv_update:

 endif

 ;end FUV build
 ;;;;;;;;;;;;;;

 
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;update extraction summary log

 ;if count gt 0 then begin ;report even if no data found

  if nuvstatus then begin
   hdr=headfits(outpath+strcompress(id[i],/rem)+'_NUV.fits.gz')
   ntime=sxpar(hdr,'MEDRR')
   nuvadded=sxpar(hdr,'NADDED')
  endif else begin
   ntime=0
   nuvadded=0
  endelse
  
  if fuvstatus then begin
   hdr=headfits(outpath+strcompress(id[i],/rem)+'_FUV.fits.gz')
   ftime=sxpar(hdr,'MEDRR')
   fuvadded=sxpar(hdr,'NADDED')
  endif else begin
   ftime=0
   fuvadded=0
  endelse


  get_date, dte
  openw, slog, summarylog, /get_lun, /append
  printf, slog, dte, id[i], fuvadded, nuvadded, $
    ftime, ntime, $
    ra[i], dec[i], majordiam[i], outpath,$
    format = '(a10,2x,a20,2x,i4,2x,i4,2x,f10.2,f10.2,f12.6,f12.6,2x,f8.2,2x,a)'
  free_lun, slog

  print, strcompress(i+1, /rem)+' of '+strcompress(n_elements(ra), /rem)
  print, id[i],  fuvadded, nuvadded,  $
    ftime, ntime, $
    ra[i], dec[i], majordiam[i], $
    format = '(a25,2x,i4,2x,i4,2x,f10.2,f10.2,f12.6,f12.6,2x,f8.2)'
  
 ;endif

 ;end extraction summary status
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 ;goto, skipdss

 ;;;;;;;;;;;;;;;
 ;get dss image
 dss:

 if (count gt 0 and (fuvstatus or nuvstatus) and keyword_set(build_dss)) $
    then begin 

  outpath = dssoutpath[i]

  if file_exist(outpath+strcompress(id[i],/rem)+'*') then goto, skipdss 

  get_dss, strcompress(id[i],/rem), ra[i], dec[i], msize[i], $
     'red',outdir=outpath, status=dssstatus

  if not dssstatus then $
  get_dss, strcompress(id[i],/rem), ra[i], dec[i], msize[i], $
     'blue',outdir=outpath, status=dssstatus

  if not dssstatus then $
  get_dss, strcompress(id[i],/rem), ra[i], dec[i], msize[i], $
     'ir',outdir=outpath, status=dssstatus

  if not dssstatus then begin
    openw, elog, errorlog, /get_lun, /append
    printf, elog, i, id[i]+' DSS download fail',format='(i5,2x,a)'
    free_lun, elog
  endif

  skipdss:

 endif 

 ;end get dss image
 ;;;;;;;;;;;;;;;;;;


skipover:
endfor 

end


