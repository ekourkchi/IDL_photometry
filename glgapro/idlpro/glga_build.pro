pro glga_build, id, ra, dec, majordiam , xd=xd, $
build_dss=build_dss, $
uvoutpath=uvoutpath,$
dssoutpath=dssoutpath,$
summarylog=summarylog,$
errorlog=errorlog,$
missionstatusfile=missionstatusfile,$
gi_age_min=gi_age_min,$
sample_id=sample_id,$
forceupdate=forceupdate,$
maxeclipse=maxeclipse

;options to be added
;fuv=fuv,nuv=nuv,
;forceupdate=forceupdate,$


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
; 3 arcmin is minimum size
msize=(majordiam*xd)>10.
bg = where(majordiam ge 50., nbg)
if nbg gt 0 then msize[bg] = majordiam[bg] * 1.5
bg = where(majordiam lt 50. and majordiam ge 30., nbg)
if nbg gt 0 then msize[bg] = majordiam[bg] * 2.0
bg = where(majordiam lt 30 and majordiam ge 20., nbg)
if nbg gt 0 then msize[bg] = majordiam[bg] * 2.5
bg = where(majordiam lt 20 and majordiam ge 10., nbg)
if nbg gt 0 then msize[bg] = majordiam[bg] * 3.0

if not keyword_set(missionstatusfile) then $
 missionstatusfile=!GLGA_ROOT+'work/mission_status_sg.fits'
;'/home/galex/fltops/logs/mission_status/mission_status.fits'

if not keyword_set(gi_age_min) then $
 gi_age_min = 8*30 ; GI data should be 8 months old to include

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
;read mission status file

hdr=headfits(missionstatusfile,ext = 0)
msdate = sxpar(hdr, 'RUNDAT')
g=mrdfits(missionstatusfile,1, /silent)
g.nuv_exptime = g.nuv_exptime > 0
g.fuv_exptime = g.fuv_exptime > 0
g.surv_type = strcompress(g.surv_type,/rem)
g.grelease = strcompress(g.grelease,/rem)
g.qa_grade = strcompress(g.qa_grade,/rem)

if not keyword_set(maxeclipse) then $
 maxeclipse = max(g.eclipse)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;filter on usable image data

a = where( (g.nuv_exptime gt 0 or g.fuv_exptime gt 0) and $
           g.surv_type ne 'SIO' and g.surv_type ne 'SSO' and $
           g.surv_type ne 'SOO' and g.ow eq 'd' and g.qa_grade ne 'FAIL' and $
           g.asp_ave_ra_rta ne -9999. and g.asp_ave_dec_rta ne -9999. and $
           g.eclipse le maxeclipse, count)

g = g[a]

;filter out GI data less thean age minimum

year='20'+strmid(g.ecl_start,0,2)
month=strmid(g.ecl_start,2,2)
day=strmid(g.ecl_start,4,2)
obsdate=julday(month,day,year)
curdate=systime(/julian)
elapsed_days = curdate-obsdate


a = where(g.surv_type eq 'GII' and elapsed_days lt gi_age_min, complement=b)
if a[0] ge 0 then g=g[b]


;if no qa_grade change to NA
a = where(g.qa_grade eq '')
if a[0] ge 0 then g[a].qa_grade =  'NA'


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;loop through the batch from input catalog

for i = 0, n_elements(ra)-1 do begin

 print, '** '+id[i]+': '+strcompress(string(i+1),/rem)+' of '+$
        strcompress(string(n_elements(ra)),/rem)


 gcirc, 1, ra[i]/15., dec[i], g.asp_ave_ra_rta/15., g.asp_ave_dec_rta,  dist



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

  ; see if updates are needed 
  if file_exist(detailfile) and not forceupdate then begin
    print,'** NUV Coadd already built. Checking data...'
    readcol, detailfile, eclipse, survey, grelease, qa_grade, tile, $
      base, fov_dist, ntime, format='l,a,a,a,a,a,f,f'
    if total(g[n].nuv_exptime) le total(ntime) then begin
      print,'** No need to update NUV.'
      nuvstatus=1
      goto, skip_nuv_update
    endif
    print,'** NUV Update needed.'
  endif 

  print,'Extracting/building NUV: '+strcompress(countnuv)+' visits'

  galex_como, strcompress(id[i],/rem), ra[i], dec[i], msize[i], $
    strcompress(g[n].base, /rem), strcompress(g[n].info_str, /rem), $  
    outdir = outpath, /nuv, status=nuvstatus


  print,'NUV status = '+strcompress(nuvstatus)

  if nuvstatus then begin
   openw, lun, detailfile, /get_lun
   !textunit = lun
   printf, lun, "# NUV visit level data suggested for coadd"
   printf, lun,$ 
    "# eclipse, survey, grelease, qa_grade, tile, base, fov_dist('), ntime,  path"
   forprint, g[n].eclipse, g[n].surv_type, g[n].grelease, g[n].qa_grade, $
        strcompress(g[n].tile, /rem), strcompress(g[n].base, /rem),  $
        dist[n]/60, g[n].nuv_exptime, $
        strcompress(g[n].info_str, /rem), $
        format = 'i6,1x,a3,1x,a12,1x,a7,1x,a31,1x,a36,1x,f6.2,f10.2, a100', $
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


 ;;;;;;;;;;
 ;build fUV

 fuvstatus=0

 if countfuv gt 0 then begin

  detailfile = outpath+strcompress(id[i],/rem)+'_detail_fuv.txt'

  ; see if updates are needed
  if file_exist(detailfile) and not forceupdate then begin
    print,'** FUV Coadd already built. Checking data...'
    readcol, detailfile, eclipse, survey, grelease, qa_grade, tile, $
      base, fov_dist, ftime, format='l,a,a,a,a,a,f,f'
    if total(g[f].fuv_exptime) le total(ftime) then begin
      print,'** No need to update FUV.'
      fuvstatus=1
      goto, skip_fuv_update
    endif
    print,'** FUV Update needed.'
  endif

  print,'Extracting/building FUV: '+strcompress(countfuv)+' visits'

  galex_como, strcompress(id[i],/rem), ra[i], dec[i], msize[i], $
    strcompress(g[f].base, /rem), strcompress(g[f].info_str, /rem), $  
    outdir = outpath, /fuv, status=fuvstatus


  print,'FUV status = '+strcompress(fuvstatus)

  if fuvstatus then begin
   openw, lun, detailfile, /get_lun
   !textunit = lun
   printf, lun, "# FUV visit level data suggested for coadd"
   printf, lun, $
    "# eclipse, survey, grelease, qa_grade, tile, base, fov_dist('), ftime,  path"
   forprint, g[f].eclipse, g[f].surv_type, g[f].grelease, g[f].qa_grade, $
        strcompress(g[f].tile, /rem), strcompress(g[f].base, /rem),  $
        dist[f]/60, g[f].fuv_exptime, $
        strcompress(g[f].info_str, /rem), $
        format = 'i6,1x,a3,1x,a12,1x,a7,1x,a31,1x,a36,1x,f6.2,f10.2, a100', $
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


