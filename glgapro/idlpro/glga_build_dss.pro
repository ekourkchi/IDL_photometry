pro glga_build_dss, id, ra, dec, majordiam , xd=xd, $
dssoutpath=dssoutpath,$
errorlog=errorlog,$
sample_id=sample_id,$
forceupdate=forceupdate

;subdirectory designation (000D to 359D)
deg = string(floor(ra), format='(i3.3)')+'D'

if not keyword_set(dssoutpath) then $
 dssoutpath = !GLGA_ROOT+'data/'+deg+'/dss/fits/'

if not keyword_set(sample_id) then get_date,sample_id

if not keyword_set(errorlog) then $
 errorlog = !GLGA_ROOT+'work/logs/'+sample_id+'_dss_errors.txt'

if not keyword_set(xd) then $
 xd = 4 ; extract data over the area xd*majordiam

if not keyword_set(forceupdate) then $
 forceupdate = 0

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

;;;;;;;;;;;;;;;;;;;;;;;;
; open log files

if not file_exist(errorlog) then begin
 openw, elog, errorlog, /get_lun, /append
 printf, elog, "# i, id, issue"
 free_lun, elog
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;loop through the batch from input catalog

for i = 0, n_elements(ra)-1 do begin

 print, '** '+id[i]+': '+strcompress(string(i+1),/rem)+' of '+$
        strcompress(string(n_elements(ra)),/rem)

 ;;;;;;;;;;;;;;;
 ;get dss image
 dss:

  outpath = dssoutpath[i]

  if file_exist(outpath+strcompress(id[i],/rem)+'*') and not forceupdate then $
	  goto, skipdss 

  get_dss, strcompress(id[i],/rem), ra[i], dec[i], msize[i],'red', $
	  outdir=outpath, status=dssstatus

  if not dssstatus then $
  	get_dss, strcompress(id[i],/rem), ra[i], dec[i], msize[i],'blue', $
  		outdir=outpath, status=dssstatus

  if not dssstatus then $
  	get_dss, strcompress(id[i],/rem), ra[i], dec[i], msize[i],'ir', $
		outdir=outpath, status=dssstatus

  if not dssstatus then begin
    openw, elog, errorlog, /get_lun, /append
    printf, elog, i, id[i]+' DSS download fail',format='(i5,2x,a)'
    free_lun, elog
  endif

  skipdss:

endfor 

end
