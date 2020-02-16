pro filestamp,lfile,verbose=verbose,archdir=archdir
;+
; filestamp - append a date stamp to the end of the current file
;
; keywords:
;	verbose - print operation
;	archdir - move stamped file into this directory
;-
fdecomp,lfile,disk,idir,file,ext
if keyword_set(archdir) then begin
	if size(archdir,/type) eq 7 then $
		arch = archdir $
	else	arch = 'arch'
	adir = idir + arch
	if not file_test(adir,/directory) then $
			file_mkdir,adir
	adir = adir + '/'
endif else adir = idir

if file_test(lfile) then begin
	finfo=file_info(lfile)
	mfile=adir+file+'.'+ext+'_'+timestr(finfo.ctime)
	if file_test(mfile) then $
		mfile=adir+file+'.'+ext+'_'+timestr(finfo.ctime,/time)
	file_move,lfile,mfile
	if keyword_set(verbose) then $
		print,'Moved '+lfile+' to '+mfile
endif
;
return
end
