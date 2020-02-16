; Copyright (c) 1998, Forschungszentrum Juelich GmbH ICG-1
; All rights reserved.
; Unauthorized reproduction prohibited.
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.
;
;+
; NAME:
;   file_exist
;
; PURPOSE:
;   The result of this function is 1 if a file exist and 0 if not
;
; CATEGORY:
;   DATAFILES
;
; CALLING SEQUENCE:
;   Result=file_exist(file_name,[valid=valid])
;
; INPUTS:
;   file_name: The name of the File or an array of file_names
;
; KEYWORD PARAMETERS:
;   valid: if arg_present the valid indices are returned in valid
;   /VERBOSE: If set a 'file not found' message is issued for any file not found
;
; OUTPUTS:
;   This function returns 1 if the file exist and 0 if not
;
; EXAMPLE:
;   result=file_exist('otto.nc')
;
; MODIFICATION HISTORY:
;   Written by: R.Bauer (ICG-1),  1998-May-18
;   1999-June-29: searching of multiple files now possible
;   2001-Jun-22 : idl bug removed if the filename is longer than 319 sign. (idl5.4 crahes)
;   27.11.2002: Added keyword parameter VERBOSE
;   10.3.2003:  TB: undefined file_name stops IF n_file EQ 0 THEN RETURN, -1
;
;-

FUNCTION file_exist,file_name,valid=valid, VERBOSE = verbose
   IF N_PARAMS() LT 1 THEN BEGIN
      ;MESSAGE,call_help(),/cont
      RETURN,-1
   ENDIF

   n_file=N_ELEMENTS(file_name)
   IF n_file EQ 0 THEN RETURN, -1
   result=BYTARR(n_file)
   FOR i=0,n_file-1 DO BEGIN
      if strlen(file_name[i]) lt 320 then begin
      file=file_name[i]

      OPENR,lun,file,err=err,/GET_LUN
      IF N_ELEMENTS(lun) GT 0 THEN FREE_LUN,lun
      IF err NE 0 THEN BEGIN
        result[i]=0
        IF KEYWORD_SET(verbose) THEN MESSAGE, 'File not found: ' + file, /INFO
      ENDIF ELSE result[i]=1
      endif
   ENDFOR

   IF ARG_PRESENT(valid) THEN valid=WHERE(result NE 0)
   IF N_ELEMENTS(result) EQ 1 THEN RETURN,(result)[0] ELSE RETURN,result


END








