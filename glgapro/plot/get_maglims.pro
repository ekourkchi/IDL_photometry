function get_maglims,mg,mge,pad_fraction=pad_fraction
;+
; get_maglims - get range of magnitudes suitable for plots
;
; INPUTS:
;	mg, mge 	- magnitudes and errors
;
; KEYWORDS:
;	pad_fraction	- fraction to pad the range with
;-
if keyword_set(pad_fraction) then $
	pf = pad_fraction $
else	pf = 0.1	; default 10% padding
;
xmax = max(mg+mge)
xmin = min(mg-mge)
rang = xmax - xmin
pad  = rang * pf
;
return,[xmax+pad,xmin-pad]
end
