pro array_clean,value,low,high,ar1,ar2,ar3,ar4,ar5,ar6,ar7,ar8,ar9,nan=nan
;+
; array_clean - replace array values outside range [low,high] with VALUE
;
;-
if n_params(0) lt 4 then begin
	print,'Usage: array_clean, value, low, high, array1 <, array2, ...>'
	return
endif
;
switch n_params(0) of
12: begin 
	if keyword_set(nan) then $
		t=where(ar9 lt low or ar9 gt high or (finite(ar9) eq 0), nt) $
	else	t=where(ar9 lt low or ar9 gt high, nt)
	if nt gt 0 then ar9[t] = value
    end
11: begin 
	if keyword_set(nan) then $
		t=where(ar8 lt low or ar8 gt high or (finite(ar8) eq 0), nt) $
	else	t=where(ar8 lt low or ar8 gt high, nt)
	if nt gt 0 then ar8[t] = value
    end
10: begin 
	if keyword_set(nan) then $
		t=where(ar7 lt low or ar7 gt high or (finite(ar7) eq 0), nt) $
	else	t=where(ar7 lt low or ar7 gt high, nt)
	if nt gt 0 then ar7[t] = value
    end
 9: begin 
	if keyword_set(nan) then $
		t=where(ar6 lt low or ar6 gt high or (finite(ar6) eq 0), nt) $
	else	t=where(ar6 lt low or ar6 gt high, nt)
	if nt gt 0 then ar6[t] = value
    end
 8: begin 
	if keyword_set(nan) then $
		t=where(ar5 lt low or ar5 gt high or (finite(ar5) eq 0), nt) $
	else	t=where(ar5 lt low or ar5 gt high, nt)
	if nt gt 0 then ar5[t] = value
    end
 7: begin 
	if keyword_set(nan) then $
		t=where(ar4 lt low or ar4 gt high or (finite(ar4) eq 0), nt) $
	else	t=where(ar4 lt low or ar4 gt high, nt)
	if nt gt 0 then ar4[t] = value
    end
 6: begin 
	if keyword_set(nan) then $
		t=where(ar3 lt low or ar3 gt high or (finite(ar3) eq 0), nt) $
	else	t=where(ar3 lt low or ar3 gt high, nt)
	if nt gt 0 then ar3[t] = value
    end
 5: begin 
	if keyword_set(nan) then $
		t=where(ar2 lt low or ar2 gt high or (finite(ar2) eq 0), nt) $
	else	t=where(ar2 lt low or ar2 gt high, nt)
	if nt gt 0 then ar2[t] = value
    end
 4: begin 
	if keyword_set(nan) then $
		t=where(ar1 lt low or ar1 gt high or (finite(ar1) eq 0), nt) $
	else	t=where(ar1 lt low or ar1 gt high, nt)
	if nt gt 0 then ar1[t] = value
    end
endswitch
;
return
end
