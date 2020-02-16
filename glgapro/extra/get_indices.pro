;=============================================================================
;
; GET_INDICES
;
; In IDL format, with indices of form (i,j,k) running from 0 to I-1,J-1,K-1,
; respectively, a single index n is broken into its components via the
; following steps :
;
;  i = n mod I
;  j = fix(n/I) mod J
;  k = fix(n/(I*J))
;
;  e.g. for I=2,J=3,K=4
;
;  n  i  j  k         n  i  j  k         n  i  j  k
; ------------       ------------       ------------
;  0  0  0  0         6  0  0  1        12  0  0  2
;  1  1  0  0         7  1  0  1        13  1  0  2
;  2  0  1  0         8  0  1  1        14  0  1  2
;  3  1  1  0         9  1  1  1        15  1  1  2
;  4  0  2  0        10  0  2  1        16  0  2  2
;  5  1  2  0        11  1  2  1        17  1  2  2
;
;  EXAMPLE:
;      Given a 2 dimensional image of dimensions 512 x 512.  Where is
;      the maximum value?
;
;        m=where(image eq max(image))
;        xy=get_indices(m,[512,512])
;        print,xy
;
;      or more succinctly: 
;
;        print,get_indices(where(image eq max(image)),[512,512])
;
;      Generates pixel coordinates in x and y
;
;===========================================================================
function get_indices,list,dimensions

dim = n_elements(dimensions)
listout=lonarr(dim,n_elements(list))

FOR i=0l,n_elements(list)-1 DO BEGIN
  ii=list[i] 
  indices = lonarr(dim)
  case dim of
    0: indices = [0]
    1: indices = [ii]
    2: begin
         indices(0) = ii mod dimensions(0)
         indices(1) = fix(ii/dimensions(0))
       end
    3: begin
         indices(0) = ii mod dimensions(0)
         indices(1) = fix(float(ii)/dimensions(0)) mod dimensions(1)
         indices(2) = fix(ii/(dimensions(0)*dimensions(1)))
       end
    4: begin
         indices(0) = ii mod dimensions(0)
         indices(1) = fix(float(ii)/dimensions(0)) mod dimensions(1)
         indices(2) = fix(ii/(dimensions(0)*dimensions(1))) mod dimensions(2)
         indices(3) = fix(ii/(dimensions(0)*dimensions(1)*dimensions(2)))
       end
    5: begin
         indices(0) = ii mod dimensions(0)
         indices(1) = fix(float(ii)/dimensions(0)) mod dimensions(1)
         indices(2) = fix(ii/(dimensions(0)*dimensions(1))) mod dimensions(2)
         indices(3) = fix(ii/(dimensions(0)*dimensions(1)*dimensions(2))) $
		      mod dimensions(3)
         indices(4) = fix(ii/(dimensions(0)*dimensions(1)*dimensions(2)* $
		      dimensions(3)))
       end
    ELSE: begin
         print,'Too many indices!! What are you trying to do, anyways?'
         indices = [0]
       end
  endcase

  listout[*,i]=indices


 ENDFOR 

  return,listout
end
