PRO TABINV,XARR,X,IEFF
;+
; NAME:
;   TABINV     
; PURPOSE:  
;   To find the effective index of a function value in
;   an ordered vector.
; CALLING SEQUENCE:
;   TABINV,XARR,X,IEFF
; INPUTS:
;        XARR - the vector array to be searched, must be monotonic
;               increasing or decreasing
;        X    - the function value(s) whose effective
;               index is sought (scalar or vector)
; OUTPUT:
;        IEFF - the effective index or indices of X in XARR
;               real or double precision, same # of elements as X
; RESTRICTIONS:
;        TABINV will abort if XARR is not monotonic.  (Equality of 
;        neighboring values in XARR is allowed but results may not be
;        unique.)  This requirement may mean that input vectors with padded
;        zeroes could cause routine to abort.
; PROCEDURE:
;        A binary search is used to find the values XARR(I)
;        and XARR(I+1) where XARR(I) LE X LT XARR(I+1).
;        IEFF is then computed using linear interpolation 
;        between I and I+1.
;                IEFF = I + (X-XARR(I)) / (XARR(I+1)-XARR(I))
;        Let N = number of elements in XARR
;                if x < XARR(0) then IEFF is set to 0
;                if x > XARR(N-1) then IEFF is set to N-1
; EXAMPLE:
;        Set all flux values of a spectrum (WAVE vs FLUX) to zero
;        for wavelengths less than 1150 Angstroms.
;           TABINV,WAVE,1150.0,I
;           FLUX( 0:FIX(I) ) = 0.                         
; FUNCTIONS CALLED:
;    ISARRAY
; REVISION HISTORY:
;   Adapted from the IUE RDAF                     January, 1988         
;   More elegant code  W. Landsman                August, 1989
;-
IF N_PARAMS(0) LT 3 THEN BEGIN
     PRINT,STRING(7B),'CALLING SEQUENCE- tabinv,xarr,x,i'
     RETURN
ENDIF
NPOINTS = N_ELEMENTS(XARR) & NPT= NPOINTS - 1
;
; Initialize binary search area and compute number of divisions needed
;
ILEFT = INTARR( N_ELEMENTS(X) ) & IRIGHT = ILEFT
NDIVISIONS = FIX( ALOG10(NPOINTS)/ALOG10(2.0)+1.0 )
;
; Test for monotonicity 
;
I = XARR - SHIFT(XARR,1)
A = WHERE(I GE 0)
IF (!ERR EQ NPT) THEN $ ; Increasing array ?
    IRIGHT = IRIGHT + NPT $
ELSE BEGIN
     A = WHERE(I LE 0)  ; Test for decreasing array
     IF (!ERR EQ NPT) THEN ILEFT = ILEFT + NPT $
     ELSE BEGIN
           PRINT,STRING(7B)+ $ 
           'TABINV: ERROR - First parameter must be a monotonic vector' 
           RETURN 
           END 
END
;
; Perform binary search by dividing search interval in
; half NDIVISIONS times
;
FOR I=1,NDIVISIONS DO BEGIN
    IDIV = (ILEFT+IRIGHT)/2      ;Split interval in half
    XVAL = XARR(IDIV)            ;Find function values at center
    GREATER= FIX(X GT XVAL)      ;Determine which side X is on
    LESS   = FIX(X LE XVAL)  
    ILEFT =  ILEFT*LESS + IDIV*GREATER ;Compute new search area
    IRIGHT = IRIGHT*GREATER + IDIV*LESS
ENDFOR
;
; Interpolate between interval of width = 1
;
XLEFT =  XARR(ILEFT)              ;Value on left side
XRIGHT = XARR(IRIGHT)             ;Value on right side
IEFF = (XRIGHT-X)*ILEFT + (X-XLEFT)*IRIGHT + ILEFT*(XRIGHT EQ XLEFT)
IEFF = IEFF/FLOAT(XRIGHT-XLEFT+(XRIGHT EQ XLEFT)) ;Interpolate
IEFF = IEFF > 0.0 < NPT        ;Do not allow extrapolation beyond ends
IF NOT ISARRAY(X) THEN IEFF = IEFF(0)  ;Make scalar if X was scalar
RETURN
END
