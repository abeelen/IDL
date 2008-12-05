PRO plot_repart, data, X,Y, $
                 XTITLE=xtitle,YTITLE=ytitle, XRANGE=xrange,$
                 OVERPLOT=OVERPLOT, LINESTYLE=linestyle,COLOR=color, $
                 REVERSED=REVERSED,INVERSED=INVERSED,_EXTRA = extra, $
                 NOPLOT=NOPLOT,SORTED=SORTED,MEDIAN=MEDIAN, $
                 NOT_NORMALIZED = not_normalized, _EXTRA = _extra
;+
; NAME:       plot_repart
; 
;
;
; PURPOSE:    plot the repartition function of a dataset
;
;
;
; CATEGORY:   plotting
;
;
;
; CALLING SEQUENCE:
;             plot_repart,data,[X,Y],$
;                         [XTITLE=],[YTITLE=],[XRANGE=],$
;                         [OVERPLOT=],[LINESTYLE=],[COLOR=],$
;                         [REVERSED=],[INVERSED=],[NOPLOT=],$
;                         [SORTED=],[MEDIAN=],[NOT_NORMALIZED=]
;
;
; INPUTS:
;             data : 1D array containing the dataset
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;             REVERSED       : reverse the X axis
;             INVERSED       : plot the decreasing function (1-X)
;             NOPLOT         : do not plot (just output to X and Y array)
;             MEDIAN         : overplot a median line
;             NOT_NORMALIZED : do not normalized the partition
;                              function to 1
;             
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;            X, Y : the repartition function coordinates
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;w
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;           2005 - AB - created
;-


max_data = MAX(data,/NAN)
min_data = MIN(data,/NAN)

IF NOT KEYWORD_SET(XRANGE) THEN BEGIN
IF KEYWORD_SET(REVERSED) THEN $
xrange = [max_data+1.0, min_data-1.0] $
ELSE $
xrange = [min_data-1.0, max_data+1.0]
ENDIF


IF KEYWORD_SET(not_normalized) THEN $
 YRANGE=[0,N_ELEMENTS(data)]

IF  NOT KEYWORD_SET(OVERPLOT) AND NOT KEYWORD_SET(NOPLOT) THEN $
PLOT,xrange,[0,1],XRANGE=xrange,yrange=YRANGE,/XSTYLE,/YSTYLE, $
  XTITLE=xtitle, YTITLE=ytitle, COLOR=color,$
      /nodata, _EXTRA=_extra


sorted = SORT(data)

IF KEYWORD_SET(REVERSED) THEN $
sorted = REVERSE(sorted)

nelements = N_ELEMENTS(WHERE(FINITE(data)) NE 0)

repartition=(FINDGEN(N_ELEMENTS(data))+1)/nelements

IF KEYWORD_SET(not_normalized) THEN $
  repartition = repartition*nelements


IF KEYWORD_SET(REVERSED) THEN $
X = [max_data,data(sorted(0)),data(sorted(0))] $
ELSE $
X = [min_data,data(sorted(0)),data(sorted(0))] 

Y = [0,0,repartition(0)]





FOR index=1L,N_ELEMENTS(data(sorted))-1 DO BEGIN
    X=[X,data(sorted(index)),data(sorted(index))]
ENDFOR
FOR index=1L,N_ELEMENTS(data(sorted))-1 DO BEGIN
    Y = [Y,repartition(index-1),repartition(index)]
ENDFOR

IF KEYWORD_SET(REVERSED) THEN $
X=[X,MIN(!X.CRANGE)] $
ELSE $
X=[X,MAX(!X.CRANGE)]


IF KEYWORD_SET(NOT_NORMALIZED) THEN $
  Y=[Y,nelements] $
ELSE $
  Y=[Y,1]



IF KEYWORD_SET(INVERSED) THEN $ 
Y = 1 - Y 


IF NOT KEYWORD_SET(NOPLOT) THEN $
OPLOT,X,Y,color=color,LINESTYLE=linestyle

IF NOT KEYWORD_SET(NOPLOT) AND KEYWORD_SET(MEDIAN) THEN $
OPLOT,!X.CRANGE,[0.5,0.5],linestyle=2

END
