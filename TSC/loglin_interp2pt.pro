FUNCTION loglin_interp2pt,x1,y1,x2,y2,x

m=((y1-y2)/(x1-x2))
b=(y1-(m*x1))
y=((m*x)+b)

x1l=alog10(x1)
y1l=alog10(y1)
x2l=alog10(x2)
y2l=alog10(y2)
xl=alog10(x)

ml=((y1l-y2l)/(x1l-x2l))
bl=(y1l-(ml*x1l))
yl=((ml*xl)+bl)

yl=10.^(yl)

;print,"Linear: ",y
;print,"Log-Linear: ",yl

return,yl

END
