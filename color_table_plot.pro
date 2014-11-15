pro color_table_plot, color_table

set_plot, 'ps'
!p.font = 0
!p.thick = 8 & !x.thick = 5 & !y.thick = 5
device, filename = '/Users/yaolun/colortable_fsc_brewer'+strtrim(string(color_table),1)+'.eps', /helvetica, /portrait, /encapsulated, isolatin = 1, font_size = 14, decomposed = 0, /color
x = fltarr(10,255)
y = fltarr(10,255)
for i = 0, 9 do begin
    for j = 0, 254 do begin
        y[i,j] = i+0.5
        x[i,j] = j
    endfor
endfor
plot, x[0,*], y[0,*], /nodata, yrange=[0,10]
colorFile = Filepath(SUBDIRECTORY=['resource','colors'], 'fsc_brewer.tbl')
loadct, color_table,file=colorfile
for k = 0, 254 do begin
    oplot, x[*,k], y[*,k], color = x[0,k]
endfor
device, /close_file, decomposed=1
set_plot,'x'
end
        