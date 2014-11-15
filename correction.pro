pro correction, scaled_data_structure

data = read_ascii('~/co_TRANS_fit.txt', data_start = 3)
data2 = read_ascii('~/co_TRANS_fit_pacs.txt', data_start = 17)

for i = 0, n_elements(data.field01[0,*])-1 do begin
	for j = 0, n_elements(data2.field01[0,*])-1 do begin
		if data2.field01[0,j] eq data.field01[0,i] then begin
			f = data.field01[9,i]/data2.field01[9,j] 
			data2.field01[3,j] = data2.field01[3,j]*f
			data2.field01[4,j] = data2.field01[4,j]*f
			data2.field01[5,j] = data2.field01[5,j]*f
			data2.field01[6,j] = data2.field01[6,j]*f
			data2.field01[9,j] = data2.field01[9,j]*f
			data2.field01[10,j] = data2.field01[10,j]*f
			data.field01[*,i] = data2.field01[*,j]
		endif
	endfor
endfor
scaled_data_structure = data

openw, lun, 'co_TRANS_fit_combined.txt', /get_lun, /append  ;you can change the filename you want
for i = 0, n_elements(data.field01[0,*])-1 do begin 
printf, lun, format='((a8,2X),(F8.4,2X),9(e10.4,2X),(F5.2,2X))',data.field01[*,i]
endfor
free_lun, lun
close, lun

end
