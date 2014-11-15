pro multi_plots
set_plot, 'ps'
device, filename = '~/sample.eps', /helvetica, /portrait, /encapsulated, font_size = 12, isolatin = 1, decomposed = 0, /color
loadct, 12
!p.thick = 1 & !x.thick = 3 & !y.thick = 3
readcol, '~/bhr71/data/pixel_spectrum/extended_correction/SLWC3.txt', format='D,D', wl, flux
plot, wl, flux, xtitle = '!9m!3m', ytitle = 'Jy'
device, /close_file, decomposed = 1
end