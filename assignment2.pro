set_plot, 'ps'
!p.font = 0
device, filename = '~/collision_profile.ps', /helvetica, /portrait, /encapsulated, isolatin = 1, font_size = 16, decomposed = 0, /color
loadct, 13
!p.thick = 4 & !x.thick = 3 & !y.thick = 3
i = indgen(20) & i = i-10
  plot, i, 1/((2*!PI)^2)*1/(1+(i-2)^2), xtitle = '!9w', ytitle = '|E(!9w!3)|!u2!n/T'


device, /close_file, decomposed = 1