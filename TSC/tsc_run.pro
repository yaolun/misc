pro tsc_run, outdir=outdir, grid=grid, time=time, c_s=c_s, omega=omega, rstar=rstar, renv_min=renv_min, renv_max=renv_max
	; Script for executing the tsc calculation
	; Constants setup
	c         = 2.998e10
	AU        = 1.49598e13     ; Astronomical Unit       [cm]
	pc        = 3.08572e18     ; Parsec                  [cm]
	MS        = 1.98892e33     ; Solar mass              [g]
	LS        = 3.8525e33      ; Solar luminosity        [erg/s]
	RS        = 6.96e10        ; Solar radius            [cm]
	G         = 6.67259e-8     ; Gravitational constant  [cm3/g/s^2]
	yr        = 60*60*24*365.0 ; Years in seconds

	; ; Grid Parameters
	; nx        = double(100.0)
	; ny        = double(400.0)
	; nz        = double(50.0)
	nx = double(grid[0])
	ny = double(grid[1])
	nz = double(grid[2])

	; variable setup for tsc.pro
	; modeltime = 32815.0   
	; c_s = 0.505957626614 * 1e5  ; 0.35 km/s from Terebey 1984
	; Omega_0 = 2.62836312414e-13   ; Terebey 1984
	modeltime = time
	c_s = c_s  ; 0.35 km/s from Terebey 1984
	Omega_0 = omega   ; Terebey 1984
	; Model Parameters
	;
	; Manually input the model Parameters
	; rstar     = 5 * RS
	; R_env_max = 1.000000e+04 * AU
	; R_env_min = 8.000000e-01 * AU
	; ; R_cen     = 1.500000e+01 * AU
	; rin       = rstar
	; rout      = R_env_max
	; rcen      = R_cen
	rstar     = rstar
	R_env_max = renv_max
	R_env_min = renv_min
	rin       = rstar
	rout      = R_env_max

	; Make the Coordinates
	;
	ri           = rin * (rout/rin)^(dindgen(nx+1)/nx)
	ri           = [0.0, ri]
	thetai       = !pi*dindgen(ny+1)/ny  ; the angle respect to the rotaional axis
	phii         = !pi*2.0*dindgen(nz+1)/nz

	; Keep the constant cell size in r-direction
	;
	ri_cellsize = ri[1:-1]-ri[0:-2]
	ind = (where(ri_cellsize/AU gt 100.0))[0]      ; The largest cell size is 100 AU
	ri = [ri[0:ind-1], ri[ind]+dindgen(ceil((rout-ri[ind])/100/AU))*100*AU]
	nx = n_elements(ri)-1
	print, nx
	print, n_elements(ri)

	; Assign the coordinates of the center of cell as its coordinates.
	;
	rc           = 0.5*( ri[0:nx-1]     + ri[1:nx] )
	thetac       = 0.5*( thetai[0:ny-1] + thetai[1:ny] )
	phic         = 0.5*( phii[0:nz-1]   + phii[1:nz] )

	indir = '~/programs/misc/TSC/'
	if not keyword_set(outdir) then outdir = indir

	; Compile
	; .r ~/programs/misc/TSC/cubesolve.pro
	; .r ~/programs/misc/TSC/loglin_interp2pt.pro
	; .r ~/programs/misc/TSC/tsc.pro
	; Run the calculation
	tsc, modeltime, c_s, Omega_0, rc, thetac, rhoenv, 1e-40, R_env_min, indir=indir, outdir=outdir

	; Print the results into file so that I can use python read them in
	openw, lun, outdir+'rhoenv.dat', /get_lun
	for i = 0, n_elements(rhoenv[0,*])-1 do printf, lun, format='('+strtrim(string(n_elements(rhoenv[*,0])),1)+'(g,2x))',rhoenv[*,i]
	free_lun, lun
	close, lun
end