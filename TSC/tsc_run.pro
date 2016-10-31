pro tsc_run, indir=indir, outdir=outdir, rc=rc, thetac=thetac, time=time, c_s=c_s, omega=omega, renv_min=renv_min;, rstar=rstar, renv_min=renv_min, renv_max=renv_max;, r_inf=r_inf
    ; Paremeters:
    ; indir - The directory contains "grid.plt34.mod"
    ; outdir - The output directory for writing out the density profile
    ; rc - The radius grid in cm.  This should be the coordinates of the center of each cell.
    ; thetac - The thetat grid in radians
    ; time - The time since the collapse began in years.
    ; c_s - The effective sound speed (including non-thermal turbulance in some cases) in cm/s.
    ; omega - The rotational speed of the envelope in s^-1
    ; renv_min - The inner radius of the envelope.
    ;
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

	modeltime = time
	c_s = c_s
	Omega_0 = omega

	if not keyword_set(outdir) then outdir = indir

	; Compile
	; .r ~/programs/misc/TSC/cubesolve.pro
	; .r ~/programs/misc/TSC/loglin_interp2pt.pro
	; .r ~/programs/misc/TSC/tsc.pro
	; Run the calculation
	tsc, modeltime, c_s, Omega_0, rc, thetac, rhoenv, 1e-40, renv_min, indir=indir, outdir=outdir

	; Print the results into file so that I can use python read them in
	openw, lun, outdir+'rhoenv.dat', /get_lun
	for i = 0, n_elements(rhoenv[0,*])-1 do printf, lun, format='('+strtrim(string(n_elements(rhoenv[*,0])),1)+'(g,2x))',rhoenv[*,i]
	free_lun, lun
	close, lun
end
