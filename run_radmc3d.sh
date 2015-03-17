#!/bin/bash
check_process() {
      #echo "$ts: checking $1"
        [ "$1" = "" ]  && return 0
	  [ `pgrep -n $1` ] &&  return 1 || return 0
}

# main=/home/bettyjo/yaolun/radmc_simulation/bhr71/
main=/home/bettyjo/yaolun/radmc_simulation/bhr71/rough_grid/
# testing the initial guess from Robitaille's online SED fitting tool
# model=(3017936 3000358 3006440 3004050 3016836 3010725 3014664 3004356 3016082 3006416 3017928 3014719)

if [ "$1" = "cavity_combo" ]
then
    model=("cavity_combo_2/cav_mod_1" "cavity_combo_2/cav_mod_2" "cavity_combo_2/cav_mod_3" "cavity_combo_2/cav_mod_4" "cavity_combo_2/cav_mod_5" "cavity_combo_2/cav_mod_6" "cavity_combo_2/cav_mod_7" "cavity_combo_2/cav_mod_8" "cavity_combo_2/cav_mod_9" "cavity_combo_2/cav_mod_10" "cavity_combo_2/cav_mod_11")
elif [ "$1" = "best_model" ]
then
    model=("best_model")
elif [ "$1" = "mstar" ]
then
    model=("mstar/mstar_1.0e-1" "mstar/mstar_5.0e-1" "mstar/mstar_8.0e-1" "mstar/mstar_1.0e+0" "mstar/mstar_1.3e+0" "mstar/mstar_1.5e+0")
elif [ "$1" = "rstar" ]
then
    model=("rstar/rstar_1.0e+0" "rstar/rstar_3.0e+0" "rstar/rstar_5.0e+0" "rstar/rstar_8.0e+0" "rstar/rstar_1.0e+1")
elif [ "$1" = "tstar" ]
then
    model=("tstar/tstar_3500K" "tstar/tstar_4000K" "tstar/tstar_4500K" "tstar/tstar_5000K" "tstar/tstar_5500K")
elif [ "$1" = "menvdot" ]
then
    model=("menvdot/menvdot_1.0e-5" "menvdot/menvdot_7.0e-6" "menvdot/menvdot_5.0e-6" "menvdot/menvdot_2.0e-6" "menvdot/menvdot_8.0e-7" "menvdot/menvdot_5.0e-7")
elif [ "$1" = "renv_max" ]
then
    model=("renv_max/renv_max_1.5e+4" "renv_max/renv_max_1.0e+4" "renv_max/renv_max_9.0e+3" "renv_max/renv_max_8.0e+3" "renv_max/renv_max_7.0e+3" "renv_max/renv_max_6.0e+3" "renv_max/renv_max_5.0e+3")
elif [ "$1" = "r_min" ]
then
    model=("r_min/r_min_0.3au" "r_min/r_min_0.5au" "r_min/r_min_1.0au" "r_min/r_min_5.0au" "r_min/r_min_8.0au" "r_min/r_min_10.au")
elif [ "$1" = "theta" ]
then
    # model=("theta/theta_5" "theta/theta_10" "theta/theta_15" "theta/theta_20" "theta/theta_25")
    model=("theta/theta_5" "theta/theta_7" "theta/theta_9" "theta/theta_15" "theta/theta_20" "theta/theta_25") # "theta/theta_15" "theta/theta_20" "theta/theta_25")
elif [ "$1" = "theta_temp" ]
then
    model=("theta_5/4000K" "theta_5/4500K" "theta_5/5000K" "theta_5/5500K")
elif [ "$1" = "mdisk" ]
then
    model=("mdisk/mdisk_1.0e-2" "mdisk/mdisk_5.0e-2" "mdisk/mdisk_1.0e-1" "mdisk/mdisk_2.5e-1" "mdisk/mdisk_5.0e-1")
elif [ "$1" = "rdisk_max" ]
then
    model=("rdisk_max/rdisk_max_30" "rdisk_max/rdisk_max_50" "rdisk_max/rdisk_max_100" "rdisk_max/rdisk_max_150" "rdisk_max/rdisk_max_200")
elif [ "$1" = "beta" ]
then
    model=("beta/beta_1.00" "beta/beta_1.05" "beta/beta_1.10" "beta/beta_1.15" "beta/beta_1.20" "beta/beta_1.50")
elif [ "$1" = "h100" ]
then
    model=("h100/h100_5au" "h100/h100_7au" "h100/h100_10au" "h100/h100_15au" "h100/h100_20au")
elif [ "$1" = "wall" ]
then
    model=("wall/wall_100au" "wall/wall_150au" "wall/wall_200au" "wall/wall_250au")
elif [ "$1" = "rho_wall" ]
then
    model=("wall/rho_1e-17" "wall/rho_1e-18" "wall/rho_1e-19" "wall/rho_1e-20")
elif [ "$1" = "no_cavity" ]
then
    model=("no_cavity/menvdot_1e-3" "no_cavity/menvdot_1e-4" "no_cavity/menvdot_1e-5" "no_cavity/menvdot_1e-6" "no_cavity/menvdot_1e-7")
fi

# model=("tstar_v2_4500K_test")
num=`ps -Al | grep radmc3d | wc -l`
echo "Check system availibility..."
echo $num
while [ "$num" -gt 35 ]; do
    sleep 60
    num=`ps -Al | grep radmc3d | wc -l`
done
date
echo "Calculating the dust temperature..."
for mod in "${model[@]}"; do
    cd "$main""$mod"
    #cp /home/bettyjo/yaolun/radmc_simulation/bhr71/model_tstar_2500K/wavelength_micron.inp .
    #cp "$main"rho_cav_test/aperture/model_4500K_1e-6/radmc3d.inp .
    cp /home/bettyjo/yaolun/radmc_simulation/dustkappa_oh5_extended.inp .
    radmc3d mctherm > log_mctherm.txt &
done
sleep 5

while [ `pgrep -nx "radmc3d"` ]; do
# while [ $num -gt 35 ]; do
#    date
#    echo "Calculating the dust temperature..."
    sleep 30
done
date
echo "Initiating the ray-tracing calculation..."

# model=("4000K_6e-6" "4000K_4e-6" "4000K_2e-6" "4000K_1e-6" "4000K_8e-7" "4000K_6e-7" "4000K_4e-7" "4000K_2e-7" "4000K_1e-7" "4500K_6e-6" "4500K_4e-6" "4500K_2e-6" "4500K_1e-6" "4500K_8e-7" "4500K_6e-7" "4500K_4e-7" "4500K_2e-7" "4500K_1e-7")
# model=("5000K_6e-6" "5000K_4e-6" "5000K_2e-6" "5000K_1e-6" "5000K_8e-7" "5000K_6e-7" "5000K_4e-7" "5000K_2e-7" "5000K_1e-7")
for mod in "${model[@]}"; do
    cd "$main""$mod"
#    mv aper_info.inp aperture_info.inp
#    cp /home/bettyjo/yaolun/radmc_simulation/bhr71/model_tstar_3000K/wavelength_micron.inp .
    radmc3d spectrum incl 82 loadlambda useapert dpc 178 > log_sed.txt &
#    radmc3d image incl 41 lambdarange 2. 600. nlam 100 > log_image_incl_41.txt & 
#    radmc3d image incl 82 lambdarange 2. 600. nlam 100 secondorder > log_image_incl_82.txt &
    radmc3d image incl 82 loadlambda npix 300 > log_image_incl_82.txt &
#    radmc3d image incl 82 lambdarange 200. 600. nlam 20 > log_image.txt &
#    radmc3d tausurf 1.0 incl 82 lambdarange 2. 669. nlam 100 > log_tau.txt &
done
sleep 5
# while [ `pgrep -nx "radmc3d"` ]; do
#     sleep 30
# done
date
echo "Ray-tracing calculation finished."
