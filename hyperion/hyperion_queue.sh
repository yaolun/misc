#!/bin/bash

while [ `pgrep -nx "hyperion_sph_mp"` ]; do
    sleep 30
done
date
echo "initiate dust RT (python script)"

time python main_hyperion.py 