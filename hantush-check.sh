#!/bin/bash

for z in {0.025,0.125,0.225,0.325,0.425,0.525,0.625,0.725,0.825,0.925,0.975}; do
    sed -e "s/%%%%z%%%%/${z}/g" < hantush-input.tpl > hantush-input.dat
    echo "${z}"
    ./unconfined ./hantush-input.dat 
done
