#!/bin/bash

for z in {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}; do
    sed -e "s/%%%%z%%%%/${z}/g" < hantush-input.tpl > hantush-input.dat
    echo "${z}"
    ./unconfined ./hantush-input.dat 
done
