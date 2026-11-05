#!/bin/bash

n_q=5
n_atoms=3
filename='qpoints.yaml'

n_modes=$((n_atoms*3))


for iq in `seq 1 $n_q ` ; do
        for im in `seq 1 $n_modes ` ; do
                nn=$(( 9+3*$im-1+($im-1)*(4*$n_atoms)+($iq-1)*(3+$n_modes*(3+4*$n_atoms)) ))
                sed -n "$nn"p $filename > temp
                awk '{printf("% 18.12f  \n",$2)}' < temp  >> freq.dat
        done
done
