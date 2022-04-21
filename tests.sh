#!/bin/bash

a=1
m=7
n=5

mpicc forest_corr_par.c -o par
gcc forest_corr_seq.c -o seq

echo Compilazione terminata, inizio...

while [ $a -lt 15 ];
do
    mpirun --allow-run-as-root -np 4 ./par $(($m * $a)) $(($n * $a))
    ./seq $(($m * $a)) $(($n * $a))

    TEXT=$('./check')

    rm correttezza_par.txt
    rm correttezza_seq.txt
    
    echo "$TEXT Righe:$(($m * $a)) Colonne: $(($n * $a))"

    ((a++))
done

echo Fine...
