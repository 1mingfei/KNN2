#!/bin/bash
read_E()
{
    grep -i entropy OUTCAR > tmp1
    tail -1 tmp1 > tmp2
    E=`awk '{print $7}' tmp2`
    rm tmp1 tmp2
}
main()
{
    mkdir -p NEB/00 NEB/06
    cp s/CONTCAR NEB/POSCAR0
    cp f/CONTCAR NEB/POSCAR1
    cp s/OUTCAR NEB/00/OUTCAR
    cp s/OUTCAR NEB/OUTCAR0
    cp f/OUTCAR NEB/06/OUTCAR
    cp f/OUTCAR NEB/OUTCAR1
    cp s/KPOINTS s/POTCAR  NEB/
    cd NEB
    ../../../vtstscripts/nebmake.pl POSCAR0 POSCAR1 5
    cd 00
    read_E
    E1=$E
    cd ..
    cd 06
    read_E
    E2=$E
    cd ..
cat << EOF > INCAR
NWRITE = 2

Electronic Relaxation 1:
PREC   = Acc
ISYM   = 2
NELM =  240
NELMIN = 4       

Ionic Relaxation:
NSW    = 10000
NBLOCK = 1
KBLOCK = 1
IBRION = 3
POTIM  = 0
IOPT = 3
ISIF   = 2

DOS Related Values:
ISMEAR = 1
SIGMA  = 0.4

Electronic Relaxation 2:
IALGO  = 48
LREAL  = AUTO
ENCUT  = 450
ENAUG  = 600.0

LCLIMB = .TRUE.
ICHAIN = 0
IMAGES = 5
SPRING = -5

LWAVE = .FALSE.
EFIRST =  lala
ELAST =   lblb

EDIFFG = -0.05
EOF

    sed  -i s/lala/${E1}/g INCAR
    sed  -i s/lblb/${E2}/g INCAR
    cd ..
}
time main;
