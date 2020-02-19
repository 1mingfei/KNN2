
getEBarrier() {
    Estart=`sed '1,1!d' neb.dat | awk '{print $3}'`
    Esaddle=`awk 'BEGIN{a=   0}{if ($3>0+a) a=$3} END{print a}' neb.dat`
    Eend=`sed '7,7!d' neb.dat | awk '{print $3}'`
    Eforward=`echo "($Esaddle)-($Estart)"| bc`
    Ebackward=`echo "($Esaddle)-($Eend)"| bc`
    Estartshow=`grep EFIRST INCAR | awk '{print $3}'`
    Eendshow=`grep ELAST INCAR | awk '{print $3}'`
    Dist=`dist.pl POSCAR0 POSCAR1`
}
for j in `seq 12 23`
do
    cd config$j/
    max=`ls | sed 's/NEB_//' | sort -n | tail -1`
    for k in `seq 0 "${max}"`
    do
        cd NEB_$k
        getEBarrier
        cd ../..
        if [ $(echo "$Dist < 3.5" |bc -l) -eq 1 ]
        then
            echo $i $j $k $Eforward $Ebackward $Estartshow $Eendshow >> ./E_neb_all_start_end.dat
        fi
        cd config$j/
    done
    cd ..
done
