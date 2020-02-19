for i in `seq 0 62`
do
    cd config$i;
    for j in $(ls -d */)
    do
        cd ${j}
        sbatch submit.stampede2
        cd ..
    done
    cd ..
done
