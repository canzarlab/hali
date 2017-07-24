if [ "$#" -ne 2 ]
then
    echo "usage: <n-leaves> <n-tree-pairs>"
    exit 1
fi
for m in {0..3}
do
    if [ ! -d "data$m" ]
    then
        mkdir "data$m"
    fi
done
./bgen $1 $(($2*2))
rm -f rand.log
for m in {0..1}
do
    b=$((m+2))
    rm -f "dists$b" "rfm$m"
    for ((n=0; n<$2; n++))
    {
        a=$((n+$2))
        for k in 1 5 10 20
        do
            ./phylo2tc -t1 "data$m/a$n" -t2 "data$m/a$a" -k $k >& /dev/null
            mv t1mod "data$b/a$n"
            mv t2mod "data$b/a$a"
            mv sim_t1_t2 "data$b/s_a"$n"_a$a"
            ./solver "data$b/a$n" "data$b/a$a" "data$b/s_a"$n"_a$a" 2>>rand.log >> "dists"$b"k$k"
        done
        cat "data$m/a$n" "data$m/a$a" >> "rfn$m"
    }
    sed -i "s/://g" "rfn$m"
    sed -i "s/$/;/g" "rfn$m"
    for d in mc rc ns tt mp
    do
        java -jar bin/TreeCmp.jar -w 2 -d $d -i "rfn$m" -o $d"n$m.out" -A
        ./filter $d"n$m.out" $d"dists$m" $d
        rm $d"n$m.out"
    done
    rm "rfn$m"
done

