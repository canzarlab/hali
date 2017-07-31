if [ "$#" -ne 2 ]
then
    echo "usage: <n-leaves> <n-tree-pairs>"
    exit 1
fi
for m in 0 1 2 3
do
    if [ ! -d "data$m" ]
    then
        mkdir "data$m"
    fi
done
make >& /dev/null
g++ generator.cpp newick.cpp -o bgen
g++ filter.cpp -o filter
( cd genPhylo2 && make >& /dev/null && mv phylo2tc .. )
./bgen $1 $(($2*2))
rm -f rand.log dists*
for m in 0 1
do
    b=$((m+2))
    rm -f "rfm$m"
    for ((n=0; n<$2; n++))
    {
        a=$((n+$2))
        for k in 1 2 3 4 5
        do
            for mt in j s
            do
                ./phylo2tc -t1 "data$m/a$n" -t2 "data$m/a$a" -k $k -d $mt >& /dev/null
                mv t1mod "data$b/a$n"
                mv t2mod "data$b/a$a"
                mv sim_t1_t2 "data$b/s_a"$n"_a"$a"_k"$k"_d"$mt
                for c in 2 1 0
                do
                    ./solver "data$b/a$n" "data$b/a$a" "data$b/s_a"$n"_a"$a"_k"$k"_d"$mt $c $mt 2>>/dev/null >> "dists"$b"k"$k"c"$c"d"$mt &
                    if [ $k != 1 ]; then
                        break
                    fi
                done
                if [ $k != 1 ]; then
                    break
                fi
            done
        done
        cat "data$m/a$n" "data$m/a$a" >> "rfn$m"
        r=$((n%8));
        if [ $r == 7 ]; then
            wait
        fi
    }
    sed -i "s/://g" "rfn$m"
    sed -i "s/$/;/g" "rfn$m"
    for d in mc rc ns tt mp
    do
        java -jar TreeCmp/bin/TreeCmp.jar -w 2 -d $d -i "rfn$m" -o $d"n$m.out" -A > /dev/null
        ./filter $d"n$m.out" $d"dists$m" $d
        rm $d"n$m.out"
    done
    rm "rfn$m"
done

