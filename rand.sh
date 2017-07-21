for m in {0..3}
do
    if [ ! -d "data$m" ]
    then
        mkdir "data$m"
    fi
done
./bgen 50 20000
for m in {0..1}
do
    for n in {0..9999}
    do
        a=$((n + 10000))
        ./phylo2tc -t1 "data$m/a$n" -t2 "data$m/a$a" -k 1 >& /dev/null
        b=$((m+2))
        mv t1mod "data$b/a$n"
        mv t2mod "data$b/a$a"
        mv sim_t1_t2 "data$b/s_a"$n"_a$a"
        rm -f "dists$b"
        ./solver "data$m/a$n" "data$m/a$a" "data$m/s_a"$n"_a$a" 2 >> rand.log >> "dists$b"
        rm -f "rfm$m"
        cat "data$m/a$n" "data$m/a$a" >> "rfn$m"
    done
    sed -i "s/$/;/g" "rfn$m"
    java -jar bin/TreeCmp.jar -w 2 -d rf -i "rfn$m" -o "rfn$m.out" -P
    rm "rfn$m"
    ./filter "rfn$m.out" "dists$m"
    rm "rfn$m.out"
done

