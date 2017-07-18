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
    done
done
exit 0
