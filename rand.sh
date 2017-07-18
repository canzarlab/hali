for m in {0..1}
do
    for n in {0..9999}
    do
        a=$((n + 10000))
        ./solver "data$m/a$n" "data$m/a$a" "data$m/s_a"$n"_a$a"
    done
done
