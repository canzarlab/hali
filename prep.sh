for m in {2..3}
do

    for n in {0..9999}
    do

        a=$((n + 10000))
        cat "data$m/a$n" "data$m/a$a" >> "rfn$m"
    done
done
