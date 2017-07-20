for m in {2..3}
do
    for n in {0..9999}
    do
        a=$((n + 10000))
        cat "data$m/a$n" "data$m/a$a" >> "rfn$m"
    done
    sed -i "s/$/;/g" "rfn$m"
    java -jar TreeCmp.jar -w 2 -d rf -i "rfn$m" -o "rfn$m.out" -P
    rm "rfn$m"
    b=$((m + 2))
    ./filter "rfn$m.out" "dists$b"
    rm "rfn$m.out"
done
