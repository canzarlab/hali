function foreach()
{
    for k in 1 2 3 4 5
    do
        for mt in j s
        do
            for c in 2 1 0
            do
                $1
                if [ $k != 1 ]; then
                    break
                fi
            done
            if [ $k != 1 ]; then
                break
            fi
        done
    done
}

function sync()
{
    cat "dists"$m"k"$k"c"$c"d"$mt"r"$l >> "dists"$m"k"$k"c"$c"d"$mt 2>/dev/null
}

function clean()
{
    rm -f "dists"$m"k"$k"c"$c"d"$mt"r"$l
}

function solve()
{
    ./solver "data$m/a$n" "data$m/a$a" $c $mt $k 2>/dev/null > "dists"$m"k"$k"c"$c"d"$mt"r"$r &
}

if [ "$#" -ne 2 ]
then
    echo "usage: <n-leaves> <n-tree-pairs>"
    exit 1
fi

for m in 0 1
do
    if [ ! -d "data$m" ]
    then
        mkdir "data$m"
    fi
done

j=32

make >& /dev/null

./bgen $1 $(($2*2))
rm -f rand.log dists*
for m in 0 1
do
    rm -f "rfm$m"
    for ((n=0; n<$2; n++))
    {
        a=$((n+$2))
        r=$((n%j))
        v=$((j-1))
        w=$(($2-1))
        foreach solve
        cat "data$m/a$n" "data$m/a$a" >> "rfn$m"
        if [ $r == $v ] || [ $n == $w ]
        then
            wait
            for ((l=0; l<$j; l++))
            {
                foreach sync
            }
        fi
    }
    for ((l=0; l<$j; l++))
    {
        foreach clean
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

