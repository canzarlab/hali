function foreach()
{
    for k in 1 2 3 4 5
    do
        for mt in j s
        do
            $1
            for c in 2 1 0
            do
                $2
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

function null()
{
    :
}

function sync()
{
    cat "dists"$b"k"$k"c"$c"d"$mt"r"$l >> "dists"$b"k"$k"c"$c"d"$mt 2>/dev/null
}

function clean()
{
    rm -f "dists"$b"k"$k"c"$c"d"$mt"r"$l
}

function prep()
{
    ./phylo2tc -t1 "data$m/a$n" -t2 "data$m/a$a" -k $k -d $mt >& /dev/null
    mv t1mod "data$b/a$n"
    mv t2mod "data$b/a$a"
    mv sim_t1_t2 "data$b/s_a"$n"_a"$a"_k"$k"_d"$mt
}

function solve()
{
    ./solver "data$b/a$n" "data$b/a$a" "data$b/s_a"$n"_a"$a"_k"$k"_d"$mt $c $mt 2>/dev/null > "dists"$b"k"$k"c"$c"d"$mt"r"$r &
}

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

j=32

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
        r=$((n%j))
        v=$((j-1))
        w=$(($2-1))
        foreach prep solve
        cat "data$m/a$n" "data$m/a$a" >> "rfn$m"
        if [ $r == $v ] || [ $n == $w ]
        then
            wait
            for ((l=0; l<$j; l++))
            {
                foreach null sync
            }
        fi
    }
    for ((l=0; l<$j; l++))
    {
        foreach null clean
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

