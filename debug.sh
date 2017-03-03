sh compile.sh -g
gdb -ex run --args ./a.out data/t1mod_d data/t2mod_d sim_t1_t2 0.1
rm a.out
