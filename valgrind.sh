sh compile.sh -g -O1
valgrind --tool=memcheck --leak-check=full ./a.out t1mod_d t2mod_d sim_t1_t2 0.1
rm a.out
