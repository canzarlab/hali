make
valgrind --tool=memcheck --leak-check=full ./solver t1mod_d t2mod_d sim_t1_t2
rm solver
