make
valgrind --tool=memcheck --leak-check=full ./solver data/t1mod_d data/t2mod_d data/sim_t1_t2
rm solver
