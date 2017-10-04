make
valgrind --tool=memcheck --leak-check=full ./solver dag1 map1 dag2 map2 align 2 j 1
rm solver
