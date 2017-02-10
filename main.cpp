#include "newick.h"
#include <iostream>

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        cout << "usage: " << argv[0] << " <filename>" << endl;
        return EXIT_FAILURE;
    }
    
    newick_node* root = load_tree(argv[1]);
    print_tree(root);
    delete root;
}
