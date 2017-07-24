#include <fstream>
#include <iostream>
using namespace std;

bool input(ifstream& ifs, string& in)
{
    for (int i = 0; i < 4; ++i)
        if (!(ifs >> in))
            return false;
    return true;
}

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        cout << "usage: " << argv[0] << " <input> <output>\n";
        return EXIT_FAILURE;
    }
    
    ifstream ifs(argv[1]);
    ofstream ofs(argv[2]);
    string in;
    getline(ifs, in);
    while (input(ifs, in))
        ofs << stod(in) * 2 << ' ';
}
