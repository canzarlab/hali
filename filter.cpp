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
    if (argc != 4)
    {
        cout << "usage: " << argv[0] << " <input> <output> <metric>\n";
        return EXIT_FAILURE;
    }

    double mult = (argv[3] == string("rc")) ? 2 : 1;
    ifstream ifs(argv[1]);
    ofstream ofs(argv[2]);
    string in;
    getline(ifs, in);
    while (input(ifs, in))
        ofs << stod(in) * mult << ' ';
}
