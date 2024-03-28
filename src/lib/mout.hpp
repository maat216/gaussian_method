#include <fstream>
#include <ctime>
#include <iostream>
#include <cstring>

using namespace std;

class __mout_class
{
private:
    ofstream fout;

public:
    __mout_class()
    {
        const int MAXLEN = 80;
        char str[MAXLEN];
        time_t t = time(0);
        strcat(str, "output_");

        strftime(str, MAXLEN, "%d_%m_%Y_%H_%M_%S", localtime(&t));

        strcat(str, ".txt");

        this->fout = ofstream(str);
    }

    template <typename T>
    __mout_class &operator<<(const T &data)
    {
        cout << data;
        this->fout << data;
        return *this;
    }

    __mout_class &operator<<(std::ostream &(*manip)(std::ostream &))
    {
        manip(cout);
        manip(this->fout);
        return *this;
    }
};

__mout_class mout;