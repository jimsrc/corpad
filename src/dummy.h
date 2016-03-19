#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
//#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>


//using namespace std;


template <class T>
class dummy{
    private:
        int a, b;
        T *v;
    public:
        double aa, bb;
        void func();
        dummy();
        explicit dummy(int);
};

typedef dummy<double> dummy_doub;
