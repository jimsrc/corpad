//#include <stacktrace/call_stack_gcc.cpp>
//#include <stacktrace/stack_exception.hpp>
#include "nr3.h"


int main(){
    VecInt a(4, 9);
    MatDoub b(3, 4, 5.6);

    try{
        printf(" b: %f ", b[10][1]);
        //printf(" a: %d \n", a[11]);
    }
    catch(NRerror s) { 
        printf(" ---> @main: LINE: %d \n", __LINE__);
        NRcatch(s); 
    }
}
