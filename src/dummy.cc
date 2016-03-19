#include "dummy.h"

template <class T>
void dummy<T>::func(){
    printf(" ---> my func!!\n ");
}

template <class T>
dummy<T>::dummy(int _a){
    a = _a;
}

template <class T>
dummy<T>::dummy(){
    a = 999;
}


//typedef dummy<double> dummy_doub;
template class dummy<double>;
