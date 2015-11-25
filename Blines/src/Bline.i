%module Bline
%include "std_vector.i"
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubVector) vector<double>;
   %template(DoubMatrix) vector<vector<double> >;
}

%include "Bline.h"
 
%{
#include "Bline.h"
%}

