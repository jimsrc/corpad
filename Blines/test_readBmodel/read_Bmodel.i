%module read_Bmodel
%include "std_vector.i"
namespace std {
   %template(IntVector) vector<int>;            // en parantesis, el nombre de la funcion en python
   %template(DoubleVector) vector<double>;      // <bis>
}

%include "read_Bmodel.h"
 
%{
#include "read_Bmodel.h"
%}

// COMENTARIOS:
// (*) estas librerias funcionan pero para usarlas, yo tendria
// q saber como traducir los punteros de c++ hacia python. En vez 
// de eso, uso "std.vector.i", xq con vector<double> SI se 
// traducirlos a python.
