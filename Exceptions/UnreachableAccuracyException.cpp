//------------------------------------------------------------------------------
// Arquivo: UnreachableAccuracyException.cpp (vers�o 1.0)
// Autor:   Ivan Ramos Pagnossin
// Data:    2006.02.20
//------------------------------------------------------------------------------
#ifndef UNREACHABLEACCURACYEXCEPTION_CPP
#define UNREACHABLEACCURACYEXCEPTION_CPP

#include "RungeKuttaFehlbergException.cpp"

class UnreachableAccuracyException : public RungeKuttaFehlbergException{
public:
    UnreachableAccuracyException() : RungeKuttaFehlbergException( "It wasn't possible to achieve required accuracy." ) {}
};

#endif
//----- FIM DA DECLARA��O E IMPLEMENTA��O DA CLASSE UNREACHABLEACCURACYEXCEPTION
 