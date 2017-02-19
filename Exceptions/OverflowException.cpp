//------------------------------------------------------------------------------
// Arquivo: OverflowException.cpp (vers�o 1.0)
// Autor:   Ivan Ramos Pagnossin
// Data:    2006.12.06
//------------------------------------------------------------------------------
#ifndef OVERFLOWEXCEPTION_CPP
#define OVERFLOWEXCEPTION_CPP

#include "RungeKuttaFehlbergException.cpp"

class OverflowException : public RungeKuttaFehlbergException{
public:
    OverflowException() : RungeKuttaFehlbergException( "OVERFLOW" ) {}
};

#endif
//--------------- FIM DA DECLARA��O E IMPLEMENTA��O DA CLASSE OVERFLOWTEXCEPTION
 