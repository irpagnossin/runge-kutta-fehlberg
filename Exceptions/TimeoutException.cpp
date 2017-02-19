//------------------------------------------------------------------------------
// Arquivo: TimeoutException.cpp (vers�o 1.0)
// Autor:   Ivan Ramos Pagnossin
// Data:    2006.02.20
//------------------------------------------------------------------------------
#ifndef TIMEOUTEXCEPTION_CPP
#define TIMEOUTEXCEPTION_CPP

#include "RungeKuttaFehlbergException.cpp"

class TimeoutException : public RungeKuttaFehlbergException{
public:
    TimeoutException() : RungeKuttaFehlbergException( "I spent too much time... it seems it will not converge.\n" ) {}
};

#endif
//----------------- FIM DA DECLARA��O E IMPLEMENTA��O DA CLASSE TIMEOUTEXCEPTION
 