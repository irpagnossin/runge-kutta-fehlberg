//------------------------------------------------------------------------------
// Arquivo: DiscontinuityException.cpp (versão 1.0)
// Autor:   Ivan Ramos Pagnossin
// Data:    2006.02.20
//------------------------------------------------------------------------------
#ifndef DISCONTINUITYEXCEPTION_CPP
#define DISCONTINUITYEXCEPTION_CPP

#include "RungeKuttaFehlbergException.cpp"

class DiscontinuityException : public RungeKuttaFehlbergException{
public:
    DiscontinuityException() : RungeKuttaFehlbergException( "It seems there is a discontinuity here." ) {}
};

#endif
//----------- FIM DA DECLARAÇÃO E IMPLEMENTAÇÃO DA CLASSE DISCONTINUITYEXCEPTION
 