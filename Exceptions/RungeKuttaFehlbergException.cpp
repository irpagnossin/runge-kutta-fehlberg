//------------------------------------------------------------------------------
// Arquivo: RungeKuttaFehlbergException.cpp (vers�o 1.0)
// Autor:   Ivan Ramos Pagnossin
// Data:    2006.02.06
//------------------------------------------------------------------------------
#ifndef RUNGEKUTTAFEHLBERGEXCEPTION_CPP
#define RUNGEKUTTAFEHLBERGEXCEPTION_CPP

class RungeKuttaFehlbergException{
public:
    RungeKuttaFehlbergException( char * message ) : message( message ){}
    const char * getMessage( void ) const{ return( this->message ); }
private:
    const char * message;
};

#endif
//------ FIM DA DECLARA��O E IMPLEMENTA��O DA CLASSE RUNGEKUTTAFEHLBERGEXCEPTION
