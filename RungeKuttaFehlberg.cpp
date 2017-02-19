//------------------------------------------------------------------------------
// Arquivo: RungeKuttaFehlberg.cpp (implementação da classe RungeKuttaFehlberg).
// autor:   Ivan Ramos Pagnossin (irpagnossin@hotmail.com).
// data:    agosto/setembro de 2005.
//------------------------------------------------------------------------------
/*
    - Troquei "double" para "long double" e deu problema! O que houve?
    - É possível restringir ainda mais o erro? Como definí-lo como um erro relativo?
    - Alterar para Cash-Karp (numerical recipes, p. 711)

    **INCRIVEL: foi só trocar de ângstrom para nanometros que passou a resolver muito mais facilmente!
 */
/* TODO
    - Escrever log no arquivo e analisar a avaliação de dx: há algo errado aqui.
    - Por que o programa está sempre caindo em dx_minimo? (mensagem "unable to reach accuracy")
 */
#include <time>
#include <math>
#include <iostream>
#include <exception>
#include <string>
#include <fstream>

#include "DiscontinuityException.cpp"
#include "TimeoutException.cpp"
#include "OverflowException.cpp"
#include "UnreachableAccuracyException.cpp"
#include "RungeKuttaFehlberg.h"

#define EPS                 (1e-30)
#define INFINITY            (1e+30)

#define DEFAULT_TIMEOUT     (5.0)
#define DEFAULT_ERROR       (1e-3)

#define SHRINK_EXPONENT     (0.25)
#define STRETCH_EXPONENT    (0.20)
#define MAX_SHRINK          (0.10)
#define MAX_STRETCH         (4.00)

#define SAFETY              (0.9)

using std::cout;
using std::endl;
//------------------------------------------------------------------------------
// Construtores.
//------------------------------------------------------------------------------
/*RungeKuttaFehlberg::RungeKuttaFehlberg( const unsigned int order, double (*F)(double, double *) ){
    this->RungeKuttaFehlberg( order, F, false );
};*/

RungeKuttaFehlberg::RungeKuttaFehlberg( const unsigned int order, double (*F)(double, double *), const bool isLinear )
    : order( order ), isLinear( isLinear ) {

    this->F = F;
    this->setError( DEFAULT_ERROR );
    this->setTimeout( DEFAULT_TIMEOUT );
};
//------------------------------------------------------------------------------
// Destrutor.
//------------------------------------------------------------------------------
RungeKuttaFehlberg::~RungeKuttaFehlberg(){
};
//------------------------------------------------------------------------------
// Informa o valor de <timeout>.
//------------------------------------------------------------------------------
double RungeKuttaFehlberg::getTimeout( void ) const{
    return( this->timeout );
}
//------------------------------------------------------------------------------
// Redefine o valor de <timeout>.
//------------------------------------------------------------------------------
void RungeKuttaFehlberg::setTimeout( double timeout ){
    if( timeout > EPS && timeout < INFINITY ){
        this->timeout    = timeout;
        this->hasTimeout = true;
    }
    else{
        this->timeout    = 0;
        this->hasTimeout = false;
    }
}
//------------------------------------------------------------------------------
// Informa o valor do erro aceitável.
//------------------------------------------------------------------------------
double RungeKuttaFehlberg::getError( void ) const{
    return( this->error );
}
//------------------------------------------------------------------------------
// Redefine o valor do erro aceitável [deve pertencer ao intervalo (EPS,INFINITY)].
//------------------------------------------------------------------------------
void RungeKuttaFehlberg::setError( double error ){
    if( error > EPS && error < INFINITY ){
        this->error = error;
    }
    else{
        this->error = DEFAULT_ERROR;
    }
}
//------------------------------------------------------------------------------
// Resolve a EDO.
// O método populará o vetor **y até que termine ou uma exceção ocorra.
//------------------------------------------------------------------------------
void RungeKuttaFehlberg::solve( unsigned int nPtos, const double *x, double **y, double *yI, char * filename ){

    std::ofstream fileOut ( filename );


    /*
     * Define e aloca memória para os ponteiros dinâmicos necessários. Lança std::bad_alloc em caso de erro.
     */
    double  *last_y,         // Última solução VÁLIDA calculada.
            *next_y,         // Solução (válida ou não) em cálculo.
            *error_y,        // Erro associado a *next_y.
            **f;             // Fatores auxiliares.

    last_y  = new double[this->order];
    next_y  = new double[this->order];
    error_y = new double[this->order];

    f = new double * [6];
    for( register int i = 0; i <= 5; i++ ) f[i] = new double[this->order];

    for( int i = 0; i <= 5; i++ )
        for( int j = 0; j <= this->order-1; j++ )
            f[i][j] = 0;


    if( nPtos < 2 ) nPtos = 2;

    double dx = (x[nPtos-1] - x[0]) / (nPtos - 1);      // Estimativa inicial para o passo.

    const bool crescente = ( dx > 0 );                  // Identifica se x[] está crescente ou decrescentemente ordenado.

    double  last_x    = 0.0,
            next_x    = 0.0,
            last_dx   = 0.0,
            soma      = 0.0,
            dx_maximo = dx,
            factor = 1,
            exponent  = SHRINK_EXPONENT;

    register unsigned int   i_ERROR = 0,
                            nCalc   = 0,
                            pto     = 1;

    bool atribui = false;                           // Identifica se aproximação calculada corresponde a um elemento de x[].

	const double a[5]    =  {       1.0/4,      3.0/8,      1.2/1.3,           1.0,    1.0/2 },
                 b[5][5] = {{       1.0/4,          0,            0,             0,        0 },
   				            {      3.0/32,     9.0/32,            0,             0,        0 },
   					        { 1.932/2.197, -7.2/2.197,  7.296/2.197,             0,        0 },
   					        {   4.39/2.16,       -8.0,   3.68/0.513,   -8.45/41.04,        0 },
   					        {     -8.0/27,        2.0, -3.544/2.565,   1.859/4.104, -11.0/40 }},
                 c[6]    =  {    1.6/13.5,          0, 6.656/12.825, 2.8561/5.6430,  -9.0/50,  2.0/55 },
                 d[6]    =  {     1.0/360,          0,  -1.28/42.75, -2.197/75.240,   1.0/50,  2.0/55 };

    //--------------------------------------------------------------------------
    //                      >> FLUXO PRINCIPAL <<
    //--------------------------------------------------------------------------

    // Zera os erros nas soluções.
    for( register int i = 0; i <= this->order-1; i++ ) error_y[i] = 0;

    // Condição inicial.
    last_x = x[0];
    for( register int i = 0; i <= this->order-1; i++ ){
        last_y[i] = yI[i];
	    y[i][0] = last_y[i];

	    next_y[i] = last_y[i];
    }

    // Inicia o cronômetro.
    time_t t = time(NULL);

    // Calcula x e *y enquanto x[0] < x < x[nPtos-1].
    do{
		
      	// Contabiliza as iterações.
      	++nCalc;
        fileOut << nCalc << "\t";

      	// Próximo valor de x.
     	next_x = last_x + dx;
        fileOut << next_x << "\t";

        // next_x é um ponto de *x? Se sim, atribuirá a solução a **y (se aprovada).
        if( fabs( next_x - x[pto] ) < EPS ) atribui = true;
        else atribui = false;

	    try{

            /*
             * ---------------------------------------------------------------------------------------
             * Uma vez definido qual será o próximo ponto a resolver, procede a integração: calcula...
             * a) ... os fatores de Fehlberg;
             * b) ... next_y[ordem-1], next_y[ordem-2], ..., next_y[0] nesta seqüência;
             * c) ... os erros cometidos.
             * ---------------------------------------------------------------------------------------
             */

            // a) ... os fatores de Fehlberg.
            if( this->order > 1 ){
                for( int i = 0; i <= this->order-2; i++ ){
                    f[0][i] = last_y[i+1];
                    checkOverflow( f[0][i] );
                }
            }

            f[0][this->order-1] = (*F)( last_x, last_y );
            checkOverflow( f[0][this->order-1] );

            for( register int n = 1; n <= 5; n++ ){
                for( register int i = this->order-1; i >= 0; i-- ){

                    soma = 0;
                    for( register int j = 0; j <= n-1; j++ ){
                        soma += b[n-1][j] * f[j][i];
                        checkOverflow( soma );
                    }

                    next_y[i] = last_y[i] + soma * dx;
                    checkOverflow( next_y[i] );
                }

                if( this->order > 1 ){
                    for( register int i = 0; i <= this->order-2; i++ ){
                        f[n][i] = next_y[i+1];
                        checkOverflow( f[n][i] );
                    }
                }

                f[n][this->order-1] = (*F)( last_x + a[n-1]*dx, next_y );
                checkOverflow( f[n][this->order-1] );
            }

            // b) ... next_y[ordem-1], next_y[ordem-2], ..., next_y[0] nesta seqüência.
            // c) ... os erros cometidos.
            for( register int i = 0; i <= this->order-1; i++ ){
	
                next_y[i] = error_y[i] = 0;
	
                for( register int j = 0; j <= 5; j++ ){
                    next_y[i]  += c[j] * f[j][i];
                    checkOverflow( next_y[i] );

                    error_y[i] += d[j] * f[j][i];
                    checkOverflow( error_y[i] );
                }

                next_y[i] = last_y[i] + next_y[i] * dx;
                checkOverflow( next_y[i] );

                // Seleciona o maior erro cometido na aproximação acima.
                if( fabs(error_y[i]) > fabs(error_y[i_ERROR]) ) i_ERROR = i;
            }
	
            error_y[i_ERROR] = fabs(error_y[i_ERROR]);
            fileOut << error_y[i_ERROR] << "\t";
            /*
             * -----------------------------------------------------
             * Calculada a solução de quinta ordem, aprova-a ou não:
             * a) compara o erro obtido com o exigido;
             * b) se aprovada, ...
             *   b.1) ... escolhe a potência apropriada para correção de dx;
             *   b.2) ... redefine o último ponto válido calculado (last_x e last_y);
             *   b.3) ... se o ponto calculado for um de satisfação, ...
             *     b.3.i)  ... atribui os valores encontrados ao vetor de satisfação (**y);
             *     b.3.ii) ... incrementa o ponto de satisfação;
             * c) se reprovada, ...
             *   c.1) ... escolhe a potência apropriada para a correção de dx.
             * -----------------------------------------------------
             */
            if( error_y[i_ERROR] <= this->error ){

                exponent = STRETCH_EXPONENT;

                // Valida a aproximação recém-calculada (move-se para o próximo ponto).
                last_x = next_x;
                for( register int i = 0; i <= this->order-1; i++ ){
                    last_y[i] = next_y[i];
                }

                // Atribui a solução recém-calculada a **y caso next_x seja um ponto de *x.
                if( atribui ){
                    for( register int i = 0; i <= this->order-1; i++ ){
                        y[i][pto] = last_y[i];
                    }

                    //atribui = false;
                    ++pto; // <---------------------- mudou o ponto de satisfação.
                }
            }
            else{
                exponent = SHRINK_EXPONENT;
            }

            /*
             * ---------------------------------
             * Define dx para o próximo ponto.
             * ---------------------------------
             */
            last_dx = dx;
            fileOut << last_dx << "\t";

            dx_maximo = x[pto] - last_x;
            if ( fabs(dx_maximo) > INFINITY ) dx_maximo = INFINITY;
            else if( fabs(dx_maximo) < EPS ) dx_maximo = EPS;

            if( error_y[i_ERROR] < EPS ){
                dx = dx_maximo;
            }
            else{
                factor = SAFETY * pow( this->error/error_y[i_ERROR], exponent );
                if( factor > MAX_STRETCH ) factor = MAX_STRETCH;
                else if( factor < MAX_SHRINK ) factor = MAX_SHRINK;

                dx = factor * last_dx;
                if( fabs(dx) > dx_maximo) dx = dx_maximo;
                else if( fabs(dx) < EPS) dx = EPS;
            }

            fileOut << dx << endl;
	    }
        catch( OverflowException e ){
            cout << "Overflow... trying to recover." << endl;

            if ( this->isLinear ){
                for( register int i = 0; i <= this->order-1; i++ ){

                    last_y[i] *= EPS;

                    for( register unsigned int j = 0; j <= pto-1; j++ ){
                        y[i][j] *= EPS;
                    }
                }
            }
            else{
                throw( e );
            }
        }

        // Impede que a integração leve mais tempo que timeout.
        if( this->hasTimeout && difftime(time(NULL),t) > this->timeout ){
            throw TimeoutException();
        }

    }while(( last_x < x[nPtos-1] == crescente ) && ( pto < nPtos ));

    fileOut.close();
    delete[] last_y, next_y, error_y, f;

    cout << "Number of iterations: " << nCalc << ".\n";
};
//------------------------------------------------------------------------------
void RungeKuttaFehlberg::checkOverflow( double number ){
  if ( fabs(number) > INFINITY ){
    throw OverflowException();
  }
}

//---------------------------- FIM DA IMPLEMENTAÇÃO DA CLASSE RUNGEKUTTAFEHLBERG
