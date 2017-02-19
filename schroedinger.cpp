#include "RungeKuttaFehlberg.h"
#include <stdio.h>
#include <math.h>
#include <iostream.h>

#include "DiscontinuityException.cpp"
#include "TimeoutException.cpp"
#include "UnreachableAccuracyException.cpp"

//#define SELF_ENERGY

using std::cout;


double energia = 0.006076;


double **psi;
//double **potential;

//------------------------------------------------------------------------------
// Equação de Schrödinger: y'' = 2m/hbar^2 * (V-E) * y.
//------------------------------------------------------------------------------
double schroedinger(double x0, double *y0){

    const double mEff  = 0.067,
                 hbar2 = 7.6199682 * mEff;

    //double V = 1/( 1 + exp((x0-40)/10) ) + 1/( 1 + exp((60-x0)/10) );
    double V = ( x0 > 30 && x0 < 70 ? 0 : 0.399 );

    return( (2 * mEff / hbar2)*(V-energia)*y0[0] );
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int main( int argc, char *argv[] ){

    int nPtos = 1000, ordem = 2;

    double xI  = 0,
           xF  = 100,
           y0[2] = { 0, 1e-5 };

    double dx = (xF-xI)/(nPtos-1), next_x;
                                   
    //------------------------------------
    // Abre o arquivo para escrita.
    //------------------------------------
    FILE *arquivo;

    if( (arquivo = fopen( argv[1], "w+" )) == NULL ){
        printf( "I could not open the output file %s\n", argv[1] );
        return( -1 );
    }

    //------------------------------------
    // Monta os vetores *xE e *xD.
    //------------------------------------
    int ptosComuns = 3;
    if( nPtos < ptosComuns ){
        if( nPtos % 2 != 0 ) ptosComuns = nPtos;
        else ptosComuns = nPtos-1;
    }

    int meio = (nPtos-1) >> 1;

    ptosComuns = ptosComuns >> 1;
    if( ptosComuns == 0 ) ptosComuns = 1;

    const int nPtos_E =     1 + meio + ptosComuns,
              nPtos_D = nPtos - meio + ptosComuns;

    double *xE = new double[nPtos_E],
           *xD = new double[nPtos_D];

    for( register int i = 0; i <= nPtos-1; i++ ){

        next_x = xI + i*dx;

	    if( i <= meio+ptosComuns ) xE[i] = next_x;
	    if( i >= meio-ptosComuns ) xD[nPtos-1-i] = next_x;

    }

    //------------------------------------
    // Resolve a eq. de Schrödinger.
    //------------------------------------
    RungeKuttaFehlberg *EDO = new RungeKuttaFehlberg( ordem, schroedinger, true );
    EDO->setError( 1.0E-1 );
    EDO->setTimeout(10.0);
    /*
     * Define e aloca o estado inicial dos vetores-solução.
     */
    double **yE = new double * [ordem],
           **yD = new double * [ordem];

    for( register int i = 0; i <= ordem-1; i++ ){
        yE[i] = new double[nPtos];
        yD[i] = new double[nPtos];
    }

    #ifdef SELF_ENERGY
    const double E_i = 1e-3, E_f = 0.4;
    const int N = 100;

    for( register int iEnergy = 0; iEnergy <= N-1; iEnergy++ ){

        energia = E_i + iEnergy * ( E_f - E_i ) / ( N - 1 );
        cout << energia << "\n";
    #endif
        for( register int i = 0; i <= ordem-1; i++ )
            for( register int j = 0; j <= nPtos-1; j++ ){
                yE[i][j] = 0;
            }

        try{
            EDO->solve( nPtos_E, xE, yE, y0, "RKFesq.log" );
            EDO->solve( nPtos_D, xD, yD, y0, "RKFdir.log" );

            #ifdef SELF_ENERGY
            double avaliacao = fabs(yE[1][meio] / yE[0][meio] - yD[1][meio] / yD[0][meio] );
            cout << "avaliacao = " << avaliacao;


            /*
             * Constrói a função de onda.

            for( register int j = 0; j <= nPtos-1; j++ ){
                if( j < meio )
                    psi[j] = yE[];
                else
                    psi[j] = yD[];
            } */

            /*
             * Normaliza a função de onda.
             */
            //normalize( psi );

            fprintf( arquivo, "%f\t%f\n", energia, avaliacao );
            #endif
        }
        catch( std::bad_alloc ){
            cout << "not enough memory to run.\n";
            exit( -1 );
        }
        catch( DiscontinuityException e ){
            cout << e.getMessage() << "\n";
        }
        catch( TimeoutException e ){
            cout << e.getMessage() << "\n";
        }
        catch( UnreachableAccuracyException e ){
            cout << e.getMessage() << "\n";
        }
        catch( RungeKuttaFehlbergException e ){
            cout << e.getMessage() << "\n";
        }

    #ifdef SELF_ENERGY
    }
    #endif

    delete EDO;

    /*
     * Imprime a solução no arquivo.
     */
    #ifndef SELF_ENERGY
    for( register int i = 0; i <= nPtos-1; i++ ){

        if( i <= meio+ptosComuns )
            fprintf( arquivo, "%f\t%f\t", xE[i], yE[0][i] );

   	    if( i >= meio-ptosComuns ){
            if( i > meio+ptosComuns ){
                fprintf( arquivo, "%f\t\t", xD[nPtos-1-i] );
            }
            fprintf( arquivo, "%f", yD[0][nPtos-1-i] );
        }
        fprintf( arquivo, "\n" );
    }
    #endif

    fclose( arquivo );

    return( 0 );
}


