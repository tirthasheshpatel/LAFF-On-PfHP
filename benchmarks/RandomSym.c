#define _XOPEN_SOURCE
#include <stdlib.h>

#define alpha( i,j ) A[ (j)*ldA + (i) ]

void RandomSym(int m, double *A, int ldA)
{
    int i,j;
    for( i=0 ; i<m ; i++ ) {
        for( j=0 ; j<=i ; j++ ) {
            alpha( i,j ) = drand48();
            alpha( j,i ) = drand48();
        }
    }
}