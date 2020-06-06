#define alpha( i,j ) A[ (j)*lda + (i) ]
#define     chi( i ) x[ (i)*incx ]
#define     psi( i ) y[ (i)*incy ]

void Axpy(int, double,
          double *, int,
          double *, int);

void MyGer(int m, int n,
           double *x, int incx,
           double *y, int incy,
           double *A, int lda)
{
    int j;
    for( j=0 ; j<n ; j++ )
        Axpy( m, psi( j ), x, incx, &alpha( 0,j ), 1 );
}
