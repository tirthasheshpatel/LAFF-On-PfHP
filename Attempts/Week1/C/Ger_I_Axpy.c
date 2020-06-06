#define alpha( i,j ) A[ (j)*lda + (i) ]
#define   chi( i ) x[ (i)*incx ]
#define   psi( i ) y[ (i)*incy ]

void Axpy(int, double,
          double *, int,
          double *, int);

void MyGer(int m, int n,
           double *x, int incx,
           double *y, int incy,
           double *A, int lda)
{
    int i;
    for( i=0 ; i<m ; i++ )
        Axpy(n, chi( i ), y, incy, &alpha( i,0 ), lda);
}