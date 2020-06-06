#define chi( i ) x[ (i)*incx ]   // map chi( i ) to array x 
#define psi( i ) y[ (i)*incy ]   // map psi( i ) to array y

void Axpy( int m, double alpha, double *x, int incx, double *y, int incy )
{
  for ( int i=0; i<m; i++ )
    psi( i ) += alpha * chi( i );
}
