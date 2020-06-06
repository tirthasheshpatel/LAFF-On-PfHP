#define chi( i ) x[ (i)*incx ]
#define psi( i ) y[ (i)*incy ]

void Dots(int n,
          const double *x, int incx,
          const double *y, int incy,
          double *gamma)
/**
 * Dot product with an update to a scalar.
 * Performs `gamma = gamma + x^T y` where
 * `gamma` is a scalar and `x` and `y` are vectors.
 */
{
    int i;
    for( i=0 ; i<n ; i++ )
        *gamma += chi( i ) * psi( i );
}
