#define alpha( i,j ) A[ (j)*lda + (i) ]
#define  beta( i,j ) B[ (j)*ldb + (i) ]
#define gamma( i,j ) C[ (j)*ldc + (i) ]

void MyGemv(int, int,
            double *, int,
            double *, int,
            double *, int);

void MyGemm(int m, int n, int k,
            double *A, int lda,
            double *B, int ldb,
            double *C, int ldc)
/**
 * General Matrix-Matrix Multiplication using Matrix Vector Multiplies
 * Computes `C = C + AB` using either Dots or Axpy as chosen by the user.
 * Use `make J_Gemv_J_Axpy` to use Axpy operation to compute `C = C + AB`
 * Use `make J_Gemv_I_Dots` to use Dots operation to compute `C = C + AB`
 */
{
    int j;
    for( j=0 ; j<n ; j++ )
        MyGemv(m, n, A, lda, &beta( 0,j ), 1, &gamma( 0,j ), 1 );
}
