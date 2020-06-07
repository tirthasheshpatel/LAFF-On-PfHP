#define alpha( i,j ) A[ (j)*lda + (i) ]
#define  beta( i,j ) B[ (j)*ldb + (i) ]
#define gamma( i,j ) C[ (j)*ldc + (i) ]

void MyGer(int, int, double *, int, double *, int, double *, int);

void MyGemm(int m, int n, int k,
            double *A, int lda,
            double *B, int ldb,
            double *C, int ldc)
{
    int p;
    for( p=0 ; p<k ; p++ )
        MyGer( m, n, &alpha( 0,p ), 1, &beta( p,0 ), ldb, C, ldc );
}
