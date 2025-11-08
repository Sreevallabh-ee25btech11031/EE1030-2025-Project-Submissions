#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define STB_IMAGE_IMPLEMENTATION
#include "../c_libs/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../c_libs/stb_image_write.h"

double norm(double vec[], int n){
    double sum = 0;
    for(int i=0; i<n; i++){
        sum += vec[i]*vec[i];
    }
    
    return sqrt(sum);
}

void copy(double **B, double **A, int m, int n){
//matrix A is being copied into matrix A
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            *(*(B+i)+j) = *(*(A+i)+j);
        }
    }
}

double Fnorm(double **A, int m, int n){

    double f = 0;

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            f += pow(*(*(A+i)+j),2);
        }
    }

    f = sqrt(f);

    return f;

}

double **transpose(double **A, int r, int c){
//B is created and made into A transpose
    double **B = (double **)malloc(c*sizeof(double *));
    for(int i=0; i<c; i++){
        B[i] = (double *)malloc(r*sizeof(double));
    }
    
    for(int i=0; i<r; i++){
        for(int j=0; j<c; j++){
            *(*(B+j)+i) = *(*(A+i)+j);
        }
    }

    return B;
    
}

double IP(double a[], double b[], int n){
    double sum = 0;
    for(int i=0; i<n; i++){
        sum += a[i]*b[i];
    }

    return sum;
}

double **MM(double **A, int ma, int na, double **B, int mb, int nb){

    if(na != mb){
        printf("invalid");
        return NULL;
    }

    double **M = (double **)malloc(ma*sizeof(double *));
    for(int i = 0; i<ma; i++){
        M[i] = (double *)malloc(nb*sizeof(double));
    }

    for(int i=0; i<ma; i++){
        for(int j=0; j<nb; j++){
            *(*(M+i)+j) = 0;
            for(int k=0; k<na; k++){
                *(*(M+i)+j) += *(*(A+i)+k)*(*(*(B+k)+j));
            }
        }
    }

    return M;
}

void QR(double **A, int m, int n, double **R){

    //logic is that A becomes the orthogonal one and R is the upper diagonal. 
    //while implementing, n = k.
    struct vectors{
        double v[m];
        double q[m];
    }; 

    struct vectors *qr = malloc(n*sizeof(struct vectors));

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            R[i][j] = 0;
        }
    }
    
    //inilitalising the vectors v
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            qr[i].v[j] = *(*(A+j)+i);   
        }
    }

    for(int i=0; i<n; i++){
        *(*(R+i)+i) = norm(qr[i].v, m);
        for(int j=0; j<m; j++){
            qr[i].q[j] = (qr[i].v[j])/(*(*(R+i)+i));
        }

        for(int j=i+1; j<n; j++){
            *(*(R+i)+j) = IP(qr[i].q, qr[j].v, m);
            for(int k=0; k<m; k++){
                qr[j].v[k] -= *(*(R+i)+j)*(qr[i].q[k]);
            }
        }

    }

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            *(*(A+i)+j) = qr[j].q[i];
        }
    }

}

void blockpow(double **A, int m, int n, int k, double **U, double **S, double **V){

// U, V, S are the U, V and Sigma in the SVD. R is the QR's R,which upon iterations goes to S. 
    double **At = transpose(A, m, n);
    double **R = (double **)malloc(k*sizeof(double *));
    for(int i = 0; i<k; i++){
        R[i] = (double *)malloc(k*sizeof(double));
    } //irrespective, R is kxk.

    //filling random entries into V, our initial matrix. 
    for(int i=0; i<n; i++){
        for(int j=0; j<k; j++){
            //useful to generate numbers between -1.0 and 1.0 for prevention of overflow in subsequent steps. 
           *(*(V+i)+j) = 2.0*((double)rand()/RAND_MAX)- 1.0;
        }
    }

    double err = 1;
    double val = 1;
    int iter = 0;

    //creating space for 3 matrices, used inside the for loop.
    double **C = malloc(m*sizeof(double *));
    for(int i=0; i<m; i++){
        C[i] = malloc(k*sizeof(double));
    }

    double **Y = malloc(m*sizeof(double *));
    for(int i=0; i<m; i++){
        Y[i] = malloc(k*sizeof(double));
    }

    double **Z = malloc(n*sizeof(double *));
    for(int i=0; i<n; i++){
        Z[i] = malloc(k*sizeof(double));
    }

    while(err/val>1e-6 && iter<20){
        Y = MM(A, m, n, V, n, k);

        QR(Y, m, k, R);

        copy(U, Y, m, k);

        Z = MM(At, n, m, U, m, k);

        QR(Z, n, k, R);

        copy(V, Z, n, k);

        err = 0; 
        val = 0;
        for(int i=0; i<k; i++){
            for(int j=0; j<k; j++){
                err += fabs(*(*(R+i)+j) - *(*(C+i)+j));
                val += fabs(*(*(R+i)+j));
            }
        }

        //Current value of R copied into a matrix C so that we can use it in the convergence checking condition next iteration. 
        //matrix C is not edited until the convergence condition is checked. 
        copy(C, R, k, k);

        iter++;

    }
    //freeing allocated matrices
    for(int i = 0; i<m; i++){
            free(*(Y+i));
        }

        free(Y);

    for(int i = 0; i<n; i++){
            free(*(Z+i));
        }

        free(Z);

    for(int i = 0; i<k; i++){
            free(*(C+i));
        }

        free(C);

    //since our loop has completed, we conclude that R is a good approximation of S, and hence copy it. 
    copy(S, R, k, k);

}

int main(){

    //width is n and height is m

    int width, height, channels;
    unsigned char *img = stbi_load("../../figs/globe.jpg", &width, &height, &channels, 1);

    if(img == NULL){
        printf("Error in loading image");
        return 1;
    }

    int  m = height;
    int n = width;
    int k = 50;

    double **A = malloc(m*sizeof(double*));
    for(int i=0; i<m; i++){
        A[i] = malloc(n*sizeof(double));
    }

    //converting the conitguous 1d array to 2d array
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            *(*(A+i)+j) = (double)(img[i*n + j]);
        }
    }

    //Creating V (nxk)
    double **V = malloc(n*sizeof(double *));
    for(int i=0; i<n; i++){
        V[i] = malloc(k*sizeof(double));
    }    

    //Creating U (mxk)
    double **U = malloc(m*sizeof(double *));
    for(int i=0; i<m; i++){
        U[i] = malloc(k*sizeof(double));
    }

    //Creating S (kxk)
    double **S = malloc(k*sizeof(double *));
    for(int i=0; i<k; i++){
        S[i] = malloc(k*sizeof(double));
    }

    blockpow(A, m, n, k, U, S, V);

    //finding Ak, the final matrix (low rank approximation)
    double **Ak = MM(MM(U, m, k, S, k , k), m, k, transpose(V, n, k), k, n);

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            if(*(*(Ak+i)+j)>255.0){
                *(*(Ak+i)+j) = 255.0;
            }

            if(*(*(Ak+i)+j)<0.0){
                *(*(Ak+i)+j) = 0.0;
            }
        }
    }

    // Converting Ak into a contiguous memory block
    unsigned char *out;
    out = malloc(m*n*(sizeof(unsigned char)));

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            out[i*n+j] = (unsigned char)*(*(Ak+i)+j);  
        }
    }

    //stbi_write_png("../../figs/greyscale_400.png", width, height, 1, out, width);
    stbi_write_jpg("../../figs/globe_50.jpg", width, height, 1, out, 90);

    //creating a matrix M = A - Ak, to compute error using frobenius norm. 
    double **M = malloc(m*sizeof(double *));
        for(int i=0; i<m; i++){
            M[i] = malloc(n*sizeof(double));
        }

        for(int i=0;i<m; i++){
            for(int j=0; j<n; j++){
                *(*(M+i)+j) = *(*(A+i)+j) - *(*(Ak+i)+j);
            }
        }

    //printing the error calculated by using frobenius norm and also percentage error. 
    printf("%lf\n", Fnorm(M,m,n));
    printf("%lf", (Fnorm(M, m, n)/Fnorm(A,m,n))*100);
    


    return 0;

    }

