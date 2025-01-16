#include <iostream>

using namespace std;

void MatrixMultiplication(double* A, int m, int n, double* B, int p, double* C);
void MatrixMultVector(double* M, int n, double* V,  double* Res);
int main()
{
    double cx = 800.0;
    double cy = 550.0;
    double focus = 0.025;
    double pixel = 9 * 1e-6;

    double shift = -50.0;

    double x_r = 650.0;
    double y_r = 650.0;
    double R[3];
    R[0] = x_r;
    R[1] = y_r;
    R[2] = 1.0;

    double M[9];
    M[0] = 1.0; M[1] = 0.0; M[2] = shift;
    M[3] = 0.0; M[4] = 1.0; M[5] = shift;
    M[6] = 0.0; M[7] = 0.0; M[8] = 1.0;





    double L[9];
    //MatrixMultiplication(M, 3, 3, R, 1, L);
    MatrixMultVector(M, 3, R, L);

//    printf("L:\n");
//    for(int i = 0; i < 9; i++)  printf("%lf ", L[i]);
//    printf("\n");

    double L_vector[3];
    L_vector[0] = L[0]/L[2];
    L_vector[1] = L[1]/L[2];
    L_vector[2] = L[2]/L[2];

    printf("L_vector:\n");
    for(int i = 0; i < 3; i++)  printf("%lf ", L_vector[i]);
    printf("\n");

    /*********************************************************/
    printf("/*********************************************************/\n");
    double  R_4d[4];
    R_4d[0] = R[0]*focus/pixel;
    R_4d[1] = R[1]*focus/pixel;
    R_4d[2] = R[2]*focus/pixel;
    R_4d[3] = focus/pixel;

    double M_4D[16];
    M_4D[0] = 1.0;
    M_4D[1] = 0.0;
    M_4D[2] = 0.0;
    M_4D[3] = shift;

    M_4D[4] = 0.0;
    M_4D[5] = 1.0;
    M_4D[6] = 0.0;
    M_4D[7] = shift;

    M_4D[8] = 0.0;
    M_4D[9] = 0.0;
    M_4D[10] = 1.0;
    M_4D[11] = shift;

    M_4D[12] = 0.0;
    M_4D[13] = 0.0;
    M_4D[14] = 0.0;
    M_4D[15] = 1.0;

    double L_4D[4];

    MatrixMultVector(M_4D, 4, R_4d, L_4D);

    double L_vector_4d[4];
    L_vector_4d[0] = L_4D[0]/L_4D[3];
    L_vector_4d[1] = L_4D[1]/L_4D[3];
    L_vector_4d[2] = L_4D[2]/L_4D[3];
    L_vector_4d[3] = L_4D[3]/L_4D[3];

    printf("L_vector_4d:\n");
    for(int i = 0; i < 4; i++)  printf("%lf ", L_vector_4d[i]);
    printf("\n");

    cout << "Hello World!" << endl;
    return 0;
}

void MatrixMultiplication(double* A, int m, int n, double* B, int p, double* C) {
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            {
                C[i*n + j] = 0.0;
                for(int k = 0; k < p; k++) {
                    C[i*n + j] += A[i*m + k] * B[k*p + j];
                }
            }
        }
    }
}

void MatrixMultVector(double* M, int n, double* V,  double* Res) {
    for(int i = 0; i < n; i++) Res[i] = 0.0;

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            Res[j] += M[j*n + i] * V[i];
        }
    }
}
