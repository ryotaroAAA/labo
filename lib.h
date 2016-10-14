//
// Created by 鈴木凌太郎 on 2016/10/05.
//


//////////////////////////////////////////////////////////////////////////////////

#ifndef CHANNEL_POLARIZATION_LIB_H
#define CHANNEL_POLARIZATION_LIB_H
#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/src/KroneckerProduct/KroneckerTensorProduct.h"
#include "time.h"
#include <chrono>
#include <iomanip>
#include <string>    // useful for reading and writing
#include <fstream>   // ifstream, ofstream


//////////////////////////////////////////////////////////////////////////////////

using namespace Eigen;
using namespace std;

//params
const int  N=1024;
const int  K=2;
const float e = 0.5f;
const int A[] ={};

#define PRINT(X) cout << #X << ":\n" << setprecision(10) << X << endl << endl;
#define ARR(array)     (sizeof(array) / sizeof(array[0]));

////////////////////////////////////////////////////////////////////////////////

int compare_int(const void *a, const void *b);
void dispArray(double *x);
void dispArray(int n,int *x);
void dispArray(int n,double *x);
void outputArray(double *x);
int vsize(VectorXi &x);
double calcW(int y, int x);
VectorXi index_o(VectorXi &x);
VectorXi index_e(VectorXi &x);
VectorXi retBinary(VectorXi &x);
VectorXi generateUi(VectorXi &x,VectorXi &u_Ac);
double calcW_i(int i, int n, VectorXi &u, VectorXi &u_i, VectorXi &y);
double calcCapacityForBec(int i, int n);
void makeArrayCapacityForBec(double *array);
double calcL_i(int i, int n ,VectorXi &y ,VectorXi &u, int u_i_est);
VectorXi encoder(int n, VectorXi &input);
VectorXi channel(VectorXi &input);
VectorXi decoder(VectorXi &input, VectorXi &u_Ac);
double calcBhatForBec(int i, int n);
void probErrBound(double *array);

//////////////////////////////////////////////////////////////////////////////////

int compare_int(const void *x, const void *y) {
    const double* a=(const double*)x;
    const double* b=(const double*)y;
    if(*a>*b)return 1;
    if(*a<*b)return -1;
    return 0;
}
void dispArray(double *x){
    for(int i = 0; i < N; i++){
        printf("x[%d] = %f\n",i, x[i]);
    }
}
void dispArray(int n, double *x){
    for(int i = 0; i < n; i++){
        printf("x[%d] = %f\n",i, x[i]);
    }
}
void dispArray(int n, int *x){
    for(int i = 0; i < n; i++){
        printf("x[%d] = %d\n",i, x[i]);
    }
}
void outputArray(double *x){
    string filename = "/Users/ryotaro/labo/symmetric_capacity";
    ofstream w_file;
    w_file.open(filename, ios::out);
    for (int i = 0; i<N; i++)
    {
        w_file << i << " " << x[i]<< endl;
    }
}
int vsize(VectorXi &x){
    return (x.array() >= 0).count();
}
VectorXi index_o(VectorXi &x){
    VectorXi ret(vsize(x)/2);
    for(int i = 0; i < vsize(x); i++){
        if ( i % 2 == 0) {
            ret[i/2] = x[i];
        }
    }
    return ret;
}
VectorXi index_e(VectorXi &x){
    VectorXi ret(vsize(x)/2);
    for(int i = 0; i < vsize(x); i++){
        if (!(i % 2 == 0)) {
            ret[(i-1)/2] = x[i];
        }
    }
    return ret;
}
VectorXi retBinary(VectorXi &x) {
    VectorXi ret(vsize(x));
    for(int i = 0; i < vsize(x) ; i++){
        ret[i] = x[i] % 2;
    }
    return ret;
}
VectorXi generateUi(int set, VectorXi &x,VectorXi &u_Ac){
    VectorXi ret(N);
    srand((int) time(NULL));
    for (int i = 0; i < N; i++) {
        int in_A = 0;
        int index = 0;
        for (int j = 0; j < K; j++) {
            if ( i+1 == A[j]) {
                in_A = 1;
                index = j;
            }
        }
        // Aに含まれるindexなら既知
        if (in_A == 1) {
            ret[i] = u_Ac[index];
        } else {
            ret[i] = rand() % 2;
        }

        if (set == 0) {
            ret[i] = 0;
        } else if (set == 1) {
            ret[i] = 1;
        }

    }
    return ret;
}
double calcW(int y, int x) {
    double retVal = 0.0f;
    if (x == y) {
        retVal = 1 - e;
    } else if (y == 2) {
        retVal = e;
    } else {
        retVal = 0;
    }
    return retVal;
}
double calcW_i(int i, int n, VectorXi &u, int u_i, VectorXi &y) {
    double W_i = 0.0;

    if ( n == 2 ){
        if ( i == 1 ) {
            W_i = 0.5 * (calcW(y[0] ,(u_i+0) % 2) * calcW(y[1],0)
                         + calcW(y[0] ,(u_i+1) % 2) * calcW(y[1],1));
        } else {
            W_i = 0.5 * calcW(y[0] ,(u[0]+u_i) % 2) * calcW(y[1],u_i);
        }
    } else {
        VectorXi tempY1(vsize(y)/2);
        VectorXi tempY2(vsize(y)/2);
        for (int j = 0; j < vsize(y) ; j++) {
            if(j < vsize(y)/2){
                tempY1[j]= y[j];
            } else {
                tempY2[j - vsize(y)/2] = y[j];
            }
        }

        VectorXi tempU;
        VectorXi tempU_bin;
        VectorXi u_e;

        tempU = index_e(u)+index_o(u);
        u_e = index_e(u);
        tempU_bin = retBinary(tempU);

        if ( i % 2 == 0 ) {
            W_i = 0.5 * (calcW_i(i/2, n/2, tempU_bin ,(1 + u_i) % 2,tempY1) * calcW_i(i/2, n/2, u_e, 1, tempY2)
                         + calcW_i(i/2, n/2, tempU_bin ,(0 + u_i) % 2,tempY1) * calcW_i(i/2, n/2, u_e, 0, tempY2));
        } else {
            W_i = 0.5 * calcW_i((i-1)/2, n/2, tempU_bin ,(u_i+u[i-1]) % 2, tempY1)
                  * calcW_i( i/2   , n/2, u_e       , u_i            , tempY2);
        }
    }
//    cout << "W_" << "n:" << n << ",i:" << i << endl;
//    PRINT(W_i);
    return W_i;
}
double calcL_i(int i, int n ,VectorXi &y ,VectorXi &u, int u_i_est) {
    double lr = 0.0;
    if ( n == 1 ) {
        lr = calcW(y[0],0) / calcW(y[0],1);
    } else {
        VectorXi tempY1(vsize(y)/2);
        VectorXi tempY2(vsize(y)/2);
        for (int j = 0; j < vsize(y) ; j++) {
            if(j < vsize(y)/2){
                tempY1[j]= y[j];
            } else {
                tempY2[j - vsize(y)/2] = y[j];
            }
        }

        VectorXi tempU;
        VectorXi tempU_bin;
        VectorXi u_e;

        tempU = index_e(u)+index_o(u);
        u_e = index_e(u);
        tempU_bin = retBinary(tempU);
        double temp1 = calcL_i(i/2, n/2, tempY1, tempU_bin, u_i_est);
        double temp2 = calcL_i(i/2, n/2, tempY2, u_e, u_i_est);
        if ( i % 2 == 0) {
            lr = ( 1 + temp1 * temp2 ) / ( temp1 + temp2 );
        } else {
            lr = pow(temp1, 1-2*u[i-1]) * temp2;
        }
    }
    if (isinf(lr) || isnan(lr)) {
        lr = 1;
    }
    //cout << i << " " << n << " " << lr << endl;
    return lr;
}
double calcCapacityForBec(int i, int n) {
    double cap =0.0;
    if ( i == 0 && n == 1 ) {
        cap = 1 - e;
    } else {
        if ( i % 2 == 0) {
            cap = pow(calcCapacityForBec(i/2, n/2),2);
        } else {
            double tempCap = calcCapacityForBec((i-1)/2, n/2);
            cap = 2 * tempCap - pow(tempCap,2);
        }
    }
    //cout << "i:" << i << ", n:" << n <<endl;
    PRINT(cap);
    return cap;
}
double calcBhatForBec(int i, int n){
    double bha =0.0;
    if ( i == 0 && n == 1 ) {
        bha = e;
    } else {
        if ( i % 2 == 0) {
            double tempBha = calcBhatForBec(i/2, n/2);
            bha = 2 * tempBha - pow(tempBha,2);
        } else {
            bha = pow(calcBhatForBec((i-1)/2, n/2),2);
        }
    }
//    cout << "i:" << i << ", n:" << n <<endl;
//    PRINT(bha);
    return bha;
}
void makeArrayCapacityForBec(double *array) {
    for(int i = 0; i < N; i++){
        array[i] = calcCapacityForBec(i,N);
    }
}
void probErrBound(double *array) {
    int N = pow(2,20);
    int count = 0;
    double tempArr[N];
    double sumArr[N];
    double sum = 0.0;

    for(int i = 0; i < N; i++) {
        tempArr[i] = calcBhatForBec(i, N);
    }
    //double bhatArr[N];
//    memcpy(bhatArr, tempArr, sizeof(double) * N);
    qsort(tempArr, N, sizeof(double), compare_int);
//    dispArray(tempArr);

//    for(int i = 0; i < j; i++) {
//        array[i] = tempArr[i];
//    }

    for(int i=0; i < N; i++){
        sum += tempArr[i];
        sumArr[i] = sum;
    }
//    for(int i=0; i < N; i++){
//        cout << (double)i/N << " " << sumArr[i] << " " << tempArr[i] << endl;
//    }

    string filename = "/Users/ryotaro/labo/err_bound3";
    ofstream w_file;
    w_file.open(filename, ios::out);
    for (int i = 0; i < N; i++)
    {
        w_file << (double)i/N << " " << sumArr[i] << " " << tempArr[i] << endl;
    }
}


//////////////////////////////////////////////////////////////////////////////////

#endif //CHANNEL_POLARIZATION_LIB_H
