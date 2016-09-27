#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/src/KroneckerProduct/KroneckerTensorProduct.h"
#include "time.h"

using namespace std;
using namespace Eigen;

//params
const int  N=4;
const int  K=2;
const float e = 0.5f;
const int A[] ={2,4};

float calcW(int y, int x) {
    float retVal = 0.0f;
    if (x == y) {
        retVal = 1 - e;
    } else if (y == 2) {
        retVal = e;
    } else {
        retVal = 0;
    }
    return retVal;
}

VectorXi index_o(int n, VectorXi &x){
    VectorXi ret(n/2);
    for(int i = 0; i < n; i++){
        if ( i % 2 == 0) {
            ret[i/2] = x[i];
        }
    }
    return ret;
}

VectorXi index_e(int n, VectorXi &x){
    VectorXi ret(n/2);
    for(int i = 0; i < n; i++){
        if (!(i % 2 == 0)) {
            ret[(i-1)/2] = x[i];
        }
    }
    return ret;
}

float calcW_i(int i, int n, VectorXi &u, VectorXi &u_i, VectorXi &y) {
    float W_i = 0.0f;
    //初期化
    if ( n == 2 ){
        if ( i == 1 ) {
            W_i = 0.5 * (calcW(y[0] ,(u[0]+0) % 2) * calcW(y[1],0)
                    + calcW(y[0] ,(u[0]+1) % 2) * calcW(y[1],1));
        } else {
            W_i = 0.5 * calcW(y[0] ,(u[0]+u[1]) % 2) * calcW(y[1],u[1]);
        }
    } else {
        VectorXi tempY1;
        VectorXi tempY2;
        for (int j = 0; j < n ; j++) {
            if(i < n/2){
                tempY1[j]= y[j];
            } else {
                tempY2[j - n/2] = y[j];
            }
        }
        if ( i % 2 == 0 ) {    //odd index
            W_i = 0.5 * calcW_i(i/2, n/2, index_e(i-1,u)+index_o(i-1,u) ,(u[i]+u[i+1]) % 2,tempY1)
                  * calcW_i(i/2, n/2, index_e(i-1,u), u[i], tempY2);
        } else {    //even index
            W_i = 0.5 * calcW_i(i/2, n/2, index_e(i-2,u)+index_o(i-2,u) ,(u[i]+u[i-1]) % 2,tempY1)
                  * calcW_i(i/2, n/2, index_e(i-2,u), u[i], tempY2);
        }
    }
    return W_i;
}


VectorXi encoder(int n, VectorXi &input){
    VectorXi x_n(n);
    VectorXi s_n(n);
    VectorXi v_n(n);

    for (int i = 0; i < n/2 ; i++) {
        s_n[2*i]   = (input[2*i+1] + input[2*i]) % 2;
        s_n[2*i+1] = input[2*i+1];
    }

    for (int i = 0; i < n/2 ; i++) {
        v_n[i] = s_n[2*i];
        v_n[n/2 + i] = s_n[2*i +1];
    }

    if(n == 2){
        x_n[0] = (v_n[0] + v_n[1]) % 2;
        x_n[1] = v_n[1];

    } else {
        VectorXi tempV_n1(n/2);
        VectorXi tempV_n2(n/2);

        for (int i = 0; i < n ; i++) {
            if(i < n/2){
                tempV_n1[i]= v_n[i];
            } else {
                tempV_n2[i - n/2] = v_n[i];
            }
        }

        VectorXi tempX_n1 = encoder(n/2, tempV_n1);
        VectorXi tempX_n2 = encoder(n/2, tempV_n2);

        for (int i = 0; i < n ; i++) {
            if(i < n/2){
                x_n[i] = tempX_n1[i];
            } else {
                x_n[i] = tempX_n2[i - n/2];
            }
        }

    }
    return x_n;
}

VectorXi channel(VectorXi &input){
    VectorXi y(N);

    for (int i = 0; i < N; i++) {
        y[i] = input[i];
        if((float)rand() / RAND_MAX < e){
            y[i] = 2;
        }
    }
    return y;
}

VectorXi decoder(VectorXi &input, VectorXi &u_Ac){
//    VectorXd w_n_i(N);
//    VectorXd w_n = VectorXd::Constant(N,e);
//    VectorXd h_i(N);
    VectorXi u_n_est(N);
//
//    //h_i計算
//    for (int i = 0; i < N; i++) {
//        int count = 0;
//        for (int j = i; j < N; j++) {
//            cout << input[j] << endl;
//            if( input[j] == 0 || input[j] == 1 ) {
//                count++;
//            }
//        }
//        w_n_i =  count * w_n / (pow(2, N-1));
//        float llr = 1; //w_n_i[i]  / w_n_i[i];
//        if (llr >= 1) {
//            h_i[i] = 0;
//        } else {
//            h_i[i] = 1;
//        }
//    }
//
//    cout << "w_n_i" << endl;
//    cout << w_n_i << endl;
//
//    cout << "h_i" << endl;
//    cout << h_i << endl;
//
////
////    //u_n_est計算
////    for (int i = 0; i < N; i++) {
////        int in_A = 0;
////        for (int j = 0; j < K; j++) {
////            if ( i == A[j]) {
////                in_A = 1;
////            }
////        }
////
////        // Aに含まれるindexなら既知
////        if (in_A == 1) {
////            u_n_est[i] = u_Ac[i];
////        } else {
////            u_n_est[i] = h_i[i];
////        }
////    }
//
    return u_n_est;
}

int main(void) {
    VectorXi u_Ac(K);
    u_Ac << 1, 0;
    VectorXi u_n(10);
    u_n << 1,2,3,4,5,6,7,8,9,0;  //input
    cout << "u_n" << endl;
    cout << u_n << endl;

    srand((int) time(NULL));

//    VectorXi x_n = encoder(N, u_n);
//    VectorXi y_n = channel(x_n);
//    float W_i = calcW_i(1, 2, u_n, y_n);

    //VectorXi u_n_est = decoder(y_n,u_Ac);
//    cout << "x_n" << endl;
//    cout << x_n << endl;
//    cout << "y_n" << endl;
//    cout << y_n << endl;
//    cout << "W_i" << endl;
//    cout << W_i << endl;
//    cout << "u_n_est" << endl;
    //cout << u_n_est << endl;

    return 0;
}
