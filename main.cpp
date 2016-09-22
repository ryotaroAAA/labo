#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/src/KroneckerProduct/KroneckerTensorProduct.h"
#include "time.h"

using namespace std;
using namespace Eigen;

VectorXi encoder(int tempN, VectorXi &input);
VectorXi channel(VectorXi &input);
VectorXi decoder(VectorXi &input, VectorXi &u_Ac);

void get_gn_coset(MatrixXi &input);
void get_uA(MatrixXi &input);
void get_g_nA(MatrixXi &input);
VectorXf W_n(VectorXf &input);

//params
const int  N=16;
const int  K=2;
const float e = 0.5f;
const int A[] ={2,4};


int main(void) {
    VectorXi u_Ac(K);
    u_Ac << 1,0;
    VectorXi u_n(N);
    u_n << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;  //input
    cout << "u_n" << endl;
    cout << u_n << endl;

    srand((int)time(NULL));

    VectorXi x_n = encoder(N,u_n);
//    VectorXi y_n = channel(x_n);
//    VectorXi u_n_est = decoder(y_n,u_Ac);
    cout << "x_n" << endl;
    cout << x_n << endl;
//    cout << "y_n" << endl;
//    cout << y_n << endl;
//    cout << "u_n_est" << endl;
    //cout << u_n_est << endl;

    return 0;
}

VectorXf W_i_n(int i, int n, int u_i, VectorXf &y){
    VectorXf ret;
    //初期化
//    if ( n == 2 ){
//
//
//    } else {
//
//
//    }
}


VectorXi encoder(int tempN, VectorXi &input){
    cout << "////////////////////////////////" << endl;

    VectorXi x_n(tempN);
    VectorXi s_n(tempN);
    VectorXi v_n(tempN);

    for (int i = 0; i < tempN/2 ; i++) {
        s_n[2*i]   = (input[2*i+1] + input[2*i]) % 2;
        s_n[2*i+1] = input[2*i+1];
    }

    cout << "s_n" << endl;
    cout << s_n << endl;

    for (int i = 0; i < tempN/2 ; i++) {
        v_n[i] = s_n[2*i];
        v_n[tempN/2 + i] = s_n[2*i +1];
    }

    cout << "v_n" << endl;
    cout << v_n << endl;

    if(tempN == 2){
        x_n[0] = v_n[0] + v_n[1];
        x_n[1] = v_n[1];

    } else {
        VectorXi tempV_n1(tempN/2);
        VectorXi tempV_n2(tempN/2);

        for (int i = 0; i < tempN ; i++) {
            if(i < tempN/2){
                tempV_n1[i]= v_n[i];
            } else {
                tempV_n2[i - tempN/2] = v_n[i];
            }
        }
        cout << "tempN" << endl;
        cout << tempN/2 << endl;
        cout << "tempV_n1" << endl;
        cout << tempV_n1 << endl;
        cout << "tempV_n2" << endl;
        cout << tempV_n2 << endl;

        VectorXi tempX_n1 = encoder(tempN/2, tempV_n1);
        VectorXi tempX_n2 = encoder(tempN/2, tempV_n2);

        for (int i = 0; i < tempN ; i++) {
            if(i < tempN/2){
                x_n[i] = tempX_n1[i];
            } else {
                x_n[i] = tempX_n2[i - tempN/2];
            }
        }

    }
    cout << "return x_n" << endl;
    cout << x_n << endl;
    return x_n;
}



VectorXi channel(VectorXi &input){
    VectorXi y(N);

    for (int i = 0; i < N; i++) {
        y[i] = input[i];
//        cout << y[i];
        if((float)rand() / RAND_MAX < e){
            y[i] = 2;
        }
    }
//    cout << rand() << endl;
//    cout << (float)rand() / RAND_MAX << endl;
    return y;
}

VectorXi decoder(VectorXi &input, VectorXi &u_Ac){
    VectorXd w_n_i(N);
    VectorXd w_n = VectorXd::Constant(N,e);
    VectorXd h_i(N);
    VectorXi u_n_est(N);

    //h_i計算
    for (int i = 0; i < N; i++) {
        int count = 0;
        for (int j = i; j < N; j++) {
            cout << input[j] << endl;
            if( input[j] == 0 || input[j] == 1 ) {
                count++;
            }
        }
        w_n_i =  count * w_n / (pow(2, N-1));
        float llr = 1; //w_n_i[i]  / w_n_i[i];
        if (llr >= 1) {
            h_i[i] = 0;
        } else {
            h_i[i] = 1;
        }
    }

    cout << "w_n_i" << endl;
    cout << w_n_i << endl;

    cout << "h_i" << endl;
    cout << h_i << endl;

//
//    //u_n_est計算
//    for (int i = 0; i < N; i++) {
//        int in_A = 0;
//        for (int j = 0; j < K; j++) {
//            if ( i == A[j]) {
//                in_A = 1;
//            }
//        }
//
//        // Aに含まれるindexなら既知
//        if (in_A == 1) {
//            u_n_est[i] = u_Ac[i];
//        } else {
//            u_n_est[i] = h_i[i];
//        }
//    }

    return u_n_est;
}
