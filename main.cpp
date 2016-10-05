#include "lib.h"

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

    if (n == 2) {
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
        if((long double)rand() / RAND_MAX < e){
            y[i] = 2;
        }
    }
    return y;
}
VectorXi decoder(VectorXi &y, VectorXi &u, VectorXi &u_Ac){

    VectorXd h_i(N);
    VectorXi u_n_est(N);

    //h_i計算
    for (int i = 0; i < N; i++) {
        int count = 0;
        for (int j = i; j < N; j++) {
            if( y[j] == 0 || y[j] == 1 ) {
                count++;
            }
        }
        long double llr = calcL_i(i, N, y, u, u[i]);
        if (llr >= 1) {
            h_i[i] = 0;
        } else {
            h_i[i] = 1;
        }
    }
    PRINT(h_i);

    //u_n_est計算
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
            u_n_est[i] = u_Ac[index];
        } else {
            u_n_est[i] = h_i[i];
        }
    }

    return u_n_est;
}

int main(void) {
    int i = 1;
    VectorXi u_Ac(K);
    u_Ac << 1, 0;
    VectorXi u_n(N);

    //calc exec time
    const auto startTime = chrono::system_clock::now();

    u_n = generateUi(u_n, u_Ac);
//    PRINT(u_n);
    VectorXi x_n = encoder(N, u_n);
    VectorXi y_n = channel(x_n);
//    PRINT(x_n);
//    PRINT(y_n);

//    PRINT(calcCapacityForBec(i-1, N));

    //long double cap[N] = {0};
    //makeArrayCapacityForBec(cap);
    //dispArray(cap);
    //outputArray(cap);
    calcL_i(i, N, y_n, u_n, x_n(i));

//    long double W_i = calcW_i(i, N, u_n, u_n[i-1], y_n);
//    PRINT(W_i);
//
    VectorXi u_n_est = decoder(y_n, x_n, u_Ac);
    PRINT(u_n_est);
    const auto endTime = chrono::system_clock::now();
    const auto timeSpan = endTime - startTime;
    cout << "処理時間:" << chrono::duration_cast<chrono::milliseconds>(timeSpan).count() << "[ms]" << endl;
    const auto astartTime = chrono::system_clock::now();
    
    return 0;
}