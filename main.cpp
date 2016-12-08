#include "main.h"


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
        v_n[n/2 + i] = s_n[2*i + 1];
    }

    if (n == 2) {
        x_n[0] = (input[0] + input[1]) % 2;
        x_n[1] = input[1];

    } else {
        VectorXi tempV_n1(n/2);
        VectorXi tempV_n2(n/2);

        for (int i = 0; i < n ; i++) {
            if (i < n/2) {
                tempV_n1[i]= v_n[i];
            } else {
                tempV_n2[i - n/2] = v_n[i];
            }
        }
        hoge++;
//        cout << "encoder size:: " << vsize(tempV_n1) << endl;
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
        if((double)rand() / RAND_MAX < e){
            y[i] = 2;
        }
    }
    return y;
}

VectorXi decoder(VectorXi &y, VectorXi &u, int *Ac, int *A){
    VectorXd h_i(N);
    VectorXi u_n_est(N);
    int size = log2(N);

    vector<vector<bool>> isCache (size, vector<bool>(N,false));
    vector<vector<double>> cache (size, vector<double>(N,0.0));

    //u_n_est計算
    for (int i = 0; i < N; i++) {
        // Acに含まれるindexなら既知
        if (containNumInArray(i, N-K, Ac)) {
            u_n_est[i] = u[i];
        } else {
            cout << i << endl;
            //処理時間計測//
            const auto startTime = chrono::system_clock::now();

            int cache_i = i;// == 7 ? 7 : (4*i) % 7;
            double lr = calcL_i(i, N, cache_i, 0, y, u, u[i], isCache, cache);

            const auto endTime = chrono::system_clock::now();
            const auto timeSpan = endTime - startTime;
            cout << "処理時間:" << chrono::duration_cast<chrono::milliseconds>(timeSpan).count() << "[ms]" << endl;
            hogetime += chrono::duration_cast<chrono::milliseconds>(timeSpan).count();
            const auto astartTime = chrono::system_clock::now();

            if (lr >= 1) {
                h_i[i] = 0;
            } else {
                h_i[i] = 1;
            }
            u_n_est[i] = h_i[i];
        }
    }
//    for (int j = 0; j < size; j++) {
//        for (int k = 0; k < N; k++) {
//            cout << "level:" << j << " ,i:" << k << " = " << cache[j][k] << "," << (isCache[j][k] ? "true" : "false")<< endl;
//        }
//    }
    return u_n_est;
}

int main(void) {
    int i = 1;
    int u_Ac[N-K] = {0};
    int u_A[K] = {0};
    int A[K] = {0};

    double temp[N] = {0.0};
    probErrBound(temp);
    defineFixedAndFree(u_Ac, u_A);
    VectorXi u_n(N);

    //処理時間計測//
    const auto startTime = chrono::system_clock::now();
    u_n = generateUi(2, u_n, u_Ac, A);

    VectorXi x_n = encoder(N, u_n);
    VectorXi y_n = channel(x_n);


//    double W_i = calcW_i(i, N, u_n, u_n[i-1], y_n);
    VectorXi u_n_est = decoder(y_n, u_n, u_Ac, u_A);
    string filename = "/Users/ryotaro/labo/log";
    ofstream log;
    log.open(filename, ios::app);

    cout << "error　probability:" << errorRate(u_n,u_n_est) << endl;
    cout << "rate:" << (double)K/N << endl;
    log << "==================================================" << endl;
    log << "(N,K) = (" << N << "," << K << ")" << endl;
    log << "error　probability:" << errorRate(u_n,u_n_est) << endl;
    log << "rate:" << (double)K/N << endl;

    //処理時間計測//
    const auto endTime = chrono::system_clock::now();
    const auto timeSpan = endTime - startTime;
    cout << "総LR計算時間:" << hogetime << "[ms]" << endl;
    cout << "処理時間:" << chrono::duration_cast<chrono::milliseconds>(timeSpan).count() << "[ms]" << endl;
    log << "総LR計算時間:" << hogetime << "[ms]" << endl;
    log << "処理時間:" << chrono::duration_cast<chrono::milliseconds>(timeSpan).count() << "[ms]" << endl;
    const auto astartTime = chrono::system_clock::now();
    cout << hoge << endl;
    cout << hoge2 << endl;
    log << hoge << endl;
    log << hoge2 << endl;
    log << "==================================================" << endl;

    return 0;
}