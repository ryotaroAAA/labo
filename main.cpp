#include "main.h"

vector<int> encoder(int n, vector<int> &input){
    vector<int> x_n(n,0);
    int s_n[n];
    int v_n[n];

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
        vector<int> tempV_n1(0);
        vector<int> tempV_n2(0);

        for (int i = 0; i < n ; i++) {
            (i < n/2) ? tempV_n1.push_back(v_n[i]) : tempV_n2.push_back(v_n[i]);
        }
        hoge++;
        vector<int> tempX_n1 = encoder(n/2, tempV_n1);
        vector<int> tempX_n2 = encoder(n/2, tempV_n2);

        for (int i = 0; i < n ; i++) {
            x_n[i] = (i < n/2) ? tempX_n1[i] : tempX_n2[i - n/2];
        }
    }
    return x_n;
}

vector<int> channel(vector<int> &input){
    vector<int> y(N);
    for (auto val : input) {
        y.push_back(val);
        if((double)rand() / RAND_MAX < e){
            y.push_back(2);
        } else {
            y.push_back(val);
        }
    }
    return y;
}

vector<int> decoder(vector<int> &y, vector<int> &u, vector<int> &Ac, vector<int> &A){
    vector<double> h_i(N);
    vector<int> u_n_est(N);
    int size = log2(N);

    vector<vector<bool> > isCache (size, vector<bool>(N,false));
    vector<vector<double> > cache (size, vector<double>(N,0.0));
    double lr = 1.0;
    int cache_i = 0;

    //u_n_est計算
    for (int i = 0; i < N; i++) {
        // Acに含まれるindexなら既知
        if (containNumInArray(i, N-K, Ac)) {
            u_n_est[i] = u[i];
        } else {
            cout << i << endl;
            //処理時間計測//
            const auto startTime = chrono::system_clock::now();
            cache_i = i;
            lr = calcL_i(i, N, cache_i, 0, y, u, u[i], isCache, cache);

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

double calcL_i(int i, int n ,int cache_i,int level ,vector<int> &y ,vector<int> &u, int u_i_est, vector<vector<bool> > &isCache , vector<vector<double> > &cache) {
    double lr = 0.0;
    if ( n == 1 ) {
        double wc = calcW(y[0],0);
        double wp = calcW(y[0],1);
        lr = wc / wp;
    } else {
        vector<int> tempY1(n/2);
        vector<int> tempY2(n/2);

        for (int j = 0; j < n ; j++) {
            (j < n/2) ? tempY1.push_back(y[j]) : tempY2.push_back(y[j]);
        }

        vector<int> tempU(n/2);
        vector<int> tempU_bin(n/2);
        vector<int> u_e = index_e(u);
        vector<int> u_o = index_o(u);

        int k = 0;
        for(auto val : u_e){
            tempU[k] = u_e[k] + u_o[k];
            k++;
        }
        tempU_bin = retBinary(tempU);
        hoge2++;
        double temp1 = 1.0;
        double temp2 = 1.0;
        temp1 = calcL_i(i/2, n/2, cache_i, level+1, tempY1, tempU_bin, u_i_est, isCache, cache);
        temp2 = calcL_i(i/2, n/2, cache_i, level+1, tempY2, u_e, u_i_est, isCache, cache);

        if(((cache_i >> (int)(log2(n)-1))%2) != 1){
            //横の辺を作る
            if (isCache[level][cache_i]) {
                temp1 = cache[level][cache_i];
            } else {
                temp1 = calcL_i(i/2, n/2, cache_i, level+1, tempY1, tempU_bin, u_i_est, isCache, cache);
                isCache[level][cache_i] = true;
                cache[level][cache_i] = temp1;
            }
            //斜め下の辺を作る
            if (isCache[level][cache_i+(n/2)]) {
                temp2 = cache[level][cache_i+(n/2)];
            } else {
                temp2 = calcL_i(i/2, n/2, cache_i+(n/2), level+1, tempY2, u_e, u_i_est, isCache, cache);
                isCache[level][cache_i+(n/2)] = true;
                cache[level][cache_i+(n/2)] = temp2;
            }
        } else {
            //斜め上の辺を作る
            if (isCache[level][cache_i-(n/2)]) {
                temp1 = cache[level][cache_i-(n/2)];
            } else {
                temp1 = calcL_i(i/2, n/2, cache_i-(n/2), level+1, tempY1, tempU_bin, u_i_est, isCache, cache);
                isCache[level][cache_i-(n/2)] = true;
                cache[level][cache_i-(n/2)] = temp1;
            }
            //横の辺を作る
            if (isCache[level][cache_i]) {
                temp2 = cache[level][cache_i];
            } else {
                temp2 = calcL_i(i/2, n/2, cache_i, level+1, tempY2, u_e, u_i_est, isCache, cache);
                isCache[level][cache_i] = true;
                cache[level][cache_i] = temp2;
            }
        }

        if ( i % 2 == 0) {
            lr = ( 1 + temp1 * temp2 ) / ( temp1 + temp2 );
        } else {
            lr = pow(temp1, 1-2*u[i-1]) * temp2;
        }
    }
    if (isinf(lr) || isnan(lr)) {
        lr = 1;
    }
//    cout << i << " " << n << " " << lr  << endl;
    return lr;
}


int main(void) {
    long int n = pow(2,8);
    calcBlockErrorRate("ord", n, 1);
//    vector<int> u{1,0,0,1};
//
//    vector<int> tempU(n);
//    vector<int> tempU_bin(n);
//    vector<int> u_e = index_e(u);
//    vector<int> u_o = index_o(u);
//
//    int k = 0;
//    for(auto val : u_e){
//        tempU[k] = u_e[k] + u_o[k];
//        k++;
//    }
//    tempU_bin = retBinary(tempU);
//
//    for(auto val : tempU){
//        cout << val;
//    }

//    vector<int> u_Ac(N,0);
//    vector<int> u_A(N,0);
//    vector<int> A(N,0);
//
//    defineFixedAndFree(N, u_Ac, u_A);
//    vector<int> u_n(N,0);
//
//    u_n = generateUi(2, u_n, u_Ac, A);
//    for(auto val: u_n){
//        cout << val << endl;
//    }
//    cout << "======================" << endl;
//    vector<int> x_n = encoder(N, u_n);
//    for(auto val: x_n){
//        cout << val << endl;
//    }
//    cout << "======================" << endl;
//    vector<int> y_n = channel(x_n);
//    for(auto val: y_n){
//        cout << val << endl;
//    }
//    cout << "======================" << endl;



//    vector<int> u_n_est = decoder(y_n, u_n, u_Ac, u_A);
//    for(auto val: u_n_est){
//        cout << val << endl;
//    }
//    cout << "======================" << endl;
    return 0;
}