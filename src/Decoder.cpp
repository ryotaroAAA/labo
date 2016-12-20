#include "../lib/Decoder.h"
vector<int> Decoder::decode(vector<int> &y, vector<int> &u, vector<int> &Ac, vector<int> &A){
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
    return u_n_est;
}

double Decoder::calcL_i(int i, int n ,int cache_i,int level ,vector<int> &y ,vector<int> &u, int u_i_est, vector<vector<bool> > &isCache , vector<vector<double> > &cache) {
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
