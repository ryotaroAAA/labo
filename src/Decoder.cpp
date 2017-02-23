#include "../lib/Decoder.h"
Decoder::Decoder(){

}

Decoder::~Decoder(){

}

vector<int> Decoder::decode(vector<int> &y, vector<int> &u, vector<int> &Ac, vector<int> &A){
    vector<double> h_i(Params::get_N());
    vector<int> u_n_est(Params::get_N());
    int size = log2(Params::get_N());

    vector<vector<bool> > isCache (size, vector<bool>(Params::get_N(),false));
    vector<vector<double> > cache (size, vector<double>(Params::get_N(),0.0));

    double lr;
    int cache_i = 0;

    //u_n_est計算
    for (int i = 0; i < Params::get_N(); i++) {
        // Acに含まれるindexなら既知
        if (Common::containNumInArray(i, Params::get_N()-Params::get_K(), Ac)) {
            u_n_est[i] = u[i];
        } else {
            cout << i << endl;
            this->startTimer();

            cache_i = i;
            lr = calcL_i(i, Params::get_N(), cache_i, 0, y, u, u[i], isCache, cache);

            this->stopTimer();
            this->outTime();

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
    this->addCount();
    if ( n == 1 ) {
        double wc = Channel::calcW(y[0],0);
        double wp = Channel::calcW(y[0],1);
        lr = wc / wp;
    } else {
        vector<int> tempY1;
        vector<int> tempY2;

        for (int j = 0; j < n ; j++) {
            (j < n/2) ? tempY1.push_back(y[j]) : tempY2.push_back(y[j]);
        }

        vector<int> tempU(Params::get_N());
        vector<int> tempU_bin(Params::get_N());
        vector<int> u_e;
        vector<int> u_o;
        u_e = Common::index_e(n/2,u);
        u_o = Common::index_o(n/2,u);

        int k = 0;

        for(auto val : u_e){
            tempU[k] = u_e[k] + u_o[k];
            k++;
        }

        tempU_bin = Common::retBinary(n/2, tempU);

        double temp1 = 0.0;
        double temp2 = 0.0;
        temp1 = calcL_i(i/2, n/2, cache_i, level+1, tempY1, tempU_bin, u_i_est, isCache, cache);
        temp2 = calcL_i(i/2, n/2, cache_i, level+1, tempY2, u_e, u_i_est, isCache, cache);

        // if(((cache_i >> (int)(log2(n)-1))%2) != 1){
        //     //横の辺を作る
        //     if (isCache[level][cache_i]) {
        //         temp1 = cache[level][cache_i];
        //     } else {
        //         temp1 = calcL_i(i/2, n/2, cache_i, level+1, tempY1, tempU_bin, u_i_est, isCache, cache);
        //         isCache[level][cache_i] = true;
        //         cache[level][cache_i] = temp1;
        //     }
        //     //斜め下の辺を作る
        //     if (isCache[level][cache_i+(n/2)]) {
        //         temp2 = cache[level][cache_i+(n/2)];
        //     } else {
        //         temp2 = calcL_i(i/2, n/2, cache_i+(n/2), level+1, tempY2, u_e, u_i_est, isCache, cache);
        //         isCache[level][cache_i+(n/2)] = true;
        //         cache[level][cache_i+(n/2)] = temp2;
        //     }
        // } else {
        //     //斜め上の辺を作る
        //     if (isCache[level][cache_i-(n/2)]) {
        //         temp1 = cache[level][cache_i-(n/2)];
        //     } else {
        //         temp1 = calcL_i(i/2, n/2, cache_i-(n/2), level+1, tempY1, tempU_bin, u_i_est, isCache, cache);
        //         isCache[level][cache_i-(n/2)] = true;
        //         cache[level][cache_i-(n/2)] = temp1;
        //     }
        //     //横の辺を作る
        //     if (isCache[level][cache_i]) {
        //         temp2 = cache[level][cache_i];
        //     } else {
        //         temp2 = calcL_i(i/2, n/2, cache_i, level+1, tempY2, u_e, u_i_est, isCache, cache);
        //         isCache[level][cache_i] = true;
        //         cache[level][cache_i] = temp2;
        //     }
        // }

        if ( i % 2 == 0) {
            lr = ( 1 + temp1 * temp2 ) / ( temp1 + temp2 );
        } else {
            lr = pow(temp1, 1-2*u[i-1]) * temp2;
        }
    }
    if (isinf(lr) || isnan(lr)) {
         lr = 1;
    }
    return lr;
}
