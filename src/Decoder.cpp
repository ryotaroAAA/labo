#include "../lib/Decoder.h"
Decoder::Decoder(){

}

Decoder::~Decoder(){

}

vector<int> Decoder::decode(vector<int> &y, vector<int> &u, vector<int> &x, vector<int> &A){
    vector<double> h_i(Params::get_N());
    vector<int> u_n_est(Params::get_N());
    int size = log2(Params::get_N());

    vector<vector<bool> > isCache (size, vector<bool>(Params::get_N(),false));
    vector<vector<double> > cache (size, vector<double>(Params::get_N(),0.0));

    double lr;
    int cache_i = 0;
    Common::bar();

    //u_n_est計算
    for (int i = 0; i < Params::get_N(); i++) {
        // Aに含まれないindexなら既知
        if (Common::containNumInArray(i, Params::get_N()-Params::get_K(), A) == false) {
            u_n_est[i] = u[i];
        } else {
//            cout << "------------------" << endl;
            cout << i+1 << endl;
            this->startTimer();

            cache_i = i;
            lr = calcL_i(i+1, Params::get_N(), cache_i, 0, y, u_n_est, isCache, cache);

            this->stopTimer();
            this->outTime();
//            cout << lr <<endl;

            if (lr >= 1.0) {
                h_i[i] = 0;
            } else {
                h_i[i] = 1;
            }
            u_n_est[i] = h_i[i];
        }
    }
    return u_n_est;
}

double Decoder::calcL_i(int i, int n ,int cache_i,int level ,vector<int> &y ,vector<int> &u, vector<vector<bool> > &isCache , vector<vector<double> > &cache) {
    double lr = 0.0;
    this->addCount();
    if ( n == 1 ) {
        double wc = Channel::calcW(y[0],0);
        double wp = Channel::calcW(y[0],1);
        lr = 1.0 * wc / wp;
    } else {
        vector<int> tempY1(n/2);
        vector<int> tempY2(n/2);
        for (int j = 0; j < n ; j++) {
            if(j < n/2){
                tempY1[j] = y[j];
            } else {
                tempY2[j-n/2] = y[j];
            }
        }

        int size_u_eo = (i % 2) == 0 ? i-2 : i-1;
        int size_aug_u = size_u_eo/2;

        vector<int> tempU(size_aug_u);
        vector<int> tempU_bin(size_aug_u);
        vector<int> u_e(size_aug_u);
        vector<int> u_o(size_aug_u);

        u_e = Common::index_e(size_u_eo, u);
        u_o = Common::index_o(size_u_eo, u);

        for(int k=0; k<size_aug_u; k++){
            tempU[k] = u_e[k] + u_o[k];
        }

        tempU_bin = Common::retBinary(size_aug_u, tempU);

        double temp1 = 0.0;
        double temp2 = 0.0;

        int temp_i = (i % 2 == 1) ? (i+1)/2 : i/2;

        if ( 1 ) {
            //通常盤
            temp1 = calcL_i(temp_i, n / 2, cache_i, level + 1, tempY1, tempU_bin, isCache, cache);
            temp2 = calcL_i(temp_i, n / 2, cache_i, level + 1, tempY2, u_e, isCache, cache);
        } else {
            //効率的な実装(要テスト)
            if (((cache_i >> (int) (log2(n) - 1)) % 2) != 1) {
                //横の辺を作る
                if (isCache[level][cache_i]) {
                    temp1 = cache[level][cache_i];
                } else {
                    temp1 = calcL_i(temp_i, n / 2, cache_i, level + 1, tempY1, tempU_bin, isCache, cache);
                    isCache[level][cache_i] = true;
                    cache[level][cache_i] = temp1;
                }
                //斜め下の辺を作る
                if (isCache[level][cache_i + (n / 2)]) {
                    temp2 = cache[level][cache_i + (n / 2)];
                } else {
                    temp2 = calcL_i(temp_i, n / 2, cache_i + (n / 2), level + 1, tempY2, u_e, isCache, cache);
                    isCache[level][cache_i + (n / 2)] = true;
                    cache[level][cache_i + (n / 2)] = temp2;
                }
            } else {
                //斜め上の辺を作る
                if (isCache[level][cache_i - (n / 2)]) {
                    temp1 = cache[level][cache_i - (n / 2)];
                } else {
                    temp1 = calcL_i(temp_i, n / 2, cache_i - (n / 2), level + 1, tempY1, tempU_bin, isCache,
                                    cache);
                    isCache[level][cache_i - (n / 2)] = true;
                    cache[level][cache_i - (n / 2)] = temp1;
                }
                //横の辺を作る
                if (isCache[level][cache_i]) {
                    temp2 = cache[level][cache_i];
                } else {
                    temp2 = calcL_i(temp_i, n / 2, cache_i, level + 1, tempY2, u_e, isCache, cache);
                    isCache[level][cache_i] = true;
                    cache[level][cache_i] = temp2;
                }
            }
        }
//        cout << "[" << pow(temp1, 1-2*u[2*(i-1)-1]) << "]" << "[" << temp2 << "]" << endl;
//        cout << ">>[" << temp1 << "][" << temp2 << "]" << endl;

        if ( i % 2 == 1) {
            //odd index
            lr = ( 1.0 + temp1 * temp2 ) / ( temp1 + temp2 );

            if (isinf(lr) || isnan(lr)) {
                lr = 10.0;
            }

//            cout << "(a*b+1)/(a+b) =" << "(" << temp1 << "*" << temp2 << "+1)/(" << temp1 << "+" << temp2 << ")" << " = " << lr << endl;
        } else {
            //even index
            lr = pow(temp1, 1-2*u[size_u_eo]) * temp2;

//            if( u[size_u_eo] == 1){
//                cout << "u[last] = 1 ::" << "b/a =" << temp2 << "/" << temp1 << " = " << lr << endl;
//            } else {
//                cout << "u[last] = 0 ::" << "b*a =" << temp1 << "*" << temp2 << " = " << lr << endl;
//            }
//            cout << "  " << size_u_eo << "  " << u[size_u_eo] << endl;
//            cout << "    " << pow(temp1, 1-2*u[size_u_eo]) << "    " << endl;
        }
//        cout << "[" << i << "][" << n << "]" <<  ">> " << lr << endl;
//        cout << "~~~~~~~~~~" << endl;
    }

    if (isinf(lr) || isnan(lr)) {
         lr = 10.0;
    }
    return lr;
}
