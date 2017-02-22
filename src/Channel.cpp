#include "../lib/Channel.h"
vector<int> Channel::channel_output(vector<int> &input){
    vector<int> y;
    for (auto val : input) {
        if((double)rand() / RAND_MAX < Params::get_e()){
            y.push_back(2);
        } else {
            y.push_back(val);
        }
    }
    return y;
}

double Channel::calcW(int y, int x) {
    double retVal = 0.0;
    if (x == y) {
        retVal = 1 - Params::get_e();
    } else if (y == 2) {
        retVal = Params::get_e();
    } else {
        retVal = 0;
    }
    return retVal;
}

double Channel::calcW_i(int i, int n, vector<int> &u, int u_i, vector<int> &y) {
    double W_i = 0.0;

    if ( n == 2 ){
        if ( i == 1 ) {
            W_i = 0.5 * (Channel::calcW(y[0] ,(u_i+0) % 2) * Channel::calcW(y[1],0)
                         + Channel::calcW(y[0] ,(u_i+1) % 2) * Channel::calcW(y[1],1));
        } else {
            W_i = 0.5 * Channel::calcW(y[0] ,(u[0]+u_i) % 2) * Channel::calcW(y[1],u_i);
        }
    } else {
        vector<int> tempY1(0);
        vector<int> tempY2(0);

        for (int j = 0; j < n ; j++) {
            (j < n/2) ? tempY1.push_back(y[j]) : tempY2.push_back(y[j]);
        }

        vector<int> tempU;
        vector<int> tempU_bin;
        vector<int> u_e = Common::index_e(u);
        vector<int> u_o = Common::index_o(u);

        int k = 0;
        for(auto val : u_e){
            tempU[k] = u_e[k] + u_o[k];
            k++;
        }
        tempU_bin = Common::retBinary(tempU);

        if ( i % 2 == 0 ) {
            W_i = 0.5 * (Channel::calcW_i(i/2, n/2, tempU_bin ,(1 + u_i) % 2,tempY1) * calcW_i(i/2, n/2, u_e, 1, tempY2)
                         + Channel::calcW_i(i/2, n/2, tempU_bin ,(0 + u_i) % 2,tempY1) * calcW_i(i/2, n/2, u_e, 0, tempY2));
        } else {
            W_i = 0.5 * Channel::calcW_i((i-1)/2, n/2, tempU_bin ,(u_i+u[i-1]) % 2, tempY1)
                  * Channel::calcW_i( i/2   , n/2, u_e       , u_i            , tempY2);
        }
    }
    return W_i;
}