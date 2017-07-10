#include "../lib/Channel.h"
vector<int> Channel::channel_output(vector<int> &input){
    vector<int> y;
    for (auto val : input) {
        if((double)rand() / RAND_MAX < Params::get_e()){
            if( Params::get_s() == BEC){
                y.push_back(2);
            } else if (Params::get_s() == BSC) {
                y.push_back(val?0:1);
            }
        } else {
            y.push_back(val);
        }
    }
    return y;
}

double Channel::calcW(int y, int x) {
    double retVal = 0.0;
    if (y == x) {
        retVal = 1 - Params::get_e();
    } else if (y == 2) {
        retVal = Params::get_e();
    } else {
        if(Params::get_s() == BEC){
            retVal = 0.00001;
        } else if(Params::get_s() == BSC){
            retVal = Params::get_e();
        }
    }
    return retVal;
}

double Channel::calcW_i(int i, int n, vector<int> &u, int u_i, vector<int> &y) {
    double W_i = 0.0;

    if ( n == 2 ){
        if ( i == 1 ) {
            W_i = 0.5 * (Channel::calcW(y[0] ,(u_i+0) % 2) * Channel::calcW(y[1],0)
                         + Channel::calcW(y[0] ,(u_i+1) % 2) * Channel::calcW(y[1],1));
            cout << "W[" << i << "][" << n << "] = 0.5*W(" << y[0] << "|" << (u_i+0) % 2 << ")W(" << y[1] << "|" << 0 << ")"
                 << " + 0.5*W(" << y[0] << "|" << (u_i+1) % 2 << ")W(" << y[1] << "|" << 1 << ") = "
                 << "0.5*" << Channel::calcW(y[0] ,(u_i+0) % 2) << "*" << Channel::calcW(y[1],0)
                 << " + 0.5*" << Channel::calcW(y[0] ,(u_i+1) % 2) << "*" << Channel::calcW(y[1],1) << " = " << W_i << endl;
        } else {
            W_i = 0.5 * Channel::calcW(y[0] ,(u[0]+u_i) % 2) * Channel::calcW(y[1],u_i);
            cout << "W[" << i << "][" << n << "] = 0.5*W(" << y[0] << "|" << (u[0]+u_i) % 2 << ")W(" << y[1] << "|" << u_i << ") = "
                 << "0.5*"<< Channel::calcW(y[0] ,(u[0]+u_i) % 2) << "*" << Channel::calcW(y[1],u_i) << " = " << W_i << endl;
        }
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

//        cout << "[" << i << "][" << n << "]" << " " << u.size()-1 << " " << size_u_eo <<endl;
//        cout << "[" << i << "][" << n << "]"  << " " << u[size_u_eo] <<endl;
//        cout << size_u_eo <<endl;

        double t1=0.0;
        double t2=0.0;
        double t3=0.0;
        double t4=0.0;
        if ( i % 2 == 1 ) {
            t1 = calcW_i(temp_i, n/2, tempU_bin ,(1 + u_i) % 2,tempY1);
            t2 = calcW_i(temp_i, n/2, u_e, 1, tempY2);
            t3 = calcW_i(temp_i, n/2, tempU_bin ,(0 + u_i) % 2,tempY1);
            t4 = calcW_i(temp_i, n/2, u_e, 0, tempY2);
            W_i = 0.5 * (t1*t2 + t3*t4);
            cout << "W[" << i << "][" << n << "] = 0.5*(" << t1 << "*" << t2 << " + " << t3 << "*" << t4 << ") = " << W_i << endl;
        } else {
            t1 = calcW_i(temp_i, n/2, tempU_bin ,(u_i+u[size_u_eo-1]) % 2, tempY1);
            t2 = calcW_i(temp_i, n/2, u_e, u_i, tempY2);
            W_i = 0.5 * t1 * t2;
            cout << "W[" << i << "][" << n << "] = 0.5*W[" << temp_i << "][" << n/2 << "]*W[" << temp_i << "][" << n/2 << "] = 0.5*" << t1 << "*" << t2 << " = " << W_i << endl;
        }
    }
    if (isinf(W_i) || isnan(W_i)) {
        W_i = 1.0;
    }
    return W_i;
}
