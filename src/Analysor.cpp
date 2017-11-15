#include "../lib/Analysor.h"

Analysor::Analysor(){

}

Analysor::~Analysor(){

}

void Analysor::errorCount(vector<int> &u, vector<int> &u_est, int* error_count){
    for(int i=0; i<Params::get_N(); i++){
        if(u[i] != u_est[i]){
            *error_count = *error_count + 1;
        }
    }
}

double Analysor::errorCalc(vector<int> &u, vector<int> &u_est, int* error_count){
    for(int i=0; i<Params::get_N(); i++){
        if(u[i] != u_est[i]){
            *error_count = *error_count + 1;
        }
    }

    return (double)(*error_count)/Params::get_N();
}

double Analysor::calcCapacity(int i, int n) {
    double cap =0.0;
    if(Params::get_s() == BEC){
        if ( i == 0 && n == 1 ) {
            cap = 1 - Params::get_e();
        } else {
            if ( i % 2 == 0) {
                cap = pow(Analysor::calcCapacity(i/2, n/2),2);
            } else {
                double tempCap = Analysor::calcCapacity((i-1)/2, n/2);
                cap = 2 * tempCap - pow(tempCap,2);
            }
        }
    } else if (Params::get_s() == BSC){

    }
    return cap;
}

double Analysor::calcBhatforBEC(int i, int n){
    double bha =0.0;

    if ( i == 0 && n == 1 ) {
        bha = Params::get_e();
    } else {
        if ( i % 2 == 0) {
            double tempBha = Analysor::calcBhatforBEC(i/2, n/2);
            bha = 2 * tempBha - pow(tempBha, 2);
        } else {
            bha = pow(Analysor::calcBhatforBEC((i-1)/2, n/2),2);
        }
    }
    return bha;
}

void Analysor::makeArrayCapacity(vector<double> &array) {
    for(int i = 0; i < Params::get_N(); i++) {
        array[i] = Analysor::calcCapacity(i, Params::get_N());
    }
}

double Analysor::calc_inv(double y){
    double x1 = 10.0;
    double x2 = 500.0;
    double x3 = 0.0;
    double fx = 0.0;
    int count = 1;
    while (abs(x1 - x2) >= 0.00001) {
        x3 = (x1+x2) / 2.0;
        fx = sqrt(M_PI/x3) * exp(-1.0 * (x3/4.0)) * (1.0 - 10.0/(7.0 * x3));
        if(fx > y){
            x1 = x3;
        } else {
            x2 = x3;
        }
        count++;
//        cout << "est:" <<  x1 << " " << x2 << " " << m_func(x3)<< " true:" << y << endl;
    }
    return x3;
}

double Analysor::inv_m_func(double y){
    double ret = 0.0;
    const double a = -0.4527;
    const double b = 0.0218;
    const double c = 0.86;
    double temp = 0.0;
    //0.03847)
    if (y > 0.03847){
        ret = pow((log(y)-b)/a, 1.0/c);
    } else if(y < 0.03847) {
        ret = calc_inv(y);
    } else {
//        ret = __nan();
    }
    return ret;
}

double Analysor::m_func(double x){
    double ret = 0.0;
    const double a = -0.4527;
    const double b = 0.0218;
    const double c = 0.86;
    double temp1 = 0.0;
    double temp2 = 0.0;
    if (x < 10.0) {
        ret = exp(a * pow(x,c) + b);
    } else if (x > 10.0) {
        temp1 = 1.0 - 3.0/x;
        temp2 = 1.0 + 1.0/(7.0 * x);
        ret = 0.5 * sqrt(M_PI/x) * exp(-1.0 * x / 4.0)  * (temp1 + temp2);
    } else {
//        ret = __nan();
    }
    return ret;
}

double Analysor::calc_m_in(int i, int n){
    double ret = 0.0;
    double temp1 = 0.0;
    double temp2 = 0.0;
    int temp_i = (i/2 == 0)? 1 : i/2;
    if (n == 1){
        ret = 2.0 / Params::get_e();
    } else {
        temp1 = Analysor::calc_m_in(temp_i, n/2);
        if (i % 2 == 0){
            ret = 2.0 * temp1;
        } else {
            temp2 = 1.0 - pow((1.0 - Analysor::m_func(temp1)), 2);
            ret = Analysor::inv_m_func(temp2);
        }
    }
//    cout << "[" << temp_i << "][" << n  << "] "  << ret << endl;
    return ret;
}

void Analysor::makeArrayBhat(vector<double> &array) {
    if ( Params::get_s() == BSC) {
        Performance performance;
        Decoder decoder;
        Encoder encoder;
        Logger logger;

        vector<int> u(Params::get_N(), 0);
        vector<int> x(Params::get_N(), 0);
        vector<double> y(Params::get_N(), 0.0);
        int u_n = 0;

        double sumbha = 0.0;
        double tempbha = 0.0;
        int size = log2(Params::get_N());

        vector<vector<bool> > isCache (size, vector<bool>(Params::get_N(),false));
        vector<vector<double> > cache (size, vector<double>(Params::get_N(),0.0));

        double lr;
        int cache_i = 0;
        int repeatNum = Params::get_monteNum();

        for (int m = 1; m <= repeatNum; m++) {
            performance.startTimer();
//            cout << "[" << m << "]" << endl;
            Preseter::preset_u(RAND, u);
            vector<int> temp_u = u;
            x = encoder.encode(Params::get_N(), u);

            for (int i = 0; i < Params::get_N(); i++) {
                temp_u = u;
                u_n = temp_u[i];
                temp_u.resize(i);

                y = Channel::channel_output(x);
                cache_i = decoder.makeTreeIndex(Params::get_N())[i] - 1;
                lr = exp(decoder.calcL_i(i+1, Params::get_N(), cache_i, 0, y, temp_u, isCache, cache));
                lr = pow(lr, -1+2 * u_n);
                tempbha = sqrt(lr);
                if (tempbha > 1.0) tempbha = 1.0/tempbha;
                if (isnan(tempbha)) tempbha = 0.0;
                array[i] += 1.0 * tempbha/repeatNum;
            }
        }
    } else {
        vector<double> m_in(Params::get_N(), 0.0);
        vector<double> sigma2_in(Params::get_N(), 0.0);
        double temp = 0.0;
        for (int i = 0; i < Params::get_N(); i++) {
            temp = Analysor::calc_m_in(i+1, Params::get_N());
//            if(isnan(temp)){
//                m_in[i] = 0.0;
//            } else {
                m_in[i] = temp;
//            }
//            cout << i << " mean " << m_in[i] <<endl;
            sigma2_in[i] = 2.0/m_in[i];
            array[i] = exp(-1.0/(2.0 * sigma2_in[i]));
        }
    }

    string filename = Params::get_rvbDir() + " Bhat";
    ofstream w_file;
    w_file.open(filename, ios::out);
    for (int i = 0; i < Params::get_N(); i++)
    {
        w_file << i << " " << array[i] << endl;
        cout << i << " " << array[i] << endl;
    }
    Common::bar();
}

void Analysor::probErrBound() {
    int count = 0;
    vector<double> tempArr;
    vector<double> sumArr;
    double sum = 0.0;

    for(int i = 0; i < Params::get_N(); i++) {
        tempArr.push_back(Analysor::calcBhatforBEC(i, Params::get_N()));
    }

    sort(tempArr.begin(), tempArr.end(), greater<int>());

    for(int i=0; i < Params::get_N(); i++){
        sum += tempArr[i];
        sumArr.push_back(sum);
    }

    string filename = "log/eb10";
    ofstream w_file;
    w_file.open(filename, ios::out);
    for (int i = 0; i < Params::get_N(); i++)
    {
        w_file << (double)i/Params::get_N() << " " << sumArr[i] << " " << tempArr[i] << endl;
        cout << (double)i/Params::get_N() << " " << sumArr[i] << " " << tempArr[i] << endl;
    }
}


string Analysor::get_itrfn(){
    time_t now = time(NULL);
    struct tm *pnow = localtime(&now);
    stringstream itrfn;
    EXP_MODE em = Params::get_exp_mode();
    MID_MODE mm = Params::get_m_mode();
    string ename;
    switch (em) {
        case NORMAL: ename = ""; break;
        case PUNC: ename = "punc"; break;
        case MID:
            switch (mm) {
                case MID_IU: ename = "mid_iu"; break;
                case MID_ID: ename = "mid_id"; break;
                case MID_BU: ename = "mid_bu"; break;
                case MID_BD: ename = "mid_bd"; break;
            }
            break;
        case QUP: ename = "qup"; break;
        case WANG: ename = "wang"; break;
        case M_QUP: ename = "m_qup"; break;
        case M_WANG: ename = "m_wang"; break;
        case VALERIO_P: ename = "valerio_p"; break;
        case VALERIO_S: ename = "valerio_s"; break;
        case M_VALERIO_P: ename = "m_valerio_p"; break;
        case M_VALERIO_S: ename = "m_valerio_s"; break;
    }
    itrfn << "log/"
       << to_string(pnow->tm_year + 1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday)
       << "/itr_N" << to_string(Params::get_N())
       << "e" << Params::get_e()
       << "BNum" << to_string(Params::get_blockNum())
       << "UPBn" << to_string(Params::get_upperBlockErrorNum())
       << "c" << (Params::get_s() ? "BSC" : "BEC") << ":" << ename;
    return itrfn.str();
}

//rate vs BERのグラフ作成用
void Analysor::calcBlockErrorRate() {
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    logger.setRvbDir(Params::get_rvbDir());
    vector<int> A(Params::get_K(), -1);
    vector<int> Ac(Params::get_N()-Params::get_K(), -1);
    vector<int> u_n(Params::get_N(), 0);
    vector<int> x_n(Params::get_N(), 0);
    vector<double> y_n(Params::get_N(), 0);

    vector<int> u_est(Params::get_N(), 0);
    vector<pair<int, double> > cap_map;
    Preseter::makeMutualInfoArray(cap_map);

    double BER = 0.0;
    double sumBER = 0.0;
    double rate = 0.0;
    int n = Params::get_N();
    int error_count = 0;
    int block_error_count = 0;
    int loopi = 0;
    int repeatNum = Params::get_blockNum();

    int tmp = Params::get_N()/Params::get_K();
    int tmpK = Params::get_K();
    for (int i = 1; i <= (tmp+1)/2; i++) {

        performance.startTimer();
        Params::set_K(i * tmpK);
        A.resize(Params::get_K(), -1);
        Preseter::represet_A(A, Ac, cap_map);
        Preseter::preset_u(RAND, u_n);
        x_n = encoder.encode(Params::get_N(), u_n);

        while (block_error_count < Params::get_upperBlockErrorNum()) {
            loopi++;

            y_n.assign(Params::get_N(), 0);
            u_est.assign(Params::get_N(), 0);
            y_n = Channel::channel_output(x_n);
//            u_est = (Params::get_decode_mode() == BP)?decoder.BP(Params::get_rp(), y_n, u_n, A):decoder.decode(y_n, u_n, A);

//            Common::pp(u_est);
            Analysor::errorCount(u_n, u_est, &error_count);
            if(error_count > 0) block_error_count++;
            if (loopi % 1 == 0 ) cout << loopi << " " << error_count << " " << block_error_count  << (double)block_error_count/loopi << endl;
            error_count = 0;
            rate = (double) Params::get_K() / Params::get_N();

//            if(loopi >= repeatNum || i * tmpK < Params::get_N()/7) break;
            if(loopi >= repeatNum) break;
        }
        performance.stopTimer();

        cout << Params::get_rvbDir() << endl;
        BER = (double)block_error_count/loopi;

        logger.outLog("=================================");
        logger.outLog(performance.outTime("処理時間"));
        logger.outLog("(N,K) = (" + to_string(Params::get_N()) + "," + to_string(Params::get_K()) + ")");
        logger.outLog("BER:" + to_string(BER));
        logger.outLog("Rate:" + to_string(rate));
        logger.outLog(encoder.outCount("encoder_count"));
        logger.outLog(decoder.outCount("decoder_count"));
        performance.outHMS();

        logger.outLogRVB(rate, BER);
        loopi = 0;
        sumBER = 0.0;
        block_error_count = 0;
//        if(mode == TEST) break;
    }
}

//rate vs BERのグラフ作成用
void Analysor::calcBlockErrorRate_BP() {
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    logger.setRvbDir(Params::get_rvbDir());

    int n = Params::get_N(), error_count = 0, block_error_count = 0, loopi = 0, itr = 0;
    int tmp = Params::get_N() / Params::get_K(), tmpK = Params::get_K();
    int size = log2(Params::get_N()) + 1;
    double BER = 0.0, sumBER = 0.0, rate = 0.0, mitr = 0.0;;

    vector<int> A;
    vector<int> Ac;
    vector<int> u(Params::get_N(), 0);
    vector<int> x(Params::get_N(), 0);
    vector<double> y(Params::get_N(), 0);
    vector<vector<int> > tmp_x(log2(Params::get_N()), vector<int>(Params::get_N(), 0));
    vector<vector<int> > xm(size, vector<int>(Params::get_N(), 0));
    vector<vector<double> > ym(size, vector<double>(Params::get_N(), 0.0));
    vector<int> u_est(Params::get_N(), 0);
    vector<pair<int, double> > cap_map;

    Preseter::makeMutualInfoArray(cap_map);
    Preseter::preset_A_Ac(A, Ac);

    vector<int> param(2, 0);
    //パンクチャorショートンパターンをp
    vector<int> p_0(Params::get_M(), -1);
    vector<int> p(Params::get_M(), -1);

    ofstream itrfile;
    itrfile.open(get_itrfn(), ios::out);

    int start = 0.32/Common::get_rate();
    int end = 0.37/Common::get_rate();
    EXP_MODE em = Params::get_exp_mode();
    if(em == VALERIO_S || em == M_VALERIO_S){
        start = (Params::get_M()+(Params::get_N()-Params::get_M())*0.32)/Params::get_K();
        end = (Params::get_M()+(Params::get_N()-Params::get_M())*0.37)/Params::get_K();
    }
    cout << "start::" << start << " end::" << end << endl;

    for (int i = start; i <= end; i++) {
        loopi = 0, sumBER = 0.0, mitr = 0.0, block_error_count = 0;
        performance.startTimer();
        //k,u,frozen設定
        Params::set_K(i * tmpK);
        Preseter::preset_u(RAND, u);
        Preseter::set_params(cap_map, A, Ac, p_0, p);

        rate = Common::get_rate();
        cout << "rate" << rate << endl;

        if(Params::get_is_outlog()) {
            cout << "A";
            Common::pp(A);
            cout << "Ac";
            Common::pp(Ac);
            cout << "punc pattern";
            Common::pp(p);
            vector<int> table = Preseter::makeTable(Params::get_N());
            cout << "BN";
            Common::pp(table);
        }

        //encode
        if( Common::is_mid_send() ){
            tmp_x[log2(Params::get_N()) - 1] = encoder.encode_m(Params::get_N(), 0, 0, u, tmp_x);
            //uとxをまとめて送る
            for (int i = 0; i < log2(Params::get_N()) + 1; i++) {
                for (int j = 0; j < Params::get_N(); j++) {
                    if (i == 0) {
                        xm[i][j] = u[j];
                    } else {
                        xm[i][j] = tmp_x[i - 1][j];
                    }
                }
            }
        } else {
            x = encoder.encode(Params::get_N(), u);
        }

        //decode
        while (block_error_count < Params::get_upperBlockErrorNum()) {
            loopi++;
            u_est.assign(Params::get_N(), 0);
            if( Common::is_mid_send() ){
                Channel::channel_output_m(xm, ym);
            } else {
                y = Channel::channel_output(x);
            }
            u_est = decoder.calcBP(p, param, u, x, y, xm, ym, A, Ac);
            itr = param[0];
            Analysor::errorCount(u, u_est, &error_count);
            if (error_count > 0) block_error_count++;
            if (loopi % 1 == 0 ) cout << loopi << " " << error_count << " " << block_error_count  << " " << itr << " " << (double)block_error_count/loopi << endl;
            error_count = 0;
            mitr += itr; // 0:itr, 1:nochecked

            if(loopi >= Params::get_blockNum()) break;
        }
        mitr = (mitr/loopi)>=70 ? 70 : mitr/loopi;
        itrfile << rate << " " << mitr << endl;
        cout << "mitr : "<< mitr << endl;
        cout << "rate : "<< rate << endl;

        performance.stopTimer();

        cout << Params::get_rvbDir() << endl;
        BER = (double)block_error_count/loopi;

        logger.outLog("=================================");
        logger.outLog(performance.outTime("処理時間"));
        logger.outLog("(N,K,M,MN) = (" + to_string(Params::get_N()) + "," + to_string(Params::get_K())  + "," + to_string(Params::get_M())  + "," + to_string(Params::get_MN())+ ")");
        logger.outLog("BER:" + to_string(BER));
        logger.outLog("Rate:" + to_string(rate));
        logger.outLog(encoder.outCount("encoder_count"));
        logger.outLog(decoder.outCount("decoder_count"));
        performance.outHMS();

        logger.outLogRVB(rate, BER);
//        if(mode == TEST) break;
    }
}