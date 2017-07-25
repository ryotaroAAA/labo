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

void Analysor::makeArrayBhat(vector<double> &array) {
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;

    vector<int> u(Params::get_N(), 0);
    vector<int> x(Params::get_N(), 0);
    vector<int> y(Params::get_N(), 0);
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
        cout << "[" << m << "]" << endl;
//        u.resize(Params::get_N());
        Preseter::preset_u(RAND, u);
        vector<int> temp_u = u;
        x = encoder.encode(Params::get_N(), u);
        y = Channel::channel_output(x);
        for (int i = 0; i < Params::get_N(); i++) {
            temp_u = u;
            if(i >= 1) {
                u_n = temp_u[i];
            }
            temp_u.resize(i);
            cache_i = decoder.makeTreeIndex(Params::get_N())[i] - 1;
            lr = exp(decoder.calcL_i(i+1, Params::get_N(), cache_i, 0, y, temp_u, isCache, cache));
            lr = pow(lr, -1+2 * u_n);
            tempbha = sqrt(lr);
            if (tempbha > 1.0) tempbha = 1.0/tempbha;
            if (isnan(tempbha)) tempbha = 0.0;
//            if (tempbha > 1.0) tempbha = 1.0;
            array[i] += 1.0 * tempbha/repeatNum;
        }
    }

    string filename = "/Users/ryotaro/labo/log/test";
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

    string filename = "/Users/ryotaro/labo/log/eb10";
    ofstream w_file;
    w_file.open(filename, ios::out);
    for (int i = 0; i < Params::get_N(); i++)
    {
        w_file << (double)i/Params::get_N() << " " << sumArr[i] << " " << tempArr[i] << endl;
        cout << (double)i/Params::get_N() << " " << sumArr[i] << " " << tempArr[i] << endl;
    }
}

//rate vs BERのグラフ作成用
void Analysor::calcBlockErrorRate(MODE mode) {
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    logger.setRvbDir(Params::get_rvbDir());
    vector<int> A(Params::get_K(), -1);
    vector<int> u_n(Params::get_N(), 0);
    vector<int> x_n(Params::get_N(), 0);
    vector<int> y_n(Params::get_N(), 0);
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
    for (int i = 1; i <= tmp/2; i++) {
        performance.startTimer();
        Params::set_K(i * tmpK);
        A.resize(Params::get_K(), -1);
        Preseter::represet_A(A, cap_map);
        Preseter::preset_u(RAND, u_n);
        x_n = encoder.encode(Params::get_N(), u_n);

        while (block_error_count < Params::get_upperBlockErrorNum()) {
            loopi++;

            y_n.assign(Params::get_N(), 0);
            u_est.assign(Params::get_N(), 0);

            y_n = Channel::channel_output(x_n);
            u_est = decoder.decode(y_n, u_n, A);

            Analysor::errorCount(u_n, u_est, &error_count);
            if(error_count > 0) block_error_count++;
            if (loopi % 100 == 0 ) cout << loopi << " " << error_count << " " << block_error_count << endl;
//            if (loopi % 1000 == 0 ) cout << "block_error_count:" << block_error_count << endl;
            error_count = 0;
            rate = (double) Params::get_K() / Params::get_N();

//            if(loopi >= repeatNum || i * tmpK < Params::get_N()/5) break;
            if(loopi >= repeatNum) break;
        }
        performance.stopTimer();

        cout << Params::get_rvbDir() << endl;
//        BER = (double)sumBER/repeatNum;
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
        if(mode == TEST) break;
    }
}