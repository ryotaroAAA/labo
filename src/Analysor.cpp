#include "../lib/Analysor.h"

double Analysor::errorRate(vector<int> &u, vector<int> &u_est){
    int error_count = 0;
    for(int i=0; i<Params::N; i++){
        if(u[i] != u_est[i]){
            error_count++;
        }
    }
    return (double)error_count/Params::N;
}
double Analysor::calcCapacityForBec(int i, int n) {
    double cap =0.0;
    if ( i == 0 && n == 1 ) {
        cap = 1 - Params::e;
    } else {
        if ( i % 2 == 0) {
            cap = pow(Analysor::calcCapacityForBec(i/2, n/2),2);
        } else {
            double tempCap = Analysor::calcCapacityForBec((i-1)/2, n/2);
            cap = 2 * tempCap - pow(tempCap,2);
        }
    }
    return cap;
}

double Analysor::calcBhatForBec(int i, int n){
    double bha =0.0;
    if ( i == 0 && n == 1 ) {
        bha = Params::e;
    } else {
        if ( i % 2 == 0) {
            double tempBha = Analysor::calcBhatForBec(i/2, n/2);
            bha = 2 * tempBha - pow(tempBha,2);
        } else {
            bha = pow(Analysor::calcBhatForBec((i-1)/2, n/2),2);
        }
    }
    return bha;
}

void Analysor::makeArrayCapacityForBec(vector<double> &array) {
    for(int i = 0; i < Params::N; i++) {
        array.push_back(Analysor::calcCapacityForBec(i, Params::N));
    }
}

void Analysor::probErrBound(vector<double> &array) {
    int count = 0;
    double tempArr[Params::N];
    double sumArr[Params::N];
    double sum = 0.0;

    for(int i = 0; i < Params::N; i++) {
        tempArr[i] = Analysor::calcBhatForBec(i, Params::N);
    }

    qsort(tempArr, Params::N, sizeof(double), Common::compare_asc);

    for(int i=0; i < Params::N; i++){
        sum += tempArr[i];
        sumArr[i] = sum;
    }

    string filename = "e18";
    ofstream w_file;
    w_file.open(filename, ios::out);
    for (int i = 0; i < Params::N; i++)
    {
        w_file << (double)i/Params::N << " " << sumArr[i] << " " << tempArr[i] << endl;
    }
}

void Analysor::calcBlockErrorRate(MODE mode, int n, double dist) {
//    int tempK = 0;
//
//    Performance performance;
//    Decoder decoder;
//    Encoder encoder;
//    Logger logger;
//
//    int eachK = Params::N*dist;

//    for (int k=0; k<Params::N ; k + eachK) {
//        vector<int> u_Ac(Params::N, 0);
//        vector<int> u_A(Params::N, 0);
//        vector<int> A(Params::N, 0);
//        vector<int> u_n(Params::N, 0);
//        vector<int> x_n(Params::N, 0);
//        vector<int> y_n(Params::N, 0);
//        vector<int> u_est(Params::N, 0);
//
//        Preseter::preset(u_n, u_Ac, u_A);
//        performance.startTimer();
//        x_n = encoder.encode(Params::N, u_n);
//        y_n = Channel::channel_output(x_n);
//        u_est = decoder.decode(y_n, u_n, u_Ac, u_A);
//
//        double BER = Analysor::errorRate(u_n, u_est);
//        double rate = (double) Params::K / Params::N;
//
//        logger.outLogTime(performance.outTime("処理時間"));
//        logger.outLog("(N,K) = (" + to_string(Params::N) + "," + to_string(Params::K) + ")");
//        logger.outLog("error　probability:" + to_string(BER));
//        logger.outLog("rate:" + to_string(rate));
//        logger.outLogCount(encoder.outCount("encoder_count"));
//        logger.outLogCount(decoder.outCount("decoder_count"));
//
//        logger.outLog("================================");
//        logger.outLog("BER:" + to_string(BER));
//        logger.outLog("Rate:" + to_string(rate));
//
//        performance.stopTimer();
//
//        logger.outLog(performance.outTime("処理時間"));
//        logger.outLog(encoder.outCount("encoder_count"));
//        logger.outLog(decoder.outCount("decoder_count"));
//        logger.outLogRVB(rate,BER);
//    }
}