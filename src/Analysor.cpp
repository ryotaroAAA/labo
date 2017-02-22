#include "../lib/Analysor.h"

Analysor::Analysor(){

}

Analysor::~Analysor(){

}

double Analysor::errorRate(vector<int> &u, vector<int> &u_est){
    int error_count = 0;
    for(int i=0; i<Params::get_N(); i++){
        if(u[i] != u_est[i]){
            error_count++;
        }
    }
    return (double)error_count/Params::get_N();
}
double Analysor::calcCapacityForBec(int i, int n) {
    double cap =0.0;
    if ( i == 0 && n == 1 ) {
        cap = 1 - Params::get_e();
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
        bha = Params::get_e();
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
    for(int i = 0; i < Params::get_N(); i++) {
        array.push_back(Analysor::calcCapacityForBec(i, Params::get_N()));
    }
}

void Analysor::probErrBound(vector<double> &array) {
    int count = 0;
    vector<double> tempArr(Params::get_N());
    vector<double> sumArr(Params::get_N());
    double sum = 0.0;

    for(int i = 0; i < Params::get_N(); i++) {
        tempArr[i] = Analysor::calcBhatForBec(i, Params::get_N());
    }

    sort(tempArr.begin(), tempArr.end(), greater<int>());

    for(int i=0; i < Params::get_N(); i++){
        sum += tempArr[i];
        sumArr[i] = sum;
    }

    string filename = "e18";
    ofstream w_file;
    w_file.open(filename, ios::out);
    for (int i = 0; i < Params::get_N(); i++)
    {
        w_file << (double)i/Params::get_N() << " " << sumArr[i] << " " << tempArr[i] << endl;
    }
}

void Analysor::calcBlockErrorRate(MODE mode, int n, double dist) {
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;

    int tmp = Params::get_N()/100;
    for (int i = 0; i < tmp; i++) {
        Params::set_K(i*5);
        vector<int> u_Ac(Params::get_N(), 0);
        vector<int> u_A(Params::get_N(), 0);
        vector<int> A(Params::get_N(), 0);
        vector<int> u_n(Params::get_N(), 0);
        vector<int> x_n(Params::get_N(), 0);
        vector<int> y_n(Params::get_N(), 0);
        vector<int> u_est(Params::get_N(), 0);

        Preseter::preset(RAND, u_n, u_Ac, u_A);
        performance.startTimer();
        x_n = encoder.encode(Params::get_N(), u_n);
        y_n = Channel::channel_output(x_n);
        u_est = decoder.decode(y_n, u_n, u_Ac, u_A);

        double BER = Analysor::errorRate(u_n, u_est);
        double rate = (double)Params::get_K() / Params::get_N();

        logger.outLogTime(performance.outTime("処理時間"));
        logger.outLog("(N,K) = (" + to_string(Params::get_N()) + "," + to_string(Params::get_K()) + ")");
        logger.outLog("error　probability:" + to_string(BER));
        logger.outLog("rate:" + to_string(rate));
        logger.outLogCount(encoder.outCount("encoder_count"));
        logger.outLogCount(decoder.outCount("decoder_count"));

        logger.outLog("================================");
        logger.outLog("BER:" + to_string(BER));
        logger.outLog("Rate:" + to_string(rate));

        performance.stopTimer();

        logger.outLog(performance.outTime("処理時間"));
        logger.outLog(encoder.outCount("encoder_count"));
        logger.outLog(decoder.outCount("decoder_count"));
        logger.outLogRVB(rate,BER);
    }
}