#include "../lib/Analysor.h"

Analysor::Analysor(){

}

Analysor::~Analysor(){

}

double Analysor::bitErrorRate(vector<int> &u, vector<int> &u_est){
    int error_count = 0;
    for(int i=0; i<Params::get_N(); i++){
        if(u[i] != u_est[i]){
            error_count++;
        }
    }
    cout << "[[error count::" << error_count << "]]" << endl;
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
    vector<double> tempArr;
    vector<double> sumArr;
    double sum = 0.0;

    for(int i = 0; i < Params::get_N(); i++) {
        tempArr.push_back(Analysor::calcBhatForBec(i, Params::get_N()));
    }

    sort(tempArr.begin(), tempArr.end(), greater<int>());

    for(int i=0; i < Params::get_N(); i++){
        sum += tempArr[i];
        sumArr.push_back(sum);
    }

    string filename = "/Users/ryotaro/labo/log/rve";
    ofstream w_file;
    w_file.open(filename, ios::out);
    for (int i = 0; i < Params::get_N(); i++)
    {
        w_file << (double)i/Params::get_N() << " " << sumArr[i] << " " << tempArr[i] << endl;
    }
}

//rate vs BERのグラフ作成用
void Analysor::calcBlockErrorRate(MODE mode, CHANNEL_TYPE channel) {
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    logger.setDir("/Users/ryotaro/labo/log");
    vector<int> A(Params::get_K(), 0);
    vector<int> u_n(Params::get_N(), 0);
    vector<int> x_n(Params::get_N(), 0);
    vector<int> y_n(Params::get_N(), 0);
    vector<int> u_est(Params::get_N(), 0);
    double BER = 0.0;
    double rate = 0.0;
    int n = Params::get_N();

    int tmp = Params::get_N()/Params::get_K();
    int tmpK = Params::get_K();
    for (int i = 1; i <= tmp; i++) {
        cout << i << endl;
        Params::set_K(i*tmpK);

        A.assign(Params::get_K(), 0);
        u_n.assign(Params::get_N(), 0);
        x_n.assign(Params::get_N(), 0);
        y_n.assign(Params::get_N(), 0);
        u_est.assign(Params::get_N(), 0);
        Preseter::preset(RAND, u_n, A);
        performance.startTimer();

        x_n = encoder.encode(Params::get_N(), u_n);
        y_n = Channel::channel_output(x_n, channel);
        u_est = decoder.decode(y_n, u_n, x_n, A);

        BER = Analysor::bitErrorRate(u_n, u_est);
        rate = (double)Params::get_K() / Params::get_N();
        performance.stopTimer();

        logger.outLog("=================================");
        logger.outLog(performance.outTime("処理時間"));
        logger.outLog("(N,K) = (" + to_string(Params::get_N()) + "," + to_string(Params::get_K()) + ")");
        logger.outLog("BitER:" + to_string(BER));
        logger.outLog("Rate:" + to_string(rate));
        logger.outLog(encoder.outCount("encoder_count"));
        logger.outLog(decoder.outCount("decoder_count"));

        logger.outLogRVB(rate, BER);
        if(mode != TEST) break;
    }
}