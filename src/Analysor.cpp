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
    //cout << "i:" << i << ", n:" << n <<endl;
//    PRINT(cap);
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
//    cout << "i:" << i << ", n:" << n <<endl;
//    PRINT(bha);
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

    string filename = "/Users/ryotaro/labo/e18";
    ofstream w_file;
    w_file.open(filename, ios::out);
    for (int i = 0; i < Params::N; i++)
    {
        w_file << (double)i/Params::N << " " << sumArr[i] << " " << tempArr[i] << endl;
    }
}

void Analysor::calcBlockErrorRate(MODE mode, int n, double dist) {
    vector<int> u_Ac(Params::N, 0);
    vector<int> u_A(Params::N, 0);
    vector<int> A(Params::N, 0);
    vector<int> u_n(Params::N, 0);
    vector<int> x_n(Params::N, 0);
    vector<int> y_n(Params::N, 0);
    vector<int> u_est(Params::N, 0);

    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;

    for (int k=0; k<Params::N ; k + Param::N*dist) {
        Preseter::preset(u_n, u_Ac, u_A);
        performance.startTimer();
        x_n = encoder.encode(Params::N, u_n);
        y_n = Channel::channel_output(x_n);
        u_est = decoder.decode(y_n, u_n, u_Ac, u_A);

        logger.outLog("================================");
        logger.outLog("error　probability:" + to_string(Analysor::errorRate(u_n, u_est)));
        logger.outLog("rate:" + to_string((double) Params::K / Params::N));

        performance.stopTimer();

        logger.outLog(performance.outTime("処理時間"));
        logger.outLog(encoder.outCount("encoder_count"));
        logger.outLog(decoder.outCount("decoder_count"));
    }
}