
#include "../lib/Common.h"
#include "../lib/Preseter.h"
#include "../lib/Decoder.h"
#include "../lib/Encoder.h"
#include "../lib/Channel.h"
#include "../lib/Analysor.h"
#include "../lib/Logger.h"
#include "../lib/Performance.h"

int main(void) {
    Params::set_e(0.5);
    Params::set_N(pow(2,10));
    Params::set_K(pow(2,8));
    Params::set_s(BEC);
//    int divNum = 7;
//    Params::set_K(Params::get_N() / (divNum * 2));

    vector<int> A(Params::get_K(), -1);
    vector<int> u_n(Params::get_N(), 0);
    vector<int> x_n(Params::get_N(), 0);
    vector<int> y_n(Params::get_N(), 0);
    vector<int> u_est(Params::get_N(), 0);

    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    Analysor analysor;

//    analysor.probErrBound();

    performance.startTimer();

//    preseter.defineFixedAndFree(A, BEC);

//    cout << "[[" << 63 << "]]" << Analysor::calcBhat(63,Params::get_N()) << endl;
//    vector<double> bha(Params::get_N());
//    for (int i = 0; i < Params::get_N(); i++) {
//        cout << i << endl;
//        bha[i] = Analysor::calcBhat(i,Params::get_N());
//    }
//    vector<double> z(Params::get_N());
//    Analysor::makeArrayBhat(z);
//    Common::pp(z);


    analysor.calcBlockErrorRate(TEST);
    performance.stopTimer();
    logger.outLog("=================================");
    logger.outLog(performance.outTime("処理時間"));
    logger.outLog("(N,K) = (" + to_string(Params::get_N()) + "," + to_string(Params::get_K()) + ")");

    performance.outHMS();

    return 0;
}




