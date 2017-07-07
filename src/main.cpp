
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
    Params::set_N(pow(2,6));
    Params::set_K(pow(2,5));
    Params::set_s(BEC);
//    int divNum = 20;
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

    cout << "[[" << 63 << "]]" << Analysor::calcBhat(63,Params::get_N()) << endl;
//    vector<double> bha(Params::get_N());
//    for (int i = 0; i < Params::get_N(); i++) {
//         bha[i] = Analysor::calcBhat(i,Params::get_N(),BEC);
//    }
//    string filename = "/Users/ryotaro/labo/log/test_bha_monte_carlo_1000";
//    string filename = "/Users/ryotaro/labo/log/test_bha_true";
//    ofstream w_file;
//    w_file.open(filename, ios::out);
//    for (int i = 0; i < Params::get_N(); i++)
//    {
//        w_file << i << " " << bha[i] << endl;
//        cout << i << " " << bha[i] << endl;
//    }
//    Common::pp(bha);

//    analysor.calcBlockErrorRate(RUN, BEC);
    performance.stopTimer();
    logger.outLog("=================================");
    logger.outLog(performance.outTime("処理時間"));
    logger.outLog("(N,K) = (" + to_string(Params::get_N()) + "," + to_string(Params::get_K()) + ")");

    performance.outHMS();

    return 0;
}




