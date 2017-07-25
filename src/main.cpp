
#include "../lib/Common.h"
#include "../lib/Preseter.h"
#include "../lib/Decoder.h"
#include "../lib/Encoder.h"
#include "../lib/Channel.h"
#include "../lib/Analysor.h"
#include "../lib/Logger.h"
#include "../lib/Performance.h"


int main(void) {
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    Analysor analysor;

    Params::set_e(0.11);
    Params::set_N(pow(2,12));
    Params::set_s(BSC);
    Params::set_monteNum(10000);
    Params::set_blockNum(100000);
    Params::set_upperBlockErrorNum(1000000);

    time_t now = time(NULL);
    struct tm *pnow = localtime(&now);
    if (Params::get_s() == BEC) {
        Params::set_rvbDir("/Users/ryotaro/labo/log/"
                           + to_string(pnow->tm_year+1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday)
                           + "/N=" + to_string(Params::get_N())
                           + " e=" + to_string(Params::get_e())
                           + " channel=" + (Params::get_s()?"BSC":"BEC")
        );
    } else {
        Params::set_rvbDir("/Users/ryotaro/labo/log/"
                           + to_string(pnow->tm_year+1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday)
                           + "/N=" + to_string(Params::get_N())
                           + " e=" + to_string(Params::get_e())
                           + " channel=" + (Params::get_s()?"BSC":"BEC")
                           + " monteNum=" + to_string(Params::get_monteNum())
                           + " blockNum=" + to_string(Params::get_blockNum())
        );
    }

    int divNum = 20;
    Params::set_K(Params::get_N() / (divNum * 2));
//    Params::set_K(22);

    vector<int> A(Params::get_K(), -1);
    vector<int> u_n(Params::get_N(), 0);
    vector<int> x_n(Params::get_N(), 0);
    vector<int> y_n(Params::get_N(), 0);
    vector<int> u_est(Params::get_N(), 0);

    performance.startTimer();

    analysor.calcBlockErrorRate(RUN);
    performance.stopTimer();
    logger.outLog("=================================");
    logger.outLog(performance.outTime("処理時間"));
    logger.outLog("(N,K) = (" + to_string(Params::get_N()) + "," + to_string(Params::get_K()) + ")");

    performance.outHMS();
    cout << Params::get_rvbDir() << endl;

    return 0;
}




