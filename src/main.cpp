
#include "../lib/Common.h"
#include "../lib/Preseter.h"
#include "../lib/Decoder.h"
#include "../lib/Encoder.h"
#include "../lib/Channel.h"
#include "../lib/Analysor.h"
#include "../lib/Logger.h"
#include "../lib/Performance.h"
//#include "../lib/ps.h"


void calcBER(){
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    Analysor analysor;

    Params::set_e(0.5);
    Params::set_N(pow(2,2));
    Params::set_s(BEC);
    Params::set_monteNum(1);
    Params::set_blockNum(1);
    Params::set_upperBlockErrorNum(10000);

    time_t now = time(NULL);
    struct tm *pnow = localtime(&now);
    if (Params::get_s() == BEC) {
        Params::set_rvbDir("/Users/ryotaro/labo/log/"
                           + to_string(pnow->tm_year+1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday)
                           + "/N=" + to_string(Params::get_N())
                           + " e=" + to_string(Params::get_e())
                           + " channel=" + (Params::get_s()?"BSC":"BEC")
        );
    } else if (Params::get_s() == BSC) {
        Params::set_rvbDir("/Users/ryotaro/labo/log/"
                           + to_string(pnow->tm_year+1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday)
                           + "/N=" + to_string(Params::get_N())
                           + " e=" + to_string(Params::get_e())
                           + " channel=" + (Params::get_s()?"BSC":"BEC")
                           + " monteNum=" + to_string(Params::get_monteNum())
                           + " blockNum=" + to_string(Params::get_blockNum())
        );
    } else if (Params::get_s() == AWGN) {
        Params::set_rvbDir("/Users/ryotaro/labo/log/"
                           + to_string(pnow->tm_year+1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday)
                           + "/N=" + to_string(Params::get_N())
                           + " sdiv=" + to_string(Params::get_e())
                           + " channel=AWGN"
                           + " monteNum=" + to_string(Params::get_monteNum())
                           + " blockNum=" + to_string(Params::get_blockNum())
        );
    }

//    int divNum = 10;
//    Params::set_K(Params::get_N() / (divNum * 2));
    Params::set_K(2);

    performance.startTimer();

    analysor.calcBlockErrorRate(RUN);
    performance.stopTimer();
    logger.outLog("=================================");
    logger.outLog(performance.outTime("処理時間"));
    logger.outLog("(N,K) = (" + to_string(Params::get_N()) + "," + to_string(Params::get_K()) + ")");

    performance.outHMS();
    cout << Params::get_rvbDir() << endl;
}


int main(void) {

//    calcBER();
    Params::set_e(0.5);
    Params::set_N(pow(2,6));
    Params::set_K(15);
    Params::set_s(BEC);
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    logger.setRvbDir(Params::get_rvbDir());
    vector<int> A(Params::get_K(), -1);
    vector<int> u_n(Params::get_N(), 0);
    vector<int> u_est(Params::get_N(), 0);
    vector<int> x_n(Params::get_N(), 0);
    vector<double> y_n(Params::get_N(), 0);

    vector<pair<int, double> > cap_map;
    Preseter::makeMutualInfoArray(cap_map);
    Preseter::represet_A(A, cap_map);

    Preseter::preset_u(ALL0, u_n);

    Common::pp(u_n);
    x_n = encoder.encode(Params::get_N(), u_n);
    Common::pp(x_n);
    Common::bar();
    y_n = Channel::channel_output(x_n);
    Common::pp(y_n);
    Common::bar();
    u_est = decoder.BP(100, y_n, u_n, A);
    Common::pp(u_est);

    int error_count=0;
    Analysor::errorCount(u_n, u_est, &error_count);
    cout << "error " << error_count << "/" << Params::get_N() <<  endl;
    return 0;
}