
#include "../lib/Common.h"
#include "../lib/Preseter.h"
#include "../lib/Decoder.h"
#include "../lib/Encoder.h"
#include "../lib/Channel.h"
#include "../lib/Analysor.h"
#include "../lib/Logger.h"
#include "../lib/Performance.h"

void calcBER(){
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    Analysor analysor;

    Params::set_e(0.5);
    Params::set_N(pow(2,10));
    Params::set_s(BEC);
    Params::set_monteNum(1);
    Params::set_blockNum(100);
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
    Params::set_N(pow(2,5));
    Params::set_s(BEC);
    Encoder encoder;
    vector<int> u_n(Params::get_N());
    vector<vector<int> > x_n(log2(Params::get_N()), vector<int>(Params::get_N()) );
    Preseter::preset_u(ALL1, u_n);
    x_n[log2(Params::get_N())- 1] = encoder.encode_m(Params::get_N(), 0, 0, u_n, x_n);
    Common::pp(x_n[log2(Params::get_N())- 1]);

    
    return 0;
}