#include "../lib/Common.h"
#include "../lib/Preseter.h"
#include "../lib/Decoder.h"
#include "../lib/Encoder.h"
#include "../lib/Channel.h"
#include "../lib/Analysor.h"
#include "../lib/Logger.h"
#include "../lib/Performance.h"

int main(void) {
//    Params params(1024, 5, 0.5);
//    Params(int _N, int _K, double _e){
//        N = _N;
//        K = _K;
//        e = _e;
//    }
    Params::set_N(1024);
    Params::set_K(5);
    Params::set_e(0.5);

    vector<int> u_Ac(Params::get_N(),0);
    vector<int> u_A(Params::get_N(),0);
    vector<int> A(Params::get_N(),0);
    vector<int> u_n(Params::get_N(),0);
    vector<int> x_n(Params::get_N(),0);
    vector<int> y_n(Params::get_N(),0);
    vector<int> u_est(Params::get_N(),0);
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    Analysor analysor;

//    analysor.calcBlockErrorRate(TEST, N, 5);
    Preseter::preset(RAND, u_n, u_Ac, u_A);
//
    performance.startTimer();
    Common::pp(u_n);
//    Common::pp(u_Ac);
//    Common::pp(u_A);
//    x_n = encoder.encode(N, u_n);
//    Common::pp(x_n);
//    y_n = Channel::channel_output(x_n);
//    Common::pp(y_n);
//    u_est = decoder.decode(y_n, u_n, u_Ac, u_A);
//    Common::pp(u_est);
//
//    cout << "error　probability:" << Analysor::errorRate(u_n, u_est) << endl;
//    cout << "rate:" << (double)K/N << endl;
//
//    performance.stopTimer();

//    logger.outLog("================================");
//    logger.outLogTime(performance.outTime("処理時間"));
//    logger.outLog("(N,K) = (" + to_string(N) + "," +to_string(K) + ")");
//    logger.outLog("error　probability:" + to_string(Analysor::errorRate(u_n, u_est)));
//    logger.outLog("rate:" + to_string((double)K/N));
//    logger.outLogCount(encoder.outCount("encoder_count"));
//    logger.outLogCount(decoder.outCount("decoder_count"));
//    logger.outLog("================================");

    return 0;
}