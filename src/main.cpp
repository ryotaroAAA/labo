#include "../lib/Common.h"
#include "../lib/Preseter.h"
#include "../lib/Decoder.h"
#include "../lib/Encoder.h"
#include "../lib/Channel.h"
#include "../lib/Analysor.h"
#include "../lib/Logger.h"
#include "../lib/Performance.h"

int main(void) {
    Params::set_N(1024);
    Params::set_K(250);
    Params::set_e(0.5);

    vector<int> u_Ac;
    vector<int> u_A;
    vector<int> A;
    vector<int> u_n;
    vector<int> x_n;
    vector<int> y_n;
    vector<int> u_est;
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    Analysor analysor;

//    analysor.calcBlockErrorRate(TEST, Params::get_N());
    Preseter::preset(RAND, u_n, u_Ac, u_A);

    performance.startTimer();
    x_n = encoder.encode(Params::get_N(), u_n);
    Common::pp(x_n);
    y_n = Channel::channel_output(x_n);
    Common::pp(y_n);
    u_est = decoder.decode(y_n, u_n, u_Ac, u_A);
    Common::pp(u_est);

    cout << "error　probability:" << Analysor::errorRate(u_n, u_est) << endl;
    cout << "rate:" << (double)Params::get_K()/Params::get_N() << endl;

    performance.stopTimer();

    logger.outLog("================================");
    logger.outLog(performance.outTime("処理時間"));
    logger.outLog("(N,K) = (" + to_string(Params::get_N()) + "," +to_string(Params::get_K()) + ")");
    logger.outLog("error　probability:" + to_string(Analysor::errorRate(u_n, u_est)));
    logger.outLog("rate:" + to_string((double)Params::get_K()/Params::get_N()));
    logger.outLog(encoder.outCount("encoder_count"));
    logger.outLog(decoder.outCount("decoder_count"));
    logger.outLog("================================");

    return 0;
}