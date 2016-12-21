#include "../lib/Common.h"
#include "../lib/Preseter.h"
#include "../lib/Decoder.h"
#include "../lib/Encoder.h"
#include "../lib/Channel.h"
#include "../lib/Analysor.h"
#include "../lib/Logger.h"
#include "../lib/Performance.h"

int main(void) {
    vector<int> u_Ac(Params::N,0);
    vector<int> u_A(Params::N,0);
    vector<int> A(Params::N,0);
    vector<int> u_n(Params::N,0);
    vector<int> x_n(Params::N,0);
    vector<int> y_n(Params::N,0);
    vector<int> u_est(Params::N,0);
    Performance performance;
    Decoder decoder;
    Encoder encoder;

    Preseter::preset(u_n, u_Ac, u_A);

    performance.startTimer();
    Common::pp(u_n);
    x_n = encoder.encode(Params::N, u_n);
    Common::pp(x_n);
    y_n = Channel::channel_output(x_n);
    Common::pp(y_n);
    u_est = decoder.decode(y_n, u_n, u_Ac, u_A);

    cout << "error　probability:" << Analysor::errorRate(u_n, u_est) << endl;
    cout << "rate:" << (double)Params::K/Params::N << endl;

    performance.stopTimer();
    performance.outTime("処理時間");
    encoder.outCount("encoder_count");
    decoder.outCount("decoder_count");
    return 0;
}

