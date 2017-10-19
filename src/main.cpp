
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

    int divNum = 15;
    Params::set_K(Params::get_N() / (divNum * 2));
//    Params::set_K(16);

    performance.startTimer();

    analysor.calcBlockErrorRate(RUN);
    performance.stopTimer();
    logger.outLog("=================================");
    logger.outLog(performance.outTime("処理時間"));
    logger.outLog("(N,K) = (" + to_string(Params::get_N()) + "," + to_string(Params::get_K()) + ")");

    performance.outHMS();
    cout << Params::get_rvbDir() << endl;
}

void testNormal(){
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    logger.setRvbDir(Params::get_rvbDir());

    int c = 0;
    vector<int> A;
    vector<int> Ac;
    Preseter::preset_A_Ac(A,Ac);
    string efn = "/Users/ryotaro/labo/log/erlog " + to_string(Params::get_N());
    ofstream efile;
    efile.open(efn, ios::out);
    efile << "i" << "\t" << "err" << "\t" << "itr" << "\t" << "no_checked" << endl;
    Common::pp(A);
    Common::pp(Ac);
//    A = {3};
    int num_1 = 1;

//    while(c < Params::get_N()){
        vector<int> u_n(Params::get_N(), 0);
        vector<int> u_est(Params::get_N(), 0);
        vector<int> x_n(Params::get_N(), 0);
        vector<double> y_n(Params::get_N(), 0.0);
        Preseter::preset_u(RAND, u_n);

        cout << "u_n" << endl;
        Common::pp(u_n);
        x_n = encoder.encode(Params::get_N(), u_n);

        Common::pp(x_n);
        cout << "x_n" << endl;
//        y_n = Channel::channel_output(x_n);

//        for (int i = 0; i < num_1; i++) {
//            y_n[A[i]] = 1.0;
//        }
        y_n[c] = 1.0;
        Common::pp(y_n);
        cout << "y_n" << endl;

        int error_count=0;
        int sc_error = 0;
        int bp_error = 0;

        //SC
        Common::bar();
        u_est = decoder.decode(y_n, u_n, A);
        Analysor::errorCount(u_n, u_est, &error_count);
        sc_error = error_count;
        cout << "error " << sc_error << "/" << Params::get_N() <<  endl;

        error_count=0;
        //BP
        Common::bar();
        vector<int> param(2,0);
        u_est = decoder.calcBP(param, y_n, u_n, A);
        Analysor::errorCount(u_n, u_est, &error_count);
        bp_error = error_count;
        cout << "error " << bp_error << "/" << Params::get_N() <<  endl;

        efile << c << "\t" << error_count << "\t" << param[0] << "\t" << param[1] << endl;
        cout << efn << endl;

        Common::bar();
        cout << "count:::" << c << endl;
        Common::bar();
        c++;

//        if(bp_error > sc_error) break;
//    }
}


void testMiddleval(){
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    logger.setRvbDir(Params::get_rvbDir());
    vector<int> A;
    vector<int> Ac;
    Preseter::preset_A_Ac(A,Ac);
    Common::pp(A);
    Common::pp(Ac);
    int size = log2(Params::get_N())+1;

    vector<int> u_n(Params::get_N(), 0);
    vector<int> u_est(Params::get_N(), 0);
    vector<vector<int> > tmp_x_n(log2(Params::get_N()), vector<int>(Params::get_N(),0) );
    vector<vector<int> > x_n(size, vector<int>(Params::get_N(),0) );
    vector<vector<double> > y_n(size, vector<double>(Params::get_N(), 0.0));

    A = {31,30,29,27,23,15,28,26,25,22,21,14,13,24,19,20,11,12,7,3};

    int count = 1;
    Preseter::preset_u(ALL0, u_n);
//        u_n = {1,1,0,1,0,0,0,1,0,0,0,0,1,1,0,0,1,1,1,0,1,0,0,1,1,1,0,0,1,1,1,1};
    Common::pp(u_n);

    //最後の列だけ返り血で入れる。それ以外はencode_m中で入る。
    tmp_x_n[log2(Params::get_N()) - 1] = encoder.encode_m(Params::get_N(), 0, 0, u_n, tmp_x_n);

    //uとxをまとめて送る
    for (int i = 0; i < log2(Params::get_N()) + 1; i++) {
        for (int j = 0; j < Params::get_N(); j++) {
            if (i == 0) {
                x_n[i][j] = u_n[j];
            } else {
                x_n[i][j] = tmp_x_n[i - 1][j];
            }
        }
    }
//        Channel::channel_output_m(x_n, y_n);

    for (int i = 0; i < Params::get_N(); i++) {
        if( i == 31){
            y_n[size-1][i] = 1.0;
        } else {
            y_n[size-1][i] = 0.0;
        }
    }
//        y_n[size-1] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1};
//        y_n[size-1][Params::get_N()/2-1] = 0;
//        y_n[size-1][Params::get_N()-1] = 0;
    Common::pp(x_n[size - 1]);
    Common::pp(y_n[size - 1]);

    int error_count = 0;
    int bp_error = 0;

    error_count = 0;
    vector<int> param(2, 0);
    u_est = decoder.calcBP_m(param, y_n, u_n, A);
    Analysor::errorCount(u_n, u_est, &error_count);
    bp_error = error_count;
    cout << "error " << bp_error << "/" << Params::get_N() << endl;

    Common::bar();
    cout << count << endl;
    count++;
    Common::bar();
//    if(error_count) break;
}

int main(void) {
    Params::set_e(0.11);
    Params::set_N(pow(2,3));
    Params::set_K(5);
    Params::set_s(BSC);
    Params::set_decode_mode(SC);
    Params::set_monteNum(1);
    Params::set_rp(100);
    Params::set_upperBlockErrorNum(10000);
//    calcBER();
    testNormal();
//    testMiddleval();
    return 0;
}