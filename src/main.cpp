
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

    EXP_MODE em = Params::get_exp_mode();
    string ename;
    switch (em) {
        case NORMAL: ename = ""; break;
        case PUNC: ename = " punc"; break;
        case MID: ename = " mid"; break;
        case QUP: ename = " qup"; break;
        case WANG: ename = " wang"; break;
        case M_QUP: ename = " m_qup"; break;
        case M_WANG: ename = " m_wang"; break;
        case VALERIO_P: ename = " valerio_p"; break;
        case VALERIO_S: ename = " valerio_s"; break;
        case M_VALERIO_P: ename = " m_valerio_p"; break;
        case M_VALERIO_S: ename = " m_valerio_s"; break;
    }

    time_t now = time(NULL);
    struct tm *pnow = localtime(&now);
    CHANNEL_TYPE channel_type = Params::get_s();
    if (channel_type == BEC || channel_type == BSC) {
        Params::set_rvbDir("log/"
                           + to_string(pnow->tm_year + 1900) + to_string(pnow->tm_mon + 1) +
                           to_string(pnow->tm_mday)
                           + "/N=" + to_string(Params::get_N())
                           + " e=" + to_string(Params::get_e())
                           + " blockNum=" + to_string(Params::get_blockNum())
                           + " upperblockNum=" + to_string(Params::get_upperBlockErrorNum())
                           + " channel=" + (Params::get_s() ? "BSC" : "BEC") + ename
        );
    } else if (channel_type == AWGN) {
        Params::set_rvbDir("log/"
                           + to_string(pnow->tm_year+1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday)
                           + "/N=" + to_string(Params::get_N())
                           + " sdiv=" + to_string(Params::get_e())
                           + " channel=AWGN"
                           + " monteNum=" + to_string(Params::get_monteNum())
                           + " blockNum=" + to_string(Params::get_blockNum())
                           + " upperblockNum=" + to_string(Params::get_upperBlockErrorNum()) + ename
        );
    }


    performance.startTimer();

//    if( mid ){
//        analysor.calcBlockErrorRate_mid_BP(RUN);
//    } else if(wang) {
//        analysor.calcBlockErrorRate_BP_wang(RUN);
//    } else if(mwang) {
//        analysor.calcBlockErrorRate_mid_wang_BP(RUN);
//    } else {
        analysor.calcBlockErrorRate_BP(TEST);
//    }


    performance.stopTimer();

    performance.outHMS();
    cout << Params::get_rvbDir() << endl;
}

void testNormal() {
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    logger.setRvbDir(Params::get_rvbDir());

    int c = 0;
    vector<int> A;
    vector<int> Ac;
    Preseter::preset_A_Ac(A, Ac);
    string efn = "/Users/ryotaro/labo/log/erlog2 " + to_string(Params::get_N());
    ofstream efile;
    efile.open(efn, ios::out);

    string itrfn = "/Users/ryotaro/labo/log/itrlog2 " + to_string(Params::get_N());
    ofstream itrfile;
    itrfile.open(itrfn, ios::out);
//    efile << "i" << "\t" << "err" << "\t" << "itr" << "\t" << "no_checked" << endl;
    vector<int> u_n(Params::get_N(), 0);
    vector<int> u_est(Params::get_N(), 0);
    vector<int> x_n(Params::get_N(), 0);
    Preseter::preset_u(ALL0, u_n);
    x_n = encoder.encode(Params::get_N(), u_n);

    for (int i = 0; i < Params::get_N(); i++) {
        for (int j = 0; j < Params::get_N(); j++) {

            vector<double> y_n(Params::get_N(), 0.0);
            y_n[i] = 1.0;
            y_n[j] = 1.0;
            int error_count = 0;

            error_count = 0;
            //BP
            Common::bar();
            vector<int> param(2, 0);//param[0]:itr, param[1]:no_checked
//            u_est = decoder.calcBP(param, y_n, u_n, A, Ac);
            Analysor::errorCount(u_n, u_est, &error_count);

            cout << "error " << error_count << "/" << Params::get_N() << endl;
//            cout << "rate :" << 1.0 * Params::get_K() / (Params::get_N()) << endl;

            efile << i+1 << " " << j+1 << " " << error_count << endl;
            itrfile << i+1 << " " << j+1 << " " << param[0] << endl;
//            efile << c << "\t" << error_count << "\t" << param[0] << "\t" << param[1] << endl;
            Common::bar();
            cout << "count:::" << c << endl;
            Common::bar();
            c++;
        }
        efile << endl;
        itrfile << endl;
    }


//    cout << efn << endl;
}

void testPunc() {
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    logger.setRvbDir(Params::get_rvbDir());

    int c = 0;
    vector<int> A;
    vector<int> Ac;
//    Preseter::preset_A_Ac(A, Ac);
    vector<pair<int, double> > cap_map;
    Preseter::makeMutualInfoArray(cap_map);
    vector<int> p_0(Params::get_M());
    vector<int> p(Params::get_M());
    Analysor::set_params(cap_map,A,Ac,p_0,p);
    cout << "A";
    Common::pp(A);
    cout << "Ac";
    Common::pp(Ac);

    vector<int> u_n(Params::get_N(), 0);
    vector<int> u_est(Params::get_N(), 0);
    vector<int> x_n(Params::get_N(), 0);
    vector<double> y_n(Params::get_N(), 0.0);
    Preseter::preset_u(RAND, u_n);
    x_n = encoder.encode(Params::get_N(), u_n);
    y_n = Channel::channel_output(x_n);

//    y_n[15] = 1.0;
    int error_count = 0;

    error_count = 0;
    //BP
    Common::bar();
    vector<int> param(2, 0);//param[0]:itr, param[1]:no_checked

//    u_est = decoder.calcBP_wang(p, param, x_n, y_n, u_n, A, Ac);
    Analysor::errorCount(u_n, u_est, &error_count);

    cout << "error " << error_count << "/" << Params::get_N() << endl;
    cout << "rate :" << 1.0 * Params::get_K() / (Params::get_N() - Params::get_M()) << endl;

    Common::bar();
    cout << "count:::" << c << endl;
    Common::bar();
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
//    Common::pp(A);
//    Common::pp(Ac);
    int size = log2(Params::get_N())+1;

    vector<int> u_n(Params::get_N(), 0);
    vector<int> u_est(Params::get_N(), 0);
    vector<vector<int> > tmp_x_n(log2(Params::get_N()), vector<int>(Params::get_N(),0) );
    vector<vector<int> > x_n(size, vector<int>(Params::get_N(),0) );
    vector<vector<double> > y_n(size, vector<double>(Params::get_N(), 0.0));

    int count = 1;
    Preseter::preset_u(RAND, u_n);
//    Common::pp(u_n);

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
    Channel::channel_output_m(x_n, y_n);


    for (int i = 0; i < Params::get_N(); i++) {
        if (i == 31) {
            y_n[size - 1][i] = 1.0;
        } else {
            y_n[size - 1][i] = 0.0;
        }
    }

    Common::pp(x_n[size - 1]);
    Common::pp(y_n[size - 1]);

    int error_count = 0;
    int bp_error = 0;

    error_count = 0;
    vector<int> param(2, 0);
//    u_est = decoder.calcBP_m(param, y_n, u_n, A);
    Analysor::errorCount(u_n, u_est, &error_count);
    bp_error = error_count;
    cout << "error :" << bp_error << "/" << Params::get_N() << endl;
    cout << "rate :" << 1.0*Params::get_K()/(Params::get_N()+Params::get_M()) << endl;

    Common::bar();
    cout << count << endl;
    count++;
    Common::bar();
//    if(error_count) break;
}

void testMidWang(){
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    logger.setRvbDir(Params::get_rvbDir());

    vector<int> A;
    vector<int> Ac;
//    Preseter::preset_A_Ac(A, Ac);
    vector<pair<int, double> > cap_map;
    Preseter::makeMutualInfoArray(cap_map);
    vector<int> p_0(Params::get_M());
    vector<int> p(Params::get_M());
    Analysor::set_params(cap_map,A,Ac,p_0,p);
    cout << "A";
    Common::pp(A);
    cout << "Ac";
    Common::pp(Ac);

    int size = log2(Params::get_N())+1;

    vector<int> u_n(Params::get_N(), 0);
    vector<int> u_est(Params::get_N(), 0);
    vector<vector<int> > tmp_x_n(log2(Params::get_N()), vector<int>(Params::get_N(),0) );
    vector<vector<int> > x_n(size, vector<int>(Params::get_N(),0) );
    vector<vector<double> > y_n(size, vector<double>(Params::get_N(), 0.0));

    int count = 1;
    Preseter::preset_u(RAND, u_n);
//    Common::pp(u_n);

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
    Channel::channel_output_m(x_n, y_n);


//    for (int i = 0; i < Params::get_N(); i++) {
//        if (i == 31) {
//            y_n[size - 1][i] = 1.0;
//        } else {
//            y_n[size - 1][i] = 0.0;
//        }
//    }

    Common::pp(x_n[size - 1]);
    Common::pp(y_n[size - 1]);

    int error_count = 0;
    int bp_error = 0;

    error_count = 0;
    vector<int> param(2, 0);
//    u_est = decoder.calcBP_m_wang(p, param, x_n, y_n, u_n, A);
    Analysor::errorCount(u_n, u_est, &error_count);
    bp_error = error_count;
    cout << "error :" << bp_error << "/" << Params::get_N() << endl;
    cout << "rate :" << 1.0 * Params::get_K()/(Params::get_N()) << endl;

    Common::bar();
    cout << count << endl;
    count++;
    Common::bar();

//    if(error_count) break;
}


inline string get_current_directory()
{
    char cwd[1024];
    getcwd(cwd, sizeof(cwd));
    return cwd;
}

bool checkFileExistence(const std::string& str)
{
    std::ifstream ifs(str);
    return ifs.is_open();
}

//mid>normal
int main(void) {
//    Params::set_e(0.5);
//    Params::set_N(512);
//    Params::set_M(32);
//    Params::set_s(BEC);
//    Params::set_is_outlog(false);
//    Params::set_decode_mode(BP);
//    Params::set_monteNum(20);
//    Params::set_rp(100);
//    Params::set_blockNum(1000);
//    Params::set_upperBlockErrorNum(1000);

    Params::set_e(0.5);
    Params::set_N(32);
    Params::set_K(1);
    Params::set_M(3);
    Params::set_s(BEC);
    Params::set_is_outlog(false);
    Params::set_decode_mode(BP);
    Params::set_monteNum(1);
    Params::set_rp(70);
    Params::set_blockNum(10000);
    Params::set_upperBlockErrorNum(100);
    Params::set_exp_mode(NORMAL);

//    vector<int> A;
//    vector<int> Ac;
//    vector<pair<int, double> > cap_map;
//    Preseter::makeMutualInfoArray(cap_map);
//    vector<int> p_0(Params::get_M());
//    vector<int> p(Params::get_M());
//    Analysor::set_params(cap_map,A,Ac,p_0,p);
//
//    Common::pp(A);
//    Common::pp(Ac);
//    Common::pp(p);
    EXP_MODE em = Params::get_exp_mode();
    if(em == MID){
//        Params::set_K(17*2); //中間ノード
        Params::set_K(40); //中間ノード
    } else if(em == QUP || em == WANG || em == PUNC) {
//        Params::set_K(15*2);
        Params::set_K(24);
    } else {
//        Params::set_K(16*2);
    }
    calcBER();

//    testNormal();
//    testPunc();
//    testMiddleval();
//    testMidWang();

    return 0;
}