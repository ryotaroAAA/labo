
#include "../lib/Common.h"
#include "../lib/Preseter.h"
#include "../lib/Decoder.h"
#include "../lib/Encoder.h"
#include "../lib/Channel.h"
#include "../lib/Analysor.h"
#include "../lib/Logger.h"
#include "../lib/Performance.h"
#include <sys/time.h>
//#include "../lib/ps.h"
//#include <bits/stdc++.h>

void calcBER(){
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    Analysor analysor;

    EXP_MODE em = Params::get_exp_mode();
    MID_MODE mm = Params::get_m_mode();
    string ename;
    switch (em) {
        case NORMAL: ename = ""; break;
        case PUNC: ename = "punc"; break;
        case MID:
            switch (mm) {
                case MID_BLUTE: {
                    int bl = 0;
                    if(Params::get_Bloop()){
                        bl = Params::get_Bloop();
                    }
                    ename = "mid_blute" + to_string(bl);
                    break;
                }
                case MID_ADOR: ename = "mid_ador"; break;
                case MID_AOR: ename = "mid_aor"; break;
                case MID_DOR: ename = "mid_dor"; break;
                case MID_AOB: ename = "mid_aob"; break;
                case MID_DOB: ename = "mid_dob"; break;
                case MID_AOV: ename = "mid_aov"; break;
                case MID_DOV: ename = "mid_dov"; break;
            }
            break;
        case QUP: ename = "qup"; break;
        case WANG: ename = "wang"; break;
        case M_QUP: ename = "m_qup"; break;
        case M_WANG: ename = "m_wang"; break;
        case VALERIO_P: ename = "valerio_p"; break;
        case VALERIO_S: ename = "valerio_s"; break;
        case M_VALERIO_P: ename = "m_valerio_p"; break;
        case M_VALERIO_S: ename = "m_valerio_s"; break;
    }
    cout << ename << endl;
//    cout << "aaaaaaaaaaa" << endl;

    time_t now = time(NULL);
    struct tm *pnow = localtime(&now);
    CHANNEL_TYPE channel_type = Params::get_s();

    stringstream fn;
    if (channel_type == BEC || channel_type == BSC) {
        fn << "log/"
             << to_string(pnow->tm_year + 1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday)
             << "/N" << to_string(Params::get_N())
             << "M" << to_string(Params::get_M())
             << "MN" << to_string(Params::get_MN())
             << "e" << Params::get_e()
             << "rp" + to_string(Params::get_rp())
             << "BNum" << to_string(Params::get_blockNum())
             << "UPBn" << to_string(Params::get_upperBlockErrorNum())
             << "_" << (Params::get_s() ? "BSC" : "BEC") << "_" << ename;
        Params::set_rvbDir(fn.str());
    } else if (channel_type == AWGN) {
        if(Params::get_is_exp_awgn()) {
            double p[3];
            Params::get_awgn_p(p);
            double rate = p[0];
            fn << "log/"
               << to_string(pnow->tm_year + 1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday)
               << "/N" + to_string(Params::get_N())
               << "M" << to_string(Params::get_M())
               << "MN" << to_string(Params::get_MN())
               << "R" << (rate)
               << "K" << (Analysor::get_eachK(rate))
               << "AWGN"
               << "rp" + to_string(Params::get_rp())
               << "MNum" + to_string(Params::get_monteNum())
               << "BNum" + to_string(Params::get_blockNum())
               << "UPBn" + to_string(Params::get_upperBlockErrorNum()) + "_" + ename;
        } else {
            fn << "log/"
               << to_string(pnow->tm_year + 1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday)
               << "/N" + to_string(Params::get_N())
               << "M" << to_string(Params::get_M())
               << "MN" << to_string(Params::get_MN())
               << "sdiv" << Params::get_e()
               << "AWGN"
               << "rp" + to_string(Params::get_rp())
               << "MNum" + to_string(Params::get_monteNum())
               << "BNum" + to_string(Params::get_blockNum())
               << "UPBn" + to_string(Params::get_upperBlockErrorNum()) + "_" + ename;
        }
        Params::set_rvbDir(fn.str());
    }
    cout << fn.str() << endl;

    performance.startTimer();
    analysor.calcBlockErrorRate_BP();
//    analysor.calcBlockErrorRate();
    performance.stopTimer();

    performance.outHMS();
    cout << Params::get_rvbDir() << endl;
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
//awgn なら e=0.8
//enum EXP_MODE{NORMAL, PUNC, QUP, WANG, MID, M_WANG, M_QUP, VALERIO_P, VALERIO_S, M_VALERIO_P, M_VALERIO_S};
int main(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    srand(tv.tv_sec + tv.tv_usec);
    init_genrand(tv.tv_sec + tv.tv_usec);

    Params::set_e(0.6);
    Params::set_N(256);
    Params::set_K(1);
    Params::set_M(0);
//    Params::set_M(96);
    Params::set_MN(40);
    Params::set_s(AWGN);
//    Params::set_is_outlog(true);
    Params::set_decode_mode(BP);
    Params::set_is_exp_awgn(true);
    Params::set_monteNum(1);
    Params::set_rp(100);
    Params::set_Bloop(500);
    Params::set_blockNum(2000);
    Params::set_upperBlockErrorNum(100);

    //point設定
//    double p[3] = {5,0.15,0.45};
//    double p[3] = {1,0.27,0.27};
    double p[3] = {1,1,1};
    Params::set_point(p);
    //awgn_p設定
    //{0.25,1.5,3.0},{0.5,1.5,3.0},{0.75,2.5,4.5}
//    double awgn_p[3] = {0.25,1.0,2.0};
//    double awgn_p[3] = {0.5,1.0,2.5};
    double awgn_p[3] = {0.75,2.5,4.5};
    Params::set_awgn_p(awgn_p);
    //calc_bloopを計算したいときにtrue, それ以外は必ずコメントアウト
//    Params::set_is_calc_bloop(true);

//    Params::set_m_mode(MID_BLUTE);
//    Params::set_m_mode(MID_DOR);
//    Params::set_m_mode(MID_DOB);
    Params::set_m_mode(MID_DOV);

//    Params::set_m_mode(MID_ADOR);
//    Params::set_m_mode(MID_AOR);
//    Params::set_m_mode(MID_AOB);
//    Params::set_m_mode(MID_AOV);

//    Params::set_exp_mode(NORMAL);
//    Params::set_exp_mode(QUP);
    Params::set_exp_mode(MID);
//    Params::set_exp_mode(WANG);
//    Params::set_exp_mode(VALERIO_P);
//    Params::set_exp_mode(VALERIO_S);
//    Params::set_exp_mode(M_WANG);
//    Params::set_exp_mode(M_QUP);
//    Params::set_exp_mode(M_VALERIO_P);
//    Params::set_exp_mode(M_VALERIO_S);

    calcBER();
//    Common::gauss_dist(, -1.0, 0.6);

    return 0;
}