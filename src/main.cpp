
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
                case MID_BLUTE: ename = "mid_blute"; break;
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
    cout << "[[" << ename << "]]" << endl;

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
             << "BNum" << to_string(Params::get_blockNum())
             << "UPBn" << to_string(Params::get_upperBlockErrorNum())
             << "c" << (Params::get_s() ? "BSC" : "BEC") << ":" << ename;
        Params::set_rvbDir(fn.str());
    } else if (channel_type == AWGN) {
        fn << "log/"
           << to_string(pnow->tm_year+1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday)
           << "/N" + to_string(Params::get_N())
           << "M" << to_string(Params::get_M())
           << "MN" << to_string(Params::get_MN())
           << "sdiv" << Params::get_e()
           << "_AWGN_"
           << "MNum" + to_string(Params::get_monteNum())
           << "BNum" + to_string(Params::get_blockNum())
           << "UPBn" + to_string(Params::get_upperBlockErrorNum()) +  ":" + ename;
        Params::set_rvbDir(fn.str());
    }

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
    Params::set_MN(64);
    Params::set_s(AWGN);
    Params::set_is_outlog(false);
    Params::set_decode_mode(BP);
    Params::set_monteNum(1);
    Params::set_rp(50);
    Params::set_Bloop(200);

    Params::set_blockNum(3000);
    Params::set_upperBlockErrorNum(100);

//    double point[3] = {1,1,1};
    double point[3] = {7,0.2,0.5};
    Params::set_is_calc_bloop(false);
    Params::set_point(point);

    //IDかBD
    Params::set_m_mode(MID_BLUTE);
//    Params::set_m_mode(MID_DOR);
//    Params::set_m_mode(MID_DOB);
//    Params::set_m_mode(MID_DOV);

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
//    Params::set_exp_mode(M_WANG);
//    Params::set_exp_mode(M_VALERIO_P);
//    Params::set_exp_mode(M_VALERIO_S);

//    calcBER();
//    vector<int> t1;
//    vector<int> t2;
//    Common::get_rate(t1,t2);
//    Analysor::calcBlockErrorRate();
    return 0;
}