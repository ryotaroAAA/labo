
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
    MID_MODE mm = Params::get_m_mode();
    string ename;
    switch (em) {
        case NORMAL: ename = ""; break;
        case PUNC: ename = "punc"; break;
        case MID:
            switch (mm) {
                case MID_IU: ename = "mid_iu"; break;
                case MID_ID: ename = "mid_id"; break;
                case MID_BU: ename = "mid_bu"; break;
                case MID_BD: ename = "mid_bd"; break;
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

    time_t now = time(NULL);
    struct tm *pnow = localtime(&now);
    CHANNEL_TYPE channel_type = Params::get_s();

    stringstream fn;
    if (channel_type == BEC || channel_type == BSC) {
        fn << "log/"
             << to_string(pnow->tm_year + 1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday)
             << "/N" << to_string(Params::get_N())
             << "e" << Params::get_e()
             << "BNum" << to_string(Params::get_blockNum())
             << "UPBn" << to_string(Params::get_upperBlockErrorNum())
             << "c" << (Params::get_s() ? "BSC" : "BEC") << ":" << ename;
        Params::set_rvbDir(fn.str());
    } else if (channel_type == AWGN) {
        fn << "log/"
           << to_string(pnow->tm_year+1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday)
           << "/N=" + to_string(Params::get_N())
           << "sdiv" << Params::get_e()
           << "c=AWGN"
           << "MNum=" + to_string(Params::get_monteNum())
           << "BNum=" + to_string(Params::get_blockNum())
           << "upperBn=" + to_string(Params::get_upperBlockErrorNum()) +  ":" + ename;
        Params::set_rvbDir(fn.str());
    }

    performance.startTimer();
    analysor.calcBlockErrorRate_BP();
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
//enum EXP_MODE{NORMAL, PUNC, QUP, WANG, MID, M_WANG, M_QUP, VALERIO_P, VALERIO_S, M_VALERIO_P, M_VALERIO_S};
int main(void) {
    Params::set_e(0.5);
    Params::set_N(256);
    Params::set_K(35);
    Params::set_M(32);
    Params::set_MN(32);
    Params::set_s(BEC);
    Params::set_is_outlog(false);
    Params::set_decode_mode(BP);
    Params::set_monteNum(1);
    Params::set_rp(50);
    Params::set_blockNum(1000);
    Params::set_upperBlockErrorNum(1000);
    //IDかBD
    Params::set_m_mode(MID_ID);
//    Params::set_m_mode(MID_IU);
//    Params::set_m_mode(MID_BU);
//    Params::set_m_mode(MID_BD);

//    Params::set_exp_mode(NORMAL);
//    Params::set_exp_mode(QUP);
//    Params::set_exp_mode(MID);
//    Params::set_exp_mode(WANG);
//    Params::set_exp_mode(VALERIO_P);
//    Params::set_exp_mode(VALERIO_S);
//    Params::set_exp_mode(M_WANG);
//    Params::set_exp_mode(M_VALERIO_P);
    Params::set_exp_mode(M_VALERIO_S);

//    vector<int> A;
//    vector<int> Ac;
//    vector<pair<int, double> > cap_map;
//    Preseter::makeMutualInfoArray(cap_map);
//    vector<int> p_0(Params::get_M(), -1);
//    vector<int> p(Params::get_M(), -1);
//    Preseter::set_params(cap_map,A,Ac,p_0,p);

    EXP_MODE em = Params::get_exp_mode();
    calcBER();

    return 0;
}