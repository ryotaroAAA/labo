
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

void calcBER(int mode){
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    Analysor analysor;

    EXP_MODE em = Params::get_exp_mode();
    MID_MODE mm = Params::get_m_mode();
    string ename;
    string mname;
    switch (mm) {
        case MID_BLUTE: {
            int bl = 0;
            if(Params::get_Bloop()){
                bl = Params::get_Bloop();
            }
            mname = "mid_blute" + to_string(bl);
            break;
        }
        case MID_ADOR: mname = "mid_ador"; break;
        case MID_AOR: mname = "mid_aor"; break;
        case MID_DOR: mname = "mid_dor"; break;
        case MID_AOB: mname = "mid_aob"; break;
        case MID_DOB: mname = "mid_dob"; break;
        case MID_AOV: mname = "mid_aov"; break;
        case MID_DOV: mname = "mid_dov"; break;
    }
    switch (em) {
        case NORMAL: ename = ""; break;
        case PUNC: ename = "punc"; break;
        case MID:
            ename = mname;
            break;
        case QUP: ename = "qup"; break;
        case WANG: ename = "wang"; break;
        case M_QUP: ename = "m_qup_"+mname; break;
        case M_WANG: ename = "m_wang_"+mname; break;
        case VALERIO_P: ename = "valerio_p"; break;
        case VALERIO_S: ename = "valerio_s"; break;
        case M_VALERIO_P: ename = "m_valerio_p_"+mname; break;
        case M_VALERIO_S: ename = "m_valerio_s_"+mname; break;
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
             << "MNum" + to_string(Params::get_monteNum())
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
    if(mode){
        analysor.calcBlockErrorRate_BP();
    } else {
        analysor.calcBlockErrorRate();
    }

    performance.stopTimer();

    performance.outHMS();
    cout << Params::get_rvbDir() << endl;
}
inline string get_current_directory() {
    char cwd[1024];
    getcwd(cwd, sizeof(cwd));
    return cwd;
}
bool checkFileExistence(const std::string& str) {
    std::ifstream ifs(str);
    return ifs.is_open();
}
void test_edge(){
    int size = 2*log2(Params::get_N())+2;
    cout << size << endl;
    int N = Params::get_N();
    int level;
    int index;
    Decoder d;

    cout << "check" << endl;
    for (int i = 1; i < size; i=i+2) {
        cout << i+1 << "::";
        //チェックノードの右側(できてそう)
        for (int j = 0; j < N; ++j) {
            level = i+1;
            index = j+1;
            vector<vector<int> > adjacent;
            vector<int> temp;
            //部分B_Nのn
            int n = (level != 1) ? Params::get_N()/pow(2,((level-2)/2)) : 1;
            temp = d.makeBPTreeIndex(n);
            int temp_i = (index == n) ? n - 1 : (index - 1) % n;
            int div = (index - 1) / n;
            cout << temp[temp_i] + div * n << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "val" << endl;
    for (int i = 2; i < size; i=i+2) {
        cout << i+1 << "::";
        for (int j = 0; j < N; ++j) {
            level = i+1;
            index = j+1;
            vector<vector<int> > adjacent;
            vector<int> temp;
            //部分B_Nのn
            int n = (level != 1) ? Params::get_N()/pow(2,((level-2)/2)) : 1;
            temp = d.makeBPTreeIndex(n);
            int temp_i = 0;
            int div = 0;
            vector<int> temp_r(Params::get_N(), 0);
            if (n < Params::get_N() && n != 1) {
                for (int i = 0; i < Params::get_N(); i++) {
                    temp_r[i] = temp[i % n];
                }
                for (int i = 0; i < Params::get_N(); i++) {
                    if (temp_r[i] == (index - 1) % n + 1) {
                        div = (index - 1) / n;
                        temp_i = i + 1 + div * n;
                        break;
                    }
                }
            } else {
                for (int i = 0; i < Params::get_N(); i++) {
                    if (temp[i] == index) {
                        temp_i = i + 1;
                        break;
                    }
                }
            }
            cout << temp_i << " ";
        }
        cout << endl;
    }
}
void test_dec(){
    Performance performance;
    Analysor analysor;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    logger.setRvbDir(Params::get_rvbDir());

    int n = Params::get_N(), error_count = 0, block_error_count = 0, loopi = 1, itr = 0;
    int size = log2(Params::get_N()) + 1;
    double BER = 0.0, sumBER = 0.0, rate = 0.0, mitr = 0.0;

    vector<int> A;
    vector<int> Ac;
    vector<vector<int> > B(2*log2(Params::get_N())+2, vector<int>(Params::get_N(),0));
    vector<int> u(Params::get_N(), 0);
    vector<int> x(Params::get_N(), 0);
    vector<double> y(Params::get_N(), 0);
    vector<vector<int> > tmp_x(log2(Params::get_N()), vector<int>(Params::get_N(), 0));
    vector<vector<int> > xm(size, vector<int>(Params::get_N(), 0));
    vector<vector<double> > ym(size, vector<double>(Params::get_N(), 0.0));
    vector<int> u_ests(Params::get_N(), 0);
    vector<int> u_estb(Params::get_N(), 0);
    vector<double> tests(Params::get_N(), 0);
    vector<double> testb(Params::get_N(), 0);

//    Preseter::makeMutualInfoArray(cap_map);
//    Preseter::preset_A_Ac(A, Ac);

    vector<int> param(2, 0);
    //右ノードのパンクチャorショートン位置を指定
    vector<int> p_0(Params::get_M(), -1);
    vector<int> p(Params::get_M(), -1);

    vector<pair<int, double> > cap_map;
    vector<vector<int> > node_error_count(2*log2(Params::get_N())+2, vector<int>(Params::get_N(),0));
    ofstream val_error_file;
    string val_error_fn = "/Users/ryotaro/Dropbox/labo/graph_js/val_error.json";
    val_error_file.open(val_error_fn, ios::out);
    val_error_file << "{" << endl;

    ofstream b_file;
    string b_fn = "/Users/ryotaro/Dropbox/labo/graph_js/b.json";
    b_file.open(b_fn, ios::out);
    b_file << "{" << endl;

    Preseter::makeMutualInfoArray(cap_map);
    Preseter::set_params(cap_map, A, Ac, p_0, p);

    rate = Common::get_rate();
//    cout << "Rate " <<  rate << ", K " <<  Analysor::get_eachK(rate) << ", Eb/N0 " << h << ", interval " << interval_x << ", s " << Params::get_s() <<   endl;

    //encode
    Preseter::preset_u(RAND, u, A);
    x = encoder.encode(Params::get_N(), u);
    y = Channel::channel_output(x);

//    u={0,0,1,1};
//    y={-1.64065,0.257906,1.84493,0.278905};

    //SC decode
    Common::bar();
    u_ests = decoder.decode(tests, y, u);
    Analysor::errorCount(u, u_ests, &error_count);
    if(error_count > 0) block_error_count++;
    Common::bar();
    if (loopi % 1 == 0 ) cout << loopi << " " << error_count << " " << block_error_count << " " << (double)block_error_count/loopi << endl;
    Common::bar();
    error_count = 0;
    //BP decode
    u_estb = decoder.calcBP(testb, loopi, param, u, x, y, xm, ym, node_error_count, val_error_file, B);
    Common::bar();
    Analysor::errorCount(u, u_estb, &error_count);
    block_error_count = 0;
    if (error_count > 0) {
        block_error_count++;
    }
    if (loopi % 1 == 0 ) cout << loopi << " " << error_count << " " << block_error_count  << " " << (double)block_error_count/loopi << endl;
    Common::bar();
    val_error_file << "\t}" << endl;
    val_error_file << "}" << endl;
    int count = 0;
    for (int i = 0; i < Params::get_N(); i++) {
        if(testb[i] != tests[i]){
            count++;
            cout << i+1 << " " << testb[i] << " " << tests[i] << endl;
        } else {
            cout << i+1 << " " << testb[i] << " " << tests[i] << endl;
        }
    }
    Common::bar();
    cout << "no match: " << count << endl;
}

//mid>normal
//awgn なら e=0.8
//enum EXP_MODE{NORMAL, PUNC, QUP, WANG, MID, M_WANG, M_QUP, VALERIO_P, VALERIO_S, M_VALERIO_P, M_VALERIO_S};

int main(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    srand(tv.tv_sec + tv.tv_usec);
    init_genrand(tv.tv_sec + tv.tv_usec);

    Params::set_e(0.9);
    Params::set_N(256);
    Params::set_K(1);
    cout << "memo::" << endl;

//    パンクチャ, ショートン
    Params::set_M(96);
//    Params::set_M(4);
//    Params::set_M(56);
//    Params::set_M(16);
    int bp = 1;
    Params::set_MN(40);
    Params::set_s(AWGN);
    //calc_bloopを計算したいときにtrue, それ以外は必ずコメントアウト
    Params::set_is_calc_bloop(true);
//    Params::set_is_outlog(true);
    Params::set_decode_mode(BP);
    Params::set_is_exp_awgn(true);
    //50で十分
    Params::set_monteNum(100);
    Params::set_rp(10);
    Params::set_Bloop(100);
    Params::set_blockNum(4100);
    Params::set_upperBlockErrorNum(4100);
    //point設定
    double p[3] = {1,0.25,0.25};
    Params::set_point(p);

    //awgn_p設定
    //{0.25,1.5,3.0},{0.5,1.5,3.0},{0.75,2.5,4.5}
//    double awgn_p[3] = {0.5,1,7};
//    double awgn_p[3] = {0.75,1.0,1.0};
    double awgn_p[3] = {0.75,2.5,4};
    Params::set_awgn_p(awgn_p);

//    Params::set_m_mode(MID_BLUTE);
//    Params::set_m_mode(MID_DOR);
//    Params::set_m_mode(MID_DOB);
    Params::set_m_mode(MID_DOV);

    Params::set_exp_mode(NORMAL);
//    Params::set_exp_mode(QUP);
//    Params::set_exp_mode(MID);
//    Params::set_exp_mode(WANG);
//    Params::set_exp_mode(VALERIO_P);
//    Params::set_exp_mode(VALERIO_S);
//    Params::set_exp_mode(M_WANG);
//    Params::set_exp_mode(M_QUP);
//    Params::set_exp_mode(M_VALERIO_P);
//    Params::set_exp_mode(M_VALERIO_S);

//    test_edge();
//    test_dec();
    calcBER(bp);

    double llr = 1.0;

    vector<double> val = {0.0 , 1.0/0};
    for (int i = 0; i < val.size(); i++) {
        llr *= tanh((double)val[i]/2.0);
    }
    llr = (double)2.0*atanh(llr);
//    llr = (double)abs(min(val[0],val[1])) * (val[0]*val[1])/(abs(val[0]*val[1]));
//    llr = (val[0]*val[1])/(abs(val[0]*val[1]));
    cout << llr << endl;

    return 0;
}