
#ifndef CHANNEL_POLARIZATION_PARAMS_H
#define CHANNEL_POLARIZATION_PARAMS_H

#include <string>
#include <vector>

using namespace std;

enum CHANNEL_TYPE{BEC, BSC, AWGN};
enum SOURCE_TYPE{ALL0, ALL1, RAND};
enum MID_MODE{MID_BLUTE, MID_ADOR, MID_AOR, MID_DOR, MID_AOB, MID_DOB, MID_AOV, MID_DOV};
enum DECODE_MODE{BP, SC};
enum EXP_MODE{NORMAL, PUNC, QUP, WANG, MID, M_WANG, M_QUP, VALERIO_P, VALERIO_S, M_VALERIO_P, M_VALERIO_S};

const static double inf_p =  50.0;
const static double inf_m = -50.0;
const static double send_0 = -1.0;

#define PRINT(X) cout << #X << ":\n" << setprecision(10) << X << endl << endl;
#define ARR(array)     (sizeof(array) / sizeof(array[0]));

/*
 * 1024,2048,4096,8192,16384,32768
 * 65536,131072,262144,524288,1048576
 */


class Params{
private:
    static int N;
    static int M;
    static int K;
    static int MN;
    static int Bloop;
    static double point[3];
    static double awgn_p[3];
    static double e;
    static int rp;
    static int monteNum;
    static int blockNum;
    static int upperBlockErrorNum;
    static bool is_outlog;
    static bool is_exp_awgn;
    static bool is_calc_bloop;
    static CHANNEL_TYPE s;
    static DECODE_MODE dm;
    static EXP_MODE exp_mode;
    static MID_MODE m_mode;
    static string rvbDir;
    static vector<vector<bool> > T;
    static vector<int> A;
    static vector<int> Ac;
    static vector<int> p;
public:
    static int get_N();
    static int get_M();
    static int get_K();
    static int get_MN();
    static int get_Bloop();
    static void get_point(double temp[3]);
    static void get_awgn_p(double temp[3]);
    static void get_T(vector<vector<bool> > &temp);
    static void get_A(vector<int> &temp);
    static void get_Ac(vector<int> &temp);
    static void get_p(vector<int> &temp);
    static double get_e();
    static int get_monteNum();
    static int get_blockNum();
    static int get_upperBlockErrorNum();
    static bool get_is_outlog();
    static bool get_is_exp_awgn();
    static bool get_is_calc_bloop();
    static string get_rvbDir();
    static CHANNEL_TYPE get_s();
    static DECODE_MODE get_decode_mode();
    static EXP_MODE get_exp_mode();
    static MID_MODE get_m_mode();
    static int get_rp();
    static void set_N(int _N);
    static void set_M(int _M);
    static void set_K(int _K);
    static void set_T(vector<vector<bool> > &temp);
    static void set_A(vector<int> &temp);
    static void set_Ac(vector<int> &temp);
    static void set_p(vector<int> &temp);
    static void set_MN(int _MN);
    static void set_Bloop(int _Bloop);
    static void set_point(double _point[3]);
    static void set_awgn_p(double _awgn_p[3]);
    static void set_is_outlog(bool _is_outlog);
    static void set_is_exp_awgn(bool _is_exp_awgn);
    static void set_is_calc_bloop(bool _is_calc_bloop);
    static void set_monteNum(int _monteNum);
    static void set_blockNum(int _blockNum);
    static void set_decode_mode(DECODE_MODE _dm);
    static void set_exp_mode(EXP_MODE _exp_m);
    static void set_m_mode(MID_MODE _m_mode);
    static void set_upperBlockErrorNum(int _upperBlockErrorNum);
    static void set_e(double _e);
    static void set_rvbDir(string rvbDir);
    static void set_s(CHANNEL_TYPE _s);
    static void set_rp(double _rp);
};

#endif //CHANNEL_POLARIZATION_PARAMS_H
