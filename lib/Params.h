
#ifndef CHANNEL_POLARIZATION_PARAMS_H
#define CHANNEL_POLARIZATION_PARAMS_H

#include <string>

using namespace std;

enum CHANNEL_TYPE{BEC, BSC, AWGN};
enum SOURCE_TYPE{ALL0, ALL1, RAND};
enum MODE{TEST, RUN};
enum DECODE_MODE{BP, SC};
enum EXP_MODE{NORMAL, PUNC, QUP, WANG, MID, M_WANG, M_QUP, VALERIO_P, VALERIO_S, M_VALERIO_P, M_VALERIO_S};

const static double inf_p =  30.0;
const static double inf_m = -30.0;

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
    static double e;
    static int rp;
    static int monteNum;
    static int blockNum;
    static int upperBlockErrorNum;
    static bool is_outlog;
    static CHANNEL_TYPE s;
    static DECODE_MODE dm;
    static EXP_MODE exp_mode;
    static string rvbDir;
public:
    static int get_N();
    static int get_M();
    static int get_K();
    static double get_e();
    static int get_monteNum();
    static int get_blockNum();
    static int get_upperBlockErrorNum();
    static bool get_is_outlog();
    static string get_rvbDir();
    static CHANNEL_TYPE get_s();
    static DECODE_MODE get_decode_mode();
    static EXP_MODE get_exp_mode();
    static int get_rp();
    static void set_N(int _N);
    static void set_M(int _M);
    static void set_K(int _K);
    static void set_is_outlog(bool _is_outlog);
    static void set_monteNum(int _monteNum);
    static void set_blockNum(int _blockNum);
    static void set_decode_mode(DECODE_MODE _dm);
    static void set_exp_mode(EXP_MODE _exp_m);
    static void set_upperBlockErrorNum(int _upperBlockErrorNum);
    static void set_e(double _e);
    static void set_rvbDir(string rvbDir);
    static void set_s(CHANNEL_TYPE _s);
    static void set_rp(double _rp);
};

#endif //CHANNEL_POLARIZATION_PARAMS_H
