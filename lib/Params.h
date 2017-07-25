
#ifndef CHANNEL_POLARIZATION_PARAMS_H
#define CHANNEL_POLARIZATION_PARAMS_H

#include <string>

using namespace std;

enum CHANNEL_TYPE{BEC, BSC};
enum SOURCE_TYPE{ALL0, ALL1, RAND};
enum MODE{TEST, RUN};

#define PRINT(X) cout << #X << ":\n" << setprecision(10) << X << endl << endl;
#define ARR(array)     (sizeof(array) / sizeof(array[0]));

/*
 * 1024,2048,4096,8192,16384,32768
 * 65536,131072,262144,524288,1048576
 */


class Params{
private:
    static int N;
    static int K;
    static double e;
    static int monteNum;
    static int blockNum;
    static int upperBlockErrorNum;
    static CHANNEL_TYPE s;
    static string rvbDir;
public:
    static int get_N();
    static int get_K();
    static double get_e();
    static int get_monteNum();
    static int get_blockNum();
    static int get_upperBlockErrorNum();
    static string get_rvbDir();
    static CHANNEL_TYPE get_s();
    static void set_N(int _N);
    static void set_K(int _K);
    static void set_monteNum(int _monteNum);
    static void set_blockNum(int _blockNum);
    static void set_upperBlockErrorNum(int _upperBlockErrorNum);
    static void set_e(double _e);
    static void set_rvbDir(string rvbDir);
    static void set_s(CHANNEL_TYPE _s);
};

#endif //CHANNEL_POLARIZATION_PARAMS_H
