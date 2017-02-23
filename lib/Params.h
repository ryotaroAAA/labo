
#ifndef CHANNEL_POLARIZATION_PARAMS_H
#define CHANNEL_POLARIZATION_PARAMS_H
enum CHANNEL_TYPE{BEC};
enum SOURCE_TYPE{ALL0, ALL1, RAND};
enum MODE{ORD,TEST};

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
public:
    static int get_N();
    static int get_K();
    static double get_e();
    static void set_N(int _N);
    static void set_K(int _K);
    static void set_e(double _e);

};

#endif //CHANNEL_POLARIZATION_PARAMS_H
