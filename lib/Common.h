#ifndef CHANNEL_POLARIZATION_COMMON_H
#define CHANNEL_POLARIZATION_COMMON_H
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

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
public:
    static int N;
    static int K;
    static double e;
    Params(int _N, int _K, double _e){
        N = _N;
        K = _K;
        e = _e;
    }
    static int get_N();
    static int get_K();
    static double get_e();
    static void set_N(int _N);
    static void set_K(int _K);
    static void set_e(double _e);

//    constexpr static int N = 1024;
//    constexpr static int K = 1000;
//    constexpr static double e = 0.5;
};
int Params::N = 0;
int Params::K = 0;
double Params::e = 0;

int Params::get_N(){
    return N;
}

int Params::get_K(){
    return K;
}

double Params::get_e(){
    return e;
}

void Params::set_N(int _N){
    N = _N;
}

void Params::set_K(int _K){
    K = _K;
}

void Params::set_e(double _e){
    e = _e;
}

class Common{
public:
    typedef pair<int, double> ass_arr;

    static void bar();
    static void outputArray(double *x);
    static void dispArray(double *x);
    static void dispArray(int n,int *x);
    static void dispArray(int n,double *x);
    static void pp(vector<int> &collection);
    static void pp(vector<double> &a);

    static int ithIndexDesc(int i, double *array, double *descArray);
    static bool containNumInArray(int i, int n, vector<int> &array);
    static bool containVal(int value, vector<int> m_array);
    static bool sort_less(const ass_arr& left,const ass_arr& right);
    static bool sort_greater(const ass_arr& left,const ass_arr& right);
    static int compare_int(const void *a, const void *b);
    static int compare_desc(const void *x, const void *y);
    static int compare_asc(const void *x, const void *y);

    static vector<int> index_o(vector<int> &x);
    static vector<int> index_e(vector<int> &x);
    static vector<int> retBinary(vector<int> &x);
};

#endif //CHANNEL_POLARIZATION_COMMON_H
