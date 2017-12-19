#ifndef CHANNEL_POLARIZATION_COMMON_H
#define CHANNEL_POLARIZATION_COMMON_H
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/stat.h>
#include <fstream>
#include <random>
#include "Params.h"
#include "sys/stat.h"
#include <unistd.h>
#include <sstream>

using namespace std;

class Common{
public:
    typedef pair<int, double> ass_arr;

    static void bar();
    static void outputArray(double *x);
    static void dispArray(double *x);
    static void dispArray(int n,int *x);
    static void dispArray(int n,double *x);
    static void pp(vector<int> &a);
    static void rpp(vector<int> &a);
    static void dpp(vector<int> &a);
    static void drpp(vector<int> &a);
    static void pp(vector<double> &a);

    static int ithIndexDesc(int i, double *array, double *descArray);
    static bool containVal(int value, vector<int> m_array);
    static bool sort_less(const ass_arr& left,const ass_arr& right);
    static bool sort_greater(const ass_arr& left,const ass_arr& right);
    static int compare_int(const void *a, const void *b);
    static int compare_desc(const void *x, const void *y);
    static int compare_asc(const void *x, const void *y);
    static double gauss_dist(double x, double mean, double sigma);

    static bool containNumInArray(int i, int n, vector<int> &array);
    static vector<int> index_o(int n, vector<int> &x);
    static vector<int> index_e(int n, vector<int> &x);
    static vector<int> retBinary(int n, vector<int> &x);

    static bool is_mid_send();
    static double get_rate();
};
#endif //CHANNEL_POLARIZATION_COMMON_H