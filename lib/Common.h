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
#include "Params.h"
using namespace std;

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
