
#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include <chrono>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

//params
const int  N=8;    //number of
const int  K=3;     //number of data bits
const double e = 0.5f;
static int hoge = 0;
static int hoge2 = 0;
static double hogetime = 0.0;
enum CHANNEL_TYPE{BEC};
enum SOURCE_TYPE{ALL0, ALL1, RAND};

#define PRINT(X) cout << #X << ":\n" << setprecision(10) << X << endl << endl;
#define ARR(array)     (sizeof(array) / sizeof(array[0]));

class Common{
public:
    typedef pair<int, double> ass_arr;

    static void outputArray(double *x);
    static void dispArray(double *x);
    static void dispArray(int n,int *x);
    static void dispArray(int n,double *x);

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

