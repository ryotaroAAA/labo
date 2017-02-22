#include "../lib/Common.h"

void Common::bar(){
    cout << "======================" << endl;
}

bool Common::containVal(int value, vector<int> m_array) {
    auto it = find(begin(m_array), m_array.end(), value);
    size_t index = distance(begin(m_array), it);
    if(index == m_array.size()) {
        return false;
    }
    return true;
}

bool Common::sort_less(const ass_arr& left,const ass_arr& right){
    return left.second < right.second;
}

bool Common::sort_greater(const ass_arr& left,const ass_arr& right){
    return left.second > right.second;
}

int Common::compare_asc(const void *x, const void *y) {
    const double* a=(const double*)x;
    const double* b=(const double*)y;
    if(*a>*b)return 1;
    if(*a<*b)return -1;
    return 0;
}

int Common::compare_desc(const void *x, const void *y) {
    const double* a=(const double*)x;
    const double* b=(const double*)y;
    if(*a>*b)return -1;
    if(*a<*b)return 1;
    return 0;
}

void Common::dispArray(double *x){
    for(int i = 0; i < Params::get_N(); i++){
        printf("x[%d] = %f\n",i, x[i]);
    }
}
void Common::dispArray(int n, double *x){
    for(int i = 0; i < n; i++){
        printf("x[%d] = %f\n",i, x[i]);
    }
}

void Common::dispArray(int n, int *x){
    for(int i = 0; i < n; i++){
        printf("x[%d] = %d\n",i, x[i]);
    }
}

void Common::pp(vector<int> &a){
    bar();
    for(auto var : a){
        cout << var << endl;
    }
}

void Common::pp(vector<double> &collection){
    for(auto var : collection){
        cout << var << endl;
    }
}

void Common::outputArray(double *x){
    string filename = "/Users/ryotaro/labo/symmetric_capacity";
    ofstream w_file;
    w_file.open(filename, ios::out);
    for (int i = 0; i<Params::get_N(); i++)
    {
        w_file << i << " " << x[i]<< endl;
    }
}

vector<int> Common::index_o(vector<int> &x){
    vector<int> ret;
    int i = 0;
    for(auto val : x){
        if (i % 2 == 0) {
            ret.push_back(val);
        }
        i++;
    }
    return ret;
}

vector<int> Common::index_e(vector<int> &x){
    vector<int> ret;
    int i = 0;
    for(auto val : x){
        if (!(i % 2 == 0)) {
            ret.push_back(val);
        }
        i++;
    }
    return ret;
}

vector<int> Common::retBinary(vector<int> &x) {
    vector<int> ret;
    for(auto val : x){
        ret.push_back(val % 2);
    }
    return ret;
}

int Common::ithIndexDesc(int i, double *array, double *descArray){
    for (int j=0; j<Params::get_N(); j++) {
        if(descArray[i] == array[j]){
            return j;
        }
    }
    return 0;
}


bool Common::containNumInArray(int i, int n, vector<int> &array){
    for(auto val: array ){
        if(i == val){
            return true;
        }
    }
    return false;
}
