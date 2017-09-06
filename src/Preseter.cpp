#include "../lib/Preseter.h"

void Preseter::preset_A(vector<int> &free){
    Preseter::defineFixedAndFree(free);
}

void Preseter::preset_u(SOURCE_TYPE mode, vector<int> &u){
    u = Preseter::generateUi(mode, u);
}

void Preseter::represet_A(vector<int> &free, vector<pair<int,double> > &cap_map){
    for(int i = 0; i < Params::get_K(); i++){
        free[i] = cap_map[i].first;
    }
}

void Preseter::defineFixedAndFree(vector<int> &free){
    vector<pair<int, double> > cap_map;
    Preseter::makeMutualInfoArray(cap_map);
    for(int i = 0; i < Params::get_K(); i++){
        free[i] = cap_map[i].first;
    }
}

void Preseter::makeMutualInfoArray(vector<pair<int, double> > &cap_map){
    vector<double> cap(Params::get_N());

    if(Params::get_s() == BEC) {
        Analysor::makeArrayCapacity(cap);
//        Analysor::makeArrayBhat(cap);
    } else {
        Analysor::makeArrayBhat(cap);
    }
    for(int i=0; i<Params::get_N(); i++){
        cap_map.push_back(pair<int, double>(i, cap[i]));
    }

    //昇順ソート
    if(Params::get_s() == BEC) {
        sort(begin(cap_map), end(cap_map), Common::sort_greater);
    } else {
        sort(begin(cap_map), end(cap_map), Common::sort_less);
    }
}

vector<int> Preseter::generateUi(SOURCE_TYPE set, vector<int> &x){
    vector<int> ret;
    init_genrand((int) time(NULL));
//    srand(0);
    for (int i = 0; i < Params::get_N(); i++) {
        if (set == ALL0) {
            ret.push_back(0);
        } else if (set == ALL1) {
            ret.push_back(1);
        } else if (set == RAND){
            ret.push_back(genrand_int32() % 2);
        }
    }
    return ret;
}