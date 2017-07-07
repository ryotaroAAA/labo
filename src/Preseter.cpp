#include "../lib/Preseter.h"

void Preseter::preset(SOURCE_TYPE mode, vector<int> &u, vector<int> &free){
    Preseter::defineFixedAndFree(free);
    u = Preseter::generateUi(mode, u);
}


void Preseter::defineFixedAndFree(vector<int> &free){
    vector<double> cap;
    vector<double> cap_desc;
    vector<pair<int, double> > cap_map;
    if(Params::get_s() == BEC) {
        Analysor::makeArrayCapacityForBec(cap);
    } else if (Params::get_s() == BSC) {
        Analysor::makeArrayCapacityForBec(cap);
    }
    for(int i=0; i<Params::get_N(); i++){
        cap_map.push_back(pair<int, double>(i, cap[i]));
    }
    //昇順ソート
    sort(begin(cap_map), end(cap_map), Common::sort_greater);
    for(int i = 0; i < Params::get_K(); i++){
        free[i] = cap_map[i].first;
    }
}

vector<int> Preseter::generateUi(SOURCE_TYPE set, vector<int> &x){
    vector<int> ret;
    srand((int) time(NULL));
//    srand(0);
    for (int i = 0; i < Params::get_N(); i++) {
        if (set == ALL0) {
            ret.push_back(0);
        } else if (set == ALL1) {
            ret.push_back(1);
        } else if (set == RAND){
            ret.push_back(rand() % 2);
        }
    }
    return ret;
}