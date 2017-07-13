#include "../lib/Preseter.h"

void Preseter::preset_A(vector<int> &free){
    Preseter::defineFixedAndFree(free);
}

void Preseter::preset_u(SOURCE_TYPE mode, vector<int> &u){
    u = Preseter::generateUi(mode, u);
}


void Preseter::defineFixedAndFree(vector<int> &free){
    vector<double> cap(Params::get_N());
    vector<pair<int, double> > cap_map;

    if(1) {
        Analysor::makeArrayCapacity(cap);
//        Analysor::makeArrayBhat(cap);
    } else if (Params::get_s() == BSC || 0) {
        Analysor::makeArrayBhat(cap);
    }
    for(int i=0; i<Params::get_N(); i++){
        cap_map.push_back(pair<int, double>(i, cap[i]));
    }

    //昇順ソート
    if(1) {
        sort(begin(cap_map), end(cap_map), Common::sort_greater);
    } else if (Params::get_s() == BSC || 0) {
        sort(begin(cap_map), end(cap_map), Common::sort_less);
    }
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