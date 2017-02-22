#include "../lib/Preseter.h"

void Preseter::preset(SOURCE_TYPE mode, vector<int> &u, vector<int> &fixed, vector<int> &free){
    Preseter::defineFixedAndFree(fixed,free);
    u = Preseter::generateUi(mode, u, fixed, free);
}

void Preseter::defineFixedAndFree(vector<int> &fixed, vector<int> &free){
    vector<double> cap(Params::get_N(),0);
    vector<double> cap_desc(Params::get_N(),0);
    cout << Params::get_N();
    Analysor::makeArrayCapacityForBec(cap);

    vector<pair<int, double> > cap_map;
    for(int i=0; i<Params::get_N(); i++){
        cap_map.push_back(pair<int, double>(i, cap[i]));
    }

    //昇順ソート
    sort(begin(cap_map), end(cap_map), Common::sort_greater);

    int i = 0;
    for(auto val : cap_map){
        i < Params::get_K() ? free.push_back(val.first) : fixed.push_back(val.first);
        i++;
    }
}

vector<int> Preseter::generateUi(SOURCE_TYPE set, vector<int> &x, vector<int> &u_Ac, vector<int> &A){
    vector<int> ret(Params::get_N(),0);
    srand((int) time(NULL));
    for (int i = 0; i < Params::get_N(); i++) {
        cout << i;
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