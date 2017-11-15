#include "../lib/Preseter.h"

void Preseter::preset_A_Ac(vector<int> &free, vector<int> &fixed){
    Preseter::defineFixedAndFree(free,fixed);
}

void Preseter::preset_u(SOURCE_TYPE mode, vector<int> &u){
    u = Preseter::generateUi(mode, u);
}

void Preseter::represet_A(vector<int> &free, vector<int> &fixed_0, vector<pair<int,double> > &cap_map){
    int count = 0;
    int i = 0;
    while(count < Params::get_K()){
        //初期frozen bitに含まれないものだけをdata bitにする
        if( !Common::containVal(cap_map[i].first, fixed_0)){
            free[count] = cap_map[i].first;
            count++;
        }
        i++;
    }
    int j = 0;
    for(int i = 0; i < Params::get_N(); i++){
        if (!Common::containVal(cap_map[Params::get_N()-1-i].first, fixed_0) && !Common::containVal(cap_map[Params::get_N()-1-i].first, free)){
            //fixedの方は相互情報量の低い順に入れていく
            fixed_0.push_back(cap_map[Params::get_N()-1-j].first);
            j++;
        }
    }
}

void Preseter::defineFixedAndFree(vector<int> &free, vector<int> &fixed){
    vector<pair<int, double> > cap_map;
    Preseter::makeMutualInfoArray(cap_map);
    vector<int> temp;
    int j = 0;
    for(int i = 0; i < Params::get_N(); i++){
        if (i < Params::get_K()){
            free.push_back(cap_map[i].first);
        } else {
            //fixedの方は相互情報量の低い順に入れていく
            fixed.push_back(cap_map[Params::get_N()-1-j].first);
            j++;
        }
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


inline vector<int>Preseter::makeTable(int n){
    vector<int> ret(n);
    if (n == 1) {
        ret[0] = 1;
    } else {
        vector<int> temp = makeTable(n/2);
        for (int i = 0; i < n ; i++) {
            if(i < n/2){
                ret[i] = 2*temp[i]-1;
            } else {
                ret[i] = 2*temp[i-n/2];
            }
        }
    }
    return ret;
}

vector<int> Preseter::get_bitReversal(vector<int> p_0){
    vector<int> p;
    vector<int> table = makeTable(Params::get_N());
    for (int i = 0; i < p_0.size(); i++) {
        p.push_back(table[p_0[i]-1]-1);
    }
    return p;
}

void Preseter::set_params(vector<pair<int, double> > &cap_map,vector<int> &A, vector<int> &Ac_0, vector<int> &p_0, vector<int> &p){
    A.resize(Params::get_K(), -1);
//    Ac_0.resize(Params::get_N()-Params::get_K(), -1);
    Ac_0 = {};
    EXP_MODE em = Params::get_exp_mode();

    //Aとして使いたくないものはAc_0にいれる
    switch (em){
        case NORMAL:
        case MID:
            Preseter::represet_A(A, Ac_0, cap_map);
            break;
        case PUNC:
            Preseter::represet_A(A, Ac_0, cap_map);
            for (int i = 0; i < Params::get_M(); i++) {
                p[i] = Ac_0[i];
            }
            break;
        case QUP:
        case M_QUP:
            for (int i = 0; i < Params::get_M(); i++) {
                p_0[i] = i+1;
            }
            p = Preseter::get_bitReversal(p_0);
            Preseter::represet_A(A, Ac_0, cap_map);
            break;
        case WANG:
        case M_WANG:
            for (int i = 0; i < Params::get_M(); i++) {
                p_0[i] = Params::get_N() - i;
                Ac_0.push_back(Params::get_N() - i - 1);
            }
            p = Preseter::get_bitReversal(p_0);
            Preseter::represet_A(A, Ac_0, cap_map);
            break;
        case VALERIO_P:
        case M_VALERIO_P:
            for (int i = 0; i < Params::get_M(); i++) {
                p_0[i] = i+1;
            }
            p = Preseter::get_bitReversal(p_0);
            Ac_0 = p;
            Preseter::represet_A(A, Ac_0, cap_map);
            break;
        case VALERIO_S:
        case M_VALERIO_S:
            for (int i = 0; i < Params::get_M(); i++) {
                p_0[i] = Params::get_N() - i;
            }
            p = Preseter::get_bitReversal(p_0);
            Ac_0 = p;
            Preseter::represet_A(A, Ac_0, cap_map);
            break;
        default:
            break;
    }
}