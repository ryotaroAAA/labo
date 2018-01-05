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
    while(count < Params::get_K() && i < cap_map.size()){
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

void Preseter::makeManyValTableAs(bool sortflg, vector<int> &table){
    vector<int> u(Params::get_N());
    vector<int> x(Params::get_N());
    Preseter::preset_u(ALL1, u);

    Encoder encoder;
    x = encoder.enc(Params::get_N(),u);
//    Common::pp(u);
//    Common::pp(x);
    vector<pair<int, int> > table_map;;

    for(int i=0; i<Params::get_N(); i++){
        table_map.push_back(pair<int, int>(i, x[i]));
    }

    //昇順ソート
    if(sortflg) {
        sort(begin(table_map), end(table_map), Common::sort_greater);
    } else {
        sort(begin(table_map), end(table_map), Common::sort_less);
    }

    for (int i = 0; i < Params::get_N(); i++) {
        table[i] = table_map[i].first;
    }
}

vector<int> Preseter::generateUi(SOURCE_TYPE set, vector<int> &x){
    vector<int> ret;
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

inline vector<int>Preseter::makeTreeIndex(int n){
    vector<int> ret(n);
    if (n == 1) {
        ret[0] = 1;
    } else {
        for (int i = 0; i < n/2 ; i++) {
            ret[2*i] = i+1;
            ret[2*i+1] = i+1 + n/2;
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

void Preseter::set_zin_all(vector<int> &A, vector<vector<int> > &sortedZn, vector<vector<bool> > &ym_isReceived){
    //1. 各ノードを一次元で解釈し, 全ノードのバタチャリア求める [] = z^i_2n
    //長さは(2log2+1)*N
    int size = log2(Params::get_N())+1;
    int tempSize = size*Params::get_N();
    int tempN = Params::get_N();
    vector<double> tempZn;
    vector<double> cap;
    int n = Params::get_N();
    for (int i = 0; i < size; i++) {
        cap.resize(n);
        Params::set_N(n);
        Analysor::makeArrayBhat(cap);
//        for(int j = 0; j < cap.size(); j++){
        int tempi = pow(2,i)-1, j = 0;
        while(j < cap.size()){
            tempZn.push_back(cap[j]);
            j++;
            if (j == cap.size()) {
                if(tempi == 0){
                    break;
                } else {
                    tempi--;
                    j = 0;
                }
            }
        }
        n = n/2;
    }
    Params::set_N(tempN);

    //2. ソートする　[] = 順位
    vector<pair<int, double> > zin_map;
    for(int i=0; i<tempZn.size(); i++){
        zin_map.push_back(pair<int, double>(i, tempZn[i]));
    }
    sort(begin(zin_map), end(zin_map), Common::sort_greater);
//    Common::pp(tempZn);

    int N = Params::get_N();
    int sorted_i = 0;
    int count = 1;
    vector<int> temp_i;
    for (int i = 0; i < zin_map.size(); i++) {
        sorted_i = zin_map[i].first;
        if( (sorted_i < Params::get_N() && Common::containVal(sorted_i, A)) || sorted_i >= Params::get_N()) {
            temp_i.push_back(sorted_i);
            if(count == Params::get_MN()) break;
            count++;
        }
    }

    //3. 一次元データを多次元に変換し, 上から送信フラグつけていく [][] = 順位
    //[i] => [i/N][i%N]として変換可能なはず
    for (int i = 0; i < Params::get_MN(); i++) {
        sorted_i = temp_i[i];
        ym_isReceived[sorted_i/N][sorted_i%N] = true;
    }
}


void Preseter::set_params(vector<pair<int, double> > &cap_map, vector<int> &A, vector<int> &Ac_0, vector<int> &p_0, vector<int> &p){
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
    Params::set_A(A);
    Params::set_Ac(Ac_0);
    Params::set_p(p);

    //中間ノード設定, yからメッセージを受信するようになる
    int size = 2*log2(Params::get_N())+2;
    int n = Params::get_N();
    vector<vector<bool> > T(size, vector<bool>(n, false));
    if( Common::is_mid_send() ){
        MID_MODE mm = Params::get_m_mode();
        vector<int> Bn_0 = Preseter::makeTable(Params::get_N());
        vector<int> v_table_0(Params::get_N());
        cout << "size::" << size << endl;

        bool sortflg = (Params::get_m_mode() == MID_DOV);
        Preseter::makeManyValTableAs(sortflg, v_table_0);
        vector<int> Bn, v_table;
        for (int i = 0; i < Params::get_N(); i++) {
            if(Common::containVal(Bn_0[i]-1, A)){
                Bn.push_back(Bn_0[i]-1);
            }
            if(Common::containVal(v_table_0[i], A)){
                v_table.push_back(v_table_0[i]);
            }
        }
        int k = Params::get_K();
        int zsize = log2(Params::get_N())+1;
        int tmpSize = log2(Params::get_N())+1;
        vector<vector<int> > sortedZn(zsize, vector<int>(Params::get_N(), 0));

        vector<int> bl;
        string bloop_fn = "save_b/N_"
                          + to_string(Params::get_N())
                          + "/Bloop_" + to_string(Params::get_Bloop());
        ifstream ifs(bloop_fn);
        string str;
        while(getline(ifs,str)) {
            cout<< str << endl;
            bl.push_back(stoi(str));
        }
        switch (mm) {
            case MID_BLUTE:{
                int t = 0;
                for (int i = 0; i < bl.size(); i++) {
                    if(t <= Params::get_MN() && bl.size()>0){
                        //if(Common::containVal(bl[t],Ac) && Common::containVal(bl[t],p)){ ?　左ノードならAに含まれるもののみ
                        if((bl[i] < n && bl[i] && Common::containVal(bl[i],A)) || n <= bl[i]){
                            T[(bl[i]/n)*2][bl[i]%n] = true;
//                            cout << "[" << bl[t]/n << "]" << "[" << bl[t]%n << "]" << endl;
                            t++;
                        }
                    }
                }
            }
                break;
            case MID_ADOR:
                Preseter::set_zin_all(A, sortedZn, T);
                break;
            case MID_DOR:
                for (int i = 0; i < Params::get_MN(); i++) {
                    T[0][A[i]] = true;
                }
                break;
            case MID_AOR:
                for (int i = 0; i < Params::get_MN(); i++) {
                    T[0][A[k-1-i]] = true;
                }
                break;
            case MID_DOB:
                for (int i = 0; i < Params::get_MN(); i++) {
                    T[0][Bn[k-1-i]-1] = true;
                }
                break;
            case MID_AOB:
                for (int i = 0; i < Params::get_MN(); i++) {
                    T[0][Bn[i]-1] = true;
                }
                break;
            case MID_DOV:
                for (int i = 0; i < Params::get_MN(); i++) {
                    T[0][v_table[i]] = true;
                }
                break;
            case MID_AOV:
                for (int i = 0; i < Params::get_MN(); i++) {
                    T[0][v_table[k-1-i]] = true;
                }
                break;
            default:
                break;
        }
    }
    Params::set_T(T);
}