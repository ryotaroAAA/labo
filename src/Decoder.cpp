#include "../lib/Decoder.h"
#include "../lib/Analysor.h"
Decoder::Decoder(){

}

Decoder::~Decoder(){

}
void Decoder::SCinit(int n, vector<double> &y, vector<int> &u, vector<int> &u_est, vector<vector<double> > &node_val, vector<vector<vector<message> > > &message_list, vector<vector<vector<message> > > &save_list, vector<vector<bool> > &node_isChecked, vector<vector<bool> > &save_isChecked) {
    double wc, wp, llr = 0.0, val = 0.0;
    int size = 2*log2(Params::get_N())+2;
    bool databit_flag = false;
    vector<int> A;
    Params::get_A(A);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < Params::get_N(); j++) {
            if(!save_isChecked[i][j]){
                node_val[i][j] = 0.0;
                node_isChecked[i][j] = false;
            }
        }
    }
    for (int i = 0; i < Params::get_N(); i++) {
        wc = Channel::calcW(y[i],send_0);
        wp = Channel::calcW(y[i],1);
        llr = 1.0 * log(wc / wp);
        node_val[size-1][i] = llr;
//        node_val[size-2][i] = llr;

        //frozen_bitのllr設定
        databit_flag = false;
        for (int j = 0; j < Params::get_K(); j++) {
            if(i == A[j]){
                databit_flag = true;
            }
        }

        //今計算してるbitよりindexが低いbitのみ設定
        if(databit_flag == false && i <= n){
            node_val[0][i] = (u[i] == 0) ? inf_p : inf_m;
            u_est[i] = u[i];
            node_isChecked[0][i] = true;
            save_isChecked[0][i] = true;
        } else if(u_est[i] != 2) {
            //一度決定したbit
            node_val[0][i] = (u_est[i] == 0) ? inf_p : inf_m;
            node_isChecked[0][i] = true;
            save_isChecked[0][i] = true;
        }
    }

    vector<vector<int> > adjacent;
    message temp = {0,0,0.0};
    int adjacent_count=0, level=0, index=0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < Params::get_N(); j++) {
            //隣接ノードの位置
            adjacent = adjacentIndex(i + 1, j + 1);

            for (int k = 0; k < adjacent.size(); k++) {
                level = adjacent[k][0];
                index = adjacent[k][1];
                val = node_val[i][j];
                temp.toLevel = level;
                temp.toIndex = index;
                temp.val = val;
                message_list[i][j].push_back(temp);
                if (n != 0 && val != 0) save_list[i][j][k].val = val;
            }
        }
    }
}
void Decoder::calcSConBP(int itr, int count, ofstream &val_file, ofstream &check_file, vector<int> &u, vector<double> &y , vector<int> &u_n_est, vector<double> &tmp_u, vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked) {
    int size = 2 * log2(Params::get_N()) + 2;

    vector<vector<vector<message> > > save_list(size, vector<vector<message> >(Params::get_N(), vector<message>()));
    vector<vector<double> > tmp_node(size, vector<double>(Params::get_N(), 0.0));
    vector<vector<double> > save_node(size, vector<double>(Params::get_N(), 0.0));
    vector<vector<bool> > save_isChecked(size, vector<bool>(Params::get_N(), false));
    vector<vector<bool> > tmp_isChecked(size, vector<bool>(Params::get_N(), false));

    vector<vector<double> > ym(size-1, vector<double>(Params::get_N(),0.0)); //dummy
    vector<int> A;
    Params::get_A(A);

    int n = Params::get_N();
    vector<vector<bool>> T(size, vector<bool>(n, false));
    Params::get_T(T);

    int tmpSize;
    int count_limit = 1;
    bool is_check_changed_left = false;
    bool is_check_changed_right = false;
    bool tmpflg = true;
    for (int i = 0; i < Params::get_N(); i++) {
        count = 0;
        tmpSize = size-2;
        is_check_changed_left = true;
        is_check_changed_right = true;
        tmpflg = true;

        vector<vector<vector<message> > > message_list(size, vector<vector<message> >(Params::get_N(), vector<message>()));
        SCinit(i, y, u, u_n_est, node_value, message_list, save_list, node_isChecked, save_isChecked);
        if( i == 0 ){
            save_list = message_list;
        }

        for (int j = 0; j < Params::get_rp(); j++) {
            //uにメッセージがきたらlogに書いてu_estに代入する
            if (node_value[0][i] != 0.0 || node_isChecked[0][i]) {
                u_n_est[i] = (node_value[0][i]>=0) ? 0 : 1;
                if(Common::containVal(i,A)){
                    // if(count < count_limit) {
                    if(is_check_changed_left && Params::get_is_outlog()) {
                        printDecodeProgress(itr, node_isChecked, check_file);
                        printDecodeProgress(itr, node_value, val_file);
                        itr++;
                    }
                }
                break;
            }

            //message passing
//                if(count < count_limit) {
            if(is_check_changed_left) {
                //左から
//                calc_mp(tmpSize, node_value, save_list, save_isChecked, ym);
                confirmIsCheck(node_value, save_isChecked);
                is_check_changed_left = isChanged(tmp_isChecked, save_isChecked);
                tmp_isChecked = save_isChecked;
                if(Params::get_is_outlog()){
                    printDecodeProgress(itr, save_isChecked, check_file);
                }
            } else {
                //右から
//                    if(count == log2(Params::get_N())){
                if(tmpflg){
                    tmpSize = size;
                    message_list = save_list;
                    node_isChecked = save_isChecked;
                    tmpflg = false;
                }
//                calc_mp(tmpSize, node_value, message_list, node_isChecked, ym);
                confirmIsCheck(node_value, node_isChecked);

                is_check_changed_right = isChanged(tmp_isChecked, node_isChecked);
                tmp_isChecked = node_isChecked;
                if(Params::get_is_outlog()){
                    printDecodeProgress(itr, node_isChecked, check_file);
                }
            }
            if(Params::get_is_outlog()){
                printDecodeProgress(itr, node_value, val_file);
            }

            itr++;

            count++;
            if((Params::get_decode_mode() == BP && j == Params::get_rp()-1) || !is_check_changed_right) {
                u_n_est[i] = (node_value[0][i]>=0.0) ? 0 : 1;
                if(Common::containVal(i,A)){
                    // if(count < count_limit) {
                    if(Params::get_is_outlog()){
                        if(is_check_changed_left) {
                            printDecodeProgress(itr, save_isChecked, check_file);
                        } else {
                            printDecodeProgress(itr, node_isChecked, check_file);
                        }
                        printDecodeProgress(itr, node_value, val_file);
                    }
                    itr++;
                }
                break;
            }
        }
        tmp_u[i] = node_value[0][i];
        node_value[0][i] = (u_n_est[i]==0) ? inf_p : inf_m;
    }

    if(Params::get_is_outlog()) {
        for (int i = 0; i < Params::get_N(); i++) {
            cout << "SConBP: " << i+1 << " " << tmp_u[i] << endl;
        }
    }
}

bool Decoder::isChanged(vector<vector<bool> > &old_node_isChecked, vector<vector<bool> > &new_node_isChecked) {
    int size = 2*log2(Params::get_N())+2;
    vector<vector<int> > adjacent;

    //check nodeが1loopで変化したかを見る
    for (int i = 1; i < size-1; i=i+2) {
        for (int j = 0; j < Params::get_N(); j++) {
            if( old_node_isChecked[i][j] != new_node_isChecked[i][j] ){
                return true;
            }
        }
    }
    return false;
}
void Decoder::printDecodeProgress(int count, vector<vector<int> > &node_value, ofstream &w_file){
    int size = 2*log2(Params::get_N())+2;
    if(count >0 ){
        w_file << "\t}," << endl;
    }
    w_file << "\t\"" << count << "\"" << ":{" << endl;
    for (int i = 0; i < size; i++) {
        w_file << "\t\t\"" << i << "\"" << ":{" << endl;
        for (int j = 0; j < Params::get_N(); j++) {
            double temp = node_value[i][j];
            if(isinf(temp) && temp > 0){
//                temp = inf_p;
//                w_file << "\t\t\t\"" << j << "\"" << ":" << temp;
                w_file << "\t\t\t\"" << j << "\"" << ":" << "\"inf\"";
            } else if(isinf(temp) && temp < 0){
//                temp = inf_m;
//                w_file << "\t\t\t\"" << j << "\"" << ":" << temp;
                w_file << "\t\t\t\"" << j << "\"" << ":" << "\"-inf\"";
            } else if(isnan(temp)){
                w_file << "\t\t\t\"" << j << "\"" << ":" << "\"NAN\"";
            } else {
                w_file << "\t\t\t\"" << j << "\"" << ":" << temp;
            }
            if (j != Params::get_N()-1) {
                w_file << ",";
            }
            w_file << endl;
        }
        if (i != size-1) {
            w_file << "\t\t}," << endl;
        } else {
            w_file << "\t\t}" << endl;
        }
    }
}
void Decoder::printDecodeProgress(int count, vector<vector<bool> > &node_value, ofstream &w_file){
    int size = 2*log2(Params::get_N())+2;
    if(count >0 ){
        w_file << "\t}," << endl;
    }
    w_file << "\t\"" << count << "\"" << ":{" << endl;
    for (int i = 0; i < size; i++) {
        w_file << "\t\t\"" << i << "\"" << ":{" << endl;
        for (int j = 0; j < Params::get_N(); j++) {
            double temp = node_value[i][j];
            if(isinf(temp) && temp > 0){
//                temp = inf_p;
//                w_file << "\t\t\t\"" << j << "\"" << ":" << temp;
                w_file << "\t\t\t\"" << j << "\"" << ":" << "\"inf\"";
            } else if(isinf(temp) && temp < 0){
//                temp = inf_m;
//                w_file << "\t\t\t\"" << j << "\"" << ":" << temp;
                w_file << "\t\t\t\"" << j << "\"" << ":" << "\"-inf\"";
            } else if(isnan(temp)){
                w_file << "\t\t\t\"" << j << "\"" << ":" << "\"NAN\"";
            } else {
                w_file << "\t\t\t\"" << j << "\"" << ":" << temp;
            }
            if (j != Params::get_N()-1) {
                w_file << ",";
            }
            w_file << endl;
        }
        if (i != size-1) {
            w_file << "\t\t}," << endl;
        } else {
            w_file << "\t\t}" << endl;
        }
    }
}
void Decoder::printDecodeProgress(int count, vector<vector<double> > &node_value, ofstream &w_file){
    int size = 2*log2(Params::get_N())+2;
    if(count >0 ){
        w_file << "\t}," << endl;
    }
    w_file << "\t\"" << count << "\"" << ":{" << endl;
    for (int i = 0; i < size; i++) {
        w_file << "\t\t\"" << i << "\"" << ":{" << endl;
        for (int j = 0; j < Params::get_N(); j++) {
            double temp = node_value[i][j];
            if(isinf(temp) && temp > 0){
//                temp = inf_p;
//                w_file << "\t\t\t\"" << j << "\"" << ":" << temp;
                w_file << "\t\t\t\"" << j << "\"" << ":" << "\"inf\"";
            } else if(isinf(temp) && temp < 0){
//                temp = inf_m;
//                w_file << "\t\t\t\"" << j << "\"" << ":" << temp;
                w_file << "\t\t\t\"" << j << "\"" << ":" << "\"-inf\"";
            } else if(isnan(temp)){
                w_file << "\t\t\t\"" << j << "\"" << ":" << "\"NAN\"";
            } else {
                w_file << "\t\t\t\"" << j << "\"" << ":" << temp;
            }
            if (j != Params::get_N()-1) {
                w_file << ",";
            }
            w_file << endl;
        }
        if (i != size-1) {
            w_file << "\t\t}," << endl;
        } else {
            w_file << "\t\t}" << endl;
        }
    }
}
//level, indexは(1,2,3,4,...)の並びで
double Decoder::take_val(vector<int> &locate, vector<vector<double> > &node_val){
    int level = locate[0];
    int index = locate[1];
    return node_val[level-1][index-1];
}
vector<vector<int> > Decoder::adjacentIndex(int level, int index){
    vector<vector<int> > adjacent;
    vector<int> temp;
    //部分B_Nのn
    int n = (level != 1) ? Params::get_N()/pow(2,((level-2)/2)) : 1;
    temp = makeBPTreeIndex(n);
    int size = 2*log2(Params::get_N())+2;

    //一番右のチャンネルファクター
    if( level ==  size){
        adjacent = {{level - 1, index}};
    }
    //一番右の変数ノード
    else if( level ==  size-1 ) {
        adjacent = {{level - 1, index}, {level + 1, index}};
//        adjacent = {{level - 1, index}};
    }

    //それ以外
    else {
        //チェックノード
        if (level % 2 == 0) {
            int temp_i = (index == n) ? n - 1 : (index - 1) % n;
            int div = (index - 1) / n;
            //index_o
            if (index % 2 == 1) {

                //一番右のval, check2つ以外
                if (level < 2 * log2(Params::get_N()))
                    adjacent = {{level - 1, index},
                                {level - 1, index + 1},
                                {level + 1, temp[temp_i] + div * n}};
                else
                    adjacent = {{level - 1, index},
                                {level - 1, index + 1},
                                {level + 1, index}};
            // /index_e
            } else {
                if (level < 2 * log2(Params::get_N()))
                    adjacent = {{level - 1, index},
                                {level + 1, temp[temp_i] + div * n}};
                else
                    adjacent = {{level - 1, index},
                                {level + 1, index}};
            }
        }//変数ノード
        else {
            int temp_i = 0;
            int div = 0;
            vector<int> temp_r(Params::get_N(), 0);
            if (n < Params::get_N() && n != 1) {
                for (int i = 0; i < Params::get_N(); i++) {
                    temp_r[i] = temp[i % n];
                }
                for (int i = 0; i < Params::get_N(); i++) {
                    if (temp_r[i] == (index - 1) % n + 1) {
                        div = (index - 1) / n;
                        temp_i = i + 1 + div * n;
                        break;
                    }
                }
            } else {
                for (int i = 0; i < Params::get_N(); i++) {
                    if (temp[i] == index) {
                        temp_i = i + 1;
                        break;
                    }
                }
            }

            //index_o
            if (index % 2 == 1) {
                //一番右
                if (level == 2 * log2(Params::get_N()) + 1) {
                    adjacent = {{level - 1, index}};
                }
                //一番左
                else if (level == 1) {
                    adjacent = {{level + 1, index}};
                }
                //それ以外
                else {
                    adjacent = {{level + 1, index},
                                {level - 1, temp_i}};
                }
            }
            //index_e
            else {
                //一番右
                if (level == 2 * log2(Params::get_N()) + 1){
                    adjacent = {{level - 1, index}};
                }
                //一番左
                else if (level == 1){
                    adjacent = {{level + 1, index},
                                {level + 1, index - 1}};
                }
                //それ以外
                else {
                    adjacent = {{level + 1, index},
                                {level + 1, index - 1},
                                {level - 1, temp_i}};
                }
            }
        }
    }
    return adjacent;
}
void Decoder::BPinit(vector<int> &u, vector<int> &x, vector<double> &y, vector<vector<int> > &xm, vector<vector<double> > &ym, vector<vector<double> > &node_val, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<int> > &B){
    vector<bool> puncFlag(Params::get_N(),false);
    vector<int> p;
    Params::get_p(p);
    //init node frozen punc etc pattern
    for (int i = 0; i < p.size(); i++) {
        puncFlag[p[i]] = true;
    }

    init_params(puncFlag, u, x, y, xm, ym, node_val, node_isChecked, B);
    init_message(puncFlag, node_val, message_list, y, ym);
}

void Decoder::init_params(vector<bool> &puncFlag, vector<int> &u, vector<int> &x, vector<double> &y, vector<vector<int> > &xm, vector<vector<double> > &ym, vector<vector<double> > &node_val, vector<vector<bool> > &node_isChecked, vector<vector<int> > &B){
    double wc, wp, llr = 0.0, val = 0.0;
    int size = 2*log2(Params::get_N())+2;
    int ysize = log2(Params::get_N())+1;
    bool databit_flag = false;
    vector<int> A;
    vector<int> Ac;
    vector<int> p;
    Params::get_A(A);
    Params::get_Ac(Ac);
    Params::get_p(p);

    //実験モード
    EXP_MODE em = Params::get_exp_mode();

    //arikan bp setting
//    vector<int> init_r_m;
//    for (int i = 0; i < Params::get_N(); i++) {
//        if(!Common::containVal(i,A)){
//            init_r_m.push_back(i);
//        }
//    }
//    init_r_m = Preseter::get_bitReversal(init_r_m);

    //punc(shorten)flg設定
    for (int i = 0; i < Params::get_N(); i++) {
        if( em == QUP || em == M_QUP || em == WANG || em == M_WANG ||
            em == VALERIO_P || em == VALERIO_S || em == M_VALERIO_P || em == M_VALERIO_S ){
            if(Common::containVal(i,p)){
                puncFlag[i] = true;
            }
        }
        else if(em == PUNC) {
            puncFlag[Ac[i]] = true;
        }
    }

    //val_node frozen shorten設定
    for (int i = 0; i < Params::get_N(); i++) {
        //channelファクターノード
        if(!puncFlag[i]){
//            if( Common::is_mid_send() ){
//                wc = Channel::calcW(ym[ysize-1][i],send_0);
//                wp = Channel::calcW(ym[ysize-1][i],1);
//            } else {
//                wc = Channel::calcW(y[i],send_0);
//                wp = Channel::calcW(y[i],1);
//            }
//            llr = 1.0 * log(wc / wp);
//            node_val[size-2][i] = llr;
//            cout << "llr" << i << " = log(" << wc << "/" << wp << ") = " << llr<< endl;
        } else {
            if( em == WANG || em == VALERIO_S){
                llr = (x[i]==0)? inf_p : inf_m;
                node_val[size-2][i] = llr;
                node_isChecked[size-2][i] = true;
            } else if( em == M_WANG || em == M_VALERIO_S){
                llr = (xm[ysize-1][i]==0)? inf_p : inf_m;
                node_val[size-2][i] = llr;
                node_isChecked[size-2][i] = true;
            }
        }
        node_isChecked[size-1][i] = true;
        //frozen_bitのllr設定
        if(!Common::containVal(i,A)){
            node_val[0][i] = (u[i] == 0) ? inf_p : inf_m;
            node_isChecked[0][i] = true;
        }
    }


    if(Params::get_is_calc_bloop()){
        //任意の位置をショートン、送信
        vector<vector<bool> > T(size, vector<bool>(Params::get_N(), false));
        Params::get_T(T);
        for (int i = 0; i < B.size(); i=i+2) {
            for (int j = 0; j < B[0].size(); j++) {
                if(B[i][j]){
                    //send_messageで送る処理はやってる
                    T[i][j] = true;
                }
            }
        }
        Params::set_T(T);
    }
}
void Decoder::init_message(vector<bool> &puncFlag, vector<vector<double> > &node_val, vector<vector<vector<message> > > &message_list , vector<double> &y, vector<vector<double> > &ym){
    double val = 0.0;
    vector<vector<int> > adjacent;
    double wc, wp, llr = 0.0;
    int ysize = log2(Params::get_N())+1;
    int size = 2*log2(Params::get_N())+2;
    message temp = {0,0,0.0};
    int adjacent_count=0, level=0, index=0;
    //グラフ構造を定義
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < Params::get_N(); j++) {
            if(i != size-1){
//                if(!puncFlag[j] || i != size-1){
                adjacent = adjacentIndex(i + 1, j + 1);
                for (int k = 0; k < adjacent.size(); k++) {
                    level = adjacent[k][0];
                    index = adjacent[k][1];
                    val = node_val[i][j];
                    temp.toLevel = level;
                    temp.toIndex = index;
                    temp.val = val;
                    message_list[i][j].push_back(temp);
                }
            }
                //右ノード受信語によるメッセージ
            else if(i == size-1) {
                adjacent = adjacentIndex(i + 1, j + 1);
                for (int k = 0; k < adjacent.size(); k++) {
                    level = adjacent[k][0];
                    index = adjacent[k][1];

                    if( Common::is_mid_send() ){
                        wc = Channel::calcW(ym[ysize-1][j],send_0);
                        wp = Channel::calcW(ym[ysize-1][j],1);
                    } else {
                        wc = Channel::calcW(y[j],send_0);
                        wp = Channel::calcW(y[j],1);
                    }
                    llr = 1.0 * log(wc / wp);
                    val = llr;
                    temp.toLevel = level;
                    temp.toIndex = index;
                    temp.val = val;
                    message_list[i][j].push_back(temp);
                }
            }
        }
    }
//    cout << " " << endl;
}

void Decoder::calc_mp(int &itr, ofstream &check_file, ofstream &val_file, int size, vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked ,vector<vector<double> > &ym){
//    calc_val_to_check(size, node_value, message_list, node_isChecked, ym);
//    calc_check_to_val(size, node_value, message_list, node_isChecked, ym);
    calc_left_to_right(itr, check_file, val_file, size, node_value, message_list, node_isChecked, ym);
    calc_right_to_left(itr, check_file, val_file, size, node_value, message_list, node_isChecked, ym);
    calc_marge(node_value, message_list, node_isChecked, ym);
}
void Decoder::calc_val_to_check(int size, vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &ym){
    vector<vector<int> > adjacent;
    bool bpflag = Params::get_decode_mode() == BP;
    for (int i = 0; i < size; i=i+2) {
        for (int j = 0; j < Params::get_N(); j++) {
            if( (!node_isChecked[i][j] || bpflag)  ) {
                send_message(i,j,message_list, node_isChecked, ym);
            }
        }
    }
}
void Decoder::calc_check_to_val(int size, vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &ym){
    vector<vector<int> > adjacent;
    bool bpflag = Params::get_decode_mode() == BP;
    for (int i = 1; i < size; i=i+2) {
        for (int j = 0; j < Params::get_N(); j++) {
            if( (!node_isChecked[i][j] || bpflag)  ) {
                send_message(i, j, message_list, node_isChecked, ym);
            }
        }
    }
}

void Decoder::calc_left_to_right(int &itr, ofstream &check_file, ofstream &val_file, int size, vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &ym){
    vector<vector<int> > adjacent;
    //left_to_right_stage
    for (int i = 0; i < size-1; i=i+2) {
        for (int j = 0; j < Params::get_N()-1; j=j+2) {
            send_message(i,   j, message_list, node_isChecked, ym);//1
            send_message(i+1, j, message_list, node_isChecked, ym);//3
            send_message(i,   j+1, message_list, node_isChecked, ym);//2
            send_message(i+1, j+1, message_list, node_isChecked, ym);//4
            send_message(i,   j+1, message_list, node_isChecked, ym);//2
            send_message(i+1, j, message_list, node_isChecked, ym);//3
            send_message(i,   j, message_list, node_isChecked, ym);//1
        }
        if(Params::get_is_outlog()) {
            calc_marge(node_value, message_list, node_isChecked, ym);
            printDecodeProgress(itr, node_value, val_file);
            printDecodeProgress(itr, node_isChecked, check_file);
            itr++;
        }
    }
}
void Decoder::calc_right_to_left(int &itr, ofstream &check_file, ofstream &val_file, int size, vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &ym){
    vector<vector<int> > adjacent;
    //right_to_left_stage
    for (int i = size-2; i >= 0; i=i-2) {
        for (int j = 0; j < Params::get_N()-1; j=j+2) {
            send_message(i,   j, message_list, node_isChecked, ym);//1
            send_message(i+1, j, message_list, node_isChecked, ym);//3
            send_message(i,   j+1, message_list, node_isChecked, ym);//2
            send_message(i+1, j+1, message_list, node_isChecked, ym);//4
            send_message(i,   j+1, message_list, node_isChecked, ym);//2
            send_message(i+1, j, message_list, node_isChecked, ym);//3
            send_message(i,   j, message_list, node_isChecked, ym);//1
        }
        if(Params::get_is_outlog()) {
            calc_marge(node_value, message_list, node_isChecked, ym);
            printDecodeProgress(itr, node_value, val_file);
            printDecodeProgress(itr, node_isChecked, check_file);
            itr++;
        }
    }
}

double Decoder::calc_message(int mode, vector<double> val) {
    double llr = 0.0;
    if (mode == 0){
        //val
        for (int i = 0; i < val.size(); i++) {
            llr += val[i];
        }
        if (isnan(llr)) {
            llr = 0.0;
//            Common::pp(val);
//            cout << "val" << endl;
        }
    } else {
        //check
        if(val.size() == 1){
            llr = val[0];
        } else {
            llr = 1.0;
            for (int i = 0; i < val.size(); i++) {
                llr *= tanh((double)val[i]/2.0);
            }
            llr = (double)2.0*atanh(llr);
//            if(val[0] != 0 && val[1] != 0){
//                llr = (double)abs(min(val[0],val[1])) * (val[0]*val[1])/(abs(val[0]*val[1]));
//            } else {
//                llr = 0.0;
//            }

        }
//        if (isnan(llr)) {
//            Common::pp(val);
//            cout << "check" << endl;
//        }
    }
    return llr;
}
//ijノードが送る値を計算
void Decoder::send_message(int i, int j, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &ym){
    int size = 2*log2(Params::get_N())+2;
    double wc,wp,llr = 0.0,val = 0.0;
    int n = Params::get_N();
    vector<vector<bool>> T(size, vector<bool>(n, false));
    Params::get_T(T);

    //中間ノードはここでおくる
    if(Common::is_mid_send()){
        if (T[i][j]) {
            wc = Channel::calcW(ym[i/2][j],send_0);
            wp = Channel::calcW(ym[i/2][j],1);
            llr = 1.0 * log(wc / wp);
        }
    }

    //次数1は無視
    if(message_list[i][j].size() == 1) return;
    //i,j番目のノードが送信するメッセージをすべて計算
    for (int k = 0; k < message_list[i][j].size(); k++) {
        vector<double> temp;

        int zero_count = 0;
        //送信元がcheckかvalか
        if(i%2 == 1){
            //check node
            for (int l = 0; l < temp.size(); l++) {
                if(temp[l] == 0.0){
                    zero_count++;
                }
            }
            //0のノードがあると計算しない
//                if(zero_count > 0) return;
        } else {
            //val node
            //frozenが出すメッセージは固定
            if(node_isChecked[i][j]) return;
        }

        //送りたい場所以外のメッセージを集めてtempに入れる
        for (int l = 0; l < message_list[i][j].size(); l++) {
            if(k != l){
                int fromLevel = message_list[i][j][l].toLevel-1;
                int fromIndex = message_list[i][j][l].toIndex-1;
                for (int m = 0; m < message_list[fromLevel][fromIndex].size(); m++) {
                    if( message_list[fromLevel][fromIndex][m].toLevel == i+1
                        && message_list[fromLevel][fromIndex][m].toIndex == j+1){
                        //今いる場所に向かってるメッセージをあつめる
                        temp.push_back(message_list[fromLevel][fromIndex][m].val);
                    }
                }
            }
        }
        int tol = message_list[i][j][k].toLevel;
        int toi = message_list[i][j][k].toIndex;

        //メッセージを送りたい場所
        val = calc_message(i%2, temp) + llr;
        message_list[i][j][k].val = val;
    }
}
//val nodeへメッセージを集める　　
void Decoder::calc_marge(vector<vector<double> > &node_value, vector<vector<vector<message> >> &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &ym){
    int size = 2*log2(Params::get_N())+2;
    int ysize = log2(Params::get_N())+1;
    int level = 0, index = 0;
    vector<vector<double> > temp(size, vector<double>(Params::get_N(),0.0));
    int n = Params::get_N();

    //check nodeが出力したメッセージをtempに集めておく
    for (int i = 1; i < size; i=i+2) {
        for (int j = 0; j < Params::get_N(); j++) {
            for (int k = 0; k < message_list[i][j].size(); k++) {
                level = message_list[i][j][k].toLevel;
                index = message_list[i][j][k].toIndex;
                temp[level-1][index-1] += (double)message_list[i][j][k].val;
            }
        }
    }

    //中間ノードからのメッセージを引き上げる
    if( Common::is_mid_send() ){
        vector<vector<bool>> T(size, vector<bool>(n, false));
        Params::get_T(T);
        for (int i = 0; i < ysize; i++) {
            for (int j = 0; j < Params::get_N(); j++) {
                if (T[i][j]) {
                    double wc = Channel::calcW(ym[i][j],send_0);
                    double wp = Channel::calcW(ym[i][j],1);
                    double llr = 1.0 * log(wc / wp);
                    //val nodeへ送るから, 2*iでいい
                    temp[2*i][j] += llr;
                }
            }
        }
    }

    //val node更新
    for (int i = 0; i < size; i=i+2) {
        for (int j = 0; j < Params::get_N(); j++) {
            //ショートンビットは変更しない
            if(!node_isChecked[i][j]){
//            if(temp[i][j] != 0.0 && !node_isChecked[i][j]){
                node_value[i][j] = temp[i][j];
            }
        }
    }
}
void Decoder::confirmIsCheck(vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked) {
    int size = 2*log2(Params::get_N())+2;
    int level = 0, index = 0, tempCheck = 0, tempVal = 0;
    double val = 0.0;
    bool zero_flag = false;
    vector<vector<int> > adjacent;

    //チェックが通るかみる
    //channel factor nodeでやる必要はない
    for (int i = 1; i < size-1; i=i+2) {
        for (int j = 0; j < Params::get_N(); j++) {
            zero_flag = false;
            tempCheck = 0;
            adjacent = adjacentIndex(i + 1, j + 1);
            for (int k = 0; k < adjacent.size(); k++) {
                val = take_val(adjacent[k], node_value);
                if(val == 0.0) {
                    zero_flag = true;
                } else {
                    tempVal = (val >= 0.0) ? 0:1;
                    tempCheck += tempVal;
                }
            }

            if(tempCheck%2 == 0 && !zero_flag) {
                //足して0ならtrue
                node_isChecked[i][j] = true;
            } else if(i != size-1) {
                //チャンネルファクター以外で上の条件満たしていない場合
                node_isChecked[i][j] = false;
            }
//            else if(zero_flag){
//                //ゼロがあると終わらない
//                node_isChecked[i][j] = false;
//            }
        }
    }
}
bool Decoder::isTerminate(int &no_checked, vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked) {
    bool flag = true;
    int size = 2*log2(Params::get_N())+2;
    vector<vector<int> > adjacent;

    no_checked = 0;
    //checkノードで回していく
    for (int i = 1; i < size; i=i+2) {
        for (int j = 0; j < Params::get_N(); j++) {
            if(!node_isChecked[i][j]){
                no_checked++;
            }
        }
    }
    flag = (no_checked == 0);
    return flag;
}
bool Decoder::isScTerminate(int &no_checked, vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked) {
    for (int i = 0; i < Params::get_N(); i++) {
        if(!node_isChecked[0][i]){
            return false;
        }
    }
    return true;
}
bool Decoder::isScTerminateG(int &no_checked, vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked) {
    Encoder e;
    int size = 2*log2(Params::get_N())+2;
    vector<int>temp_u(Params::get_N(),0);
    vector<int>temp_x(Params::get_N(),0);
    vector<int>calc_x(Params::get_N(),0);
    for (int i = 0; i < Params::get_N(); i++) {
        temp_u[i] = (node_value[0][i] >= 0.0)?0:1;
        temp_x[i] = (node_value[size-2][i] >= 0.0)?0:1;
    }
    calc_x = e.encode(Params::get_N(), temp_u);
    for (int i = 0; i < Params::get_N(); i++) {
        if(calc_x[i] != temp_x[i]){
            return false;
        }
    }
    return true;
}
void Decoder::reset_message(vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked, vector<vector<vector<message> >> &message_list){
    int level, index;
    int size = 2*log2(Params::get_N())+2;
    vector<vector<int> > adjacent;
    for (int j = 1; j < size; j=j+2) {
        for (int k = 0; k < Params::get_N(); k++) {
            //チャンネルファクター以外のチェックノードをfalse
            if(j != size-1){
                node_isChecked[j][k] = false;
            }
        }
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < Params::get_N(); j++) {
            if(!node_isChecked[i][j]){
                node_value[i][j] = 0.0;
                adjacent = adjacentIndex(i+1, j+1);
                for (int k = 0; k < adjacent.size(); k++) {
                    level = adjacent[k][0];
                    index = adjacent[k][1];
                    message_list[i][j][k].val = 0.0;
                }
            }
        }
    }
}
void Decoder::sc_decision(vector<double> &test_u, vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked, vector<vector<vector<message> >> &message_list){
    vector<int> A;
    vector<int> Ac;
    vector<int> p;
    Params::get_A(A);
    Params::get_Ac(Ac);
    Params::get_p(p);
    for (int i = 0; i < Params::get_N(); i++) {
        if(!Common::containVal(i,A)){
            test_u[i] = (node_value[0][i] >= 0) ? inf_p : inf_m;
        }
        if(Common::containVal(i,A) && !node_isChecked[0][i]){
            if(node_value[0][i] != 0.0) {
                //無限大で確定させる, 対応するメッセージを送信する
                test_u[i] = node_value[0][i];
                node_value[0][i] = (node_value[0][i] >= 0) ? inf_p : inf_m;
                node_isChecked[0][i] = true;

                vector<vector<int> > adjacent;
                double val = 0.0;
                int size = 2 * log2(Params::get_N()) + 2;
                int level = 0, index = 0;
                adjacent = adjacentIndex(1, i + 1);
                for (int k = 0; k < adjacent.size(); k++) {
                    level = adjacent[k][0];
                    index = adjacent[k][1];
                    val = node_value[0][i];
                    message_list[0][i][k].val = val;
                }
                //2n番目にあるノードを2n-1番目のノードを元に確定させる
                if(i%2 == 0){
                    vector<vector<double> > ym;
                    send_message(1, i, message_list, node_isChecked, ym);
                    calc_marge(node_value, message_list, node_isChecked, ym);
                    test_u[i+1] = node_value[0][i+1];
                    node_value[0][i+1] = (node_value[0][i+1] >= 0) ? inf_p : inf_m;
                    node_isChecked[0][i+1] = true;

                    adjacent = adjacentIndex(1, i + 2);
                    for (int k = 0; k < adjacent.size(); k++) {
                        level = adjacent[k][0];
                        index = adjacent[k][1];
                        val = node_value[0][i+1];
                        message_list[0][i+1][k].val = val;
                    }
                }
                reset_message(node_value, node_isChecked, message_list);
            }
            return;
        }
    }
}
vector<int> Decoder::calcBP(vector<double> &test_u, int loop_num, vector<int> &param, vector<int> &u, vector<int> &x, vector<double> &y, vector<vector<int> > &xm, vector<vector<double> > &ym, vector<vector<int> > &node_error_count, ofstream &val_error_file, vector<vector<int> > &B) {
    vector<double> tmp_u(Params::get_N(),0.0);
    vector<int> u_n_est(Params::get_N());
    vector<vector<int> > adjacent;

    //node初期化、奇数が変数ノード（1,3,5...）、偶数がチェックノード（2,4,6...）
    int size = 2*log2(Params::get_N())+2;
    int ysize = log2(Params::get_N())+1;

//    for (int l = 0; l < Params::get_N(); l++) {
//        ym[ysize-1][l] = y[l];
//    }

    vector<vector<vector<message> >> message_list(size, vector<vector<message> >(Params::get_N(), vector<message>()));
    vector<vector<double> > node_value(size, vector<double>(Params::get_N(),0.0));
    vector<vector<bool> > node_isChecked(size, vector<bool>(Params::get_N(),false));

    //BP
    int count = 0;
    int itr = 0;
    int no_checked = 0;
    ofstream val_file;
    ofstream check_file;

    init_outLog(val_file, check_file, val_error_file);
    BPinit(u, x, y, xm, ym, node_value, message_list, node_isChecked, B);

    for (int i = 0; i < Params::get_rp(); i++) {
        if(Params::get_decode_mode() == SC){
            if(isScTerminate(no_checked, node_value, node_isChecked)) break;
        } else {
            if(isTerminate(no_checked, node_value, node_isChecked)) break;
//            if(itr > 0){
//                if(isScTerminateG(no_checked, node_value, node_isChecked)) break;
//            }
        }
        if(Params::get_is_outlog()) {
            printDecodeProgress(itr, node_value, val_file);
            printDecodeProgress(itr, node_isChecked, check_file);
        }
        itr++;
        for (int j = 1; j < size; j=j+2) {
            for (int k = 0; k < Params::get_N(); k++) {
                //チャンネルファクター以外のチェックノードをfalse
                if(j != size-1){
                    node_isChecked[j][k] = false;
                }
            }
        }
        calc_mp(itr, check_file, val_file, size, node_value, message_list, node_isChecked, ym);
        confirmIsCheck(node_value, node_isChecked);

        if(Params::get_decode_mode() == SC){
            sc_decision(test_u, node_value, node_isChecked, message_list);
        }
        count++;
    }
//    else {
//        calcSConBP(itr, count, val_file, check_file, u, y, u_n_est, tmp_u, node_value, node_isChecked);
//    }

    int error_count = 0;
    Analysor::errorCount(u, u_n_est, &error_count);
    param[0] = itr;
    param[1] = no_checked;

    for (int i = 0; i < Params::get_N(); i++) {
        if(Params::get_is_outlog()) {
//            cout << "BP: " << i + 1 << " " << (double)test_u[i] << endl;
        }
        test_u[i] = node_value[0][i];
        u_n_est[i] = (node_value[0][i] >= 0.0)?0:1;
    }

    if(Common::is_mid_send() && Params::get_is_calc_bloop()){
        for (int i = 0; i < size; i=i+2) {
            for (int j = 0; j < Params::get_N(); j++) {
                int temp = (node_value[i][j] >= 0.0)?0:1;
//                int temp_x = xm[i/2][j];
                if(xm[i/2][j] != temp){
                    node_error_count[i][j]++;
                }
            }
        }
        printDecodeProgress(loop_num-1, node_error_count, val_error_file);
    }

    outLog(itr, no_checked, u, x, xm, y, ym ,u_n_est, val_file, check_file , node_value, node_isChecked);

    return u_n_est;
}

void Decoder::init_outLog(ofstream &val_file, ofstream &check_file, ofstream &val_error_file){
    if(Params::get_is_outlog()) {
        string val_fn = "/Users/ryotaro/Dropbox/labo/graph_js/val.json";
        string check_fn = "/Users/ryotaro/Dropbox/labo/graph_js/check.json";
        val_file.open(val_fn, ios::out);
        check_file.open(check_fn, ios::out);
        val_file << "{" << endl;
        check_file << "{" << endl;
    }
}
void Decoder::outLog(int itr, int no_checked, vector<int> &u, vector<int> &x, vector<vector<int>> &xm, vector<double> &y, vector<vector<double>> &ym, vector<int> &u_n_est, ofstream &val_file, ofstream &check_file ,vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked){
    int ysize = log2(Params::get_N())+1;
    if(Params::get_is_outlog()) {
        printDecodeProgress(itr, node_value, val_file);
        printDecodeProgress(itr, node_isChecked, check_file);
        val_file << "\t}" << endl;
        val_file << "}" << endl;
        check_file << "\t}" << endl;
        check_file << "}" << endl;

        string channel = "";
        if (Params::get_s() == BEC) {
            channel = "BEC";
        } else if (Params::get_s() == BSC) {
            channel = "BSC";
        } else {
            channel = "AWGN";
        }

        string correct_fn = "/Users/ryotaro/Dropbox/labo/graph_js/correct.json";
        ofstream correct_file;
        correct_file.open(correct_fn, ios::out);
        correct_file << "{" << endl;

        string correcty_fn = "/Users/ryotaro/Dropbox/labo/graph_js/correcty.json";
        ofstream correcty_file;
        correcty_file.open(correcty_fn, ios::out);
        correcty_file << "{" << endl;

        string y_fn = "/Users/ryotaro/Dropbox/labo/graph_js/y.json";
        ofstream y_file;
        y_file.open(y_fn, ios::out);
        y_file << "{" << endl;

        for (int i = 0; i < Params::get_N(); i++) {
            if (i == Params::get_N() - 1) {
                correct_file << "\t\"" << i << "\" : \"" << u[i] << "\"" << endl;
            } else {
                correct_file << "\t\"" << i << "\" : \"" << u[i] << "\"," << endl;
            }
            if (Common::is_mid_send()) {
                if (i == Params::get_N() - 1) {
                    correcty_file << "\t\"" << i << "\" : \"" << xm[ysize-1][i] << "\"" << endl;
                    y_file << "\t\"" << i << "\" : \"" << ym[ysize-1][i] << "\"" << endl;
                } else {
                    correcty_file << "\t\"" << i << "\" : \"" << xm[ysize-1][i] << "\"," << endl;
                    y_file << "\t\"" << i << "\" : \"" << ym[ysize-1][i] << "\"," << endl;
                }
            } else {
                if (i == Params::get_N() - 1) {
                    correcty_file << "\t\"" << i << "\" : \"" << x[i] << "\"" << endl;
                    y_file << "\t\"" << i << "\" : \"" << y[i] << "\"" << endl;
                } else {
                    correcty_file << "\t\"" << i << "\" : \"" << x[i] << "\"," << endl;
                    y_file << "\t\"" << i << "\" : \"" << y[i] << "\"," << endl;
                }
            }
        }
        correcty_file << "}" << endl;
        correct_file << "}" << endl;
        y_file << "}" << endl;

        string params_fn = "/Users/ryotaro/Dropbox/labo/graph_js/params.json";
        ofstream params_file;
        params_file.open(params_fn, ios::out);
        params_file << "{" << endl;
        params_file << "\t\"N\" : \"" << Params::get_N() << "\"," << endl;
        params_file << "\t\"Channel\" : \"" << channel << "\"," << endl;
        params_file << "\t\"e\" : \"" << Params::get_e() << "\"," << endl;
        params_file << "\t\"ITR\" : \"" << itr << "\"," << endl;
        params_file << "\t\"no_checked\" : \"" << no_checked << "\"," << endl;

        int error_count = 0;
        Analysor::errorCount(u, u_n_est, &error_count);
        params_file << "\t\"error_count\" : \"" << error_count << "\"" << endl;

        params_file << "}" << endl;
    }
}
inline vector<int>Decoder::makeBPTreeIndex(int n){
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
inline vector<int>Decoder::makeTreeIndex(int n){
    vector<int> ret(n);
    if (n == 1) {
        ret[0] = 1;
    } else {
        vector<int> temp = makeTreeIndex(n/2);
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
vector<int> Decoder::decode(vector<double> &test_u, vector<double> &y, vector<int> &u){
    vector<double> h_i(Params::get_N());
    vector<int> u_n_est(Params::get_N());
    int size = log2(Params::get_N());

    vector<vector<bool> > isCache (size, vector<bool>(Params::get_N(),false));
    vector<vector<double> > cache (size, vector<double>(Params::get_N(),0.0));
    vector<int> A;
    Params::get_A(A);

    double llr = 0.0;
    int cache_i = 0;

    //u_n_est計算
    for (int i = 0; i < Params::get_N(); i++) {
        // Aに含まれないindexなら既知
        if (Common::containNumInArray(i, Params::get_N()-Params::get_K(), A) == false) {
            u_n_est[i] = u[i];
            test_u[i] = (u[i] == 1)?(inf_m):(inf_p);
            if(Params::get_is_outlog()) {
                cout << "SC: " << i+1 << " " << ((u[i] == 1)?(inf_m):(inf_p)) << endl;
            }
        } else {
            this->startTimer();

            cache_i = makeTreeIndex(Params::get_N())[i] - 1;
            llr = calcL_i(i+1, Params::get_N(), cache_i, 0, y, u_n_est, isCache, cache);

            this->stopTimer();
            this->outTime();
//            cout << i+1 << " " << llr <<endl;

            if (llr >= 0.0) {
                h_i[i] = 0;
            } else {
                h_i[i] = 1;
            }
            u_n_est[i] = h_i[i];
            test_u[i] = llr;
            if(Params::get_is_outlog()) {
                cout << "SC: " << i + 1 << " " << llr << endl;
            }
//            printf("%05f\n", llr);
        }
    }
    return u_n_est;
}
double Decoder::calcL_i(int i, int n, int cache_i, int level, vector<double> &y, vector<int> &u, vector<vector<bool> > &isCache, vector<vector<double> > &cache) {
    double llr = 0.0;
    this->addCount();
    if ( n == 1 ) {
        double wc = Channel::calcW(y[0],send_0);
        double wp = Channel::calcW(y[0],1);
        llr = 1.0 * log(wc / wp);

    } else {
        vector<double> tempY1(n/2);
        vector<double> tempY2(n/2);
        for (int j = 0; j < n ; j++) {
            if(j < n/2){
                tempY1[j] = y[j];
            } else {
                tempY2[j-n/2] = y[j];
            }
        }
        int size_u_eo = (i % 2) == 0 ? i-2 : i-1;
        int size_aug_u = size_u_eo/2;

        vector<int> tempU(size_aug_u);
        vector<int> tempU_bin(size_aug_u);
        vector<int> u_e(size_aug_u);
        vector<int> u_o(size_aug_u);

        u_e = Common::index_e(size_u_eo, u);
        u_o = Common::index_o(size_u_eo, u);

        for(int k=0; k<size_aug_u; k++){
            tempU[k] = u_e[k] + u_o[k];
        }

        tempU_bin = Common::retBinary(size_aug_u, tempU);

        double temp1 = 0.0;
        double temp2 = 0.0;

        int temp_i = (i % 2 == 1) ? (i+1)/2 : i/2;

        int ci = cache_i;
        if (((cache_i >> (int) (log2(n) - 1)) % 2) != 1) {
            //横の辺を作る
            if (isCache[level][cache_i]) {
                temp1 = cache[level][cache_i];
            } else {
                temp1 = calcL_i(temp_i, n / 2, cache_i, level + 1, tempY1, tempU_bin, isCache, cache);
                isCache[level][cache_i] = true;
                cache[level][cache_i] = temp1;
            }
            //斜め下の辺を作る
            if (isCache[level][cache_i + (n / 2)]) {
                temp2 = cache[level][cache_i + (n / 2)];
            } else {
                temp2 = calcL_i(temp_i, n / 2, cache_i + (n / 2), level + 1, tempY2, u_e, isCache, cache);
                isCache[level][cache_i + (n / 2)] = true;
                cache[level][cache_i + (n / 2)] = temp2;
            }
        } else {
            //斜め上の辺を作る
            if (isCache[level][cache_i - (n / 2)]) {
                temp1 = cache[level][cache_i - (n / 2)];
            } else {
                temp1 = calcL_i(temp_i, n / 2, cache_i - (n / 2), level + 1, tempY1, tempU_bin, isCache, cache);
                isCache[level][cache_i - (n / 2)] = true;
                cache[level][cache_i - (n / 2)] = temp1;
            }
            //横の辺を作る
            if (isCache[level][cache_i]) {
                temp2 = cache[level][cache_i];
            } else {
                temp2 = calcL_i(temp_i, n / 2, cache_i, level + 1, tempY2, u_e, isCache, cache);
                isCache[level][cache_i] = true;
                cache[level][cache_i] = temp2;
            }
        }
        if ( i % 2 == 1) {
            llr = 2.0 * atanh(tanh(temp1/2.0) * tanh(temp2/2.0) );
        } else {
            llr = (1-2*u[size_u_eo]) * temp1 + temp2;
        }
    }

    if (isinf(llr) && llr > 0) {
//        llr = 10.0;
    } else if (isinf(llr) && llr < 0) {
//        llr = -10.0;
    }
    return llr;

}
