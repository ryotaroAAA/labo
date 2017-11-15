#include "../lib/Decoder.h"
#include "../lib/Analysor.h"
Decoder::Decoder(){

}

Decoder::~Decoder(){

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
                temp = inf_p;
                w_file << "\t\t\t\"" << j << "\"" << ":" << temp;
            } else if(isinf(temp) && temp < 0){
                temp = inf_m;
                w_file << "\t\t\t\"" << j << "\"" << ":" << temp;
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
                temp = inf_p;
                w_file << "\t\t\t\"" << j << "\"" << ":" << temp;
            } else if(isinf(temp) && temp < 0){
                temp = inf_m;
                w_file << "\t\t\t\"" << j << "\"" << ":" << temp;
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

vector<vector<int> > Decoder::adjacentIndex(int level, int index){
    vector<vector<int> > adjacent;
    vector<int> temp;
    int n = (level != 1) ? Params::get_N()/pow(2,((level-2)/2)) : 1;
    temp = makeBPTreeIndex(n);
    int size = 2*log2(Params::get_N())+2;
    if( level ==  size ){
        adjacent = {{level - 1, index}};
    } else if( level ==  size-1 ) {
        adjacent = {{level - 1, index}, {level + 1, index}};
    } else {
        if (level % 2 == 0) {
            //check_node
            int temp_i = (index == n) ? n - 1 : (index - 1) % n;
            int div = (index - 1) / n;
            if (index % 2 == 1) {
                //index_o
                if (level < 2 * log2(Params::get_N()))
                    adjacent = {{level - 1, index},
                                {level - 1, index + 1},
                                {level + 1, temp[temp_i] + div * n}};
                else
                    adjacent = {{level - 1, index},
                                {level - 1, index + 1},
                                {level + 1, index}};
            } else {
                // /index_e
                if (level < 2 * log2(Params::get_N()))
                    adjacent = {{level - 1, index},
                                {level + 1, temp[temp_i] + div * n}};
                else
                    adjacent = {{level - 1, index},
                                {level + 1, index}};
            }
        } else {
            //val_node
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

            if (index % 2 == 1) {
                //index_o
                //一番右
                if (level == 2 * log2(Params::get_N()) + 1) adjacent = {{level - 1, index}};
                    //一番左
                else if (level == 1) adjacent = {{level + 1, index}};
                    //それ以外
                else {
                    adjacent = {{level + 1, index},
                                {level - 1, temp_i}};
                }
            } else {
                //index_e
                //一番右
                if (level == 2 * log2(Params::get_N()) + 1) adjacent = {{level - 1, index}};
                    //一番左
                else if (level == 1)
                    adjacent = {{level + 1, index},
                                {level + 1, index - 1}};
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

double Decoder::take_val(vector<int> &locate, vector<vector<double> > &node_val){
    int level = locate[0];
    int index = locate[1];
    return node_val[level-1][index-1];
}

void Decoder::send_message(int i, int j, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &y, vector<vector<bool> > &ym_isReceived){
    int size = 2*log2(Params::get_N())+2;
    double wc,wp,llr = 0.0,val = 0.0;
    bool is_send = false;
    if(Common::is_mid_send()){
        if (ym_isReceived[i][j]) {
            wc = Channel::calcW(y[i][j],0);
            wp = Channel::calcW(y[i][j],1);
            llr = 1.0 * log(wc / wp);
        }
    }
    for (int k = 0; k < message_list[i][j].size(); ++k) {
        vector<double> temp;
        if(message_list[i][j].size() == 1){
            //次数1のval nodeはメッセージをそのまま返す, channel factor nodeはほっとく
            if (!node_isChecked[i][j] && i < size-1 ) {
                //frozen bitはメッセージ固定なのでそれ以外を計算
                int fromLevel = message_list[i][j][0].toLevel-1;
                int fromIndex = message_list[i][j][0].toIndex-1;
                for (int m = 0; m < message_list[fromLevel][fromIndex].size(); m++) {
                    if( message_list[fromLevel][fromIndex][m].toLevel == i+1 && message_list[fromLevel][fromIndex][m].toIndex == j+1){
                        //メッセージを送りたい場所以外のメッセージ
                        val = message_list[fromLevel][fromIndex][m].val + llr;
                        if(val != 0.0 || isinf(val)) {
                            message_list[i][j][k].val = val;
                            is_send = true;
                        }
                        break;
                    }
                }
            }
        } else {
            for (int l = 0; l < message_list[i][j].size(); l++) {
                if(k != l){
                    int fromLevel = message_list[i][j][l].toLevel-1;
                    int fromIndex = message_list[i][j][l].toIndex-1;
                    for (int m = 0; m < message_list[fromLevel][fromIndex].size(); m++) {
                        if( message_list[fromLevel][fromIndex][m].toLevel == i+1 && message_list[fromLevel][fromIndex][m].toIndex == j+1){
                            //メッセージを送りたい場所以外のメッセージ
                            temp.push_back(message_list[fromLevel][fromIndex][m].val);
                        }
                    }
                }
            }
            int tol = message_list[i][j][k].toLevel;
            int toi = message_list[i][j][k].toIndex;

            int zero_count = 0;
            bool isCalcable = true;
            if(i%2 == 1){
                //check node
                for (int l = 0; l < temp.size(); ++l) {
                    if(temp[l] == 0.0){
                        zero_count++;
                    }
                }
                if(zero_count > 0) isCalcable = false;
            } else {
                //val node
                if(node_isChecked[i][j]) isCalcable = false;
            }

            //メッセージを送りたい場所
            if(isCalcable){
                val = calc_message(i%2, temp) + llr;
                if(val != 0.0){
                    message_list[i][j][k].val = val;
                    is_send = true;
                }
            }
        }
    }
    if (Params::get_decode_mode() == SC && is_send) node_isChecked[i][j] = true;
}

double Decoder::calc_message(int mode, vector<double> val) {
    double llr = 0.0;
    if (mode == 0){
        //val
        llr = 0.0;
        for (int i = 0; i < val.size(); i++) {
            llr += val[i];
        }
    } else {
        //check
        if(val.size() == 1){
            llr = val[0];
        } else {
            llr = 1.0;
            for (int i = 0; i < val.size(); i++) {
                llr *= tanh(val[i]/2.0);
            }
            llr = 2*atanh(llr);
        }
    }
    return llr;
}

void Decoder::confirmIsCheck(vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked) {
    int size = 2*log2(Params::get_N())+2;
    int level = 0, index = 0, tempCheck = 0, tempVal = 0;
    double val = 0.0;
    bool zero_flag = false;
    vector<vector<int> > adjacent;

    //val nodeからのメッセージをあつめる
    //check nodeを回して隣接のval nodeをみる
    //channel check nodeは見ない
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
            if(tempCheck%2 == 0 && !zero_flag){
                node_isChecked[i][j] = true;
            } else if(i != size-1) {
//                node_isChecked[i][j] = false; //なんのため？
            }
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

bool Decoder::isTerminate(int &no_checked, vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked) {
    bool flag = true;
    int size = 2*log2(Params::get_N())+2;
    vector<vector<int> > adjacent;

    no_checked = 0;
    for (int i = 1; i < size; i=i+2) {
        for (int j = 0; j < Params::get_N(); j++) {
            if(node_isChecked[i][j]){
                flag &= true;
            } else {
                no_checked++;
                flag &= false;
            }
        }
    }
    return flag;
}

void Decoder::init_message(vector<bool> &puncFlag, vector<vector<double> > &node_val, vector<vector<vector<message> > > &message_list){
    double val = 0.0;
    vector<vector<int> > adjacent;
    int size = 2*log2(Params::get_N())+2;
    message temp = {0,0,0.0};
    int adjacent_count=0, level=0, index=0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < Params::get_N(); j++) {
            //隣接ノードの位置
            if(!puncFlag[j] || i != size-1){
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
        }
    }
}

void Decoder::SCinit(int n, vector<double> &y, vector<int> &u, vector<int> &u_est, vector<int> &A, vector<vector<double> > &node_val, vector<vector<vector<message> > > &message_list, vector<vector<vector<message> > > &save_list, vector<vector<bool> > &node_isChecked, vector<vector<bool> > &save_isChecked){
    double wc, wp, llr = 0.0, val = 0.0;
    int size = 2*log2(Params::get_N())+2;
    bool databit_flag = false;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < Params::get_N(); j++) {
            if(!save_isChecked[i][j]){
                node_val[i][j] = 0.0;
                node_isChecked[i][j] = false;
            }
        }
    }
    for (int i = 0; i < Params::get_N(); i++) {
        wc = Channel::calcW(y[i],0);
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

void Decoder::init_params(vector<bool> &puncFlag, vector<int> &p, vector<int> &u, vector<int> &x, vector<double> &y, vector<vector<int> > &xm, vector<vector<double> > &ym, vector<int> &A, vector<int> &Ac, vector<vector<double> > &node_val, vector<vector<bool> > &node_isChecked, vector<vector<bool> > &ym_isReceived){
    double wc, wp, llr = 0.0, val = 0.0;
    int size = 2*log2(Params::get_N())+2;
    int ysize = log2(Params::get_N())+1;
    bool databit_flag = false;

    //実験モード
    EXP_MODE em = Params::get_exp_mode();

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
        if(!puncFlag[i]){
            //puncなし
            if( Common::is_mid_send() ){
                wc = Channel::calcW(ym[ysize-1][i],0);
                wp = Channel::calcW(ym[ysize-1][i],1);
            } else {
                wc = Channel::calcW(y[i],0);
                wp = Channel::calcW(y[i],1);
            }
            llr = 1.0 * log(wc / wp);
            node_val[size-2][i] = llr;

        } else {
            //punc(shorten)あり
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
        //channelファクターノード
        node_isChecked[size-1][i] = true;

        //frozen_bitのllr設定
        databit_flag = false;
        for (int j = 0; j < Params::get_K(); j++) {
            if(i == A[j]){
                databit_flag = true;
            }
        }
        if(databit_flag == false){
            node_val[0][i] = (u[i] == 0) ? inf_p : inf_m;
            node_isChecked[0][i] = true;
        }
    }

    //中間ノード設定, yからメッセージを受信するようになる
    if( Common::is_mid_send() ){
        for (int i = 0; i < Params::get_M(); i++) {
            ym_isReceived[0][A[i]] = true;
        }
    }

}

void Decoder::BPinit(vector<int> &p, vector<int> &u, vector<int> &x, vector<double> &y, vector<vector<int> > &xm, vector<vector<double> > &ym, vector<int> &A, vector<int> &Ac, vector<vector<double> > &node_val, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<bool> > &ym_isReceived){
    vector<bool> puncFlag(Params::get_N(),false);
    //init node frozen punc etc pattern
    init_params(puncFlag, p, u, x, y, xm, ym, A, Ac, node_val, node_isChecked, ym_isReceived);
    init_message(puncFlag, node_val, message_list);
}

void Decoder::calc_mp(int size, vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked ,vector<vector<double> > &ym, vector<vector<bool> > &ym_isReceived){
    calc_val_to_check(size, node_value, message_list, node_isChecked, ym, ym_isReceived);
    calc_check_to_val(size, node_value, message_list, node_isChecked, ym, ym_isReceived);
    calc_marge(node_value, message_list, node_isChecked, ym, ym_isReceived);
}

void Decoder::calc_val_to_check(int size, vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &ym, vector<vector<bool> > &ym_isReceived){
    vector<vector<int> > adjacent;
    bool bpflag = Params::get_decode_mode() == BP;
    for (int i = 0; i < size; i=i+2) {
        for (int j = 0; j < Params::get_N(); j++) {
            if( (!node_isChecked[i][j] || bpflag)  ) {
                send_message(i,j,message_list, node_isChecked, ym, ym_isReceived);
            }
        }
    }
}

void Decoder::calc_check_to_val(int size, vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &y, vector<vector<bool> > &ym_isReceived){
    vector<vector<int> > adjacent;
    bool bpflag = Params::get_decode_mode() == BP;
    for (int i = 1; i < size; i=i+2) {
        for (int j = 0; j < Params::get_N(); j++) {
            if( (!node_isChecked[i][j] || bpflag)  ) {
                send_message(i,j,message_list, node_isChecked, y, ym_isReceived);
            }
        }
    }
}

void Decoder::calc_marge(vector<vector<double> > &node_value, vector<vector<vector<message> >> &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &y, vector<vector<bool> > &ym_isReceived){
    int size = 2*log2(Params::get_N())+2;
    int ysize = log2(Params::get_N())+1;
    int level = 0, index = 0;
    vector<vector<double> > temp(size, vector<double>(Params::get_N(),0.0));
    //check nodeからのメッセージをあつめる
    //check nodeが出力したメッセージをtempに集めておく
    for (int i = 1; i < size; i=i+2) {
        for (int j = 0; j < Params::get_N(); j++) {
            for (int k = 0; k < message_list[i][j].size(); k++) {
                level = message_list[i][j][k].toLevel-1;
                index = message_list[i][j][k].toIndex-1;
                temp[level][index] += message_list[i][j][k].val;
            }
        }
    }

    if( Common::is_mid_send() ){
        for (int i = 0; i < ysize; i++) {
            for (int j = 0; j < Params::get_N(); j++) {
                if (ym_isReceived[i][j]) {
                    double wc = Channel::calcW(y[i][j],0);
                    double wp = Channel::calcW(y[i][j],1);
                    double llr = 1.0 * log(wc / wp);
                    //val nodeへ送るから, 2*iでいい
                    temp[2*i][j] += llr;
                }
            }
        }
    }

    //更新
    for (int i = 0; i < size; i=i+2) {
        for (int j = 0; j < Params::get_N(); j++) {
            if(temp[i][j] != 0.0 && !node_isChecked[i][j]){
                node_value[i][j] = temp[i][j];
            }
        }
    }
}

vector<int> Decoder::calcBP(vector<int> p, vector<int> &param, vector<int> &u, vector<int> &x, vector<double> &y, vector<vector<int> > &xm, vector<vector<double> > &ym, vector<int> &A, vector<int> &Ac) {
    vector<double> tmp_u(Params::get_N(),0.0);
    vector<int> u_n_est(Params::get_N());
    vector<vector<int> > adjacent;

    //node初期化、奇数が変数ノード（1,3,5...）、偶数がチェックノード（2,4,6...）
    int size = 2*log2(Params::get_N())+2;

    vector<vector<vector<message> >> message_list(size, vector<vector<message> >(Params::get_N(), vector<message>()));
    vector<vector<double> > node_value(size, vector<double>(Params::get_N(),0.0));
    vector<vector<bool> > node_isChecked(size, vector<bool>(Params::get_N(),false));
    vector<vector<bool> > ym_isReceived(size, vector<bool>(Params::get_N(),false));

    //BP
    int count = 0;
    int itr = 0;
    int no_checked = 0;

    ofstream val_file;
    ofstream check_file;
    init_outLog(val_file, check_file);
    if (Params::get_decode_mode() == BP) {
        vector<vector<vector<message> > > message_list(size, vector<vector<message> >(Params::get_N(), vector<message>()));
        BPinit(p, u, x, y, xm, ym, A, Ac, node_value, message_list, node_isChecked, ym_isReceived);

        for (int i = 0; i < Params::get_rp(); i++) {
            if(isTerminate(no_checked, node_value, node_isChecked)) break;
            if(Params::get_is_outlog()) {
                printDecodeProgress(itr, node_value, val_file);
                printDecodeProgress(itr, node_isChecked, check_file);
            }
            itr++;
            calc_mp(size, node_value, message_list, node_isChecked, ym, ym_isReceived);

            confirmIsCheck(node_value, node_isChecked);
            count++;
        }

        for (int i = 0; i < Params::get_N(); i++) {
            if(Params::get_is_outlog()) {
                cout << "BP: " << i+1 << " " << node_value[0][i] << endl;
            }
            u_n_est[i] = (node_value[0][i]>=0.0) ? 0 : 1;
        }
    } else {
        calcSConBP(itr, count, val_file, check_file, u, y, u_n_est, tmp_u, A, node_value, node_isChecked, ym_isReceived);
    }

    int error_count = 0;
    Analysor::errorCount(u, u_n_est, &error_count);
    param[0] = itr;
    param[1] = no_checked;

    for (int i = 0; i < Params::get_N(); i++) {
        if(Params::get_is_outlog()) {
            cout << "BP: " << i + 1 << " " << node_value[0][i] << endl;
        }
        u_n_est[i] = (node_value[0][i] > 0) ? 0 : 1;
    }

    outLog(itr, no_checked, u, u_n_est, val_file, check_file , node_value, node_isChecked);

    return u_n_est;
}

void Decoder::init_outLog(ofstream &val_file, ofstream &check_file){
    if(Params::get_is_outlog()) {
        string val_fn = "/Users/ryotaro/Dropbox/labo/graph_js/val.json";
        string check_fn = "/Users/ryotaro/Dropbox/labo/graph_js/check.json";
        val_file.open(val_fn, ios::out);
        check_file.open(check_fn, ios::out);
        val_file << "{" << endl;
        check_file << "{" << endl;
    }
}

void Decoder::outLog(int itr, int no_checked, vector<int> &u, vector<int> &u_n_est, ofstream &val_file, ofstream &check_file ,vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked){
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
        for (int i = 0; i < Params::get_N(); i++) {
            if (i == Params::get_N() - 1) {
                correct_file << "\t\"" << i << "\" : \"" << u[i] << "\"" << endl;
            } else {
                correct_file << "\t\"" << i << "\" : \"" << u[i] << "\"," << endl;
            }
        }
        correct_file << "}" << endl;

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

vector<int> Decoder::decode(vector<double> &y, vector<int> &u, vector<int> &A){
    vector<double> h_i(Params::get_N());
    vector<int> u_n_est(Params::get_N());
    int size = log2(Params::get_N());

    vector<vector<bool> > isCache (size, vector<bool>(Params::get_N(),false));
    vector<vector<double> > cache (size, vector<double>(Params::get_N(),0.0));

    double llr = 0.0;
    int cache_i = 0;

    //u_n_est計算
    for (int i = 0; i < Params::get_N(); i++) {
        // Aに含まれないindexなら既知
        if (Common::containNumInArray(i, Params::get_N()-Params::get_K(), A) == false) {
            u_n_est[i] = u[i];
            cout << "SC: " << i+1 << " " << ((u[i] == 1)?(inf_m):(inf_p)) << endl;
        } else {
            this->startTimer();

            cache_i = makeTreeIndex(Params::get_N())[i] - 1;
            llr = calcL_i(i+1, Params::get_N(), cache_i, 0, y, u_n_est, isCache, cache);

            this->stopTimer();
            this->outTime();
//            cout << i+1 << " " << llr <<endl;

            if (exp(llr) >= 1.0) {
                h_i[i] = 0;
            } else {
                h_i[i] = 1;
            }
            u_n_est[i] = h_i[i];
            cout << "SC: " << i+1 << " " << llr << endl;
//            printf("%05f\n", llr);
        }
    }
    return u_n_est;
}

double Decoder::calcL_i(int i, int n, int cache_i, int level, vector<double> &y, vector<int> &u, vector<vector<bool> > &isCache, vector<vector<double> > &cache) {
    double llr = 0.0;
    this->addCount();
    if ( n == 1 ) {
        double wc = Channel::calcW(y[0],0);
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

void Decoder::calcSConBP(int itr, int count, ofstream &val_file, ofstream &check_file, vector<int> &u, vector<double> &y , vector<int> &u_n_est, vector<double> &tmp_u, vector<int> &A, vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked, vector<vector<bool> > &ym_isReceived) {
    int size = 2 * log2(Params::get_N()) + 2;

    vector<vector<vector<message> > > save_list(size, vector<vector<message> >(Params::get_N(), vector<message>()));
    vector<vector<double> > tmp_node(size, vector<double>(Params::get_N(), 0.0));
    vector<vector<double> > save_node(size, vector<double>(Params::get_N(), 0.0));
    vector<vector<bool> > save_isChecked(size, vector<bool>(Params::get_N(), false));
    vector<vector<bool> > tmp_isChecked(size, vector<bool>(Params::get_N(), false));

    vector<vector<double> > ym(size-1, vector<double>(Params::get_N(),0.0)); //dummy

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
        SCinit(i, y, u, u_n_est, A, node_value, message_list, save_list, node_isChecked, save_isChecked);
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
                calc_mp(tmpSize, node_value, save_list, save_isChecked, ym, ym_isReceived);
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
                calc_mp(tmpSize, node_value, message_list, node_isChecked, ym, ym_isReceived);
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
            cout << "BP: " << i+1 << " " << tmp_u[i] << endl;
        }
    }
}