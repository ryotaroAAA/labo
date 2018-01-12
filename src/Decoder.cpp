#include "../lib/Decoder.h"
#include "../lib/Analysor.h"
Decoder::Decoder(){

}

Decoder::~Decoder(){

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

void Decoder::send_message(int i, int j, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &y){
    int size = 2*log2(Params::get_N())+2;
    double wc,wp,llr = 0.0,val = 0.0;
    bool is_send = false;
    int n = Params::get_N();
    vector<vector<bool>> T(size, vector<bool>(n, false));
    Params::get_T(T);

    //中間ノードはここでおくる
    if(Common::is_mid_send()){
        if (T[i][j]) {
            wc = Channel::calcW(y[i/2][j],send_0);
            wp = Channel::calcW(y[i/2][j],1);
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
        }
        if (isnan(llr)) {
//            Common::pp(val);
//            cout << "check" << endl;
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
                //足して0ならtrue
                node_isChecked[i][j] = true;
            } else if(i != size-1) {
                //チャンネルファクター
                node_isChecked[i][j] = false;
            } else if(zero_flag){
                //ゼロがあると終わらない
                node_isChecked[i][j] = false;
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
    if(!Params::get_is_calc_bloop()){
        for (int i = 0; i < Params::get_N(); i++) {
            if(!puncFlag[i]){
                //puncなし
                if( Common::is_mid_send() ){
                    wc = Channel::calcW(ym[ysize-1][i],send_0);
                    wp = Channel::calcW(ym[ysize-1][i],1);
                } else {
                    wc = Channel::calcW(y[i],send_0);
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
    } else {
        //channelファクターノード
        for (int i = 0; i < Params::get_N(); i++) {
            node_isChecked[size - 1][i] = true;
        }
        //任意の位置をショートン、送信
        vector<vector<bool> > T(size, vector<bool>(Params::get_N(), false));
        Params::get_T(T);
        for (int i = 0; i < B.size(); i=i+2) {
            for (int j = 0; j < B[0].size(); j++) {
//                 if(B[i][j]){
//                     llr = (xm[i/2][j]==0)? inf_p : inf_m;
//                     node_val[i][j] = llr;
//                     node_isChecked[i][j] = true;
//                 }
                 if(B[i][j]){
                     //send_messageで送る処理はやってる
                     T[i][j] = true;
                 }
            }
        }
        Params::set_T(T);
    }
}

void Decoder::BPinit(vector<int> &u, vector<int> &x, vector<double> &y, vector<vector<int> > &xm, vector<vector<double> > &ym, vector<vector<double> > &node_val, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<int> > &B){
    vector<bool> puncFlag(Params::get_N(),false);
    //init node frozen punc etc pattern
    init_params(puncFlag, u, x, y, xm, ym, node_val, node_isChecked, B);
    init_message(puncFlag, node_val, message_list);
}

void Decoder::calc_mp(int size, vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked ,vector<vector<double> > &ym){
    calc_val_to_check(size, node_value, message_list, node_isChecked, ym);
    calc_check_to_val(size, node_value, message_list, node_isChecked, ym);
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

void Decoder::calc_check_to_val(int size, vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &y){
    vector<vector<int> > adjacent;
    bool bpflag = Params::get_decode_mode() == BP;
    for (int i = 1; i < size; i=i+2) {
        for (int j = 0; j < Params::get_N(); j++) {
            if( (!node_isChecked[i][j] || bpflag)  ) {
                send_message(i,j,message_list, node_isChecked, y);
            }
        }
    }
}

//val nodeへメッセージを集める　　
void Decoder::calc_marge(vector<vector<double> > &node_value, vector<vector<vector<message> >> &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &y){
    int size = 2*log2(Params::get_N())+2;
    int ysize = log2(Params::get_N())+1;
    int level = 0, index = 0;
    vector<vector<double> > temp(size, vector<double>(Params::get_N(),0.0));

    int n = Params::get_N();

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

    //中間ノードからのメッセージを引き上げる
    if( Common::is_mid_send() ){
        vector<vector<bool>> T(size, vector<bool>(n, false));
        Params::get_T(T);
        for (int i = 0; i < ysize; i++) {
            for (int j = 0; j < Params::get_N(); j++) {
                if (T[i][j]) {
                    double wc = Channel::calcW(y[i][j],send_0);
                    double wp = Channel::calcW(y[i][j],1);
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
            if(temp[i][j] != 0.0 && !node_isChecked[i][j]){
                node_value[i][j] = temp[i][j];
            }
        }
    }
}

vector<int> Decoder::calcBP(int loop_num, vector<int> &param, vector<int> &u, vector<int> &x, vector<double> &y, vector<vector<int> > &xm, vector<vector<double> > &ym, vector<vector<int> > &node_error_count, ofstream &val_error_file, vector<vector<int> > &B) {
    vector<double> tmp_u(Params::get_N(),0.0);
    vector<int> u_n_est(Params::get_N());
    vector<vector<int> > adjacent;

    //node初期化、奇数が変数ノード（1,3,5...）、偶数がチェックノード（2,4,6...）
    int size = 2*log2(Params::get_N())+2;

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
    if (Params::get_decode_mode() == BP) {
        BPinit(u, x, y, xm, ym, node_value, message_list, node_isChecked, B);

        for (int i = 0; i < Params::get_rp(); i++) {
            if(isTerminate(no_checked, node_value, node_isChecked)) break;
            if(Params::get_is_outlog()) {
                printDecodeProgress(itr, node_value, val_file);
                printDecodeProgress(itr, node_isChecked, check_file);
            }
            itr++;
            calc_mp(size, node_value, message_list, node_isChecked, ym);

            confirmIsCheck(node_value, node_isChecked);
            count++;
        }
    } else {

    }

    int error_count = 0;
    Analysor::errorCount(u, u_n_est, &error_count);
    param[0] = itr;
    param[1] = no_checked;

    for (int i = 0; i < Params::get_N(); i++) {
        if(Params::get_is_outlog()) {
//            cout << "BP: " << i + 1 << " " << node_value[0][i] << endl;
        }
        u_n_est[i] = (node_value[0][i] > 0) ? 0 : 1;
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

    outLog(itr, no_checked, u, u_n_est, val_file, check_file , node_value, node_isChecked);

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

vector<int> Decoder::decode(vector<double> &y, vector<int> &u){
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

            if (exp(llr) >= 1.0) {
                h_i[i] = 0;
            } else {
                h_i[i] = 1;
            }
            u_n_est[i] = h_i[i];
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
