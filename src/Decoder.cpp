#include "../lib/Decoder.h"
Decoder::Decoder(){

}

Decoder::~Decoder(){

}

vector<vector<int> > Decoder::adjacentIndex(int level, int index){
    vector<vector<int> > adjacent;
    vector<int> temp;
    int n = (level != 1) ? Params::get_N()/pow(2,((level-2)/2)) : 1;
    temp = makeBPTreeIndex(n);
    if(level%2 == 0) {
        //check_node
        int temp_i = (index==n) ? n-1 : (index-1)%n;
        int div = (index-1)/n;
        if(index%2 == 1){
            //index_o
            if(level < 2*log2(Params::get_N())) adjacent = {{level-1, index}, {level-1, index+1}, {level+1, temp[temp_i]+div*n}};
            else adjacent = {{level-1, index}, {level-1, index+1}, {level+1, index}};
        } else {
            // /index_e
            if(level < 2*log2(Params::get_N())) adjacent = {{level-1, index}, {level+1, temp[temp_i]+div*n}};
            else adjacent = {{level-1, index}, {level+1, index}};
        }
    } else {
        //val_node
        int temp_i = 0;
        int div = 0;
        vector<int> temp_r(Params::get_N(),0);
        if(n < Params::get_N() && n != 1){
            for (int i = 0; i < Params::get_N(); i++) {
                temp_r[i] = temp[i%n];
            }
            for (int i = 0; i < Params::get_N(); i++) {
                if (temp_r[i] == (index-1)%n+1) {
                    div = (index-1)/n;
                    temp_i = i+1+div*n;
                    break;
                }
            }
        } else {
            for (int i = 0; i < Params::get_N(); i++) {
                if (temp[i] == index) {
                    temp_i = i+1;
                    break;
                }
            }
        }

        if(index%2 == 1){
            //index_o
            //一番右
            if(level == 2*log2(Params::get_N())+1) adjacent = {{level-1, index}};
            //一番左
            else if(level == 1) adjacent = {{level+1, index}};
            //それ以外
            else {
                adjacent = {{level + 1, index}, {level - 1, temp_i}};
            }
        } else {
            //index_e
            //一番右
            if(level == 2*log2(Params::get_N())+1) adjacent = {{level-1, index}};
            //一番左
            else if(level == 1) adjacent = {{level+1, index}, {level+1, index-1}};
            //それ以外
            else {
                adjacent = {{level+1, index}, {level+1, index-1}, {level-1, temp_i}};
            }
        }
    }
    return adjacent;
}

int Decoder::take_val(vector<int> locate, vector<vector<int> > &node_val){
    int level = locate[0];
    int index = locate[1];
    return node_val[level-1][index-1];
}

double Decoder::take_val(vector<int> locate, vector<vector<double> > &node_val){
    int level = locate[0];
    int index = locate[1];
    return node_val[level-1][index-1];
}

bool Decoder::take_val(vector<int> locate, vector<vector<bool> > &node_val){
    int level = locate[0];
    int index = locate[1];
    return node_val[level-1][index-1];
}

double Decoder::calc_node(int mode, double a, double b) {
    double llr = 0.0;
    if (mode == 0){
        //val
        llr = a + b;
    } else{
        //check
        llr = 2.0 * atanh(tanh(a/2.0) * tanh(b/2.0) );
    }
    return llr;
}

bool Decoder::isCalculable(int level, int index, vector<vector<int> > &adjacent, vector<vector<double> > &node_val, vector<vector<bool> > &node_isChecked){
    bool flag = true;
    double val=0.0,val1=0.0,val2=0.0,val3=0.0;
    bool ch=false, ch1=false, ch2=false, ch3=false;
    int adjacent_count = adjacent.size();
    int zero_count = 0;
    val = node_val[level-1][index-1];
    ch = node_isChecked[level-1][index-1];
    flag = (ch == false) ? true : false;
    val1 = take_val(adjacent[0], node_val);
    ch1 = take_val(adjacent[0], node_isChecked);
    if(val1 == 0) zero_count++;
    if(adjacent_count > 1){
        val2 = take_val(adjacent[1], node_val);
        ch2 = take_val(adjacent[1], node_isChecked);
        if(val2 == 0) zero_count++;
        if(adjacent_count == 3) {
            val3 = take_val(adjacent[2], node_val);
            ch3 = take_val(adjacent[2], node_isChecked);
            if (val3 == 0) zero_count++;
        }
    }

    if (level%2 == 0 && zero_count > 1){
        //checkが通ってるならやらない
        flag = false;
    }
    if(zero_count == adjacent_count) {
        flag = false;
    }

    return flag;
}

bool Decoder::isChecked(int level, int index, vector<vector<int> > &adjacent, vector<vector<double> > &node_val, vector<vector<bool> > &node_isChecked){
    bool flag = false;
    double val,val1,val2,val3,tmp;
    bool ch, ch1, ch2, ch3;
    int check_val = 0;
    int adjacent_count = adjacent.size();

    val = node_val[level-1][index-1];
//    ch = node_isChecked[level-1][index-1];
//    flag = (ch == false) ? true : false;

    val1 = take_val(adjacent[0], node_val);
    ch1 = take_val(adjacent[0], node_isChecked);
    val2 = take_val(adjacent[1], node_val);
    ch2 = take_val(adjacent[1], node_isChecked);

    if(val1 < 0.0) check_val++;
    if(val2 < 0.0) check_val++;

    if(adjacent_count == 3){
        val3 = take_val(adjacent[2], node_val);
        ch3 = take_val(adjacent[2], node_isChecked);
        if(val3 < 0.0) check_val++;
    }
    flag = (check_val%2 == 0) ? true : false;
    if(val1 == 0 || val2 == 0 || val3 == 0) flag = false;
    return flag;
}

bool Decoder::isTerminate(vector<vector<double> > &node_val, vector<vector<bool> > &node_isChecked) {
    vector<vector<int> > adjacent;
    bool flag = true;
    for (int i = 0; i < node_val.size(); i++) {
        for (int j = 0; j < Params::get_N(); j++) {
            //check_nodeのみ調べる
            if( (i+1)%2 == 0 ){
                //隣接ノードのチェックを調べる
                adjacent = adjacentIndex(i + 1, j + 1);
                if(isChecked(i+1, j+1, adjacent, node_val, node_isChecked) == false) {
                    flag = false;
                    break;
                } else {
//                    node_isChecked[i][j] = true;
                }
            }
        }
    }
    return flag;
}

void Decoder::BPinit(vector<double> &y, vector<int> &u, vector<int> &A, vector<vector<double> > &node_val, vector<vector<int> > &update_count, vector<vector<bool> > &node_isChecked){
    double wc, wp, llr;
    int size = 2*log2(Params::get_N())+1;
    bool databit_flag = false;
    for (int i = 0; i < Params::get_N(); i++) {
        wc = Channel::calcW(y[i],0);
        wp = Channel::calcW(y[i],1);
        llr = 1.0 * log(wc / wp);
        node_val[size-1][i] = llr;
        update_count[size-1][i] = 1;
        //frozen_bitのllr設定
        databit_flag = false;
        for (int j = 0; j < Params::get_K(); j++) {
            if(i == A[j]){
                databit_flag = true;
            }
        }
        if(databit_flag == false){
            node_val[0][i] = (u[i] == 0) ? 10.0 : -10.0;
            update_count[0][i] = 10000000;
            node_isChecked[0][i] = true;
        }
    }
}

vector<int> Decoder::BP(int limit, vector<double> &y, vector<int> &u, vector<int> &A){
    vector<int> u_n_est(Params::get_N());
    vector<vector<int> > adjacent;

    //node初期化、奇数が変数ノード（1,3,5...）、偶数がチェックノード（2,4,6...）
    int size = 2*log2(Params::get_N())+1;
    vector<vector<double> > node_value(size, vector<double>(Params::get_N(),0.0));
    vector<vector<int> > update_count(size, vector<int>(Params::get_N(),0));
    vector<vector<bool> > node_isChecked(size, vector<bool>(Params::get_N(),false));
    //yとかの値をセット
    BPinit(y, u, A, node_value, update_count, node_isChecked);

    //BP
    int count = 0;
    while(count <= limit) {
        for (int i = count%2; i < size; i=i+2){
            for (int j = 0; j < Params::get_N(); j++) {
                int adjacent_count = 0, up1 = 0, up2 = 0, up3 = 0, min_up = 0;
                double val1 = 0.0, val2 = 0.0, val3 = 0.0;
                bool ch1 = false, ch2 = false, ch3 = false;

                //隣接ノードの位置
                adjacent = adjacentIndex(i + 1, j + 1);
                adjacent_count = adjacent.size();
                if (isCalculable(i + 1, j + 1, adjacent, node_value, node_isChecked) == true &&
                    node_isChecked[i][j] == false) {
                    if (i % 2 == 0) {
                        //val_node
                        val1 = take_val(adjacent[0], node_value);
                        up1 = take_val(adjacent[0], update_count);
//                    ch1 = take_val(adjacent[0], node_isChecked);
                        switch (adjacent_count) {
                            case 1 :
                                if (val1 != 0.0) {
                                    node_value[i][j] = val1;
                                    update_count[i][j]++;
                                }
                                break;
                            case 2 :
                                val2 = take_val(adjacent[1], node_value);
                                up2 = take_val(adjacent[1], update_count);
                                if (val1 != 0.0 && val2 != 0.0) {
                                    if(i == 0){
                                        node_value[i][j] = calc_node(0, val1, val2);
                                        update_count[i][j]++;
                                    } else {
                                        //違った?
                                        if (up1 > up2) {
                                            node_value[i][j] = val1;
                                            update_count[i][j]++;
                                        } else if (up1 < up2) {
                                            node_value[i][j] = val2;
                                            update_count[i][j]++;
                                        } else {
                                            //たぶんここにはこないはず後で乱数を絡ませる
                                            node_value[i][j] = val1;
                                            update_count[i][j]++;
                                        }
                                    }
                                } else if (val1 == 0.0 && val2 != 0.0) {
                                    node_value[i][j] = val2;
                                    update_count[i][j]++;
                                } else if (val1 != 0.0 && val2 == 0.0) {
                                    node_value[i][j] = val1;
                                    update_count[i][j]++;
                                }
                                break;
                            case 3 :
                                val2 = take_val(adjacent[1], node_value);
                                up2 = take_val(adjacent[1], update_count);
                                val3 = take_val(adjacent[2], node_value);
                                up3 = take_val(adjacent[2], update_count);
                                if (val1 != 0.0 && val2 != 0.0 && val3 != 0.0) {
                                    //全部同じ更新回数だった場合は適当になる
                                    min_up = min(min(up1, up2), up3);
                                    //全部違う時
                                    if (min_up == up1) {
                                        node_value[i][j] = calc_node(0, val2, val3);
                                    } else if (min_up == up2) {
                                        node_value[i][j] = calc_node(0, val1, val3);
                                    } else {
                                        node_value[i][j] = calc_node(0, val1, val2);
                                    }
                                    update_count[i][j]++;
                                } else if (val1 == 0.0 && val2 != 0.0 && val3 != 0.0) {
                                    node_value[i][j] = calc_node(0, val2, val3);
                                    update_count[i][j]++;
                                } else if (val1 != 0.0 && val2 == 0.0 && val3 != 0.0) {
                                    node_value[i][j] = calc_node(0, val1, val3);
                                    update_count[i][j]++;
                                } else if (val1 != 0.0 && val2 != 0.0 && val3 == 0.0) {
                                    node_value[i][j] = calc_node(0, val1, val2);
                                    update_count[i][j]++;
                                } else if (val1 == 0.0 && val2 == 0.0 && val3 != 0.0) {
                                    node_value[i][j] = val3;
                                    update_count[i][j]++;
                                } else if (val1 != 0.0 && val2 == 0.0 && val3 == 0.0) {
                                    node_value[i][j] = val1;
                                    update_count[i][j]++;
                                } else if (val1 == 0.0 && val2 != 0.0 && val3 == 0.0) {
                                    node_value[i][j] = val2;
                                    update_count[i][j]++;
                                } else {
                                    //何もしない
                                    cout << "aaaaaaaaaaaaaaa" << endl;
                                }
                                break;
                            default :
                                break;
                        }
                        //一回決まったら値を確定させる
//                        if (i == 0 && node_value[i][j] != 0) {
//                            node_isChecked[i][j] = true;
//                            node_value[i][j] = (node_value[i][j] > 0) ? 10.0 : -10.0;
//                        }
                    } else {
                        //check_node
                        val1 = take_val(adjacent[0], node_value);
                        up1 = take_val(adjacent[0], update_count);
                        switch (adjacent_count) {
//                    ch2 = take_val(adjacent[1], node_isChecked);
                            case 2 :
                                val2 = take_val(adjacent[1], node_value);
                                up2 = take_val(adjacent[1], update_count);
                                if (val1 != 0.0 && val2 != 0.0) {
//                                    node_value[i][j] = calc_node(1, val1, val2);
//                                    update_count[i][j]++;
                                    //違った?
                                    if (up1 > up2) {
                                        node_value[i][j] = val1;
                                        update_count[i][j]++;
                                    } else if (up1 < up2) {
                                        node_value[i][j] = val2;
                                        update_count[i][j]++;
                                    } else {
                                        //たぶんここにはこないはず
                                        node_value[i][j] = val1;
                                        update_count[i][j]++;
                                    }
                                } else if (val1 == 0.0 && val2 != 0.0) {
                                    node_value[i][j] = val2;
                                    update_count[i][j]++;
                                } else if (val1 != 0.0 && val2 == 0.0) {
                                    node_value[i][j] = val1;
                                    update_count[i][j]++;
                                }
                                break;
                            case 3 :
                                val2 = take_val(adjacent[1], node_value);
                                up2 = take_val(adjacent[1], update_count);
                                val3 = take_val(adjacent[2], node_value);
                                up3 = take_val(adjacent[2], update_count);

                                if (val1 != 0.0 && val2 != 0.0 && val3 != 0.0) {
                                    //全部同じ更新回数だった場合は適当になる?
                                    min_up = min(min(up1, up2), up3);
                                    if (min_up == up1) {
                                        node_value[i][j] = calc_node(1, val2, val3);
                                    } else if (min_up == up2) {
                                        node_value[i][j] = calc_node(1, val1, val3);
                                    } else {
                                        node_value[i][j] = calc_node(1, val1, val2);
                                    }
                                    update_count[i][j]++;
                                } else if (val1 == 0.0 && val2 != 0.0 && val3 != 0.0) {
                                    node_value[i][j] = calc_node(1, val2, val3);
                                    update_count[i][j]++;
                                } else if (val1 != 0.0 && val2 == 0.0 && val3 != 0.0) {
                                    node_value[i][j] = calc_node(1, val1, val3);
                                    update_count[i][j]++;
                                } else if (val1 != 0.0 && val2 != 0.0 && val3 == 0.0) {
                                    node_value[i][j] = calc_node(1, val1, val2);
                                    update_count[i][j]++;
                                } else if (val1 == 0.0 && val2 == 0.0 && val3 != 0.0) {
                                    node_value[i][j] = val3;
                                    update_count[i][j]++;
                                } else if (val1 != 0.0 && val2 == 0.0 && val3 == 0.0) {
                                    node_value[i][j] = val1;
                                    update_count[i][j]++;
                                } else if (val1 == 0.0 && val2 != 0.0 && val3 == 0.0) {
                                    node_value[i][j] = val2;
                                    update_count[i][j]++;
                                } else {
                                    //何もしない
                                    cout << "aaaaaaaaaaaaaaa" << endl;
                                }
                                break;
                            default :
                                break;
                        }
                    }
                }
            }
        }
        if(isTerminate(node_value, node_isChecked)) break;
        count++;
    }
    for (int i = 0; i < Params::get_N(); i++) {
        u_n_est[i] = (node_value[0][i]>0) ? 0 : 1;
    }
    return u_n_est;
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

    double llr;
    int cache_i = 0;

    //u_n_est計算
    for (int i = 0; i < Params::get_N(); i++) {
        // Aに含まれないindexなら既知
        if (Common::containNumInArray(i, Params::get_N()-Params::get_K(), A) == false) {
            u_n_est[i] = u[i];
        } else {
            this->startTimer();

            cache_i = makeTreeIndex(Params::get_N())[i] - 1;
            llr = calcL_i(i+1, Params::get_N(), cache_i, 0, y, u_n_est, isCache, cache);

            this->stopTimer();
            this->outTime();
//            cout << i+1 << "  " << llr <<endl;

            if (exp(llr) >= 1.0) {
                h_i[i] = 0;
            } else {
                h_i[i] = 1;
            }
            u_n_est[i] = h_i[i];
        }
    }
    return u_n_est;
}

vector<int> Decoder::decode_m(vector<vector<double> > &y, vector<int> &u, vector<int> &A){
    vector<double> h_i(Params::get_N());
    vector<int> u_n_est(Params::get_N());
    int size = log2(Params::get_N());

    vector<vector<bool> > isCache (size, vector<bool>(Params::get_N(),false));
    vector<vector<double> > cache (size, vector<double>(Params::get_N(),0.0));

    int cache_i = 0;

    vector<int> tree = makeTreeIndex(Params::get_N());
    for (int i = 0; i<log2(Params::get_N())-1; i++) {
        for (int j = 0; j < Params::get_N(); j++) {
            if (y[i][j] == 0){
                isCache[i][ tree[j]-1 ] = true;
                cache[i][ tree[j]-1 ] = 100.0;
            } else if(y[i][j] == 1) {
                isCache[i][ tree[j]-1 ] = true;
                cache[i][ tree[j]-1 ] = -100.0;
            }
        }
    }

    double llr;


    //u_n_est計算
    for (int i = 0; i < Params::get_N(); i++) {
        // Aに含まれないindexなら既知
        if (Common::containNumInArray(i, Params::get_N()-Params::get_K(), A) == false) {
            u_n_est[i] = u[i];
        } else {
            this->startTimer();

            cache_i = makeTreeIndex(Params::get_N())[i] - 1;
            llr = calcL_i(i+1, Params::get_N(), cache_i, 0, y[log2(Params::get_N())- 1], u_n_est, isCache, cache);

            this->stopTimer();
            this->outTime();
//            cout << i+1 << "  " << llr <<endl;
//            if (exp(llr) >= 1.0) {
            if (llr >= 0.0) {
                h_i[i] = 0;
            } else {
                h_i[i] = 1;
            }
            u_n_est[i] = h_i[i];
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
