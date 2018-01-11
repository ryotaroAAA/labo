#include "../lib/Analysor.h"

Analysor::Analysor(){

}

Analysor::~Analysor(){

}

void Analysor::errorCount(vector<int> &u, vector<int> &u_est, int* error_count){
    for(int i=0; i<Params::get_N(); i++){
        if(u[i] != u_est[i]){
            *error_count = *error_count + 1;
        }
    }
}

double Analysor::errorCalc(vector<int> &u, vector<int> &u_est, int* error_count){
    for(int i=0; i<Params::get_N(); i++){
        if(u[i] != u_est[i]){
            *error_count = *error_count + 1;
        }
    }

    return (double)(*error_count)/Params::get_N();
}

double Analysor::calcCapacity(int i, int n) {
    double cap =0.0;
    if(Params::get_s() == BEC){
        if ( i == 0 && n == 1 ) {
            cap = 1 - Params::get_e();
        } else {
            if ( i % 2 == 0) {
                cap = pow(Analysor::calcCapacity(i/2, n/2),2);
            } else {
                double tempCap = Analysor::calcCapacity((i-1)/2, n/2);
                cap = 2 * tempCap - pow(tempCap,2);
            }
        }
    } else if (Params::get_s() == BSC){

    }
    return cap;
}

double Analysor::calcBhatforBEC(int i, int n){
    double bha =0.0;

    if ( i == 0 && n == 1 ) {
        bha = Params::get_e();
    } else {
        if ( i % 2 == 0) {
            double tempBha = Analysor::calcBhatforBEC(i/2, n/2);
            bha = 2 * tempBha - pow(tempBha, 2);
        } else {
            bha = pow(Analysor::calcBhatforBEC((i-1)/2, n/2),2);
        }
    }
    return bha;
}

void Analysor::makeArrayCapacity(vector<double> &array) {
    for(int i = 0; i < Params::get_N(); i++) {
        array[i] = Analysor::calcCapacity(i, Params::get_N());
    }
}

double Analysor::calc_inv(double y){
    double x1 = 10.0;
    double x2 = 500.0;
    double x3 = 0.0;
    double fx = 0.0;
    int count = 1;
    while (abs(x1 - x2) >= 0.00001) {
        x3 = (x1+x2) / 2.0;
        fx = sqrt(M_PI/x3) * exp(-1.0 * (x3/4.0)) * (1.0 - 10.0/(7.0 * x3));
        if(fx > y){
            x1 = x3;
        } else {
            x2 = x3;
        }
        count++;
//        cout << "est:" <<  x1 << " " << x2 << " " << m_func(x3)<< " true:" << y << endl;
    }
    return x3;
}

double Analysor::inv_m_func(double y){
    double ret = 0.0;
    const double a = -0.4527;
    const double b = 0.0218;
    const double c = 0.86;
    double temp = 0.0;
    //0.03847)
    if (y > 0.03847){
        ret = pow((log(y)-b)/a, 1.0/c);
    } else if(y < 0.03847) {
        ret = calc_inv(y);
    } else {
//        ret = __nan();
    }
    return ret;
}

double Analysor::m_func(double x){
    double ret = 0.0;
    const double a = -0.4527;
    const double b = 0.0218;
    const double c = 0.86;
    double temp1 = 0.0;
    double temp2 = 0.0;
    if (x < 10.0) {
        ret = exp(a * pow(x,c) + b);
    } else if (x > 10.0) {
        temp1 = 1.0 - 3.0/x;
        temp2 = 1.0 + 1.0/(7.0 * x);
        ret = 0.5 * sqrt(M_PI/x) * exp(-1.0 * x / 4.0)  * (temp1 + temp2);
    } else {
//        ret = __nan();
    }
    return ret;
}

double Analysor::calc_m_in(int i, int n){
    double ret = 0.0;
    double temp1 = 0.0;
    double temp2 = 0.0;
    int temp_i = (i/2 == 0)? 1 : i/2;
    if (n == 1){
        ret = 2.0 / Params::get_e();
    } else {
        temp1 = Analysor::calc_m_in(temp_i, n/2);
        if (i % 2 == 0){
            ret = 2.0 * temp1;
        } else {
            temp2 = 1.0 - pow((1.0 - Analysor::m_func(temp1)), 2);
            ret = Analysor::inv_m_func(temp2);
        }
    }
//    cout << "[" << temp_i << "][" << n  << "] "  << ret << endl;
    return ret;
}

void Analysor::makeArrayBhat(vector<double> &array) {
//    if ( Params::get_s() == BSC) {
    if ( 1 ) {
        Performance performance;
        Decoder decoder;
        Encoder encoder;
        Logger logger;

        vector<int> u(Params::get_N(), 0);
        vector<int> x(Params::get_N(), 0);
        vector<double> y(Params::get_N(), 0.0);
        int u_n = 0;

        double sumbha = 0.0;
        double tempbha = 0.0;
        int size = log2(Params::get_N());

        vector<vector<bool> > isCache (size, vector<bool>(Params::get_N(),false));
        vector<vector<double> > cache (size, vector<double>(Params::get_N(),0.0));

        double lr;
        int cache_i = 0;
        int repeatNum = Params::get_monteNum();

        for (int m = 1; m <= repeatNum; m++) {
            performance.startTimer();
//            cout << "[" << m << "]" << endl;
            Preseter::preset_u(RAND, u);
            vector<int> temp_u = u;
            x = encoder.encode(Params::get_N(), u);

            for (int i = 0; i < Params::get_N(); i++) {
                temp_u = u;
                u_n = temp_u[i];
                temp_u.resize(i);

                y = Channel::channel_output(x);
                cache_i = decoder.makeTreeIndex(Params::get_N())[i] - 1;
                lr = exp(decoder.calcL_i(i+1, Params::get_N(), cache_i, 0, y, temp_u, isCache, cache));
                lr = pow(lr, -1+2 * u_n);
                tempbha = sqrt(lr);
                if (tempbha > 1.0) tempbha = 1.0/tempbha;
                if (isnan(tempbha)) tempbha = 0.0;
                array[i] += 1.0 * tempbha/repeatNum;
            }
        }
    } else {
        vector<double> m_in(Params::get_N(), 0.0);
        vector<double> sigma2_in(Params::get_N(), 0.0);
        double temp = 0.0;
        for (int i = 0; i < Params::get_N(); i++) {
            temp = Analysor::calc_m_in(i+1, Params::get_N());
//            if(isnan(temp)){
//                m_in[i] = 0.0;
//            } else {
                m_in[i] = temp;
//            }
//            cout << i << " mean " << m_in[i] <<endl;
            sigma2_in[i] = 2.0/m_in[i];
            array[i] = exp(-1.0/(2.0 * sigma2_in[i]));
        }
    }

    if(Params::get_m_mode() != MID_ADOR){
//        string filename = Params::get_rvbDir() + " Bhat";
//        ofstream w_file;
//        w_file.open(filename, ios::out);
//        for (int i = 0; i < Params::get_N(); i++)
//        {
//            w_file << i << " " << array[i] << endl;
//            cout << i << " " << array[i] << endl;
//        }
//        Common::bar();
    }
}

void Analysor::probErrBound() {
    int count = 0;
    vector<double> tempArr;
    vector<double> sumArr;
    double sum = 0.0;

    for(int i = 0; i < Params::get_N(); i++) {
        tempArr.push_back(Analysor::calcBhatforBEC(i, Params::get_N()));
    }

    sort(tempArr.begin(), tempArr.end(), greater<int>());

    for(int i=0; i < Params::get_N(); i++){
        sum += tempArr[i];
        sumArr.push_back(sum);
    }

    string filename = "log/eb10";
    ofstream w_file;
    w_file.open(filename, ios::out);
    for (int i = 0; i < Params::get_N(); i++)
    {
        w_file << (double)i/Params::get_N() << " " << sumArr[i] << " " << tempArr[i] << endl;
        cout << (double)i/Params::get_N() << " " << sumArr[i] << " " << tempArr[i] << endl;
    }
}


string Analysor::get_itrfn(){
    time_t now = time(NULL);
    struct tm *pnow = localtime(&now);
    stringstream itrfn;
    EXP_MODE em = Params::get_exp_mode();
    MID_MODE mm = Params::get_m_mode();
    string ename;
    switch (em) {
        case NORMAL: ename = ""; break;
        case PUNC: ename = "punc"; break;
        case MID:
            switch (mm) {
                case MID_BLUTE: ename = "mid_blute"; break;
                case MID_ADOR: ename = "mid_ador"; break;
                case MID_AOR: ename = "mid_aor"; break;
                case MID_DOR: ename = "mid_dor"; break;
                case MID_AOB: ename = "mid_aob"; break;
                case MID_DOB: ename = "mid_dob"; break;
                case MID_AOV: ename = "mid_aov"; break;
                case MID_DOV: ename = "mid_dov"; break;
            }
            break;
        case QUP: ename = "qup"; break;
        case WANG: ename = "wang"; break;
        case M_QUP: ename = "m_qup"; break;
        case M_WANG: ename = "m_wang"; break;
        case VALERIO_P: ename = "valerio_p"; break;
        case VALERIO_S: ename = "valerio_s"; break;
        case M_VALERIO_P: ename = "m_valerio_p"; break;
        case M_VALERIO_S: ename = "m_valerio_s"; break;
    }
    itrfn << "log/"
       << to_string(pnow->tm_year + 1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday)
       << "/itr_N" << to_string(Params::get_N())
       << "e" << Params::get_e()
       << "BNum" << to_string(Params::get_blockNum())
       << "UPBn" << to_string(Params::get_upperBlockErrorNum())
       << "c" << (Params::get_s() ? "BSC" : "BEC") << ":" << ename;
    return itrfn.str();
}

void Analysor::printDecodeProgress(int count, vector<vector<int> > &node_value, ofstream &w_file){
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

//rate vs BERのグラフ作成用
void Analysor::calcBlockErrorRate() {
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    logger.setRvbDir(Params::get_rvbDir());
    vector<int> u_n(Params::get_N(), 0);
    vector<int> x_n(Params::get_N(), 0);
    vector<double> y_n(Params::get_N(), 0);

    vector<int> u_est(Params::get_N(), 0);
    vector<pair<int, double> > cap_map;
    Preseter::makeMutualInfoArray(cap_map);

    double BER = 0.0;
    double sumBER = 0.0;
    double rate = 0.0;
    int n = Params::get_N();
    int error_count = 0;
    int block_error_count = 0;
    int loopi = 0;

    int pointNum = 10;
    double from = 0.0;
    double to = 0.7;
    int startK = get_eachK(from);
    int endK = get_eachK(to);
    double interval = (double)(endK-startK)/pointNum;

    EXP_MODE em = Params::get_exp_mode();
    cout << "startK::" << startK << " interval::" << interval << endl;

    for (int i = 0; i < pointNum; i++) {
        loopi = 0;
        sumBER = 0.0;
        block_error_count = 0;

        Params::set_K((int)(startK + i * interval));
        performance.startTimer();
        vector<int> A(Params::get_K(), -1);
        vector<int> Ac;

        Preseter::represet_A(A, Ac, cap_map);
        cout << "k : " << Params::get_K() << endl;

        while (block_error_count < Params::get_upperBlockErrorNum()) {
            Preseter::preset_u(RAND, u_n);
            x_n = encoder.encode(Params::get_N(), u_n);
            loopi++;

            y_n.assign(Params::get_N(), 0);
            u_est.assign(Params::get_N(), 0);
            y_n = Channel::channel_output(x_n);
            u_est = decoder.decode(y_n, u_n);

//            Common::pp(u_est);
            Analysor::errorCount(u_n, u_est, &error_count);
            if(error_count > 0) block_error_count++;
            if (loopi % 50 == 0 ) cout << loopi << " " << error_count << " " << block_error_count << " " << (double)block_error_count/loopi << endl;
            error_count = 0;

//            if(loopi >= repeatNum || i * tmpK < Params::get_N()/7) break;
            if(loopi >= Params::get_blockNum()) break;
        }
        performance.stopTimer();
        rate = (double) Params::get_K() / Params::get_N();
        cout << Params::get_rvbDir() << endl;
        BER = (double)block_error_count/loopi;

        logger.outLog("=================================");
        logger.outLog(performance.outTime("処理時間"));
        logger.outLog("(N,K) = (" + to_string(Params::get_N()) + "," + to_string(Params::get_K()) + ")");
        logger.outLog("BER:" + to_string(BER));
        logger.outLog("Rate:" + to_string(rate));
        logger.outLog(encoder.outCount("encoder_count"));
        logger.outLog(decoder.outCount("decoder_count"));
        performance.outHMS();

        logger.outLogRVB(rate, BER);
//        if(mode == TEST) break;
    }
}

//fromで得たいレートを渡すと対応するKを返す
int Analysor::get_eachK(double from){
    EXP_MODE em = Params::get_exp_mode();
    int retK = 0;
    int n = Params::get_N();
    int k = Params::get_K();
    int m = Params::get_M();
    int mn = Params::get_MN();

    switch (em){
        case NORMAL:
            retK = n*from;
            break;
        case PUNC:
        case QUP:
        case WANG:
        case VALERIO_P:
            retK = (n-m)*from;
            break;
        case VALERIO_S:
            retK = m+(n-m)*from;
            break;
        case MID:
            retK = (n+mn)*from;
            break;
        case M_WANG:
        case M_QUP:
        case M_VALERIO_P:
            retK = (n-m+mn)*from;
            break;
        case M_VALERIO_S:
            retK = m+(n-m+mn)*from;
            break;
        default:
            break;
    }
    return retK;
}

//rate vs BERのグラフ作成用
void Analysor::calcBlockErrorRate_BP() {
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    logger.setRvbDir(Params::get_rvbDir());

    int n = Params::get_N(), error_count = 0, block_error_count = 0, loopi = 0, itr = 0;
    int tmp = Params::get_N() / Params::get_K(), tmpK = Params::get_K();
    int size = log2(Params::get_N()) + 1;
    double BER = 0.0, sumBER = 0.0, rate = 0.0, mitr = 0.0;;

    vector<int> A;
    vector<int> Ac;
    vector<vector<int> > B(2*log2(Params::get_N())+2, vector<int>(Params::get_N(),0));
    vector<int> u(Params::get_N(), 0);
    vector<int> x(Params::get_N(), 0);
    vector<double> y(Params::get_N(), 0);
    vector<vector<int> > tmp_x(log2(Params::get_N()), vector<int>(Params::get_N(), 0));
    vector<vector<int> > xm(size, vector<int>(Params::get_N(), 0));
    vector<vector<double> > ym(size, vector<double>(Params::get_N(), 0.0));
    vector<int> u_est(Params::get_N(), 0);
    vector<pair<int, double> > cap_map;

    Preseter::makeMutualInfoArray(cap_map);
    Preseter::preset_A_Ac(A, Ac);

    vector<int> param(2, 0);
    //右ノードのパンクチャorショートン位置を指定
    vector<int> p_0(Params::get_M(), -1);
    vector<int> p(Params::get_M(), -1);

    ofstream itrfile;
    itrfile.open(get_itrfn(), ios::out);

    int bloop = Params::get_Bloop();
    double point[3];
    Params::get_point(point);
    int pointNum = point[0];
    double from = point[1];
    double to = point[2];
    int startK = get_eachK(from);
    int endK = get_eachK(to);
    double interval = (double)(endK-startK)/pointNum;
    vector<vector<int> > node_error_count(2*log2(Params::get_N())+2, vector<int>(Params::get_N(),0));

    ofstream val_error_file;
    string val_error_fn = "/Users/ryotaro/Dropbox/labo/graph_js/val_error.json";
    val_error_file.open(val_error_fn, ios::out);
    val_error_file << "{" << endl;

    ofstream b_file;
    string b_fn = "/Users/ryotaro/Dropbox/labo/graph_js/b.json";
    b_file.open(b_fn, ios::out);
    b_file << "{" << endl;

    ofstream saveb_file;
    struct stat st;
    if(Params::get_is_calc_bloop()){
        string dir = "save_b/N_"+to_string(Params::get_N());
        int ret = stat(dir.data(), &st);
        if (ret == -1) {
            mkdir(dir.data(), S_IRUSR | S_IWUSR | S_IXUSR);
            cout <<dir.data()<<endl;
        }
        string saveb_fn = "save_b/N_"
                          + to_string(Params::get_N())
                          + "/Bloop_" + to_string(bloop);
        saveb_file.open(saveb_fn, ios::out);
    }

    EXP_MODE em = Params::get_exp_mode();
    vector<int> saveb;

    //AWWGNの場合
    double awgn_p[3];
    if(Params::get_is_exp_awgn()){
        Params::get_awgn_p(awgn_p);
        pointNum = (awgn_p[2] - awgn_p[1])/0.5 + 1.0;
    }

    for (int i = 0; i < pointNum; i++) {
        loopi = 0, sumBER = 0.0, mitr = 0.0, block_error_count = 0.0;
        performance.startTimer();

        double h = 0.0;
        //AWWGNの場合
        if(Params::get_is_exp_awgn()){
            double sig = 0.0;
            double rate = awgn_p[0];
            h = awgn_p[1] + 0.5*i;
            sig = sqrt(2.0/(rate*pow(10,h)));
            Params::set_e(sig);
            Params::set_K(get_eachK(rate));
            cout << "Rate " <<  rate << ", K " <<  get_eachK(rate) << ", Eb/N0 " << h << ", interval " << 0.5 << ", sig " << sig <<   endl;
        } else{
            Params::set_K(startK + i * interval);
        }

        //実質TPS設定
        Preseter::set_params(cap_map, A, Ac, p_0, p);

        rate = Common::get_rate();
        Common::bar();

        if(Params::get_is_outlog()) {
            cout << "T";
            int size = log2(Params::get_N())+2;
            int n = Params::get_N();
            vector<vector<bool>> T(2*log2(Params::get_N())+2, vector<bool>(n, false));
            Params::get_T(T);
            for (int j = 0; j < T.size(); ++j) {
                Common::pp(T[j]);
            }
            cout << "A";
            Common::pp(A);
            cout << "Ac";
            Common::pp(Ac);
            cout << "punc pattern";
            Common::pp(p);
            vector<int> table = Preseter::makeTable(Params::get_N());
            cout << "dBN";
            Common::dpp(table);
            cout << "dRBN";
            Common::drpp(table);
            cout << "BN";
            Common::pp(table);
            cout << "RBN";
            Common::rpp(table);
        }

        //decode
        while (block_error_count < Params::get_upperBlockErrorNum()) {
            //encode
            Preseter::preset_u(RAND, u);
            if( Common::is_mid_send() ){
                tmp_x[log2(Params::get_N()) - 1] = encoder.encode_m(Params::get_N(), 0, 0, u, tmp_x);
                //uとxをまとめて送る
                for (int i = 0; i < log2(Params::get_N()) + 1; i++) {
                    for (int j = 0; j < Params::get_N(); j++) {
                        if (i == 0) {
                            xm[i][j] = u[j];
                        } else {
                            xm[i][j] = tmp_x[i - 1][j];
                        }
                    }
                }
            } else {
                x = encoder.encode(Params::get_N(), u);
            }
            loopi++;
            u_est.assign(Params::get_N(), 0);

            if( Common::is_mid_send() ){
                Channel::channel_output_m(xm, ym);
            } else {
                y = Channel::channel_output(x);
            }

            if(loopi % bloop == 0 && Params::get_is_calc_bloop()){
                //1. error_countを一次元化
                int size = log2(Params::get_N())+1;
                int tempSize = size*Params::get_N();
                vector<int> temp_er;

                //node_error_count.size()-2だと右ノードは見ない
                for (int i = 0; i < node_error_count.size()-2; i=i+2) {
                    for (int j = 0; j < node_error_count[0].size(); j++) {
                        temp_er.push_back(node_error_count[i][j]);
                    }
                }
                //2. ソートする　[] = 順位
                vector<pair<int, int> > error_map;
                for(int i=0; i < temp_er.size(); i++){
                    error_map.push_back(pair<int, int>(i, temp_er[i]));
                }
                sort(begin(error_map), end(error_map), Common::sort_greater);
                //    Common::pp(tempZn);

                int N = Params::get_N();
                int sorted_i = 0;
                int count = 1;
                vector<int> temp_i;
                for (int i = 0; i < error_map.size(); i++) {
                    sorted_i = error_map[i].first;
                    if( !Common::containVal(sorted_i, saveb) ) {
                        temp_i.push_back(sorted_i);
                        if(count == Params::get_MN()) break;
                        count++;
                    }
                }

                //3. 一次元データを多次元に変換し, 上から送信フラグつけていく [][] = 順位
                //[i] => [i/N][i%N]として変換可能なはず
                sorted_i = temp_i[0];
                saveb.push_back(sorted_i);
                //Bに高いindexを入れていく
                B[(sorted_i/N)*2][sorted_i%N] = 1;
                vector<vector<int> > temp(2*log2(Params::get_N())+2, vector<int>(Params::get_N(),0));
                node_error_count = temp;
                Analysor::printDecodeProgress(loopi/bloop, B, b_file);
            }

            u_est = decoder.calcBP(loopi, param, u, x, y, xm, ym, node_error_count, val_error_file, B);
            itr = param[0];
            Analysor::errorCount(u, u_est, &error_count);
            if (error_count > 0) {
                block_error_count++;
            }
            if (loopi % 1 == 0 ) cout << loopi << " " << error_count << " " << block_error_count  << " " << itr << " " << (double)block_error_count/loopi << endl;

            error_count = 0;
            mitr += itr; // 0:itr, 1:nochecked

            if(loopi >= Params::get_blockNum()) break;
//            break;
        }
        val_error_file << "\t}" << endl;
        val_error_file << "}" << endl;
        b_file << "\t}" << endl;
        b_file << "}" << endl;

        mitr = (mitr/loopi) >= 70 ? 70 : mitr/loopi;
        itrfile << rate << " " << mitr << endl;

        if(Params::get_is_exp_awgn()) {
            cout << "mitr : " << mitr << endl;
            cout << "Eb/N0 : " << h << endl;
        } else {
            cout << "mitr : " << mitr << endl;
            cout << "rate : " << rate << endl;
        }

        performance.stopTimer();

        cout << Params::get_rvbDir() << endl;
        BER = (double)block_error_count/loopi;

        logger.outLog("=================================");
        logger.outLog(performance.outTime("処理時間"));
        logger.outLog("(N,K,M,MN) = (" + to_string(Params::get_N()) + "," + to_string(Params::get_K())  + "," + to_string(Params::get_M())  + "," + to_string(Params::get_MN())+ ")");
        logger.outLog("BER:" + to_string(BER));
        if(Params::get_is_exp_awgn()){
            logger.outLog("Eb/N0:" + to_string(h));
            logger.outLogRVB(h, BER);
        } else {
            logger.outLog("Rate:" + to_string(rate));
            logger.outLogRVB(rate, BER);
        }
        logger.outLog(encoder.outCount("encoder_count"));
        logger.outLog(decoder.outCount("decoder_count"));
        performance.outHMS();
//        break;
    }
    for (int i = 0; i < saveb.size(); i++) {
        saveb_file << saveb[i];
        if(i != saveb.size()-1 ){
            saveb_file << endl;
        }
    }
}