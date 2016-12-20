#include "../lib/Analysor.h"

double Analysor::errorRate(vector<int> &u, vector<int> &u_est){
    int error_count = 0;
    for(int i=0; i<N; i++){
        if(u[i] != u_est[i]){
            error_count++;
        }
    }
    return (double)error_count/N;
}
double Analysor::errorRate(vector<int> &u, vector<int> &u_est){
    int error_count = 0;
    for(int i=0; i<N; i++){
        if(u[i] != u_est[i]){
            error_count++;
        }
    }
    return (double)error_count/N;
}


double Analysor::calcCapacityForBec(int i, int n) {
    double cap =0.0;
    if ( i == 0 && n == 1 ) {
        cap = 1 - e;
    } else {
        if ( i % 2 == 0) {
            cap = pow(calcCapacityForBec(i/2, n/2),2);
        } else {
            double tempCap = calcCapacityForBec((i-1)/2, n/2);
            cap = 2 * tempCap - pow(tempCap,2);
        }
    }
    //cout << "i:" << i << ", n:" << n <<endl;
//    PRINT(cap);
    return cap;
}

double Analysor::calcBhatForBec(int i, int n){
    double bha =0.0;
    if ( i == 0 && n == 1 ) {
        bha = e;
    } else {
        if ( i % 2 == 0) {
            double tempBha = calcBhatForBec(i/2, n/2);
            bha = 2 * tempBha - pow(tempBha,2);
        } else {
            bha = pow(calcBhatForBec((i-1)/2, n/2),2);
        }
    }
//    cout << "i:" << i << ", n:" << n <<endl;
//    PRINT(bha);
    return bha;
}

void Analysor::makeArrayCapacityForBec(vector<double> &array) {
    for(int i = 0; i < N; i++) {
        array.push_back(calcCapacityForBec(i, N));
    }
}

void Analysor::probErrBound(vector<double> &array) {
    int N = pow(2,18);
    int count = 0;
    double tempArr[N];
    double sumArr[N];
    double sum = 0.0;

    for(int i = 0; i < N; i++) {
        tempArr[i] = calcBhatForBec(i, N);
    }

    qsort(tempArr, N, sizeof(double), compare_asc);

    for(int i=0; i < N; i++){
        sum += tempArr[i];
        sumArr[i] = sum;
    }

    string filename = "/Users/ryotaro/labo/e18";
    ofstream w_file;
    w_file.open(filename, ios::out);
    for (int i = 0; i < N; i++)
    {
        w_file << (double)i/N << " " << sumArr[i] << " " << tempArr[i] << endl;
    }
}

void Analysor::calcBlockErrorRate(string mode, int n, int a){
    vector<int> u_Ac(0);
    vector<int> u_A(0);
    vector<int> A(0);
    vector<int> x_n;
    vector<int> y_n;
    vector<int> u_n_est;

    for (int k = 0; k < n; k += a) {
        hoge =0;
        hoge2 =0;
        if(mode == "eval"){
            N = n;
            K = k;
        } else {
            n=N;
            k=K;
        }
        int i = 1;
        vector<int> u_Ac(0);
        vector<int> u_A(0);
        vector<int> A(0);

//        double temp[n] = {0.0};
//        probErrBound(temp);
        defineFixedAndFree(n, u_Ac, u_A);
        vector<int> u_n(n);

        //処理時間計測//
        const auto startTime = chrono::system_clock::now();
        u_n = generateUi(2, u_n, u_Ac, A);
        x_n = encoder(n,u_n);
        y_n = channel(x_n);
        u_n_est = decoder(y_n, u_n, u_Ac, u_A);

        string filename = "/Users/ryotaro/labo/log";
        ofstream log;
        log.open(filename, ios::app);

        string rate_vs_error = "/Users/ryotaro/labo/log_blockerr_vs_rate";
        ofstream rve;
        rve.open(rate_vs_error, ios::app);

        cout << "error　probability:" << errorRate(u_n,u_n_est) << endl;
        cout << "rate:" << (double)k/n << endl;
//        log << "==================================================" << endl;
//        log << "(N,K) = (" << n << "," << k << ")" << endl;
//        log << "error　probability:" << errorRate(u_n,u_n_est) << endl;
//        log << "rate:" << (double)k/n << endl;
//        rve << (double)k/n << " " << errorRate(u_n,u_n_est) << endl;
        //処理時間計測//
        const auto endTime = chrono::system_clock::now();
        const auto timeSpan = endTime - startTime;
        cout << "総LR計算時間:" << hogetime << "[ms]" << endl;
        cout << "処理時間:" << chrono::duration_cast<chrono::milliseconds>(timeSpan).count() << "[ms]" << endl;
//        log << "総LR計算時間:" << hogetime << "[ms]" << endl;
//        log << "処理時間:" << chrono::duration_cast<chrono::milliseconds>(timeSpan).count() << "[ms]" << endl;
        const auto astartTime = chrono::system_clock::now();
        cout << hoge << endl;
        cout << hoge2 << endl;
//        log << hoge << endl;
//        log << hoge2 << endl;
//        log << "==================================================" << endl;

        if(mode != "eval") break;
    }
}