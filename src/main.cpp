#include "../lib/Common.h"
#include "../lib/Preseter.h"
#include "../lib/Decoder.h"
#include "../lib/Encoder.h"
#include "../lib/Channel.h"
#include "../lib/Analysor.h"
#include "../lib/Logger.h"

int main(void) {
    long int n = pow(2,8);

    vector<int> u_Ac(N,0);
    vector<int> u_A(N,0);
    vector<int> A(N,0);
    vector<int> u_n(N,0);

    Preseter *preset = new Preseter(N, u_n, u_Ac, u_A);

    for(auto val: u_n){
        cout << val << endl;
    }
//    cout << "======================" << endl;
//    vector<int> x_n = encoder(N, u_n);
//    for(auto val: x_n){
//        cout << val << endl;
//    }
//    cout << "======================" << endl;
//    vector<int> y_n = channel(x_n);
//    for(auto val: y_n){
//        cout << val << endl;
//    }
//    cout << "======================" << endl;


//    vector<int> u_n_est = decoder(y_n, u_n, u_Ac, u_A);
//    for(auto val: u_n_est){
//        cout << val << endl;
//    }
//    cout << "======================" << endl;
    return 0;
}