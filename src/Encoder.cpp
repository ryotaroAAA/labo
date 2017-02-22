#include "../lib/Encoder.h"
Encoder::Encoder(){

}

Encoder::~Encoder(){

}

vector<int> Encoder::encode(int n, vector<int> &input){
    vector<int> x_n(n,0);
    vector<int> s_n(n,0);
    vector<int> v_n(n,0);
//    int s_n[n];
//    int v_n[n];

    for (int i = 0; i < n/2 ; i++) {
        s_n[2*i]   = (input[2*i+1] + input[2*i]) % 2;
        s_n[2*i+1] = input[2*i+1];
    }

    for (int i = 0; i < n/2 ; i++) {
        v_n[i] = s_n[2*i];
        v_n[n/2 + i] = s_n[2*i + 1];
    }

    if (n == 2) {
        x_n[0] = (input[0] + input[1]) % 2;
        x_n[1] = input[1];
    } else {
        vector<int> tempV_n1(0);
        vector<int> tempV_n2(0);

        for (int i = 0; i < n ; i++) {
            (i < n/2) ? tempV_n1.push_back(v_n[i]) : tempV_n2.push_back(v_n[i]);
        }

        vector<int> tempX_n1 = encode(n/2, tempV_n1);
        vector<int> tempX_n2 = encode(n/2, tempV_n2);

        for (int i = 0; i < n ; i++) {
            x_n[i] = (i < n/2) ? tempX_n1[i] : tempX_n2[i - n/2];
        }
    }
    this->addCount();
    return x_n;
}