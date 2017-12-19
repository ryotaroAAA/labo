#ifndef CHANNEL_POLARIZATION_ENCODER_H
#define CHANNEL_POLARIZATION_ENCODER_H

#include "Common.h"
#include "Performance.h"

class Encoder : public Performance{
public:
    explicit Encoder();
    ~Encoder();

    vector<int> encode(int n, vector<int> &input);
    vector<int> enc(int n, vector<int> &input);
    vector<int> encode_m(int n, int m, int level, vector<int>  &input, vector<vector<int> > &x);
};

#endif //CHANNEL_POLARIZATION_ENCODER_H