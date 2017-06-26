#ifndef CHANNEL_POLARIZATION_DECODER_H
#define CHANNEL_POLARIZATION_DECODER_H

#include "Common.h"
#include "Channel.h"
#include "Performance.h"

class Decoder : public Performance{
public:
    explicit Decoder();
    ~Decoder();
    vector<int> makeTreeIndex(int n);
    vector<int> decode(vector<int> &y, vector<int> &u, vector<int> &x, vector<int> &A);
    double calcL_i(int i, int n ,int cache_i, int lefvel ,vector<int> &y ,vector<int> &u, vector<vector<bool> > &isCache , vector<vector<double> > &cache);
};

#endif //CHANNEL_POLARIZATION_DECDER_H
