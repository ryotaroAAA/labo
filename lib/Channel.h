#ifndef CHANNEL_POLARIZATION_CHANNEL_H
#define CHANNEL_POLARIZATION_CHANNEL_H

#include "Common.h"
#include "MT.h"

class Channel{
public:
    explicit Channel();
    ~Channel();
    static double calcW(double y, int x);
    static double calcW_i(int i, int n, vector<int> &u, int u_i, vector<int> &y);
    static vector<double> channel_output(vector<int> &input);
};

#endif //CHANNEL_POLARIZATION_CHANNEL_H

