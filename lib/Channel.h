#ifndef CHANNEL_POLARIZATION_CHANNEL_H
#define CHANNEL_POLARIZATION_CHANNEL_H

#include "Common.h"
#include "MT.h"

class Channel{
public:
    explicit Channel();
    ~Channel();
    static double calcW(int y, int x);
    static double calcW_i(int i, int n, vector<int> &u, int u_i, vector<int> &y);
    static vector<int> channel_output(vector<int> &input);
};

#endif //CHANNEL_POLARIZATION_CHANNEL_H

