#include "Common.h"

class Encoder{
public:
    explicit Encoder();
    ~Encoder();
    vector<int> encode(int n, vector<int> &input); 
};
