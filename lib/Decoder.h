#include "Common.h"

class Decoder{
public:
    explicit Decoder();
    ~Decoder();
    vector<int> decode(vector<int> &y, vector<int> &u, vector<int> &Ac, vector<int> &A);
    double calcL_i(int i, int n ,int cache_i,int level ,vector<int> &y ,vector<int> &u, int u_i_est, vector<vector<bool> > &isCache , vector<vector<double> > &cache); 
};
