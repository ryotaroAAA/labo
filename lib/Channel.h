#include "Common.h"
class Channel{
public:
    explicit Channel(double error_prob);
    ~Channel();
    double calcW(int y, int x);
    double calcW_i(int i, int n, vector<int> &u, vector<int> &u_i, vector<int> &y);
    vector<int> output(vector<int> &input);
private:
    double error_prob;
};
