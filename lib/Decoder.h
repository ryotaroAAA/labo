#ifndef CHANNEL_POLARIZATION_DECODER_H
#define CHANNEL_POLARIZATION_DECODER_H

#include "Common.h"
#include "Channel.h"
#include "Performance.h"

class Decoder : public Performance{
public:
    explicit Decoder();
    ~Decoder();
    double calc_node(int mode, double a, double b);
    void BPinit(vector<double> &y, vector<int> &u, vector<int> &A, vector<vector<double> > &node_val, vector<vector<int> > &update_count, vector<vector<bool> > &node_isChecked);
    bool isChecked(int level, int index, vector<vector<int> > &adjacent, vector<vector<double> > &node_val, vector<vector<bool> > &node_isChecked);
    bool isTerminate(vector<vector<double> > &node_val, vector<vector<bool> > &node_isChecked);
    bool isCalculable(int level, int index, vector<vector<int> > &adjacent, vector<vector<double> > &node_val, vector<vector<bool> > &node_isChecked);
    int take_val(vector<int> index, vector<vector<int> > &node_val);
    double take_val(vector<int> index, vector<vector<double> > &node_val);
    bool take_val(vector<int> index, vector<vector<bool> > &node_val);
    vector<vector<int> > adjacentIndex(int level, int index);
    vector<int> makeTreeIndex(int n);
    vector<int> makeBPTreeIndex(int n);
    vector<int> decode(vector<double> &y, vector<int> &u, vector<int> &A);
    vector<int> decode_m(vector<vector<double> > &y, vector<int> &u, vector<int> &A);
    double calcL_i(int i, int n ,int cache_i, int lefvel ,vector<double> &y ,vector<int> &u, vector<vector<bool> > &isCache , vector<vector<double> > &cache);
    vector<int> BP(int limit, vector<double> &y, vector<int> &u, vector<int> &A);
};

#endif //CHANNEL_POLARIZATION_DECDER_H
