#ifndef CHANNEL_POLARIZATION_DECODER_H
#define CHANNEL_POLARIZATION_DECODER_H

#include "Common.h"
#include "Channel.h"
#include "Performance.h"

typedef struct {
    int toLevel;
    int toIndex;
    double val;
} message;

class Decoder : public Performance{
public:
    explicit Decoder();
    ~Decoder();
    double calc_message(int mode, vector<double> val);
    void calc_node_val(vector<vector<double>> &node_value, vector<vector<vector<message>>> &message_list, vector<vector<bool>> &node_isChecked);
    void calc_BP(int count, vector<vector<double>> &node_value, vector<vector<vector<message>>> &message_list, vector<vector<bool> > &node_isChecked);
    void BPinit(vector<double> &y, vector<int> &u, vector<int> &A, vector<vector<double> > &node_val, vector<vector<vector<message>>> &message_list, vector<vector<bool> > &node_isChecked);
    bool isTerminate(vector<vector<double> > &node_value, vector<vector<vector<message>>> &message_list, vector<vector<bool> > &node_isChecked);
    int take_val(vector<int> index, vector<vector<int> > &node_val);
    double take_val(vector<int> index, vector<vector<double> > &node_val);
    bool take_val(vector<int> index, vector<vector<bool> > &node_val);
    void set_val(bool val, vector<int> locate, vector<vector<bool> > &node_val);
    vector<vector<int> > adjacentIndex(int level, int index);
    vector<int> makeTreeIndex(int n);
    vector<int> makeBPTreeIndex(int n);
    vector<int> decode(vector<double> &y, vector<int> &u, vector<int> &A);
    vector<int> decode_m(vector<vector<double> > &y, vector<int> &u, vector<int> &A);
    double calcL_i(int i, int n ,int cache_i, int lefvel ,vector<double> &y ,vector<int> &u, vector<vector<bool> > &isCache , vector<vector<double> > &cache);
    vector<int> BP(int limit, vector<double> &y, vector<int> &u, vector<int> &A);
};

#endif //CHANNEL_POLARIZATION_DECDER_H
