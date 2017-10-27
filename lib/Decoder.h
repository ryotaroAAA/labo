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

    void calc_marge(vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked);
    void calc_marge_m(vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &y, vector<vector<bool> > &ym_isReceived);
    void calc_check_to_val(int size, vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked);
    void calc_check_to_val_m(vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &y, vector<vector<bool> > &ym_isReceived);
    void calc_val_to_check(int size, vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked);
    void calc_val_to_check_m(vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &y, vector<vector<bool> > &ym_isReceived);
    void send_message(int i, int j, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked);
    void send_message_m(int i, int j, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &y, vector<vector<bool> > &ym_isReceived);

    void BPinit(vector<double> &y, vector<int> &u, vector<int> &A, vector<vector<double> > &node_val, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked);
    void SCinit(int n, vector<double> &y, vector<int> &u, vector<int> &u_est, vector<int> &A, vector<vector<double> > &node_val, vector<vector<vector<message> > > &message_list, vector<vector<vector<message> > > &save_list, vector<vector<bool> > &node_isChecked, vector<vector<bool> > &save_isChecked);
    void BPinit_m(vector<vector<double> > &y, vector<int> &u, vector<int> &A, vector<vector<double> > &node_val, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<bool> > &ym_isReceived);

    bool isTerminate(int &no_checked, vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked);
    void confirmIsCheck(vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked);
    bool isChanged(vector<vector<bool> > &old_node_isChecked, vector<vector<bool> > &new_node_isChecked2);

    void printDecodeProgress(int count, vector<vector<bool> > &node_value, ofstream &w_file);
    void printDecodeProgress(int count, vector<vector<double> > &node_value, ofstream &w_file);

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
    vector<int> calcBP(vector<int> &param, vector<double> &y, vector<int> &u, vector<int> &A);
    vector<int> calcBP_m(vector<int> &param, vector<vector<double> > &y, vector<int> &u, vector<int> &A);

    void calc_mp(int size, vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked);
};

#endif //CHANNEL_POLARIZATION_DECDER_H
