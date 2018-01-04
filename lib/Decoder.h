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

    //main BP
    vector<int> calcBP(int loop_num, vector<int> &param, vector<int> &u, vector<int> &x, vector<double> &y, vector<vector<int> > &xm, vector<vector<double> > &ym, vector<vector<int> > &node_error_count, ofstream &val_error_file, vector<vector<int> > &B);
    void calcSConBP(int itr, int count, ofstream &val_file, ofstream &check_file, vector<int> &u, vector<double> &y , vector<int> &u_n_est, vector<double> &tmp_u, vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked);

    //BP init
    void init_message(vector<bool> &puncFlag, vector<vector<double> > &node_val, vector<vector<vector<message> > > &message_list);
    void init_params(vector<bool> &puncFlag, vector<int> &u, vector<int> &x, vector<double> &y, vector<vector<int> > &xm, vector<vector<double> > &ym, vector<vector<double> > &node_val, vector<vector<bool> > &node_isChecked, vector<vector<int> > &B);
    void SCinit(int n, vector<double> &y, vector<int> &u, vector<int> &u_est, vector<vector<double> > &node_val, vector<vector<vector<message> > > &message_list, vector<vector<vector<message> > > &save_list, vector<vector<bool> > &node_isChecked, vector<vector<bool> > &save_isChecked);
    void BPinit(vector<int> &u, vector<int> &x, vector<double> &y, vector<vector<int> > &xm, vector<vector<double> > &ym, vector<vector<double> > &node_val, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<int> > &B);

    //calc message
    double calc_message(int mode, vector<double> val);
    void calc_marge(vector<vector<double> > &node_value, vector<vector<vector<message> >> &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &y);
    void calc_check_to_val(int size, vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &y);
    void calc_val_to_check(int size, vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &y);
    void send_message(int i, int j, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &y);
    void calc_mp(int size, vector<vector<double> > &node_value, vector<vector<vector<message> > > &message_list, vector<vector<bool> > &node_isChecked, vector<vector<double> > &ym);

    //calc terminal condition
    bool isTerminate(int &no_checked, vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked);
    void confirmIsCheck(vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked);
    bool isChanged(vector<vector<bool> > &old_node_isChecked, vector<vector<bool> > &new_node_isChecked2);

    //log
    void init_outLog(ofstream &val_file, ofstream &check_file, ofstream &val_error_file);
    void outLog(int itr, int no_checked, vector<int> &u, vector<int> &u_est, ofstream &val_file, ofstream &check_file ,vector<vector<double> > &node_value, vector<vector<bool> > &node_isChecked);
    void printDecodeProgress(int count, vector<vector<bool> > &node_value, ofstream &w_file);
    void printDecodeProgress(int count, vector<vector<int> > &node_value, ofstream &w_file);
    void printDecodeProgress(int count, vector<vector<double> > &node_value, ofstream &w_file);

    //gear
    bool is_mid_send();
    vector<vector<int> > adjacentIndex(int level, int index);
    vector<int> makeTreeIndex(int n);
    vector<int> makeBPTreeIndex(int n);
    double take_val(vector<int> &locate, vector<vector<double> > &node_val);
//    void set_val(bool val, vector<int> &locate, vector<vector<bool> > &node_val);

    //sc
    vector<int> decode(vector<double> &y, vector<int> &u);
    double calcL_i(int i, int n ,int cache_i, int lefvel ,vector<double> &y ,vector<int> &u, vector<vector<bool> > &isCache , vector<vector<double> > &cache);
};

#endif //CHANNEL_POLARIZATION_DECDER_H
