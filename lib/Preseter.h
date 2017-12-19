#ifndef CHANNEL_POLARIZATION_PRESETER_H
#define CHANNEL_POLARIZATION_PRESETER_H
#include "Common.h"
#include "Analysor.h"
#include "Params.h"
#include "MT.h"

class Preseter {
private:
public:
    explicit Preseter();
    ~Preseter();
    static void preset_A_Ac(vector<int> &free, vector<int> &fixed);
    static void preset_u(SOURCE_TYPE mode, vector<int> &u);
    static void defineFixedAndFree(vector<int> &free, vector<int> &fixed);
    static void represet_A(vector<int> &free, vector<int> &fixed_0, vector<pair<int,double> > &cap_map);
    static void makeMutualInfoArray(vector<pair<int, double> > &cap_map);
    static vector<int> generateUi(SOURCE_TYPE set, vector<int> &x);
    static void set_params(vector<pair<int, double> > &cap_map,vector<int> &A, vector<int> &Ac, vector<int> &p_0, vector<int> &p);
    static vector<int> makeTable(int n);
    static void makeManyValTableAs(bool sortflg, vector<int> &table_map);
    static vector<int> get_bitReversal(vector<int> p_0);
    static vector<int> makeTreeIndex(int n);
    static void set_zin_all(vector<int> &A, vector<vector<int> > &tempZn, vector<vector<bool> > &ym_isReceived);

};
#endif //CHANNEL_POLARIZATION_PRESETER_H
