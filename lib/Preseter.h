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
    static void represet_A(vector<int> &free, vector<pair<int,double> > &cap_map);
    static void represet_A_wang(vector<int> &free, vector<int> &fixed_0, vector<pair<int,double> > &cap_map);
    static void makeMutualInfoArray(vector<pair<int, double> > &cap_map);
    static vector<int> generateUi(SOURCE_TYPE set, vector<int> &x);
};
#endif //CHANNEL_POLARIZATION_PRESETER_H
