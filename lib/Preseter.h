#ifndef CHANNEL_POLARIZATION_PRESETER_H
#define CHANNEL_POLARIZATION_PRESETER_H
#include "Common.h"
#include "Analysor.h"
class Preseter {
private:
public:
    explicit Preseter();
    ~Preseter();
    static void preset(vector<int> &u, vector<int> &fixed, vector<int> &free);
    static void defineFixedAndFree(vector<int> &fixed, vector<int> &free);
    static vector<int> generateUi(SOURCE_TYPE set, vector<int> &x, vector<int> &u_Ac, vector<int> &A);
};
#endif //CHANNEL_POLARIZATION_PRESETER_H
