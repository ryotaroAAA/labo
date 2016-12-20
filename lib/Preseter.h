#ifndef CHANNEL_POLARIZATION_PRESETER_H
#define CHANNEL_POLARIZATION_PRESETER_H
#include "Common.h"

class Preseter {
private:
public:
    explicit Preseter::Preseter(int n, vector<int> &u, vector<int> &fixed, vector<int> &free);
    ~Preseter();
    static void defineFixedAndFree(int n, vector<int> &fixed, vector<int> &free);
    static vector<int> generateUi(SOURCE_TYPE set, vector<int> &x, vector<int> &u_Ac, vector<int> &A);
};



#endif //CHANNEL_POLARIZATION_PRESETER_H
