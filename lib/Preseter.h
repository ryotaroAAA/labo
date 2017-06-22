#ifndef CHANNEL_POLARIZATION_PRESETER_H
#define CHANNEL_POLARIZATION_PRESETER_H
#include "Common.h"
#include "Analysor.h"
class Preseter {
private:
public:
    explicit Preseter();
    ~Preseter();
    static void preset(SOURCE_TYPE mode, vector<int> &u, vector<int> &free);
    static void defineFixedAndFree(vector<int> &free);
    static vector<int> generateUi(SOURCE_TYPE set, vector<int> &x);
};
#endif //CHANNEL_POLARIZATION_PRESETER_H
