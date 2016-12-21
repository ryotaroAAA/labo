#ifndef CHANNEL_POLARIZATION_ANALYSOR_H
#define CHANNEL_POLARIZATION_ANALYSOR_H
#include "Common.h"
class Analysor {
public:
    explicit Analysor();
    ~Analysor();
    double calcBhatForBec(int i, int n);
    void probErrBound(vector<double> &array);

    static double calcCapacityForBec(int i, int n);
    static double errorRate(vector<int> &u, vector<int> &u_est);
    static void makeArrayCapacityForBec(vector<double> &array);
    static void calcBlockErrorRate(string mode, int n, int a);
    static void defineFixedAndFree(int n, vector<int> &fixed, vector<int> &free);
};

#endif //CHANNEL_POLARIZATION_ANALYSOR_H
