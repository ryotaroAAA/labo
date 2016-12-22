#ifndef CHANNEL_POLARIZATION_ANALYSOR_H
#define CHANNEL_POLARIZATION_ANALYSOR_H
#include "Common.h"
#include "Preseter.h"
#include "Decoder.h"
#include "Encoder.h"
#include "Channel.h"
#include "Logger.h"
#include "Performance.h"

class Analysor {
private:
    int tempK;
    vector<int> u_Ac(Params::N);
    vector<int> u_A(Params::N;
    vector<int> A(Params::N);
    vector<int> u_n(Params::N);
    vector<int> x_n(Params::N);
    vector<int> y_n(Params::N);
    vector<int> u_est(Params::N);
    Performance performance;
    Decoder decoder;
    Encoder encoder;
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
