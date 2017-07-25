#ifndef CHANNEL_POLARIZATION_ANALYSOR_H
#define CHANNEL_POLARIZATION_ANALYSOR_H
#include "Common.h"
#include "Preseter.h"
#include "Decoder.h"
#include "Encoder.h"
#include "Channel.h"
#include "Logger.h"
#include "Performance.h"
#include "Params.h"

class Analysor {
private:

public:
    explicit Analysor();
    ~Analysor();

    static double calcCapacity(int i, int n);
    static double calcBhatforBEC(int i, int n);

    static void errorCount(vector<int> &u, vector<int> &u_est, int* error_count);

    static void probErrBound();
    static double errorCalc(vector<int> &u, vector<int> &u_est, int* error_count);
    static void makeArrayCapacity(vector<double> &array);
    static void makeArrayBhat(vector<double> &array);
    static void calcBlockErrorRate(MODE mode);
    static void defineFixedAndFree(int n, vector<int> &fixed, vector<int> &free);
};

#endif //CHANNEL_POLARIZATION_ANALYSOR_H
