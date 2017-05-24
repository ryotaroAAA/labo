#ifndef CHANNEL_POLARIZATION_LOGGER_H
#define CHANNEL_POLARIZATION_LOGGER_H
#include "Common.h"
class Logger {
private:
    string dir = "/Users/ryotaro/labo/log/all_log";
    string rvb_dir = "/Users/ryotaro/labo/log/rve";
    ofstream log;
    ofstream rvb;
public:
    explicit Logger();
    ~Logger();
    void setDir(string dir);
    void outLog(string content);
    void outLogRVB(double rate, double ber);
};
#endif //CHANNEL_POLARIZATION_LOGGER_H