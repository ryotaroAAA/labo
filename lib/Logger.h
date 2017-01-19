#ifndef CHANNEL_POLARIZATION_LOGGER_H
#define CHANNEL_POLARIZATION_LOGGER_H
#include "Common.h"
class Logger {
private:
    string dir = "log/all_log";
    string rvb_dir = "log/";
    ofstream log;
    ofstream rvb;
public:
    explicit Logger();
    ~Logger();
    void setDir(string dir);
    void outLog(string content);
    void outLogRVB(int rate, int ber);
    void outLogTime(string content);
    void outLogCount(string content);
};

#endif //CHANNEL_POLARIZATION_LOGGER_H