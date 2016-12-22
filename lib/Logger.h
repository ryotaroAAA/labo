#ifndef CHANNEL_POLARIZATION_LOGGER_H
#define CHANNEL_POLARIZATION_LOGGER_H
#include "Common.h"
class Logger {
private:
    string dir = "/Users/ryotaro/labo/log/all_log";
    string rvb_dir = "/Users/ryotaro/labo/log/";
    ofstream log;
    ofstream rvb;
public:
    explicit Logger();
    ~Logger();
    void setDir(string dir);
    void outLog(string content);
    void outLogTime(string content, string title = "実行回数");
    void outLogCount(string content, string title = "処理時間");
};

#endif //CHANNEL_POLARIZATION_LOGGER_H