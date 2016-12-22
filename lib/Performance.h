#ifndef CHANNEL_POLARIZATION_PERFORMANCE_H
#define CHANNEL_POLARIZATION_PERFORMANCE_H

#include "Common.h"
class Performance{
private:
    int count = 0;
    double sumTime;
    double timeSpan;
    std::chrono::system_clock::time_point startTime;
public:
    explicit Performance();
    ~Performance();
    void addCount();
    void startTimer();
    void stopTimer();
    string outTime(string content = "実行回数");
    string outSumTime(string content = "処理時間");
    string outCount(string content = "処理時間");
};

#endif //CHANNEL_POLARIZATION_PERFORMANCE_H
