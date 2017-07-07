#include "../lib/Performance.h"

Performance::Performance(){

}

Performance::~Performance(){

}

void Performance::addCount(){
    this->count++;
}

void Performance::startTimer(){
    this->startTime = chrono::system_clock::now();
}

void Performance::stopTimer(){
    this->timeSpan = std::chrono::duration_cast<std::chrono::milliseconds>(chrono::system_clock::now() - this->startTime).count();
    this->sumTime += this->timeSpan;
}

void Performance::outHMS (){
    long int in = this->timeSpan/1000;
    const long int h { in / 3600 };
    const long int m { (in / 60) % 60 };
    const long int s { in % 60 };
    printf("%02ld:%02ld:%02ld\n",h,m,s);
}

string Performance::outCount(string content){
    string str = content + "::" + to_string(this->count) + "回";
    return str;
}

string Performance::outTime(string content){
    string str = content + "::" + to_string(this->timeSpan) + "[ms]";
    return str;
}

string Performance::outSumTime(string content){
    string str = content + "::" + to_string(this->sumTime) + "[ms]";
    return str;
}
