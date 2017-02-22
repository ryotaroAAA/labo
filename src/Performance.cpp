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
