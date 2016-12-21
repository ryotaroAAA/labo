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

void Performance::outCount(string content){
    cout << content << "::" << this->count << "回" << endl;
}

void Performance::outTime(string content){
    cout << content << "::" << this->timeSpan << "[ms]" << endl;
}

void Performance::outSumTime(string content){
    cout << content << "::" << this->sumTime << "[ms]" << endl;
}
