#include "../lib/Logger.h"
Logger::Logger(){
    this->log.open(this->dir, ios::app);
    this->rvb.open(this->rvb_dir, ios::app);
}

Logger::~Logger(){

}

void Logger::setDir(string dir){
    this->dir = dir;
}

void Logger::outLog(string content){
    this->log << content << endl;
}

void Logger::outLogRVB(int rate, int ber){
    this->rvb << rate << " " << ber << endl;
}

void Logger::outLogCount(string content, string title) {
    this->log << title << "::" << content << "回" << endl;
}

void Logger::outLogTime(string content, string title){
    this->log << title << "::" << content << "[ms]" << endl;
}