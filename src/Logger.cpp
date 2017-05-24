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
    cout << content << endl;
}

void Logger::outLogRVB(double rate, double ber){
    this->rvb << rate << " " << ber << endl;
}
