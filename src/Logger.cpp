#include "../lib/Logger.h"
Logger::Logger(){
    time_t now = time(NULL);
    struct tm *pnow = localtime(&now);

    struct stat st;
    string dir = "/Users/ryotaro/labo/log/" + to_string(pnow->tm_year+1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday);
    int ret = stat(dir.data(), &st);
    if (ret == -1) {
        mkdir(dir.data(), S_IRUSR | S_IWUSR | S_IXUSR);
        cout <<dir.data()<<endl;
    }
    cout << ret <<endl;
    this->log.open(this->dir, ios::app);
    this->rvb.open(Params::get_rvbDir(), ios::app);
}

Logger::~Logger(){

}

void Logger::setDir(string dir){
    this->dir = dir;
}

void Logger::setRvbDir(string rvb_dir){
    this->rvb_dir = rvb_dir;
}

void Logger::outLog(string content){
    this->log << content << endl;
    cout << content << endl;
}

void Logger::outLogRVB(double rate, double ber){
    this->rvb << rate << " " << ber << endl;
}
