#include "../lib/Params.h"

int Params::N = 0;
int Params::K = 0;
int Params::monteNum = 0;
int Params::blockNum = 0;
int Params::upperBlockErrorNum = 0;
int Params::rp = 0;
double Params::e = 0;
string Params::rvbDir = "";
CHANNEL_TYPE Params::s = BEC;
DECODE_MODE Params::dm = BP;

int Params::get_N(){
    return N;
}

int Params::get_K(){
    return K;
}

int Params::get_monteNum(){
    return monteNum;
}

int Params::get_blockNum(){
    return blockNum;
}

int Params::get_upperBlockErrorNum(){
    return upperBlockErrorNum;
}

double Params::get_e(){
    return e;
}

int Params::get_rp(){
    return rp;
}

string Params::get_rvbDir() {
    return rvbDir;
}

CHANNEL_TYPE Params::get_s(){
    return s;
}

DECODE_MODE Params::get_decode_mode(){
    return dm;
}

void Params::set_N(int _N){
    N = _N;
}

void Params::set_K(int _K){
    K = _K;
}

void Params::set_monteNum(int _monteNum){
    monteNum = _monteNum;
}

void Params::set_blockNum(int _blockNum){
    blockNum = _blockNum;
}

void Params::set_upperBlockErrorNum(int _upperBlockErrorNum){
    upperBlockErrorNum = _upperBlockErrorNum;
}

void Params::set_e(double _e){
    e = _e;
}

void Params::set_rp(double _rp){
    rp = _rp;
}

void Params::set_rvbDir(string _rvbDir) {
    rvbDir = _rvbDir;
}

void Params::set_s(CHANNEL_TYPE _s){
    s = _s;
}

void Params::set_decode_mode(DECODE_MODE _dm){
    dm = _dm;
}

