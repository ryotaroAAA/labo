#include "../lib/Params.h"

int Params::N = 0;
int Params::M = 0;
int Params::K = 0;
int Params::monteNum = 0;
int Params::blockNum = 0;
int Params::upperBlockErrorNum = 0;
int Params::rp = 0;
bool Params::is_outlog = false;
double Params::e = 0;
string Params::rvbDir = "";
CHANNEL_TYPE Params::s = BEC;
DECODE_MODE Params::dm = BP;
EXP_MODE Params::exp_mode = NORMAL;

int Params::get_N(){
    return N;
}

int Params::get_M(){
    return M;
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

bool Params::get_is_outlog(){
    return is_outlog;
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

EXP_MODE Params::get_exp_mode(){
    return exp_mode;
}

void Params::set_N(int _N){
    N = _N;
}

void Params::set_M(int _M){
    M = _M;
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

void Params::set_exp_mode(EXP_MODE _exp_mode){
    exp_mode = _exp_mode;
}

void Params::set_is_outlog(bool _is_outlog){
    is_outlog = _is_outlog;
}

