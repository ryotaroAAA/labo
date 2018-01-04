#include "../lib/Params.h"
#include <cmath>
#include <iostream>
//int size = 2*log2(Params::get_N())+2;
//int n = Params::get_N();
int Params::N = 0;
int Params::M = 0;
int Params::K = 0;
int Params::MN = 0;
int Params::Bloop = 0;
int Params::monteNum = 0;
int Params::blockNum = 0;
int Params::upperBlockErrorNum = 0;
int Params::rp = 0;
//vector<vector<bool> > Params::T(ysize, vector<bool>(n, false));
vector<vector<bool> > Params::T;
vector<int> Params::A;
vector<int> Params::Ac;
vector<int> Params::p;
double Params::point[] = {1,1,1};
bool Params::is_outlog = false;
bool Params::is_calc_bloop = false;
double Params::e = 0;
string Params::rvbDir = "";
CHANNEL_TYPE Params::s = BEC;
DECODE_MODE Params::dm = BP;
EXP_MODE Params::exp_mode = NORMAL;
MID_MODE Params::m_mode = MID_DOR;

int Params::get_N(){
    return N;
}

int Params::get_M(){
    return M;
}

int Params::get_K(){
    return K;
}

int Params::get_MN(){
    return MN;
}

int Params::get_Bloop(){
    return Bloop;
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

void Params::get_T(vector<vector<bool>> &temp){
//    for (int i = 0; i < temp.size(); i++) {
//        for (int j = 0; j < temp[0].size(); j++) {
//            temp[i][j] = T[i][j];
////            std::cout << "a";
//        }
//    }
    temp = T;
}

void Params::get_A(vector<int> &temp){
    temp = A;
}

void Params::get_Ac(vector<int> &temp){
    temp = Ac;
}

void Params::get_p(vector<int> &temp){
    temp = p;
}

void Params::get_point(double temp[3]){
    temp[0] = point[0];
    temp[1] = point[1];
    temp[2] = point[2];
}

int Params::get_rp(){
    return rp;
}

bool Params::get_is_outlog(){
    return is_outlog;
}

bool Params::get_is_calc_bloop(){
    return is_calc_bloop;
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

MID_MODE Params::get_m_mode(){
    return m_mode;
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

void Params::set_MN(int _MN){
    MN = _MN;
}

void Params::set_point(double _point[3]){
    point[0] = _point[0];
    point[1] = _point[1];
    point[2] = _point[2];
}

void Params::set_T(vector<vector<bool> > &temp){
    T=temp;
//    for (int i = 0; i < temp.size(); i++) {
//        T.push_back(temp[i]);
//    }
//    for (int i = 0; i < temp.size(); i++) {
//        for (int j = 0; j < temp[0].size(); j++) {
//            T[i][j] = temp[i][j];
//        }
//    }

}

void Params::set_A(vector<int> &temp){
    A = temp;
}

void Params::set_Ac(vector<int> &temp){
    Ac = temp;
}

void Params::set_p(vector<int> &temp){
    p = temp;
}

void Params::set_Bloop(int _Bloop){
    Bloop = _Bloop;
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

void Params::set_m_mode(MID_MODE _m_mode){
    m_mode = _m_mode;
}

void Params::set_is_outlog(bool _is_outlog){
    is_outlog = _is_outlog;
}

void Params::set_is_calc_bloop(bool _is_calc_bloop){
    is_calc_bloop = _is_calc_bloop;
}

