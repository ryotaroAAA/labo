#include "../lib/Params.h"

int Params::N = 0;
int Params::K = 0;
double Params::e = 0;

int Params::get_N(){
    return N;
}

int Params::get_K(){
    return K;
}

double Params::get_e(){
    return e;
}

void Params::set_N(int _N){
    N = _N;
}

void Params::set_K(int _K){
    K = _K;
}

void Params::set_e(double _e){
    e = _e;
}

