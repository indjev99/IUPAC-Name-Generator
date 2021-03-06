#ifndef ATOM_H_INCLUDED
#define ATOM_H_INCLUDED
#include<vector>
#include<stack>
#include<string>
#include "bond.h"
using namespace std;
struct atom
{
    string symbol;
    double x,y;
    vector<bond> bonds;
    stack<int> free_bonds;
    atom();
    atom(string new_symbol, int new_valance, double new_x, double new_y);
    atom(string new_symbol, int new_valance, double new_x, double new_y, int new_bond);
    int connect(int bond);
    bool canConnect(int bond);
    int isConnected(int bond);
    void removeBond(int bond);
    void changeXY(double new_x, double new_y);
};
#endif //ATOM_H_INCLUDED
