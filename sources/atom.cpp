#include "../headers/atom.h"
atom::atom() {}
atom::atom(string new_symbol, int new_valance, double new_x, double new_y)
{
    symbol=new_symbol;
    x=new_x;
    y=new_y;
    bonds.resize(new_valance);
    for (int i=new_valance-1; i>=0; --i)
    {
        bonds[i].to=-1;
        free_bonds.push(i);
    }
}
atom::atom(string new_symbol, int new_valance, double new_x, double new_y, int new_bond)
{
    symbol=new_symbol;
    x=new_x;
    y=new_y;
    bonds.resize(new_valance);
    for (int i=new_valance-1; i>=0; --i)
    {
        bonds[i].to=-1;
        free_bonds.push(i);
    }
    connect(new_bond);
}
int atom::connect(int bond)
{
    if (!free_bonds.empty())
    {
        int nc=isConnected(bond);
        if (nc==-1)
        {
            nc=free_bonds.top();
        }
        else if (symbol=="C" && bonds[nc].spots_taken.size()==bonds.size()-1) return -1;
        bonds[nc].to=bond;
        bonds[nc].spots_taken.push_back(free_bonds.top());
        free_bonds.pop();
        return 1;
    }
    return -1;
}
bool atom::canConnect(int bond)
{
    if (!free_bonds.empty())
    {
        int nc=isConnected(bond);
        if (symbol=="C" && nc>=0 && bonds[nc].spots_taken.size()==bonds.size()-1) return 0;
        return 1;
    }
    return 0;
}
int atom::isConnected(int bond)
{
    for (int i=0; i<bonds.size(); ++i)
    {
        if (bonds[i].to==bond) return i;
    }
    return -1;
}
void atom::removeBond(int bond)
{
    for (int i=0; i<bonds.size(); ++i)
    {
        if (bonds[i].to==bond)
        {
            bonds[i].to=-1;
            for (int j=0; j<bonds[i].spots_taken.size(); ++j)
            {
                free_bonds.push(bonds[i].spots_taken[j]);
            }
            bonds[i].spots_taken.resize(0);
        }
    }
}
void atom::changeXY(double new_x, double new_y)
{
    x=new_x;
    y=new_y;
}
