#ifndef COMPOUND_H_INCLUDED
#define COMPOUND_H_INCLUDED
#include<vector>
#include<stack>
#include<string>
#include "dictionary.h"
#include "atom.h"
using namespace std;

extern const int halogen_N;
extern const string halogen_symbols[];
extern unordered_map<string, int> substituent_priorities;

struct compound
{
    vector<atom> atoms;
    stack<int> free_positions;
    string name;
    bool changed;
    dictionary curr_dict;

    vector<int> atoms_unions;
    vector<int> atom_in_union;
    vector<bool> unions_cycle;

    compound();
    compound(dictionary new_dict);
    compound(dictionary new_dict, string new_symbol, int new_valance, double new_x, double new_y);

    void setDictionary (dictionary new_dict);

    int addAtom(atom& a);
    int addAtom(string new_symbol, int new_valance, double new_x, double new_y);
    int addAtom(string new_symbol, int new_valance, double new_x, double new_y, int new_bond);
    int connectAtoms(int a1, int a2);
    int removeAtom(int a);
    int moveAtom(int a, double x, double y);
    int findAtom(double x, double y);

    static string addSuffix(string name, string suffix);
    static string onlyLetters(string a);
    static bool cmpBySubName(tuple<int, string, bool> a, tuple<int, string, bool> b);
    static string intToString(int num);
    static int cmpVectors(vector<int> a, vector<int> b);
    static vector<int> convertVector(vector<tuple<int, string, bool> > a);
    static vector<int> convertVector(vector<pair<int,int> > a);

    static bool isHalogen(string s);

    string findGN(vector<int> pos, string name, bool carbon, int parent_chain_length, bool cyclic, bool benzene, bool most_important); //find group's name
    string findGNs(vector<tuple<int, string, bool> > subs, int parent_chain_length, bool cyclic, bool benzene, bool most_important); //find groups' names

    bool isConnected();
    int findAtomInCycle(int in, int out);
    vector<string> findFunctionalGroups(atom a, int out);
    pair<vector<int>, vector<int> > findFarthest(int s, int out);
    vector<int> findPathFrom(int s, int out, bool cycle);
    vector<pair<int, int> > findComplexBonds(vector<int> parent_chain);
    int findHighestPriority(vector<int> parent_chain, int out);
    vector<tuple<int, string, bool> > findSubstituents(vector<int> parent_chain, int out, bool prefix, bool want_names);
    bool isParentChainCyclic(vector<int> parent_chain);
    vector<int> directAcyclicParentChain(vector<int> parent_chain, int out);
    vector<int> directCyclicParentChain(vector<int> parent_chain, int out);
    vector<int> directParentChain(vector<int> parent_chain, int out);
    void findCyclicUnions(int in);
    void findAcyclicUnions(int in);
    void findUnions();
    vector<int> findParentChain(int in, int out);
    string generateName(int in, int out);
    string getName();
    void setName(string name, double x, double y, double distx, double disty);
};
#endif //COMPOUND_H_INCLUDED
