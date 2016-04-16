#ifndef DICTIONARY_H_INCLUDED
#define DICTIONARY_H_INCLUDED
#include<vector>
#include<stack>
#include<unordered_map>
using namespace std;
struct dictionary
{
    vector<string> CNP; //carbon_number_prefixes
    vector<string> NP; //number_prefixes
    vector<string> CBI; //complex_bond_infixes
    vector<string> SS; //substituent_suffix
    string CP; //cycle_prefix
    string NC; //not_connected
    string error;
    string help;
    string PACTC; //Press any key to continue
    unordered_map<string, string> HP; //halogen_prefixes
    unordered_map<string, string> FGP; //functional_group_prefixes
    unordered_map<string, string> FGS; //functional_group_suffixes
    string benzene;
    string phen;

    string getCNP(int CN);
    string getNP(int ST);
};
dictionary curr_dict;
vector<dictionary> dictionaries;
#endif //DICTIONARY_H_INCLUDED
