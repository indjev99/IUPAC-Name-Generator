#ifndef RUN_H_INCLUDED
#define RUN_H_INCLUDED
#include<deque>
using namespace std;

extern const string element_symbol[];
extern const int element_valence[];
extern int curr_dict_N;

void help();
void snap(double& x, double& y);
void run(GLFWwindow* w);
#endif //RUN_H_INCLUDED
