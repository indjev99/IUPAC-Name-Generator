#ifndef DRAW_H_INCLUDED
#define DRAW_H_INCLUDED
#include<GLFW/glfw3.h>
#include<string>
using namespace std;

extern const double DEG2RAD;
extern const double TEXT_COLOUR_R;
extern const double TEXT_COLOUR_G;
extern const double TEXT_COLOUR_B;
extern const double BACKGROUND_COLOUR_R;
extern const double BACKGROUND_COLOUR_R2;
extern const double BACKGROUND_COLOUR_G;
extern const double BACKGROUND_COLOUR_B;

void drawPartEllipse(float x, float y, float radiusX, float radiusY, double alpha, double beta);
double drawSymbol1(char symbol, double x, double y, bool centered);
double drawSymbol(string symbol, double x, double y, bool centered);
double drawIndex(int index, double x, double y);
void drawBond(double x1, double y1, double x2, double y2, int num);
void drawCompound(GLFWwindow* w, compound& c);
#endif //DRAW_H_INCLUDED
