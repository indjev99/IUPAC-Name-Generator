#ifndef WINDOW_FUNCTIONS_H_INCLUDED
#define WINDOW_FUNCTIONS_H_INCLUDED
#include<GLFW/glfw3.h>
#include<string>
#include "compound.h"

extern double mxpos;
extern double mypos;
extern int pressed;
extern bool selected_element[16];
extern bool snappingEnabled;

void error_callback(int error, const char* description);
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
void window_size_callback(GLFWwindow* window, int width, int height);
string setCallbacks(GLFWwindow* w);
string initializeGLFW();
string createWindow(GLFWwindow*& w);
void stopGraphics(GLFWwindow*& w);
void drawWindow(GLFWwindow*& w,compound& c);
#endif //WINDOW_FUNCTIONS_H_INCLUDED
