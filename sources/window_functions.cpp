#include<iostream>
#include<GLFW/glfw3.h>
#include "../headers/compound.h"
#include "../headers/window_size.h"
#include "../headers/draw.h"
#include "../headers/window_functions.h"

double mxpos=0;
double mypos=0;
int pressed;
bool selected_element[16];
bool snappingEnabled=1;

void error_callback(int error, const char* description)
{
    cerr<<error<<": "<<description<<endl;
}
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key==GLFW_KEY_ESCAPE && action==GLFW_PRESS) glfwSetWindowShouldClose(window,1);
    if (key==GLFW_KEY_BACKSPACE && action==GLFW_PRESS) pressed=-3;
    if ((key==GLFW_KEY_RIGHT_SHIFT || key==GLFW_KEY_LEFT_SHIFT) && action==GLFW_PRESS) snappingEnabled=!snappingEnabled;
    if (key==GLFW_KEY_R && action==GLFW_PRESS) pressed=-5;
    if (key==GLFW_KEY_L && action==GLFW_PRESS) pressed=-6;
    if (key==GLFW_KEY_N && action==GLFW_PRESS) pressed=-7;
    if (key==GLFW_KEY_H && action==GLFW_PRESS) pressed=-10;

    if (key==GLFW_KEY_O && action==GLFW_PRESS) selected_element[2]=1;
    if (key==GLFW_KEY_F && action==GLFW_PRESS) selected_element[3]=1;
    if (key==GLFW_KEY_C && action==GLFW_PRESS) selected_element[4]=1;
    if (key==GLFW_KEY_B && action==GLFW_PRESS) selected_element[5]=1;
    if (key==GLFW_KEY_I && action==GLFW_PRESS) selected_element[6]=1;

    if (key==GLFW_KEY_O && action==GLFW_RELEASE) selected_element[2]=0;
    if (key==GLFW_KEY_F && action==GLFW_RELEASE) selected_element[3]=0;
    if (key==GLFW_KEY_C && action==GLFW_RELEASE) selected_element[4]=0;
    if (key==GLFW_KEY_B && action==GLFW_RELEASE) selected_element[5]=0;
    if (key==GLFW_KEY_I && action==GLFW_RELEASE) selected_element[6]=0;
}
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (action==GLFW_PRESS)
    {
        pressed=button;
        glfwGetCursorPos(window,&mxpos,&mypos);
        mxpos=mxpos*2-WINDOWS_WIDTH;
        mypos=-mypos*2+WINDOWS_HEIGHT;
        mxpos/=ORIGINAL_WINDOWS_HEIGHT;
        mypos/=ORIGINAL_WINDOWS_HEIGHT;
    }
}
void window_size_callback(GLFWwindow* window, int width, int height)
{
    WINDOWS_WIDTH=width;
    WINDOWS_HEIGHT=height;
    pressed=-2;
}
string setCallbacks(GLFWwindow* w)
{
    glfwSetErrorCallback(error_callback);
    glfwSetKeyCallback(w,key_callback);
    glfwSetMouseButtonCallback(w,mouse_button_callback);
    glfwSetWindowSizeCallback(w,window_size_callback);
    return "Callbacks set successfully.";
}
string initializeGLFW()
{
    if (!glfwInit())
        return "Unable to initialize GLFW.";
    return "GLFW initialized successfully.";
}
string createWindow(GLFWwindow*& w)
{
    w=glfwCreateWindow(WINDOWS_WIDTH,WINDOWS_HEIGHT,"IUPAC Name Generator",NULL,NULL);
    if (!w)
    {
        glfwTerminate();
        return "Unable to open window.";
    }
    return "Window created successfully.";
}
void stopGraphics(GLFWwindow*& w)
{
    glfwDestroyWindow(w);
    glfwTerminate();
}
void drawWindow(GLFWwindow*& w,compound& c)
{
    glfwSetWindowShouldClose(w,0);
    pressed=-1;

    float ratio;
    int width, height;

    glfwMakeContextCurrent(w);
    glfwGetFramebufferSize(w,&width,&height);
    ratio=width/(float)height;
    glViewport(0,0,width,height);
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-ratio,ratio,-1.f,1.f,1.f,-1.f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    drawCompound(w,c);

    glfwSwapBuffers(w);

    while (!glfwWindowShouldClose(w) && pressed==-1)
    {
        glfwPollEvents();
    }
}
