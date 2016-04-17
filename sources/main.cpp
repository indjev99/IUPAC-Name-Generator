//IUPAC Name Generator
#include<iostream>
#include<conio.h>
#include<stdlib.h>
#include "../headers/dictionary.h"
#include "../headers/atom.h"
#include "../headers/compound.h"
#include "../headers/window_functions.h"
#include "../headers/run.h"
using namespace std;

GLFWwindow *window;

int main()
{
    system("chcp 1251");
    system("cls");

    setDictionaries();

    help();
    cout<<dictionaries[curr_dict_N].PACTC<<endl;
    getch();

    system("cls");

    string message;
    message=initializeGLFW();
    cerr<<message<<endl;
    if (message!="GLFW initialized successfully.") return -1;

    message=createWindow(window);
    cerr<<message<<endl;
    if (message!="Window created successfully.") return -1;

    message=setCallbacks(window);
    cerr<<message<<endl;
    if (message!="Callbacks set successfully.") return -1;

    system("cls");
    run(window);

    stopGraphics(window);
    return 0;
}
