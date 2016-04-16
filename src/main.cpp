//IUPAC Name Generator
#include<iostream>
#include<conio.h>
#include "dictionary.cpp"
#include "atom.cpp"
#include "compound.cpp"
#include "window_functions.cpp"
#include "run.cpp"
using namespace std;
GLFWwindow *window;



int main()
{
    system("chcp 1251");
    system("cls");

    setDictionaries();

    help();
    cout<<curr_dict.PACTC<<endl;
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
