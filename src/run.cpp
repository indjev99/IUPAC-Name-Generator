#include<deque>
#include "compound.h"
#include "window_functions.cpp"
using namespace std;
int curr_dict_N=0;
const vector<string> element_symbol= {"H","C","O","F","Cl","Br","I"};
const vector<int> element_valence= {1,4,2,1,1,1,1};
void snap(double& x, double& y)
{
    double hx,lx;
    //cerr<<x<<" "<<y<<endl;
    hx=(round(x*3)+0.2)/3;
    lx=(round(x*3)-1+0.2)/3;
    //cerr<<hx<<" > "<<x<<" > "<<lx<<endl;
    if (hx-x<x-lx) x=hx;
    else x=lx;
    y=round(y*3.7)/3.7;
    //cerr<<x<<" "<<y<<endl;
}
void run(GLFWwindow* w)
{
    string curr_symbol;
    string name;
    int curr_valence;
    double sx,sy;
    sx=0;
    sy=0;
    if (snappingEnabled) snap(sx,sy);
    compound c(element_symbol[1],element_valence[1],sx,sy);
    deque<compound> history= {c};
    int last2=-1;
    int last=-1;
    bool toMove=0;
    int result;
    while (!glfwWindowShouldClose(w))
    {
        if (history.size()>50) history.pop_front();
        name=c.getName();
        cout<<name<<endl;
        drawWindow(w,c);
        system("cls");
        curr_symbol=element_symbol[1];
        curr_valence=element_valence[1];
        for (int i=0; i<7; ++i)
        {
            if (selected_element[i])
            {
                curr_symbol=element_symbol[i];
                curr_valence=element_valence[i];
            }
        }
        if (pressed==-3)
        {
            last2=-1;
            last=-1;
            toMove=0;
            if (history.size()>1) history.pop_back();
            c=history.back();
        }
        if (pressed==-5)
        {
            last2=-1;
            last=-1;
            toMove=0;
            sx=0;
            sy=0;
            if (snappingEnabled) snap(sx,sy);
            c=*(new compound(element_symbol[1],element_valence[1],sx,sy));
            history.push_back(c);
            //BACKGROUND_COLOUR_R2=!BACKGROUND_COLOUR_R2;
        }
        if (pressed==-6)
        {
            ++curr_dict_N;
            curr_dict_N%=dictionaries.size();
            curr_dict=dictionaries[curr_dict_N];
        }
        if (pressed==-7)
        {
            system("cls");
            cin>>name;
            sx=0;
            sy=0;
            if (snappingEnabled) snap(sx,sy);
            c.setName(name,sx,sy,1.0/3,1.0/3.7);
            history.push_back(c);
        }
        if (pressed==-10)
        {
            help();
        }
        last=c.findAtom(mxpos,mypos);
        if (snappingEnabled) snap(mxpos,mypos);
        if (last==-1 && c.findAtom(mxpos,mypos)!=-1) continue;
        if (pressed==GLFW_MOUSE_BUTTON_LEFT || pressed==GLFW_MOUSE_BUTTON_MIDDLE)
        {
            //cerr<<"Connect: "<<last<<" "<<last2<<'\n'<<endl;
            if (last2!=-1)
            {
                if (toMove)
                {
                    result=c.moveAtom(last2,mxpos,mypos);
                    history.push_back(c);
                    if (pressed==GLFW_MOUSE_BUTTON_LEFT || result==-1)
                    {
                        last2=-1;
                    }
                    toMove=0;
                }
                else if (last!=-1)
                {
                    result=c.connectAtoms(last2,last);
                    if (result!=-1)
                    {
                        history.push_back(c);
                        if (pressed==GLFW_MOUSE_BUTTON_LEFT || !c.atoms[last].canConnect(-1)) last2=-1;
                        else last2=last;
                    }
                    last=-1;

                }
                else
                {
                    result=c.addAtom(curr_symbol,curr_valence,mxpos,mypos,last2);
                    if (result!=-1)
                    {
                        history.push_back(c);
                        if (pressed==GLFW_MOUSE_BUTTON_LEFT || !c.atoms[result].canConnect(-1)) last2=-1;
                        else last2=result;
                    }
                    last=-1;
                }
            }
            else
            {
                if (last==-1)
                {
                    result=c.addAtom(curr_symbol,curr_valence,mxpos,mypos);
                    history.push_back(c);
                    if (result!=-1)
                    {
                        if (pressed==GLFW_MOUSE_BUTTON_LEFT || !c.atoms[result].canConnect(-1)) last2=-1;
                        else last2=result;
                    }
                }
                else
                {
                    if (pressed==GLFW_MOUSE_BUTTON_LEFT) toMove=1;
                    else toMove=0;
                    if (toMove || c.atoms[last].canConnect(-1))
                    {
                        last2=last;
                    }
                }
            }
        }
        if (pressed==GLFW_MOUSE_BUTTON_RIGHT)
        {
            last2=-1;
            toMove=0;
            if (last!=-1)
            {
                result=c.removeAtom(last);
                history.push_back(c);
            }
        }
    }
}
