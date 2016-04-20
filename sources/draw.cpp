#include<math.h>
#include<vector>
#include "../headers/atom.h"
#include "../headers/compound.h"
#include "../headers/window_size.h"
#include "../headers/draw.h"
using namespace std;

const double DEG2RAD=3.14159/180.0;
double TEXT_COLOUR_R=0;
double TEXT_COLOUR_G=0;
double TEXT_COLOUR_B=0;
const double BACKGROUND_COLOUR_R=1;
const double BACKGROUND_COLOUR_R2=1;
const double BACKGROUND_COLOUR_G=1;
const double BACKGROUND_COLOUR_B=1;

void drawPartEllipse(float x, float y, float radiusX, float radiusY, double alpha, double beta)
{
    alpha=round(alpha*2);
    beta=round(beta*2);
    glBegin(GL_TRIANGLES);
    for(int i=alpha; i<beta; ++i)
    {
        float rad=i*0.5*DEG2RAD;
        float rad2=(i+1)*0.5;
        while (rad2>=360) rad2-=360;
        rad2*=DEG2RAD;
        glVertex2f(cos(rad)*radiusX+x,sin(rad)*radiusY+y);
        glVertex2f(cos(rad2)*radiusX+x,sin(rad2)*radiusY+y);
        glVertex2f(x,y);
    }
    glEnd();
}
double drawSymbol1(char symbol, double x, double y, bool centered)
{
    double scale=1.0*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;
    if (symbol=='C')
    {
        if (!centered) x+=0.055*scale;

        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);
        drawPartEllipse(x,y,0.055*scale,0.085*scale,0,360);

        glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);
        drawPartEllipse(x,y,0.05*scale,0.08*scale,50,310);

        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);
        drawPartEllipse(x,y,0.035*scale,0.065*scale,0,360);
        return x+0.035*scale;
    }
    if (symbol=='O')
    {
        if (!centered) x+=0.055*scale;

        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);
        drawPartEllipse(x,y,0.055*scale,0.085*scale,0,360);

        glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);
        drawPartEllipse(x,y,0.05*scale,0.08*scale,0,360);
        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

        drawPartEllipse(x,y,0.035*scale,0.065*scale,0,360);
        return x+0.035*scale;
    }
    if (symbol=='H')
    {
        if (!centered) x+=0.055*scale;

        glBegin(GL_QUADS);

        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

        glVertex2f(x-0.055*scale,y-0.085*scale);
        glVertex2f(x+0.055*scale,y-0.085*scale);
        glVertex2f(x+0.055*scale,y+0.085*scale);
        glVertex2f(x-0.055*scale,y+0.085*scale);

        glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

        glVertex2f(x-0.05*scale,y-0.08*scale);
        glVertex2f(x-0.035*scale,y-0.08*scale);
        glVertex2f(x-0.035*scale,y+0.08*scale);
        glVertex2f(x-0.05*scale,y+0.08*scale);

        glVertex2f(x+0.05*scale,y-0.08*scale);
        glVertex2f(x+0.035*scale,y-0.08*scale);
        glVertex2f(x+0.035*scale,y+0.08*scale);
        glVertex2f(x+0.05*scale,y+0.08*scale);

        glVertex2f(x-0.035*scale,y-0.0075*scale);
        glVertex2f(x+0.035*scale,y-0.0075*scale);
        glVertex2f(x+0.035*scale,y+0.0075*scale);
        glVertex2f(x-0.035*scale,y+0.0075*scale);

        glEnd();
        return x+0.05*scale;
    }
    if (symbol=='F')
    {
        if (!centered) x+=0.055*scale;

        glBegin(GL_QUADS);

        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

        glVertex2f(x-0.055*scale,y-0.085*scale);
        glVertex2f(x+0.04*scale,y-0.085*scale);
        glVertex2f(x+0.04*scale,y+0.085*scale);
        glVertex2f(x-0.055*scale,y+0.085*scale);

        glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

        glVertex2f(x-0.05*scale,y-0.08*scale);
        glVertex2f(x-0.035*scale,y-0.08*scale);
        glVertex2f(x-0.035*scale,y+0.08*scale);
        glVertex2f(x-0.05*scale,y+0.08*scale);

        glVertex2f(x-0.035*scale,y+0.08*scale);
        glVertex2f(x+0.035*scale,y+0.08*scale);
        glVertex2f(x+0.035*scale,y+0.065*scale);
        glVertex2f(x-0.035*scale,y+0.065*scale);

        glVertex2f(x-0.035*scale,y-0.0075*scale);
        glVertex2f(x+0.035*scale,y-0.0075*scale);
        glVertex2f(x+0.035*scale,y+0.0075*scale);
        glVertex2f(x-0.035*scale,y+0.0075*scale);

        glEnd();
        return x+0.05*scale;
    }
    if (symbol=='E')
    {
        if (!centered) x+=0.055*scale;

        glBegin(GL_QUADS);

        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

        glVertex2f(x-0.055*scale,y-0.085*scale);
        glVertex2f(x+0.04*scale,y-0.085*scale);
        glVertex2f(x+0.04*scale,y+0.085*scale);
        glVertex2f(x-0.055*scale,y+0.085*scale);

        glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

        glVertex2f(x-0.05*scale,y-0.08*scale);
        glVertex2f(x-0.035*scale,y-0.08*scale);
        glVertex2f(x-0.035*scale,y+0.08*scale);
        glVertex2f(x-0.05*scale,y+0.08*scale);

        glVertex2f(x-0.035*scale,y+0.08*scale);
        glVertex2f(x+0.035*scale,y+0.08*scale);
        glVertex2f(x+0.035*scale,y+0.065*scale);
        glVertex2f(x-0.035*scale,y+0.065*scale);

        glVertex2f(x-0.035*scale,y-0.0075*scale);
        glVertex2f(x+0.035*scale,y-0.0075*scale);
        glVertex2f(x+0.035*scale,y+0.0075*scale);
        glVertex2f(x-0.035*scale,y+0.0075*scale);

        glVertex2f(x-0.035*scale,y-0.08*scale);
        glVertex2f(x+0.035*scale,y-0.08*scale);
        glVertex2f(x+0.035*scale,y-0.065*scale);
        glVertex2f(x-0.035*scale,y-0.065*scale);

        glEnd();
        return x+0.05*scale;
    }
    if (symbol=='I')
    {
        if (!centered) x+=0.0275*scale;

        glBegin(GL_QUADS);

        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

        glVertex2f(x-0.0275*scale,y-0.085*scale);
        glVertex2f(x+0.0275*scale,y-0.085*scale);
        glVertex2f(x+0.0275*scale,y+0.085*scale);
        glVertex2f(x-0.0275*scale,y+0.085*scale);

        glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

        glVertex2f(x-0.0225*scale,y+0.08*scale);
        glVertex2f(x+0.0225*scale,y+0.08*scale);
        glVertex2f(x+0.0225*scale,y+0.065*scale);
        glVertex2f(x-0.0225*scale,y+0.065*scale);

        glVertex2f(x-0.0075*scale,y-0.08*scale);
        glVertex2f(x+0.0075*scale,y-0.08*scale);
        glVertex2f(x+0.0075*scale,y+0.08*scale);
        glVertex2f(x-0.0075*scale,y+0.08*scale);


        glVertex2f(x-0.0225*scale,y-0.08*scale);
        glVertex2f(x+0.0225*scale,y-0.08*scale);
        glVertex2f(x+0.0225*scale,y-0.065*scale);
        glVertex2f(x-0.0225*scale,y-0.065*scale);

        glEnd();
        return x+0.0225*scale;
    }
    if (symbol=='B')
    {
        if (!centered) x+=0.05*scale;

        glBegin(GL_QUADS);

        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

        glVertex2f(x-0.05*scale,y-0.085*scale);
        glVertex2f(x-0.025*scale,y-0.085*scale);
        glVertex2f(x-0.025*scale,y+0.085*scale);
        glVertex2f(x-0.05*scale,y+0.085*scale);

        glEnd();

        drawPartEllipse(x-0.03*scale,y+0.04*scale,0.05*scale,0.045*scale,270,360);
        drawPartEllipse(x-0.03*scale,y+0.04*scale,0.05*scale,0.045*scale,0,90);

        drawPartEllipse(x-0.03*scale,y-0.035*scale,0.055*scale,0.05*scale,270,360);
        drawPartEllipse(x-0.03*scale,y-0.035*scale,0.055*scale,0.05*scale,0,90);

        drawPartEllipse(x-0.03*scale,y,0.055*scale,0.085*scale,270,360);
        drawPartEllipse(x-0.03*scale,y,0.055*scale,0.085*scale,0,90);

        glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

        drawPartEllipse(x-0.03*scale,y+0.04*scale,0.045*scale,0.04*scale,270,360);
        drawPartEllipse(x-0.03*scale,y+0.04*scale,0.045*scale,0.04*scale,0,90);

        drawPartEllipse(x-0.03*scale,y-0.035*scale,0.05*scale,0.045*scale,270,360);
        drawPartEllipse(x-0.03*scale,y-0.035*scale,0.05*scale,0.045*scale,0,90);

        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

        drawPartEllipse(x-0.03*scale,y+0.04*scale,0.03*scale,0.025*scale,270,360);
        drawPartEllipse(x-0.03*scale,y+0.04*scale,0.03*scale,0.025*scale,0,90);

        drawPartEllipse(x-0.03*scale,y-0.035*scale,0.035*scale,0.03*scale,270,360);
        drawPartEllipse(x-0.03*scale,y-0.035*scale,0.035*scale,0.03*scale,0,90);

        glBegin(GL_QUADS);

        glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);
        glVertex2f(x-0.045*scale,y-0.08*scale);
        glVertex2f(x-0.03*scale,y-0.08*scale);
        glVertex2f(x-0.03*scale,y+0.08*scale);
        glVertex2f(x-0.045*scale,y+0.08*scale);

        glEnd();
        return x+0.02*scale;
    }
    if (symbol=='l')
    {
        if (!centered) x+=0.0125*scale;

        glBegin(GL_QUADS);

        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

        glVertex2f(x-0.0125*scale,y-0.085*scale);
        glVertex2f(x+0.0125*scale,y-0.085*scale);
        glVertex2f(x+0.0125*scale,y+0.085*scale);
        glVertex2f(x-0.0125*scale,y+0.085*scale);

        glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

        glVertex2f(x-0.0075*scale,y-0.08*scale);
        glVertex2f(x+0.0075*scale,y-0.08*scale);
        glVertex2f(x+0.0075*scale,y+0.08*scale);
        glVertex2f(x-0.0075*scale,y+0.08*scale);

        glEnd();
        return x+0.0075*scale;
    }
    if (symbol=='r')
    {
        if (!centered) x+=0.04125*scale;

        glBegin(GL_QUADS);

        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

        glVertex2f(x-0.04125*scale,y-0.085*scale);
        glVertex2f(x-0.01625*scale,y-0.085*scale);
        glVertex2f(x-0.01625*scale,y+0.01*scale);
        glVertex2f(x-0.04125*scale,y+0.01*scale);

        glVertex2f(x-0.04125*scale,y-0.085*scale);
        glVertex2f(x-0.01625*scale,y-0.085*scale);
        glVertex2f(x+0.046*scale,y-0.01*scale);
        glVertex2f(x-0.04125*scale,y-0.01*scale);

        glEnd();

        drawPartEllipse(x,y-0.06*scale,0.07*scale,0.065*scale,50,125);

        glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

        drawPartEllipse(x,y-0.06*scale,0.065*scale,0.06*scale,55,125);

        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

        drawPartEllipse(x,y-0.06*scale,0.05*scale,0.045*scale,55,125);

        glBegin(GL_QUADS);

        glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

        glVertex2f(x-0.03625*scale,y-0.08*scale);
        glVertex2f(x-0.02125*scale,y-0.08*scale);
        glVertex2f(x-0.02125*scale,y+0.005*scale);
        glVertex2f(x-0.03625*scale,y+0.005*scale);

        glEnd();
        return x+0.03625*scale;
    }
    return x;
}
double drawSymbol(string symbol, double x, double y, bool centered)
{
    double scale=1.0*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;
    //cerr<<symbol<<endl;
    double nextpos=x;
    for (int i=0; i<symbol.size(); ++i)
    {
        //cerr<<" "<<symbol[i];
        if (!i) nextpos=drawSymbol1(symbol[i],nextpos,y,centered);
        else nextpos=drawSymbol1(symbol[i],nextpos,y,0);
        //cerr<<" "<<nextpos<<endl;
    }
    return nextpos+0.02*scale;
}
double drawIndex(int index, double x, double y)
{
    double scale=1.0*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;
    if (index<=1) return x;
    vector<int> digits;
    do
    {
        digits.push_back(index%10);
        index/=10;
    }
    while (index!=0);
    for (int i=digits.size()-1; i>=0; --i)
    {
        if (i<=digits.size()-1)
        {
            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glBegin(GL_QUADS);

            glVertex2f(x,y-0.005*scale);
            glVertex2f(x+0.015*scale,y-0.005*scale);
            glVertex2f(x+0.015*scale,y+0.085*scale);
            glVertex2f(x,y+0.085*scale);

            glEnd();

            x+=0.01;
        }
        if (digits[i]==9)
        {
            glBegin(GL_QUADS);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x,y-0.005*scale);
            glVertex2f(x+0.045*scale,y-0.005*scale);
            glVertex2f(x+0.045*scale,y+0.085*scale);
            glVertex2f(x,y+0.085*scale);

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);;

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

            glVertex2f(x,y);
            glVertex2f(x+0.04*scale,y);
            glVertex2f(x+0.04*scale,y+0.08*scale);
            glVertex2f(x,y+0.08*scale);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x+0.0125*scale,y+0.045*scale);
            glVertex2f(x+0.0275*scale,y+0.045*scale);
            glVertex2f(x+0.0275*scale,y+0.07*scale);
            glVertex2f(x+0.0125*scale,y+0.07*scale);

            glVertex2f(x,y+0.01*scale);
            glVertex2f(x+0.0275*scale,y+0.01*scale);
            glVertex2f(x+0.0275*scale,y+0.035*scale);
            glVertex2f(x,y+0.035*scale);

            glEnd();

            x+=0.04*scale;
        }
        if (digits[i]==4)
        {
            glBegin(GL_QUADS);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x,y-0.005*scale);
            glVertex2f(x+0.045*scale,y-0.005*scale);
            glVertex2f(x+0.045*scale,y+0.085*scale);
            glVertex2f(x,y+0.085*scale);

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);;

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

            glVertex2f(x,y);
            glVertex2f(x+0.04*scale,y);
            glVertex2f(x+0.04*scale,y+0.08*scale);
            glVertex2f(x,y+0.08*scale);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x+0.0125*scale,y+0.045*scale);
            glVertex2f(x+0.0275*scale,y+0.045*scale);
            glVertex2f(x+0.0275*scale,y+0.08*scale);
            glVertex2f(x+0.0125*scale,y+0.08*scale);

            glVertex2f(x,y);
            glVertex2f(x+0.0275*scale,y);
            glVertex2f(x+0.0275*scale,y+0.035*scale);
            glVertex2f(x,y+0.035*scale);

            glEnd();

            x+=0.04*scale;
        }
        if (digits[i]==3)
        {
            glBegin(GL_QUADS);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x,y-0.005*scale);
            glVertex2f(x+0.0375*scale,y-0.005*scale);
            glVertex2f(x+0.0375*scale,y+0.085*scale);
            glVertex2f(x,y+0.085*scale);

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);;

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

            glVertex2f(x,y);
            glVertex2f(x+0.0325*scale,y);
            glVertex2f(x+0.0325*scale,y+0.08*scale);
            glVertex2f(x,y+0.08*scale);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x,y+0.045*scale);
            glVertex2f(x+0.02*scale,y+0.045*scale);
            glVertex2f(x+0.02*scale,y+0.07*scale);
            glVertex2f(x,y+0.07*scale);

            glVertex2f(x,y+0.01*scale);
            glVertex2f(x+0.02*scale,y+0.01*scale);
            glVertex2f(x+0.02*scale,y+0.035*scale);
            glVertex2f(x,y+0.035*scale);

            glEnd();

            x+=0.0325*scale;
        }
        if (digits[i]==2)
        {
            glBegin(GL_QUADS);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x,y-0.005*scale);
            glVertex2f(x+0.0375*scale,y-0.005*scale);
            glVertex2f(x+0.0375*scale,y+0.085*scale);
            glVertex2f(x,y+0.085*scale);

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);;

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

            glVertex2f(x,y);
            glVertex2f(x+0.0325*scale,y);
            glVertex2f(x+0.0325*scale,y+0.08*scale);
            glVertex2f(x,y+0.08*scale);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x,y+0.045*scale);
            glVertex2f(x+0.02*scale,y+0.045*scale);
            glVertex2f(x+0.02*scale,y+0.07*scale);
            glVertex2f(x,y+0.07*scale);

            glVertex2f(x+0.0125*scale,y+0.01*scale);
            glVertex2f(x+0.0325*scale,y+0.01*scale);
            glVertex2f(x+0.0325*scale,y+0.035*scale);
            glVertex2f(x+0.0125*scale,y+0.035*scale);

            glEnd();

            x+=0.0325*scale;
        }
        if (digits[i]==1)
        {
            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glBegin(GL_QUADS);

            glVertex2f(x,y-0.005*scale);
            glVertex2f(x+0.0175*scale,y-0.005*scale);
            glVertex2f(x+0.0175*scale,y+0.085*scale);
            glVertex2f(x,y+0.085*scale);

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);;

            glVertex2f(x,y);
            glVertex2f(x+0.0125*scale,y);
            glVertex2f(x+0.0125*scale,y+0.08*scale);
            glVertex2f(x,y+0.08*scale);

            glEnd();

            x+=0.0125*scale;
        }
        if (digits[i]==0)
        {
            glBegin(GL_QUADS);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x,y-0.005*scale);
            glVertex2f(x+0.04*scale,y-0.005*scale);
            glVertex2f(x+0.04*scale,y+0.085*scale);
            glVertex2f(x,y+0.085*scale);

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);;

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

            glVertex2f(x,y);
            glVertex2f(x+0.035*scale,y);
            glVertex2f(x+0.035*scale,y+0.08*scale);
            glVertex2f(x,y+0.08*scale);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x+0.0125*scale,y+0.01*scale);
            glVertex2f(x+0.0225*scale,y+0.01*scale);
            glVertex2f(x+0.0225*scale,y+0.07*scale);
            glVertex2f(x+0.0125*scale,y+0.07*scale);

            glEnd();

            x+=0.035*scale;
        }
    }
    return x;
}
void drawAtom(atom& a)
{
    double scale=1.0*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;
    double nextpos;
    if (a.free_bonds.empty())
    {
        drawSymbol(a.symbol,a.x*scale,a.y*scale,1);
    }
    else
    {
        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

        glBegin(GL_QUADS);

        glVertex2f(a.x*scale,a.y*scale-0.08*scale);
        glVertex2f(a.x*scale+0.055*scale,a.y*scale-0.08*scale);
        glVertex2f(a.x*scale+0.055*scale,a.y*scale+0.08*scale);
        glVertex2f(a.x*scale,a.y*scale+0.08*scale);

        glEnd();

        if (a.symbol=="C" || (a.symbol=="O" && a.free_bonds.size()<2))
        {
            nextpos=drawSymbol(a.symbol,a.x*scale,a.y*scale,1);
            nextpos=drawSymbol("H",nextpos,a.y*scale,0);
            drawIndex(a.free_bonds.size(),nextpos-0.02*scale,a.y*scale-0.08*scale);
        }
        else
        {
            nextpos=drawSymbol("H",a.x*scale,a.y*scale,1);
            nextpos=drawIndex(a.free_bonds.size(),nextpos-0.02*scale,a.y*scale-0.08*scale);
            nextpos+=0.02*scale;
            nextpos=drawSymbol(a.symbol,nextpos,a.y*scale,0);
        }
    }
}
void drawBond(double x1, double y1, double x2, double y2, int num)
{
    double scale=1.0*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;
    double offset,x,y,deltax,deltay;
    glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);
    glLineWidth(4.0);
    glBegin(GL_LINES);
    double alpha,beta;
    deltax=x2-x1;
    deltay=y2-y1;
    if (deltax==0)
    {
        x=1;
        y=0;
    }
    else
    {
        alpha=atan(deltay/deltax);
        beta=90*DEG2RAD-alpha;
        x=cos(beta);
        y=sin(beta);
    }
    x=-x;
    for (int i=0; i<num; ++i)
    {
        offset=((num-1)*0.5-i)*0.02;
        glVertex2f((x1+offset*x)*scale,(y1+offset*y)*scale);
        glVertex2f((x2+offset*x)*scale,(y2+offset*y)*scale);
    }
    glEnd();
}
void setTextColour(int n)
{
    int d;
    ++n;
    n%=7;
    d=n%2;
    n/=2;
    TEXT_COLOUR_R=d;
    d=n%2;
    n/=2;
    TEXT_COLOUR_G=d;
    d=n%2;
    n/=2;
    TEXT_COLOUR_B=d;
}
void drawCompound(GLFWwindow* w, compound& c)
{
    int p;
    double proportion=1.0*WINDOWS_WIDTH/WINDOWS_HEIGHT;
    atom a,a2;

    //background
    //glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);
    glColor3f(BACKGROUND_COLOUR_R2,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

    glBegin(GL_QUADS);

    glVertex2f(-1.0*proportion,-1.0);
    glVertex2f(1.0*proportion,-1.0);
    glVertex2f(1.0*proportion,1.0);
    glVertex2f(-1.0*proportion,1.0);

    glEnd();

    for (int i=0; i<c.atoms.size(); ++i)
    {
        a=c.atoms[i];
        if (a.symbol=="") continue;
        setTextColour(-1);
        for (int i=0; i<a.bonds.size(); ++i)
        {
            if (a.bonds[i].to!=-1)
            {
                a2=c.atoms[a.bonds[i].to];
                drawBond(a.x,a.y,a2.x,a2.y,a.bonds[i].spots_taken.size());
            }
        }
    }

    for (int i=0; i<c.atoms.size(); ++i)
    {
        a=c.atoms[i];
        if (a.symbol=="") continue;
        p=-1;
        if (i<c.atoms_unions.size()) p=c.atoms_unions[i];
        setTextColour(p);
        drawAtom(a);
    }
}
