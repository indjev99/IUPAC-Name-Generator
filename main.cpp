//IUPAC
#include<iostream>
#include<algorithm>
#include<vector>
#include<stack>
#include<stdlib.h>
#include<conio.h>
#include<GLFW/glfw3.h>
using namespace std;

const int ORIGINAL_WINDOWS_WIDTH=1100;
const int ORIGINAL_WINDOWS_HEIGHT=1050;
int WINDOWS_WIDTH=ORIGINAL_WINDOWS_WIDTH;
int WINDOWS_HEIGHT=ORIGINAL_WINDOWS_HEIGHT;
double DEG2RAD=3.14159/180.0;
double TEXT_COLOUR_R=0;
double TEXT_COLOUR_G=0;
double TEXT_COLOUR_B=0;
double BACKGROUND_COLOUR_R=1;
double BACKGROUND_COLOUR_G=1;
double BACKGROUND_COLOUR_B=1;

string Hsymbol="H";
int pressed;
double mxpos=0,mypos=0;

bool snappingEnabled=1;

struct connection
{
    vector<int> spots_taken;
    int to;
};
struct dictionary
{
    vector<string> prefix;
};
struct atom
{
    string symbol;
    double x,y;
    vector<connection> connections;
    stack<int> free_connections;

    void connect(int connection)
    {
        if (!free_connections.empty())
        {
            int nc=isConnected(connection);
            if (nc==-1)
            {
                nc=free_connections.top();
            }
            connections[nc].to=connection;
            connections[nc].spots_taken.push_back(free_connections.top());
            free_connections.pop();
        }
    }
    int isConnected(int connection)
    {
        for (int i=0;i<connections.size();++i)
        {
            if (connections[i].to==connection) return i;
        }
        return -1;
    }
    void removeConnection(int connection)
    {
        for (int i=0;i<connections.size();++i)
        {
            if (connections[i].to==connection)
            {
                connections[i].to=-1;
                for (int j=0;j<connections[i].spots_taken.size();++j)
                {
                    free_connections.push(connections[i].spots_taken[j]);
                }
                connections[i].spots_taken.resize(0);
            }
        }
    }
    void changeXY(double new_x, double new_y)
    {
        x=new_x;
        y=new_y;
    }
    atom() {}
    atom(string new_symbol, int new_valance, double new_x, double new_y)
    {
        symbol=new_symbol;
        x=new_x;
        y=new_y;
        connections.resize(new_valance);
        for (int i=new_valance-1;i>=0;--i)
        {
            connections[i].to=-1;
            free_connections.push(i);
        }
    }
    atom(string new_symbol, int new_valance, double new_x, double new_y, int new_connection)
    {
        symbol=new_symbol;
        x=new_x;
        y=new_y;
        connections.resize(new_valance);
        for (int i=new_valance-1;i>=0;--i)
        {
            connections[i].to=-1;
            free_connections.push(i);
        }
        connect(new_connection);
    }
};
struct compound
{
    vector<atom> atoms;
    stack<int> free_positions;
    int addAtom(atom& a)
    {
        int ind;
        if (free_positions.empty())
        {
            ind=atoms.size();
            atoms.push_back(a);
        }
        else
        {
            ind=free_positions.top();
            free_positions.pop();
            atoms[ind]=a;
        }
        return ind;
    }
    int addAtom(string new_symbol, int new_valance, double new_x, double new_y)
    {
        if (new_valance<0) return -1;
        atom a(new_symbol,new_valance,new_x,new_y);
        return addAtom(a);
    }
    int addAtom(string new_symbol, int new_valance, double new_x, double new_y, int new_connection)
    {
        int ind=-1;
        if (new_valance<=0) return -1;
        if (new_connection>=atoms.size() || atoms[new_connection].symbol=="") return -1;
        atom a(new_symbol,new_valance,new_x,new_y,new_connection);
        if (!atoms[new_connection].free_connections.empty())
        {
            ind=addAtom(a);
            if (ind!=-1) atoms[new_connection].connect(ind);
        }
        return ind;
    }
    int connectAtoms(int a1, int a2)
    {
        if (a1<atoms.size() && a2<atoms.size() && a1!=a2 && atoms[a1].symbol!="" && atoms[a2].symbol!="" && !atoms[a1].free_connections.empty() && !atoms[a2].free_connections.empty())
        {
            atoms[a1].connect(a2);
            atoms[a2].connect(a1);
            return 1;
        }
        return -1;
    }
    int removeAtom(int a)
    {
        atom a1;
        if (a<atoms.size())
        {
            a1=atoms[a];
            if (a1.symbol!="")// && a1.connections.size()-a1.free_connections.size()<=1)
            {
                for (int i=0;i<a1.connections.size();++i)
                {
                    if (a1.connections[i].to!=-1) atoms[a1.connections[i].to].removeConnection(a);
                }
                atoms[a].symbol="";
                free_positions.push(a);
                return 1;
            }
        }
        return -1;
    }
    int moveAtom(int a, double x, double y)
    {
        if (a<atoms.size() && atoms[a].symbol!="")
        {
            atoms[a].changeXY(x,y);
            return 1;
        }
        return -1;
    }
    int findAtom(double x, double y)
    {
        atom a;
        int minDistInd=-1;
        double minDist=-1;
        double dist;
        for (int i=0;i<atoms.size();++i)
        {
            a=atoms[i];
            if (a.symbol=="") continue;
            dist=(a.x-x)*(a.x-x)+(a.y-y)*(a.y-y);
            if (dist<minDist || minDist==-1)
            {
                minDist=dist;
                minDistInd=i;
            }
        }
        if (minDistInd!=-1)
        {
            a=atoms[minDistInd];
            if (x<a.x-0.05 || (x>a.x+0.05 && a.free_connections.empty()) || x>a.x+0.15 || y<a.y-0.08 || y>a.y+0.08) minDistInd=-1;
        }
        return minDistInd;
    }
    compound() {}
    compound(string new_symbol, int new_valance, double new_x, double new_y)
    {
        addAtom(new_symbol,new_valance,new_x,new_y);
    }
};

//Graphics related functions
GLFWwindow* window,*window2;
void error_callback(int error, const char* description)
{
    cout<<error<<": "<<description<<endl;
}
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key==GLFW_KEY_ESCAPE && action==GLFW_PRESS) glfwSetWindowShouldClose(window,1);
    if (key==GLFW_KEY_BACKSPACE && action==GLFW_PRESS) pressed=-3;
    if (key==GLFW_KEY_RIGHT_SHIFT || key==GLFW_KEY_LEFT_SHIFT && action==GLFW_PRESS) snappingEnabled=!snappingEnabled;
    if (key==GLFW_KEY_R && action==GLFW_PRESS) pressed=-5;
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
    w=glfwCreateWindow(WINDOWS_WIDTH,WINDOWS_HEIGHT,"IUPAC naming",NULL,NULL);
    if (!w)
    {
        glfwTerminate();
        return "Unable to open window.";
    }
    return "Window created successfully.";
}
void stopGraphics()
{
    glfwDestroyWindow(window);
    glfwTerminate();
}

//Draw functions
void drawPartEllipse(float x, float y, float radiusX, float radiusY, double alpha, double beta)
{
    alpha=round(alpha*2);
    beta=round(beta*2);
    glBegin(GL_TRIANGLES);
    for(int i=alpha;i<beta;++i)
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
    if (symbol=='C')
    {
        if (!centered) x+=0.05*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;

        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);
        drawPartEllipse(x,y,0.055*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,0,360);

        glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);
        drawPartEllipse(x,y,0.05*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,50,310);

        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);
        drawPartEllipse(x,y,0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,0.065*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,0,360);
        return x+0.055*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;
    }
    if (symbol=='O')
    {
        if (!centered) x+=0.05*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;

        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);
        drawPartEllipse(x,y,0.055*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,0,360);

        glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);
        drawPartEllipse(x,y,0.05*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,0,360);
        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

        drawPartEllipse(x,y,0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,0.065*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,0,360);
        return x+0.055*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;
    }
    if (symbol=='H')
    {
        if (!centered) x+=0.05*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;

        glBegin(GL_QUADS);

        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

        /*glVertex2f(x-0.055*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x-0.03*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x-0.03*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x-0.055*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

        glVertex2f(x+0.055*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x+0.03*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x+0.03*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x+0.055*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

        glVertex2f(x-0.04*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.0125*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x+0.04*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.0125*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x+0.04*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.0125*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x-0.04*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.0125*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);*/

        glVertex2f(x-0.055*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x+0.055*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x+0.055*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x-0.055*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

        glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

        glVertex2f(x-0.05*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x-0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x-0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x-0.05*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

        glVertex2f(x+0.05*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x+0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x+0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x+0.05*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

        glVertex2f(x-0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.0075*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x+0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.0075*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x+0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.0075*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(x-0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.0075*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

        glEnd();
        return x+0.07*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;
    }
}
double drawSymbol(string symbol, double x, double y, bool centered)
{
    double nextpos=x;
    for (int i=0;i<symbol.size();++i)
    {
        if (!i) nextpos=drawSymbol1(symbol[i],nextpos,y,centered);
        else nextpos=drawSymbol1(symbol[i],nextpos,y,0);
    }
    return nextpos;
}
void drawIndex(int index, double x, double y)
{
    vector<int> digits;
    do
    {
        digits.push_back(index%2);
        index/=2;
    }
    while (index!=0);
    for (int i=digits.size()-1;i>=0;--i)
    {
        if (i<=digits.size()-1)
        {
            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glBegin(GL_QUADS);

            glVertex2f(x,y-0.005*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.015*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.005*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.015*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glEnd();

            x+=0.01;
        }
        if (digits[i]==1)
        {
            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glBegin(GL_QUADS);

            glVertex2f(x,y-0.005*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0175*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.005*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0175*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);;

            glVertex2f(x,y);
            glVertex2f(x+0.0125*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y);
            glVertex2f(x+0.0125*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glEnd();

            x+=0.0125*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;
        }
        if (digits[i]==0)
        {
            glBegin(GL_QUADS);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x,y-0.005*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.04*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.005*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.04*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);;

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

            glVertex2f(x,y);
            glVertex2f(x+0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y);
            glVertex2f(x+0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x+0.0125*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.01*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0225*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.01*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0225*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.07*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0125*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.07*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glEnd();

            x+=0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;
        }
    }
}
void drawAtom(atom& a)
{
    //cout<<a.x<<" "<<a.y<<" "<<a.symbol<<Hsymbol<<a.free_connections.size()<<'\n';
    double nextpos;
    if (a.free_connections.empty())
    {
        drawSymbol(a.symbol,a.x*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,a.y*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,1);
    }
    else
    {
        glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

        glBegin(GL_QUADS);

        glVertex2f(a.x,a.y-0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(a.x+0.055*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,a.y-0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(a.x+0.055*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,a.y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(a.x,a.y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

        glEnd();

        nextpos=drawSymbol(a.symbol,a.x*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,a.y*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,1);
        nextpos=drawSymbol(Hsymbol,nextpos,a.y*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,0);
        drawIndex(a.free_connections.size(),nextpos-0.02*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,a.y*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT-0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
    }
}
void drawConnection(double x1, double y1, double x2, double y2, int num)
{
    glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);
    glLineWidth(4.0);
    glBegin(GL_LINES);
    for (int i=0;i<num;++i)
    {
        glVertex2f((x1+i*0.02-num*0.5*0.02)*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,(y1+i*0.02-num*0.5*0.02)*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f((x2+i*0.02-num*0.5*0.02)*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,(y2+i*0.02-num*0.5*0.02)*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
    }
    glEnd();
}
void drawCompound(GLFWwindow* w, compound& c)
{
    atom a,a2;

    //background
    glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

    glBegin(GL_QUADS);

    glVertex2f(-1.0*WINDOWS_WIDTH/WINDOWS_HEIGHT,-1.0);
    glVertex2f(1.0*WINDOWS_WIDTH/WINDOWS_HEIGHT,-1.0);
    glVertex2f(1.0*WINDOWS_WIDTH/WINDOWS_HEIGHT,1.0);
    glVertex2f(-1.0*WINDOWS_WIDTH/WINDOWS_HEIGHT,1.0);

    glEnd();

    for (int i=0;i<c.atoms.size();++i)
    {
        a=c.atoms[i];
        if (a.symbol=="") continue;
        for (int i=0;i<a.connections.size();++i)
        {
            if (a.connections[i].to!=-1)
            {
                a2=c.atoms[a.connections[i].to];
                drawConnection(a.x,a.y,a2.x,a2.y,a.connections[i].spots_taken.size());
            }
        }
    }

    for (int i=0;i<c.atoms.size();++i)
    {
        a=c.atoms[i];
        if (a.symbol=="") continue;
        drawAtom(a);
    }
}

void drawWindow(GLFWwindow* w,compound& c)
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
    //cout<<pressed<<" "<<mxpos<<" "<<mypos<<'\n';
}
void snap(double& x, double& y)
{
    double hx,lx;
    //cout<<x<<" "<<y<<endl;
    hx=(round(x*3)+0.2)/3;
    lx=(round(x*3)-1+0.2)/3;
    //cout<<hx<<" > "<<x<<" > "<<lx<<endl;
    if (hx-x<x-lx) x=hx;
    else x=lx;
    y=round(y*3.7)/3.7;
    //cout<<x<<" "<<y<<endl;
}
void run(GLFWwindow* w)
{
    string carbonSymbol="C";
    double sx,sy;
    sx=0;
    sy=0;
    if (snappingEnabled) snap(sx,sy);
    compound c(carbonSymbol,4,sx,sy);
    compound c_old=c;
    int last2=-1;
    int last=-1;
    bool toMove=0;
    int result;
    while (!glfwWindowShouldClose(w))
    {
        drawWindow(w,c);
        if (pressed==-3)
        {
            last2=-1;
            last=-1;
            toMove=0;
            swap(c,c_old);
        }
        if (pressed==-5)
        {
            c_old=c;
            last2=-1;
            last=-1;
            toMove=0;
            sx=0;
            sy=0;
            if (snappingEnabled) snap(sx,sy);
            c=*(new compound(carbonSymbol,4,sx,sy));
        }
        last=c.findAtom(mxpos,mypos);
        if (snappingEnabled) snap(mxpos,mypos);
        if (last==-1 && c.findAtom(mxpos,mypos)!=-1) continue;
        if (pressed==GLFW_MOUSE_BUTTON_LEFT || pressed==GLFW_MOUSE_BUTTON_MIDDLE)
        {
            //cout<<"Connect: "<<last<<" "<<last2<<'\n'<<endl;
            if (last2!=-1)
            {
                c_old=c;
                if (toMove)
                {
                    result=c.moveAtom(last2,mxpos,mypos);
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
                        if (pressed==GLFW_MOUSE_BUTTON_LEFT || c.atoms[last].free_connections.empty()) last2=-1;
                        else last2=last;
                    }
                    last=-1;

                }
                else
                {
                    result=c.addAtom(carbonSymbol,4,mxpos,mypos,last2);
                    if (result!=-1)
                    {
                        if (pressed==GLFW_MOUSE_BUTTON_LEFT || c.atoms[result].free_connections.empty()) last2=-1;
                        else last2=result;
                    }
                    last=-1;
                }
            }
            else
            {
                if (last==-1)
                {
                    c_old=c;
                    result=c.addAtom(carbonSymbol,4,mxpos,mypos);
                    if (result!=-1)
                    {
                        if (pressed==GLFW_MOUSE_BUTTON_LEFT || c.atoms[result].free_connections.empty()) last2=-1;
                        else last2=result;
                    }
                }
                else
                {
                    if (pressed==GLFW_MOUSE_BUTTON_LEFT) toMove=1;
                    else toMove=0;
                    if (!c.atoms[last].free_connections.empty())
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
                c_old=c;
                result=c.removeAtom(last);
            }
        }
    }
}
int main()
{
    system("cls");

    cout<<"Use the Middle Mouse Button to start or continue chains.\nUse the Left Mouse Button to end chains and to move atoms,\nUse the right mouse button to cancel chains and remove atoms.\nPress Backspace to undo.\nPress Shift to toggle snapping on and off.\nPress R to reset.\nPress Escape to end the program.\nPress any key to continue."<<endl;
    getch();

    system("cls");

    string message;
    message=initializeGLFW();
    cout<<message<<endl;
    if (message!="GLFW initialized successfully.") return -1;

    message=createWindow(window);
    cout<<message<<endl;
    if (message!="Window created successfully.") return -1;

    message=setCallbacks(window);
    cout<<message<<endl;
    if (message!="Callbacks set successfully.") return -1;

    //system("cls");
    run(window);

    stopGraphics();
    return 0;
}
