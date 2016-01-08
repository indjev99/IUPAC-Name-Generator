/*
    Copyright (C) 2015 Emil Indzhev
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version. 
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//IUPAC Name Generator
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
double BACKGROUND_COLOUR_R2=1;
double BACKGROUND_COLOUR_G=1;
double BACKGROUND_COLOUR_B=1;

string Hsymbol="H";
int pressed;
double mxpos=0,mypos=0;

bool snappingEnabled=1;

struct dictionary
{
    vector<string> CNB; //carbon_number_bases
    vector<string> SNP; //substituent_number_prefixes
    vector<string> FGTS; //functional_group_type_suffixes
    string SS; //substituent_suffix
    string CP; //cycle_prefix
    string NC; //not_connected

    string getCNB(int CN)
    {
        if (CN>=CNB.size()) CN=0;
        return CNB[CN];
    }
    string getSNP(int ST)
    {
        if (ST>=SNP.size()) ST=0;
        return SNP[ST];
    }
};

dictionary English,curr_dict;

bool cmpBySubName(pair<int, int> a, pair<int, int> b)
{
    return curr_dict.getSNP(a.second)<curr_dict.getSNP(b.second);
}

struct connection
{
    vector<int> spots_taken;
    int to;
};
struct atom
{
    string symbol;
    double x,y;
    vector<connection> connections;
    stack<int> free_connections;

    int connect(int connection)
    {
        if (!free_connections.empty())
        {
            int nc=isConnected(connection);
            if (nc==-1)
            {
                nc=free_connections.top();
            }
            else if (symbol=="C" && connections[nc].spots_taken.size()==connections.size()-1) return -1;
            connections[nc].to=connection;
            connections[nc].spots_taken.push_back(free_connections.top());
            free_connections.pop();
        }
        return 1;
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
            return atoms[a2].connect(a1);
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

    int isConnected()
    {
        int cycle=0,vc;
        atom a;
        stack<int> st;
        vector<bool> vis;
        int s=-1;
        vis.resize(atoms.size());
        for (int i=0;i<atoms.size();++i)
        {
            vis[i]=0;
            if (s==-1 && atoms[i].symbol!="") s=i;
        }
        if (s==-1) return 0;
        st.push(s);
        vis[s]=1;
        while (!st.empty())
        {
            s=st.top();
            st.pop();
            a=atoms[s];
            vc=0;
            for (int i=0;i<a.connections.size();++i)
            {
                s=a.connections[i].to;
                if (s!=-1 && atoms[s].symbol!="")
                {
                    if (!vis[s])
                    {
                        st.push(s);
                        vis[s]=1;
                    }
                    else ++vc;
                }
            }
            if (vc>1) cycle=1;
        }
        for (int i=0;i<atoms.size();++i)
        {
            if (atoms[i].symbol!="" && !vis[i]) return 0;
        }
        return 1+cycle;
    }
    int findMaxCon()
    {
        atom a;
        int maxcon=1;
        for (int i=0;i<atoms.size();++i)
        {
            a=atoms[i];
            if (a.symbol=="C") for (int j=0;j<a.connections.size();++j)
            {
                if (a.connections[j].to!=-1 && atoms[a.connections[j].to].symbol=="C" && a.connections[j].spots_taken.size()>maxcon) maxcon=a.connections[j].spots_taken.size();
            }
        }
        return maxcon;
    }
    string getNameStupid(int CS)
    {
        string name="";
        if (CS==2) name="cyclo";
        name+=curr_dict.getCNB(atoms.size()-free_positions.size());
        name+=curr_dict.FGTS[findMaxCon()];
        return name;
    }
    pair<vector<int>, int> findFarthest(int s)
    {
        vector<int> maxDistAtoms;
        pair<vector<int>, int> ans;
        int maxDist=-1;
        atom a;
        stack<pair<int, int> > st;
        vector<bool> vis;
        pair<int, int> curr;
        vis.resize(atoms.size());
        for (int i=0;i<vis.size();++i) vis[i]=0;
        curr.first=s;
        curr.second=0;
        st.push(curr);
        vis[s]=1;
        while (!st.empty())
        {
            curr=st.top();
            st.pop();
            if (curr.second>maxDist)
            {
                maxDist=curr.second;
                maxDistAtoms.resize(0);
            }
            if (curr.second==maxDist) maxDistAtoms.push_back(curr.first);
            a=atoms[curr.first];
            ++curr.second;
            for (int i=0;i<a.connections.size();++i)
            {
                s=a.connections[i].to;
                if (s!=-1 && atoms[s].symbol!="" && !vis[s])
                {
                    curr.first=s;
                    vis[s]=1;
                    st.push(curr);
                }
            }
        }
        ans.first=maxDistAtoms;
        ans.second=maxDist;
        return ans;
    }
    vector<int> findPathFrom(int s)
    {
        atom a;
        int s2;
        stack<int> st;
        vector<int> prev;
        prev.resize(atoms.size());
        for (int i=0;i<prev.size();++i) prev[i]=-2;
        st.push(s);
        prev[s]=-1;
        while (!st.empty())
        {
            s=st.top();
            st.pop();
            a=atoms[s];
            for (int i=0;i<a.connections.size();++i)
            {
                s2=a.connections[i].to;
                if (s2!=-1 && atoms[s2].symbol!="" && prev[s2]==-2)
                {
                    prev[s2]=s;
                    st.push(s2);
                }
            }
        }
        return prev;
    }
    int findSubType(vector<int> parent_chain, int s)
    {
        int subSize=0;
        atom a;
        int s2;
        stack<int> st;
        vector<bool> vis;
        vis.resize(atoms.size());
        for (int i=0;i<vis.size();++i) vis[i]=0;
        for (int i=0;i<parent_chain.size();++i) vis[parent_chain[i]]=1;
        st.push(s);
        vis[s]=1;
        while (!st.empty())
        {
            s=st.top();
            st.pop();
            ++subSize;
            a=atoms[s];
            for (int i=0;i<a.connections.size();++i)
            {
                s=a.connections[i].to;
                if (s!=-1 && atoms[s].symbol!="" && !vis[s])
                {
                    vis[s]=1;
                    st.push(s);
                }
            }
        }
        return subSize;
    }
    vector<pair<int, int> > findSubstituents(vector<int> parent_chain)
    {
        atom a;
        vector<pair<int, int> > subs;
        pair<int, int> sub;
        for (int i=1;i<parent_chain.size()-1;++i)
        {
            a=atoms[parent_chain[i]];
            for (int j=0;j<a.connections.size();++j)
            {
                if (a.connections[j].to!=-1 && a.connections[j].to!=parent_chain[i-1] && a.connections[j].to!=parent_chain[i+1])
                {
                    sub.first=i+1;
                    sub.second=findSubType(parent_chain,a.connections[j].to);
                    subs.push_back(sub);
                }
            }
        }
        return subs;
    }
    string intToString(int num)
    {
        string ans2="",ans;
        do
        {
            ans2+=num%10+'0';
            num/=10;
        }
        while (num!=0);
        ans=ans2;
        for (int i=0;i<ans.size();++i)
        {
            ans[i]=ans2[ans2.size()-1-i];
        }
        return ans;
    }
    string findSP(vector<int> pos, int t) //find substituents_prefix
    {
        string prefix="";
        if (t==-1 || pos.size()==0) return prefix;
        for (int i=0;i<pos.size();++i)
        {
            if (i>0) prefix+=',';
            prefix+=intToString(pos[i]);
        }
        prefix+='-';
        if (pos.size()>1) prefix+=curr_dict.getSNP(pos.size());
        prefix+=curr_dict.getCNB(t);
        prefix+=curr_dict.SS;
        return prefix;
    }
    string findSPs(vector<pair<int, int> > subs)//find substituents_prefixes
    {
        string prefixes="";
        vector<int> currPos;
        if (subs.size()==0) return prefixes;
        sort(subs.begin(),subs.end(),cmpBySubName);
        int curr=subs[0].second;
        for (int i=0;i<subs.size();++i)
        {
            if (subs[i].second!=curr)
            {
                prefixes+=findSP(currPos,curr);
                prefixes+='-';
                currPos.resize(0);
                curr=subs[i].second;
            }
            currPos.push_back(subs[i].first);
        }
        prefixes+=findSP(currPos,curr);
        return prefixes;
    }
    vector<int> directParentChain(vector<int> parent_chain)
    {
        atom a;
        vector<int> parent_chain2;
        vector<pair<int, int> > subs1;
        vector<pair<int, int> > subs2;
        parent_chain2.resize(parent_chain.size());
        for (int i=0;i<parent_chain.size();++i)
        {
            parent_chain2[parent_chain.size()-1-i]=parent_chain[i];
        }
        subs1=findSubstituents(parent_chain);
        subs2=findSubstituents(parent_chain2);
        for (int i=0;i<subs1.size();++i)
        {
            if (subs1[i].first<subs2[i].first)
            {
                break;
            }
            if (subs1[i].first>subs2[i].first)
            {
                parent_chain=parent_chain2;
                break;
            }
        }
        return parent_chain;
    }
    vector<int> findParentChain()
    {
        atom a;
        pair<vector<int>, int> candidates,candidates2;
        vector<pair<int, int> > finalCandidates;
        pair<int, int> currCandidate;
        vector<int> parent_chain={0};
        vector<int> prev;
        vector<vector<int> > candidateParrentChains;
        vector<vector<int> > candidateParrentChains2;
        int maxSideChains=-1;
        int curr;
        int maxDist=-1;
        int starting_atom;
        for (int i=0;i<atoms.size();++i)
        {
            if (atoms[i].symbol!="")
            {
                starting_atom=i;
                break;
            }
        }
        candidates=findFarthest(starting_atom);
        for (int i=0;i<candidates.first.size();++i)
        {
            candidates2=findFarthest(candidates.first[i]);
            if (candidates2.second>maxDist)
            {
                maxDist=candidates2.second;
                finalCandidates.resize(0);
            }
            if (candidates2.second==maxDist)
            {
                currCandidate.first=candidates.first[i];
                for (int i=0;i<candidates2.first.size();++i)
                {
                    currCandidate.second=candidates2.first[i];
                    finalCandidates.push_back(currCandidate);
                }
            }
        }
        for (int i=0;i<finalCandidates.size();++i)
        {
            if (finalCandidates[i].first>finalCandidates[i].second) swap(finalCandidates[i].first,finalCandidates[i].second);
        }
        sort(finalCandidates.begin(),finalCandidates.end());
        for (int i=0;i<finalCandidates.size();++i)
        {
            if (i==0 || finalCandidates[i].first!=finalCandidates[i-1].first)
            {
                prev=findPathFrom(finalCandidates[i].first);
            }
            if (i==0 || finalCandidates[i].first!=finalCandidates[i-1].first || finalCandidates[i].second!=finalCandidates[i-1].second)
            {
                parent_chain.resize(0);
                curr=finalCandidates[i].second;
                while (curr!=-1)
                {
                    parent_chain.push_back(curr);
                    curr=prev[curr];
                }
                candidateParrentChains.push_back(parent_chain);
            }
        }
        int sideChains;
        for (int i=0;i<candidateParrentChains.size();++i)
        {
            sideChains=0;
            parent_chain=candidateParrentChains[i];
            for (int j=1;j<parent_chain.size()-1;++j)
            {
                a=atoms[parent_chain[j]];
                for (int o=0;o<a.connections.size();++o)
                {
                    if (a.connections[o].to!=-1 && a.connections[o].to!=parent_chain[j-1] && a.connections[o].to!=parent_chain[j+1]) ++sideChains;
                }
            }
            if (sideChains>maxSideChains)
            {
                maxSideChains=sideChains;
                candidateParrentChains2.resize(0);
            }
            if (sideChains==maxSideChains)
            {
                candidateParrentChains2.push_back(parent_chain);
            }
        }
        /*cout<<"Candidate parent chains: "<<candidateParrentChains2.size()<<endl;
        for (int i=0;i<candidateParrentChains2.size();++i)
        {
            parent_chain=candidateParrentChains2[i];
            for (int j=0;j<parent_chain.size();++j)
            {
                cout<<parent_chain[j]<<" ";
            }
            cout<<endl;
        }*/
        parent_chain=candidateParrentChains2[0];
        return parent_chain;
    }
    string getName1()
    {
        string name;
        vector<int> parent_chain;
        vector<pair<int, int> > subs;
        parent_chain=findParentChain();
        parent_chain=directParentChain(parent_chain);
        subs=findSubstituents(parent_chain);
        name=findSPs(subs);
        name+=curr_dict.getCNB(parent_chain.size());
        name+=curr_dict.FGTS[findMaxCon()];
        return name;
    }
    string getName()
    {
        string name;
        int CS=isConnected(); //connection_status
        if (!CS) return curr_dict.NC;
        if (CS==2) name=getNameStupid(CS);
        else name=getName1();

        for (int i=0;i<name.size();++i)
        {
            if (name[i]>='a' && name[i]<='z')
            {
                name[i]+='A'-'a';
                break;
            }
        }
        return name;
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
    w=glfwCreateWindow(WINDOWS_WIDTH,WINDOWS_HEIGHT,"IUPAC Name Generator",NULL,NULL);
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
    if (index<=1) return;
    vector<int> digits;
    do
    {
        digits.push_back(index%10);
        index/=10;
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
        if (digits[i]==9)
        {
            glBegin(GL_QUADS);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x,y-0.005*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.045*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.005*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.045*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);;

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

            glVertex2f(x,y);
            glVertex2f(x+0.04*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y);
            glVertex2f(x+0.04*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x+0.0125*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.045*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0275*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.045*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0275*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.07*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0125*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.07*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glVertex2f(x,y+0.01*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0275*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.01*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0275*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glEnd();

            x+=0.04*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;
        }
        if (digits[i]==4)
        {
            glBegin(GL_QUADS);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x,y-0.005*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.045*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.005*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.045*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);;

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

            glVertex2f(x,y);
            glVertex2f(x+0.04*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y);
            glVertex2f(x+0.04*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x+0.0125*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.045*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0275*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.045*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0275*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0125*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glVertex2f(x,y);
            glVertex2f(x+0.0275*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y);
            glVertex2f(x+0.0275*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glEnd();

            x+=0.04*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;
        }
        if (digits[i]==3)
        {
            glBegin(GL_QUADS);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x,y-0.005*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0375*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.005*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0375*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);;

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

            glVertex2f(x,y);
            glVertex2f(x+0.0325*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y);
            glVertex2f(x+0.0325*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x,y+0.045*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.02*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.045*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.02*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.07*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.07*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glVertex2f(x,y+0.01*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.02*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.01*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.02*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glEnd();

            x+=0.0325*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;
        }
        if (digits[i]==2)
        {
            glBegin(GL_QUADS);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x,y-0.005*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0375*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y-0.005*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0375*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.085*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);;

            glColor3f(TEXT_COLOUR_R,TEXT_COLOUR_G,TEXT_COLOUR_B);

            glVertex2f(x,y);
            glVertex2f(x+0.0325*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y);
            glVertex2f(x+0.0325*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

            glVertex2f(x,y+0.045*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.02*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.045*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.02*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.07*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x,y+0.07*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glVertex2f(x+0.0125*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.01*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0325*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.01*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0325*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
            glVertex2f(x+0.0125*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,y+0.035*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

            glEnd();

            x+=0.0325*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;
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

        glVertex2f(a.x*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,a.y*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT-0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(a.x*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT+0.055*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,a.y*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT-0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(a.x*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT+0.055*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,a.y*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
        glVertex2f(a.x*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,a.y*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT+0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);

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
    //glColor3f(BACKGROUND_COLOUR_R,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);
    glColor3f(BACKGROUND_COLOUR_R2,BACKGROUND_COLOUR_G,BACKGROUND_COLOUR_B);

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
        system("cls");
        cout<<c.getName()<<endl;
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
            //BACKGROUND_COLOUR_R2=!BACKGROUND_COLOUR_R2;
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
void setDictionaries()
{
    English.CNB={"alk","meth","eth","prop","but","pent","hex","hept","oct","non",
    "dec","undec","dodec","tridec","tetradec","pentadec","hexadec","heptadec","octadec","nonadec",
    "icosan","henicos","docos","tricos","tetracos","pentacos","hexacos","heptacos","octacos","nonacos",
    "triacont","hentriacont","hentriacont","tritriacont"};
    English.FGTS={"ERROR","ane","ene","yne"};
    English.SNP={"ERROR","mono","bi","tri","tetra","penta","hexa","hepta","octa","nona",
    "deca","undeca","dodeca","trideca","tetradeca","pentadeca","hexadeca","heptadeca","octadeca","nonadeca",
    "icosa","henicosa","docosa","tricosa"};
    English.SS="yl";
    English.CP="cyclo";
    English.NC="Not Connected";

    curr_dict=English;
}
int main()
{
    system("cls");

    setDictionaries();

    cout<<"Use the Middle Mouse Button to start or continue chains.\nUse the Left Mouse Button to end chains and to move atoms,"
    <<"\nUse the Right mouse button to cancel chains and remove atoms.\nPress Backspace to undo.\nPress Shift to toggle snapping on and off."
    <<"\nPress R to reset.\nPress Escape to end the program.\nPress any key to continue."<<endl;
    getch();

    system("cls");

    string message;
    message=initializeGLFW();
    //cout<<message<<endl;
    if (message!="GLFW initialized successfully.") return -1;

    message=createWindow(window);
    //cout<<message<<endl;
    if (message!="Window created successfully.") return -1;

    message=setCallbacks(window);
    //cout<<message<<endl;
    if (message!="Callbacks set successfully.") return -1;

    //system("cls");
    run(window);

    stopGraphics();
    return 0;
}
