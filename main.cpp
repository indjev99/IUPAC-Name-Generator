//IUPAC Name Generator
#include<iostream>
#include<algorithm>
#include<vector>
#include<deque>
#include<stack>
#include<stdlib.h>
#include<conio.h>
#include<math.h>
#include<unordered_map>
#include<GLFW/glfw3.h>
using namespace std;

const int ORIGINAL_WINDOWS_WIDTH=1100;
const int ORIGINAL_WINDOWS_HEIGHT=1050;
int WINDOWS_WIDTH=ORIGINAL_WINDOWS_WIDTH;
int WINDOWS_HEIGHT=ORIGINAL_WINDOWS_HEIGHT;
const double DEG2RAD=3.14159/180.0;
double TEXT_COLOUR_R=0;
double TEXT_COLOUR_G=0;
double TEXT_COLOUR_B=0;
double BACKGROUND_COLOUR_R=1;
double BACKGROUND_COLOUR_R2=1;
double BACKGROUND_COLOUR_G=1;
double BACKGROUND_COLOUR_B=1;

int pressed;
double mxpos=0,mypos=0;

bool snappingEnabled=1;

struct dictionary
{
    vector<string> CNP; //carbon_number_prefixes
    vector<string> NP; //number_prefixes
    vector<string> CBI; //complex_bond_infixes
    vector<string> SS; //substituent_suffix
    string CP; //cycle_prefix
    string NC; //not_connected
    string error;
    string help;
    string PACTC; //Press any key to continue
    unordered_map<string, string> HP; //halogen_prefixes
    unordered_map<string, string> FGP; //functional_group_prefixes
    unordered_map<string, string> FGS; //functional_group_suffixes
    string benzene;
    string phen;

    string getCNP(int CN)
    {
        if (CN>=CNP.size())
        {
            if (CN>=NP.size()) CN=0;
            else return NP[CN];
        }
        return CNP[CN];
    }
    string getNP(int ST)
    {
        if (ST>=NP.size()) ST=0;
        return NP[ST];
    }
};

vector<string> element_symbol={"H","C","O","F","Cl","Br","I"};
vector<int> element_valence={1,4,2,1,1,1,1};
vector<string> halogen_symbol={"F","Cl","Br","I"};
unordered_map<string, int> substituent_priorities={{"alkyl",-1},{"alkyl halide",-1},{"alcohol",1},{"ketone",2},{"aldehyde",3},{"carboxylic acid",4},{"attachment",100}};

dictionary Bulgarian,English,curr_dict;
vector<dictionary> dictionaries;
int curr_dict_N;

string onlyLetters(string a)
{
    string a2;
    for (int i=0;i<a.size();++i)
    {
        if ((a[i]>='a' && a[i]<='z') || (a[i]>='а' && a[i]<='я'))
        {
            a2+=a[i];
        }
    }
    //cerr<<"OL: "<<a<<" "<<a2<<endl;
    return a2;
}
bool cmpBySubName(tuple<int, string, bool> a, tuple<int, string, bool> b)
{
    string a2,b2;
    a2=onlyLetters(get<1>(a));
    b2=onlyLetters(get<1>(b));
    return a2<b2;
}

struct bond
{
    vector<int> spots_taken;
    int to;
};
struct atom
{
    string symbol;
    double x,y;
    vector<bond> bonds;
    stack<int> free_bonds;

    int connect(int bond)
    {
        if (!free_bonds.empty())
        {
            int nc=isConnected(bond);
            if (nc==-1)
            {
                nc=free_bonds.top();
            }
            else if (symbol=="C" && bonds[nc].spots_taken.size()==bonds.size()-1) return -1;
            bonds[nc].to=bond;
            bonds[nc].spots_taken.push_back(free_bonds.top());
            free_bonds.pop();
            return 1;
        }
        return -1;
    }
    bool canConnect(int bond)
    {
        if (!free_bonds.empty())
        {
            int nc=isConnected(bond);
            if (symbol=="C" && nc>=0 && bonds[nc].spots_taken.size()==bonds.size()-1) return 0;
            if (symbol=="O" && bond!=-1 && nc==-1 && free_bonds.size()<2) return 0;
            return 1;
        }
        return 0;
    }
    int isConnected(int bond)
    {
        for (int i=0;i<bonds.size();++i)
        {
            if (bonds[i].to==bond) return i;
        }
        return -1;
    }
    void removeBond(int bond)
    {
        for (int i=0;i<bonds.size();++i)
        {
            if (bonds[i].to==bond)
            {
                bonds[i].to=-1;
                for (int j=0;j<bonds[i].spots_taken.size();++j)
                {
                    free_bonds.push(bonds[i].spots_taken[j]);
                }
                bonds[i].spots_taken.resize(0);
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
        bonds.resize(new_valance);
        for (int i=new_valance-1;i>=0;--i)
        {
            bonds[i].to=-1;
            free_bonds.push(i);
        }
    }
    atom(string new_symbol, int new_valance, double new_x, double new_y, int new_bond)
    {
        symbol=new_symbol;
        x=new_x;
        y=new_y;
        bonds.resize(new_valance);
        for (int i=new_valance-1;i>=0;--i)
        {
            bonds[i].to=-1;
            free_bonds.push(i);
        }
        connect(new_bond);
    }
};
struct compound
{
    vector<atom> atoms;
    stack<int> free_positions;
    string name;
    bool changed;
    int addAtom(atom& a)
    {
        changed=1;
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
        changed=1;
        if (new_valance<0) return -1;
        atom a(new_symbol,new_valance,new_x,new_y);
        return addAtom(a);
    }
    int addAtom(string new_symbol, int new_valance, double new_x, double new_y, int new_bond)
    {
        changed=1;
        int ind=-1;
        if (new_valance<=0) return -1;
        if (new_bond>=atoms.size() || atoms[new_bond].symbol=="") return -1;
        atom a(new_symbol,new_valance,new_x,new_y,new_bond);
        if (atoms[new_bond].canConnect(-2))
        {
            ind=addAtom(a);
            if (ind!=-1) atoms[new_bond].connect(ind);
        }
        return ind;
    }
    int connectAtoms(int a1, int a2)
    {
        int a=findAtomInCycle(-1,-1);
        vector<int> prev;
        changed=1;
        if (a1<atoms.size() && a2<atoms.size() && a1!=a2 && atoms[a1].symbol!="" && atoms[a2].symbol!="" && atoms[a1].canConnect(a2) && atoms[a2].canConnect(a1))
        {
            if (a!=-1)
            {
                prev=findPathFrom(a1,-1,0);
                if (prev[a2]!=-2 && atoms[a1].isConnected(a2)==-1)
                {
                    return -1;
                }
            }
            atoms[a1].connect(a2);
            atoms[a2].connect(a1);
            return 1;
        }
        return -1;
    }
    int removeAtom(int a)
    {

        changed=1;
        atom a1;
        if (a<atoms.size())
        {
            a1=atoms[a];
            if (a1.symbol!="")// && a1.bonds.size()-a1.free_bonds.size()<=1)
            {
                for (int i=0;i<a1.bonds.size();++i)
                {
                    if (a1.bonds[i].to!=-1) atoms[a1.bonds[i].to].removeBond(a);
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
            if (x<a.x-0.05 || (x>a.x+0.05 && a.free_bonds.empty()) || x>a.x+0.175 || y<a.y-0.08 || y>a.y+0.08) minDistInd=-1;
        }
        return minDistInd;
    }
    compound()
    {
        changed=1;
    }
    compound(string new_symbol, int new_valance, double new_x, double new_y)
    {
        changed=1;
        addAtom(new_symbol,new_valance,new_x,new_y);
    }


    string addSuffix(string name, string suffix)
    {
        bool first_letter_vowel=0;
        for (int i=0;i<suffix.size();++i)
        {
            if ((suffix[i]>='a' && suffix[i]<='z') || (suffix[i]>='а' && suffix[i]<='я'))
            {
                if (suffix[i]=='a' || suffix[i]=='e' || suffix[i]=='i' || suffix[i]=='o'
                    || suffix[i]=='u' || suffix[i]=='y' || suffix[i]=='а' || suffix[i]=='ъ'
                    || suffix[i]=='о' || suffix[i]=='у' || suffix[i]=='е' || suffix[i]=='и'
                    || suffix[i]=='й' || suffix[i]=='ю' || suffix[i]=='я' || suffix[i]=='ь')
                {
                    first_letter_vowel=1;
                }
                break;
            }
        }
        int ns=name.size()-1;
        if (first_letter_vowel)
        {
            if (name[ns]=='a' || name[ns]=='e' || name[ns]=='i' || name[ns]=='o'
                || name[ns]=='u' || name[ns]=='y' || name[ns]=='а' || name[ns]=='ъ'
                || name[ns]=='о' || name[ns]=='у' || name[ns]=='е' || name[ns]=='и'
                || name[ns]=='й' || name[ns]=='ю' || name[ns]=='я' || name[ns]=='ь')
            {
                name.resize(ns);
            }
        }
        if (suffix[0]>='0' && suffix[0]<='9' && name[name.size()-1]!='-') name+='-';
        name+=suffix;
        return name;
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

    bool isConnected()
    {
        atom a;
        stack<int> st;
        vector<bool> vis;
        int s=-1;
        vis.resize(atoms.size());
        for (int i=0;i<atoms.size();++i)
        {
            vis[i]=0;
            if (s==-1 && atoms[i].symbol=="C")
            {
                s=i;
                break;
            }
        }
        if (s==-1) return 0;
        st.push(s);
        vis[s]=1;
        while (!st.empty())
        {
            s=st.top();
            st.pop();
            a=atoms[s];
            for (int i=0;i<a.bonds.size();++i)
            {
                s=a.bonds[i].to;
                if (s!=-1 && !vis[s])
                {
                    if (atoms[s].symbol=="C") st.push(s);
                    vis[s]=1;
                }
            }
        }
        for (int i=0;i<atoms.size();++i)
        {
            if (atoms[i].symbol!="" && !vis[i]) return 0;
        }
        return 1;
    }
    int findAtomInCycle(int in, int out)
    {
        int vc;
        atom a;
        stack<int> st;
        vector<bool> vis;
        int s=-1,s2;
        vis.resize(atoms.size());
        if (in==-1)
        {
            for (int i=0;i<atoms.size();++i)
            {
                vis[i]=0;
                if (s==-1 && atoms[i].symbol=="C")
                {
                    s=i;
                    break;
                }
            }
        }
        if (s==-1) return -1;
        st.push(s);
        vis[s]=1;
        while (!st.empty())
        {
            s=st.top();
            st.pop();
            a=atoms[s];
            vc=0;
            for (int i=0;i<a.bonds.size();++i)
            {
                s2=a.bonds[i].to;
                if (s2!=-1 && atoms[s2].symbol=="C" && s2!=out)
                {
                    if (!vis[s2])
                    {
                        st.push(s2);
                        vis[s2]=1;
                    }
                    else ++vc;
                }
            }
            if (vc>1) return s;
        }
        return -1;
    }
    bool isHalogen(string s)
    {
        for (int i=0;i<halogen_symbol.size();++i)
        {
            if (s==halogen_symbol[i]) return 1;
        }
        return 0;
    }
    vector<string> findFunctionalGroups(atom a, int out)
    {
        vector<string> ans;
        if (out!=-1 && a.isConnected(out)!=-1) ans.push_back("attachment");
        int carbonyl=0;
        int hydroxyl=0;
        int carbon=0;
        for (int i=0;i<a.bonds.size();++i)
        {
            if (a.bonds[i].to==-1 || a.bonds[i].to==out) continue;
            if (atoms[a.bonds[i].to].symbol=="O")
            {
                if (atoms[a.bonds[i].to].free_bonds.empty())
                {
                    ++carbonyl;
                }
                else
                {
                    ++hydroxyl;
                }
            }
            if (atoms[a.bonds[i].to].symbol=="C")
            {
                ++carbon;
            }
        }

        //cerr<<"Carbonyl and hydroxyl: "<<carbonyl<<" "<<hydroxyl<<endl;

        while (carbonyl || hydroxyl)
        {
            if (carbonyl)
            {
                --carbonyl;
                if (hydroxyl) ans.push_back("carboxylic acid");
                //else if (!a.free_bonds.empty()) ans.push_back("aldehyde");
                else if (carbon<2) ans.push_back("aldehyde");
                else ans.push_back("ketone");
            }
            else if (hydroxyl) ans.push_back("alcohol");
            if (hydroxyl) --hydroxyl;
        }
        return ans;
    }
    pair<vector<int>, vector<int> > findFarthest(int s, int out)
    {
        vector<int> max_dist_atoms;
        pair<vector<int>, vector<int> > ans;
        vector<int> max_dist;
        atom a;
        stack<pair<int, vector<int> > > st;
        vector<bool> vis;
        pair<int, vector<int>> curr;
        int f;
        int priority;
        vector<string> groups;
        vis.resize(atoms.size());
        for (int i=0;i<vis.size();++i) vis[i]=0;
        curr.first=s;
        curr.second.resize(5);
        max_dist.resize(5);
        for (int i=0;i<5;++i)
        {
            curr.second[i]=0;
            max_dist[i]=0;
        }
        st.push(curr);
        vis[s]=1;
        while (!st.empty())
        {
            curr=st.top();
            st.pop();
            a=atoms[curr.first];
            --curr.second[2];
            priority=1;
            groups=findFunctionalGroups(a,out);
            if (groups.size()) priority=-substituent_priorities[groups[0]];
            if (priority<curr.second[0])
            {
                curr.second[0]=priority;
                curr.second[1]=0;
            }
            if (priority==curr.second[0])
            {
                --curr.second[1];
                for (int i=1;i<groups.size();++i)
                {
                    if (groups[i]==groups[i-1]) --curr.second[1];
                    else break;
                }
            }
            f=cmpVectors(curr.second,max_dist);
            if (f==1)
            {
                f=0;
                max_dist=curr.second;
                max_dist_atoms.resize(0);
            }
            if (f==0) max_dist_atoms.push_back(curr.first);
            for (int i=0;i<a.bonds.size();++i)
            {
                s=a.bonds[i].to;
                if (s!=-1 && atoms[s].symbol=="C" && !vis[s] && s!=out)
                {
                    curr.first=s;
                    vis[s]=1;
                    if (a.bonds[i].spots_taken.size()>1)
                    {
                        --curr.second[3];
                        if (a.bonds[i].spots_taken.size()==2) --curr.second[4];
                        st.push(curr);
                        ++curr.second[3];
                        if (a.bonds[i].spots_taken.size()==2) ++curr.second[4];
                    }
                    else
                    {
                        st.push(curr);
                    }
                }
            }
        }
        ans.first=max_dist_atoms;
        ans.second=max_dist;
        return ans;
    }
    vector<int> findPathFrom(int s, int out, bool cycle)
    {
        atom a;
        int s2;
        stack<pair<int, int> > st;
        pair<int, int> curr;
        vector<int> prev;
        prev.resize(atoms.size());
        for (int i=0;i<prev.size();++i) prev[i]=-2;
        st.push(make_pair(s,-1));
        while (!st.empty())
        {
            curr=st.top();
            st.pop();
            s=curr.first;
            if (prev[s]!=-2) continue;
            prev[s]=curr.second;
            a=atoms[s];
            for (int i=0;i<a.bonds.size();++i)
            {
                s2=a.bonds[i].to;
                if (s2!=-1 && atoms[s2].symbol=="C"&& s2!=out)
                {
                    if (prev[s2]==-2)
                    {
                        st.push(make_pair(s2,s));
                    }
                    else if (prev[s2]==-1 && cycle && prev[s]!=s2)
                    {
                        prev[s2]=s;
                    }
                }
            }
        }
        return prev;
    }
    vector<pair<int, int> > findComplexBonds(vector<int> parent_chain)
    {
        atom a;
        vector<pair<int, int> > comp_bonds;
        pair<int, int> comp_bond;
        //cerr<<"PC & CB: ";
        for (int i=0;i<parent_chain.size();++i)
        {
            //cerr<<" "<<parent_chain[i];
            a=atoms[parent_chain[i]];
            for (int j=0;j<a.bonds.size();++j)
            {
                if (a.bonds[j].to==parent_chain[(i+1)%parent_chain.size()] && a.bonds[j].spots_taken.size()>1 && (i<parent_chain.size()-1 || parent_chain.size()>2))
                {
                    //cerr<<"X"<<a.bonds[j].spots_taken.size();
                    comp_bond.first=i+1;
                    comp_bond.second=a.bonds[j].spots_taken.size();
                    comp_bonds.push_back(comp_bond);
                }
            }
        }
        //cerr<<endl;
        return comp_bonds;
    }
    int findHighestPriority(vector<int> parent_chain, int out)
    {
        //cerr<<"TO CALC MAX P"<<endl;
        atom a;
        int max_priority=0;
        int curr_p;
        int PCS=parent_chain.size();
        vector<string> groups;
        for (int i=0;i<PCS;++i)
        {
            a=atoms[parent_chain[i]];
            groups=findFunctionalGroups(a,out);
            for (int j=0;j<groups.size();++j)
            {
                curr_p=substituent_priorities[groups[j]];
                if (curr_p>max_priority) max_priority=curr_p;
            }
        }
        return max_priority;
    }
    vector<tuple<int, string, bool> > findSubstituents(vector<int> parent_chain, int out, bool prefix, bool want_names)
    {
        atom a;
        vector<tuple<int, string, bool> > subs;
        tuple<int, string, bool> sub;
        vector<string> groups;
        int PCS=parent_chain.size();
        int max_priority=findHighestPriority(parent_chain,out);
        //cerr<<"MAX P: "<<max_priority<<endl;
        int curr_p;
        for (int i=0;i<PCS;++i)
        {
            a=atoms[parent_chain[i]];
            for (int j=0;j<a.bonds.size();++j)
            {
                if (a.bonds[j].to!=-1 && a.bonds[j].to!=parent_chain[(i-1+PCS)%PCS] && a.bonds[j].to!=parent_chain[(i+1)%PCS])
                {
                    if (a.bonds[j].to==out);
                    else if (isHalogen(atoms[a.bonds[j].to].symbol))
                    {
                        curr_p=substituent_priorities["alkyl halide"];
                        if ((prefix && curr_p<max_priority) || (!prefix && curr_p>=max_priority))
                        {
                            get<0>(sub)=i+1;
                            if (want_names)
                            {
                                if (prefix) get<1>(sub)=curr_dict.HP[atoms[a.bonds[j].to].symbol];
                                else get<1>(sub)=curr_dict.error;
                            }
                            else get<1>(sub)="";
                            get<2>(sub)=0;
                            subs.push_back(sub);
                        }
                    }
                    else if (atoms[a.bonds[j].to].symbol=="C")
                    {
                        curr_p=substituent_priorities["alkyl"];
                        if ((prefix && curr_p<max_priority) || (!prefix && curr_p>=max_priority))
                        {
                            get<0>(sub)=i+1;
                            if (want_names)
                            {
                                if (prefix) get<1>(sub)=generateName(a.bonds[j].to,parent_chain[i]);
                                else get<1>(sub)=curr_dict.error;
                            }
                            else get<1>(sub)="";
                            get<2>(sub)=1;
                            subs.push_back(sub);
                        }
                    }
                }
            }
            groups=findFunctionalGroups(a,out);
            for (int j=0;j<groups.size();++j)
            {
                //cerr<<groups[j]<<" ";
                curr_p=substituent_priorities[groups[j]];
                if ((prefix && curr_p<max_priority) || (!prefix && curr_p>=max_priority))
                {
                    get<0>(sub)=i+1;
                    if (want_names)
                    {
                        if (groups[j]!="attachment")
                        {
                            if (prefix) get<1>(sub)=curr_dict.FGP[groups[j]];
                            else get<1>(sub)=curr_dict.FGS[groups[j]];
                        }
                        else
                        {
                            if (prefix) get<1>(sub)=curr_dict.error;
                            else get<1>(sub)=curr_dict.SS[a.bonds[a.isConnected(out)].spots_taken.size()];
                        }
                    }
                    else get<1>(sub)="";
                    get<2>(sub)=0;
                    subs.push_back(sub);
                }
            }
            //cerr<<"; ";
        }
        //cerr<<endl;
        return subs;
    }
    bool isParentChainCyclic(vector<int> parent_chain)
    {
        if (parent_chain.size()<3) return 0;
        atom a;
        a=atoms[parent_chain[parent_chain.size()-1]];
        if (a.isConnected(parent_chain[0])==-1) return 0;
        return 1;
    }
    string findGN(vector<int> pos, string name, bool carbon, int parent_chain_length, bool cyclic, bool most_important) //find group's name
    {
        string prefix="";
        bool show_locants=1;
        if (pos.size()==1 && pos[0]==1)
        {
            for (int i=1;i<=3;++i)
            {
                if (name==curr_dict.SS[i])
                {
                    show_locants=0;
                    break;
                }
            }
        }
        if (parent_chain_length==1)
        {
            show_locants=0;
        }
        if (most_important && pos.size()==1 && (cyclic || parent_chain_length<3))
        {
            show_locants=0;
        }
        if (most_important && name==curr_dict.FGS["carboxylic acid"] || name==curr_dict.FGS["aldehyde"])
        {
            show_locants=0;
        }
        if ((!cyclic && pos.size()==parent_chain_length-1) || pos.size()==parent_chain_length || (most_important && pos.size()==1 && parent_chain_length==3))
        {
            for (int i=2;i<=3;++i)
            {
                if (name==curr_dict.CBI[i])
                {
                    show_locants=0;
                    break;
                }
            }
        }
        if (parent_chain_length<4 && carbon)
        {
            show_locants=0;
        }
        if (show_locants)
        {
            prefix+='-';
            for (int i=0;i<pos.size();++i)
            {
                if (i>0) prefix+=',';
                prefix+=intToString(pos[i]);
            }
            prefix+='-';
        }
        if (pos.size()>1) prefix+=curr_dict.getNP(pos.size());
        prefix+=name;
        return prefix;
    }
    string findGNs(vector<tuple<int, string, bool> > subs, int parent_chain_length, bool cyclic, bool most_important)//find groups' names
    {
        string prefixes="";
        vector<int> curr_pos;
        if (subs.size()==0) return prefixes;
        sort(subs.begin(),subs.end(),cmpBySubName);
        string curr=get<1>(subs[0]);
        bool curr_carbon=get<2>(subs[0]);
        for (int i=0;i<subs.size();++i)
        {
            if (get<1>(subs[i])!=curr)
            {
                most_important=0;
                prefixes+=findGN(curr_pos,curr,curr_carbon,parent_chain_length,cyclic,most_important);
                curr_pos.resize(0);
                curr_carbon=get<2>(subs[i]);
                curr=get<1>(subs[i]);
            }
            curr_pos.push_back(get<0>(subs[i]));
        }
        prefixes+=findGN(curr_pos,curr,curr_carbon,parent_chain_length,cyclic,most_important);
        return prefixes;
    }
    int cmpVectors(vector<int> a, vector<int> b)
    {
        if (a.size()>b.size()) return 1;
        if (a.size()<b.size()) return -1;
        for (int i=0;i<a.size();++i)
        {
            if (a[i]<b[i])
            {
                return 1;
            }
            if (a[i]>b[i])
            {
                return -1;
            }
        }
        return 0;
    }
    vector<int> convertVector (vector<tuple<int, string, bool> > a)
    {
        vector<int> ans(a.size());
        for (int i=0;i<a.size();++i)
        {
            ans[i]=get<0>(a[i]);
        }
        return ans;
    }
    vector<int> convertVector (vector<pair<int,int> > a)
    {
        vector<int> ans(a.size());
        for (int i=0;i<a.size();++i)
        {
            ans[i]=a[i].first;
        }
        return ans;
    }
    vector<int> directAcyclicParentChain(vector<int> parent_chain, int out)
    {
        vector<int> parent_chain2;
        vector<pair<int, int> > complex_bonds1;
        vector<pair<int, int> > complex_bonds2;
        vector<tuple<int, string, bool> > subs1;
        vector<tuple<int, string, bool> > subs2;
        vector<int> a1,a2;
        int f;
        parent_chain2.resize(parent_chain.size());
        for (int i=0;i<parent_chain.size();++i)
        {
            parent_chain2[parent_chain.size()-1-i]=parent_chain[i];
        }

        subs1=findSubstituents(parent_chain,out,0,0);
        subs2=findSubstituents(parent_chain2,out,0,0);
        a1=convertVector(subs1);
        a2=convertVector(subs2);
        f=cmpVectors(a1,a2);
        if (f==1) return parent_chain;
        if (f==-1) return parent_chain2;

        complex_bonds1=findComplexBonds(parent_chain);
        complex_bonds2=findComplexBonds(parent_chain2);
        a1=convertVector(complex_bonds1);
        a2=convertVector(complex_bonds2);
        f=cmpVectors(a1,a2);
        if (f==1) return parent_chain;
        if (f==-1) return parent_chain2;

        for (int i=0;i<complex_bonds1.size();++i)
        {
            if (complex_bonds1[i].second==2) a1.push_back(complex_bonds1[i].first);
        }
        for (int i=0;i<complex_bonds2.size();++i)
        {
            if (complex_bonds2[i].second==2) a2.push_back(complex_bonds2[i].first);
        }
        f=cmpVectors(a1,a2);
        if (f==1) return parent_chain;
        if (f==-1) return parent_chain2;

        subs1=findSubstituents(parent_chain,out,1,0);
        subs2=findSubstituents(parent_chain2,out,1,0);
        a1=convertVector(subs1);
        a2=convertVector(subs2);
        f=cmpVectors(a1,a2);
        if (f==1) return parent_chain;
        if (f==-1) return parent_chain2;

        subs1=findSubstituents(parent_chain,out,1,1);
        subs2=findSubstituents(parent_chain2,out,1,1);
        sort(subs1.begin(),subs1.end(),cmpBySubName);
        sort(subs2.begin(),subs2.end(),cmpBySubName);
        a1=convertVector(subs1);
        a2=convertVector(subs2);
        f=cmpVectors(a1,a2);
        if (f==1) return parent_chain;
        if (f==-1) return parent_chain2;
        return parent_chain;
    }
    vector<int> directCyclicParentChain(vector<int> parent_chain, int out)
    {
        vector<vector<int> > parent_chains;
        vector<vector<int> > parent_chains2;
        vector<pair<int, int> > complex_bonds;
        vector<tuple<int, string, bool> > subs;
        vector<int> maxx,a1;
        int f;

        int PCS=parent_chain.size();
        parent_chains.resize(PCS*2);
        parent_chains[0]=parent_chain;
        parent_chains[1].resize(PCS);
        for (int i=0;i<PCS;++i) parent_chains[1][i]=parent_chain[PCS-1-i];
        for (int i=1;i<PCS;++i)
        {
            parent_chains[i*2].resize(PCS);
            parent_chains[i*2+1].resize(PCS);
            parent_chains[i*2][0]=parent_chains[i*2-2][PCS-1];
            parent_chains[i*2+1][0]=parent_chains[i*2-1][PCS-1];
            for (int j=1;j<PCS;++j)
            {
                parent_chains[i*2][j]=parent_chains[i*2-2][j-1];
                parent_chains[i*2+1][j]=parent_chains[i*2-1][j-1];
            }
        }

        parent_chains2.resize(0);
        maxx.resize(0);
        for (int i=0;i<parent_chains.size();++i)
        {
            parent_chain=parent_chains[i];
            subs=findSubstituents(parent_chain,out,0,0);
            a1=convertVector(subs);
            f=cmpVectors(a1,maxx);
            if (f==1)
            {
                f=0;
                parent_chains2.resize(0);
                maxx=a1;
            }
            if (f==0)
            {
                parent_chains2.push_back(parent_chain);
            }
        }

        parent_chains.resize(0);
        maxx.resize(0);
        for (int i=0;i<parent_chains2.size();++i)
        {
            parent_chain=parent_chains2[i];
            complex_bonds=findComplexBonds(parent_chain);
            a1=convertVector(complex_bonds);
            f=cmpVectors(a1,maxx);
            if (f==1)
            {
                f=0;
                parent_chains.resize(0);
                maxx=a1;
            }
            if (f==0)
            {
                parent_chains.push_back(parent_chain);
            }
        }

        parent_chains2.resize(0);
        maxx.resize(0);
        for (int i=0;i<parent_chains.size();++i)
        {
            parent_chain=parent_chains[i];
            complex_bonds=findComplexBonds(parent_chain);
            a1.resize(0);
            for (int i=0;i<complex_bonds.size();++i)
            {
                if (complex_bonds[i].second==2) a1.push_back(complex_bonds[i].first);
            }
            f=cmpVectors(a1,maxx);
            if (f==1)
            {
                f=0;
                parent_chains2.resize(0);
                maxx=a1;
            }
            if (f==0)
            {
                parent_chains2.push_back(parent_chain);
            }
        }

        parent_chains.resize(0);
        maxx.resize(0);
        for (int i=0;i<parent_chains2.size();++i)
        {
            parent_chain=parent_chains2[i];
            subs=findSubstituents(parent_chain,out,1,0);
            a1=convertVector(subs);
            f=cmpVectors(a1,maxx);
            if (f==1)
            {
                f=0;
                parent_chains.resize(0);
                maxx=a1;
            }
            if (f==0)
            {
                parent_chains.push_back(parent_chain);
            }
        }

        parent_chains2.resize(0);
        maxx.resize(0);
        for (int i=0;i<parent_chains.size();++i)
        {
            parent_chain=parent_chains[i];
            subs=findSubstituents(parent_chain,out,1,0);
            sort(subs.begin(),subs.end(),cmpBySubName);
            a1=convertVector(subs);
            f=cmpVectors(a1,maxx);
            if (f==1)
            {
                f=0;
                parent_chains2.resize(0);
                maxx=a1;
            }
            if (f==0)
            {
                parent_chains2.push_back(parent_chain);
            }
        }

        parent_chain=parent_chains2[0];
        return parent_chain;
    }
    vector<int> directParentChain(vector<int> parent_chain, int out)
    {
        if (isParentChainCyclic(parent_chain))
        {
            parent_chain=directCyclicParentChain(parent_chain,out);
        }
        else parent_chain=directAcyclicParentChain(parent_chain,out);
        return parent_chain;
    }
    vector<int> findParentChain(int in, int out)
    {
        pair<vector<int>, vector<int>> candidates,candidates2;
        vector<pair<int, int> > finalCandidates;
        pair<int, int> currCandidate;
        vector<int> parent_chain={0};
        vector<int> prev;
        vector<vector<int> > candidateParentChains;
        vector<vector<int> > candidateParentChains2;
        vector<tuple<int, string, bool> > subs;
        vector<pair<int, int> > complex_bonds;
        vector<int> double_bonds;
        vector<int> max_double_bonds;
        vector<int> maxx;
        vector<int> a1;
        int curr;
        int f;
        vector<int> maxDist;
        maxDist.resize(5);
        for (int i=0;i<5;++i)
        {
            maxDist[i]=0;
        }
        int starting_atom;
        starting_atom=findAtomInCycle(in,out);
        if (starting_atom==-1)
        {
            if (in==-1)
            {
                for (int i=0;i<atoms.size();++i)
                {
                    if (atoms[i].symbol=="C")
                    {
                        starting_atom=i;
                        break;
                    }
                }
            }
            else
            {
                starting_atom=in;
            }
            candidates=findFarthest(starting_atom,out);
            for (int i=0;i<candidates.first.size();++i)
            {
                candidates2=findFarthest(candidates.first[i],out);
                f=cmpVectors(candidates2.second,maxDist);
                if (f==1)
                {
                    f=0;
                    maxDist=candidates2.second;
                    finalCandidates.resize(0);
                }
                if (f==0)
                {
                    currCandidate.first=candidates.first[i];
                    for (int j=0;j<candidates2.first.size();++j)
                    {
                        //cerr<<"FC1: "<<candidates.first[i]<<" "<<candidates2.first[j]<<endl;
                        currCandidate.second=candidates2.first[j];
                        finalCandidates.push_back(currCandidate);
                    }
                }
            }
            for (int i=0;i<finalCandidates.size();++i)
            {
                if (finalCandidates[i].first>finalCandidates[i].second) swap(finalCandidates[i].first,finalCandidates[i].second);
            }
            sort(finalCandidates.begin(),finalCandidates.end());
            candidateParentChains2.resize(0);
            for (int i=0;i<finalCandidates.size();++i)
            {
                if (i==0 || finalCandidates[i].first!=finalCandidates[i-1].first)
                {
                    prev=findPathFrom(finalCandidates[i].first,out,0);
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
                    candidateParentChains2.push_back(directParentChain(parent_chain,out));
                }
            }

            candidateParentChains.resize(0);
            maxx.resize(0);
            for (int i=0;i<candidateParentChains2.size();++i)
            {
                parent_chain=candidateParentChains2[i];
                subs=findSubstituents(parent_chain,out,0,0);
                a1=convertVector(subs);
                f=cmpVectors(a1,maxx);
                if (f==1)
                {
                    f=0;
                    candidateParentChains.resize(0);
                    maxx=a1;
                }
                if (f==0)
                {
                    candidateParentChains.push_back(parent_chain);
                }
            }

            candidateParentChains2.resize(0);
            maxx.resize(0);
            for (int i=0;i<candidateParentChains.size();++i)
            {
                parent_chain=candidateParentChains[i];
                complex_bonds=findComplexBonds(parent_chain);
                a1=convertVector(subs);
                f=cmpVectors(a1,maxx);
                if (f==1)
                {
                    f=0;
                    candidateParentChains2.resize(0);
                    maxx=a1;
                }
                if (f==0)
                {
                    candidateParentChains2.push_back(parent_chain);
                }
            }

            candidateParentChains.resize(0);
            maxx.resize(0);
            for (int i=0;i<candidateParentChains2.size();++i)
            {
                parent_chain=candidateParentChains2[i];
                complex_bonds=findComplexBonds(parent_chain);
                a1.resize(0);
                for (int i=0;i<complex_bonds.size();++i)
                {
                    if (complex_bonds[i].second==2) a1.push_back(complex_bonds[i].first);
                }
                f=cmpVectors(a1,maxx);
                if (f==1)
                {
                    f=0;
                    candidateParentChains.resize(0);
                    maxx=a1;
                }
                if (f==0)
                {
                    candidateParentChains.push_back(parent_chain);
                }
            }

            candidateParentChains2.resize(0);
            maxx.resize(0);
            for (int i=0;i<candidateParentChains.size();++i)
            {
                parent_chain=candidateParentChains[i];
                subs=findSubstituents(parent_chain,out,1,0);
                a1=convertVector(subs);
                f=cmpVectors(a1,maxx);
                if (f==1)
                {
                    f=0;
                    candidateParentChains2.resize(0);
                    maxx=a1;
                }
                if (f==0)
                {
                    candidateParentChains2.push_back(parent_chain);
                }
            }

            candidateParentChains.resize(0);
            maxx.resize(0);
            for (int i=0;i<candidateParentChains2.size();++i)
            {
                parent_chain=candidateParentChains2[i];
                subs=findSubstituents(parent_chain,out,1,1);
                sort(subs.begin(),subs.end(),cmpBySubName);
                a1=convertVector(subs);
                f=cmpVectors(a1,maxx);
                if (f==1)
                {
                    f=0;
                    candidateParentChains.resize(0);
                    maxx=a1;
                }
                if (f==0)
                {
                    candidateParentChains.push_back(parent_chain);
                }
            }
            parent_chain=candidateParentChains[0];
        }
        else
        {
            parent_chain.resize(0);
            if (in!=-1)
            {
                parent_chain.push_back(in);
                return parent_chain;
            }
            prev=findPathFrom(starting_atom,out,1);
            curr=starting_atom;
            do
            {
                parent_chain.push_back(curr);
                curr=prev[curr];
            }
            while (curr!=starting_atom);
        }
        return parent_chain;
    }
    string generateName(int in, int out)
    {
        //cerr<<"IN/OUT: "<<in<<" "<<out<<endl;
        string name="";
        vector<string> suffixes;
        string suffix;
        string prefix;
        vector<int> parent_chain;
        vector<pair<int, int> > complex_bonds;
        vector<int> double_bonds;
        vector<int> triple_bonds;
        vector<tuple<int, string, bool> > prefix_subs;
        vector<tuple<int, string, bool> > suffix_subs;
        bool cyclic;
        bool most_important;
        bool benzene;
        parent_chain=findParentChain(in,out);
        parent_chain=directParentChain(parent_chain,out);
        prefix_subs=findSubstituents(parent_chain,out,1,1);
        suffix_subs=findSubstituents(parent_chain,out,0,1);
        complex_bonds=findComplexBonds(parent_chain);
        cyclic=isParentChainCyclic(parent_chain);
        for (int i=0;i<complex_bonds.size();++i)
        {
            //cerr<<"CB: "<<complex_bonds[i].first<<" "<<complex_bonds[i].second<<endl;
            if (complex_bonds[i].second==2) double_bonds.push_back(complex_bonds[i].first);
            else if (complex_bonds[i].second==3) triple_bonds.push_back(complex_bonds[i].first);
        }

        benzene=0;
        if (parent_chain.size()==6 && double_bonds.size()==3 && triple_bonds.size()==0 && cyclic)
        {
            benzene=1;
            for (int i=1;i<3;++i)
            {
                if (double_bonds[i]-double_bonds[i-1]!=2)
                {
                    benzene=0;
                    break;
                }
            }
        }

        most_important=0;
        if (suffix_subs.empty() && (benzene || complex_bonds.empty() || complex_bonds.size()==parent_chain.size() || (!cyclic && complex_bonds.size()==parent_chain.size()-1))) most_important=1;
        prefix=findGNs(prefix_subs,parent_chain.size(),cyclic,most_important);
        if (!prefix.empty() && prefix[0]=='-') prefix=prefix.substr(1,prefix.size()-1);

        name=prefix+name;

        if (benzene)
        {
            if (out==-1) name+=curr_dict.benzene;
            else name+=curr_dict.phen;
            goto isBenzene;
        }

        if (!double_bonds.empty())
        {
            most_important=0;
            if (suffix_subs.empty() && triple_bonds.empty()) most_important=1;
            suffix=findGN(double_bonds,curr_dict.CBI[2],0,parent_chain.size(),cyclic,most_important);
            suffixes.push_back(suffix);
        }

        if (!triple_bonds.empty())
        {
            most_important=0;
            if (suffix_subs.empty() && double_bonds.empty()) most_important=1;
            suffix=findGN(triple_bonds,curr_dict.CBI[3],0,parent_chain.size(),cyclic,most_important);
            suffixes.push_back(suffix);
        }

        if (complex_bonds.empty() && (out==-1 || get<0>(suffix_subs[0])!=1))
        {
            suffix=curr_dict.CBI[1];
            suffixes.push_back(suffix);
        }

        if (cyclic) name+=curr_dict.CP;

        name+=curr_dict.getCNP(parent_chain.size());

        isBenzene:

        suffix=findGNs(suffix_subs,parent_chain.size(),cyclic,1);
        suffixes.push_back(suffix);

        for (int i=0;i<suffixes.size();++i)
        {
            name=addSuffix(name,suffixes[i]);
        }

        if (out!=-1)
        {
            if (!prefix_subs.empty()) name="("+name+")";
            else
            {
                for (int i=0;i<name.size();++i)
                {
                    if (name[i]=='-')
                    {
                        name="("+name+")";
                        break;
                    }
                }
            }
        }
        //cerr<<"FOR IN/OUT: "<<in<<" "<<out<<" name: "<<name<<endl;
        return name;
    }
    string getName()
    {
        if (!changed) return name;
        int CS=isConnected(); //bond_status
        if (!CS) return curr_dict.NC;
        name=generateName(-1,-1);

        /*for (int i=0;i<name.size();++i)
        {
            if (name[i]>='a' && name[i]<='z')
            {
                name[i]+='A'-'a';
                break;
            }
            if (name[i]>='а' && name[i]<='я')
            {
                name[i]+='А'-'а';
                break;
            }
        }*/
        return name;
    }
    void setName(string name, double x, double y, double distx, double disty)
    {
        for (int i=0;i<atoms.size();++i) removeAtom(i);
        vector<tuple<vector<int>, int, int> > subs;
        vector<pair<vector<int>, int> > hal_subs;
        vector<tuple<int, int, int> > subs2;
        vector<pair<int, int> > hal_subs2;
        tuple<vector<int>, int, int> sub;
        pair<vector<int>, int> hal_sub;
        int parent_chain=1;
        vector<int> double_bonds;
        vector<int> triple_bonds;
        string curr,p;
        int curr_num;
        vector<int> nums;
        int last_num;
        bool was_last_num=0;
        int i2;
        for (int i=0;i<name.size();++i)
        {
            //cerr<<"i: "<<i<<" name[i]: "<<name[i]<<endl;
            if (name[i]>='0' && name[i]<='9')
            {
                //cerr<<"CHISLO E"<<endl;
                curr_num=0;
                while (i<name.size() && name[i]>='0' && name[i]<='9')
                {
                    curr_num*=10;
                    curr_num+=name[i]-'0';
                    ++i;
                }
                nums.push_back(curr_num);
                //cerr<<"curr_num: "<<curr_num<<endl;
                was_last_num=1;
            }
            else
            {
                //cerr<<"BUKVA E"<<endl;
                if (i<=name.size()-2)
                {
                    //cerr<<"IMA POVECHE OT 2 BUKVI"<<endl;
                    if (name.substr(i,2)==curr_dict.CBI[1].substr(0,2))
                    {
                        //cerr<<"ALKANE"<<endl;
                        parent_chain=last_num;
                        ++i;
                        if (i<name.size()-1 && curr_dict.CBI[1].size()>2 && name[i+1]==curr_dict.CBI[1][2]) ++i;
                    }
                    else if (name.substr(i,2)==curr_dict.CBI[2].substr(0,2))
                    {
                        //cerr<<"ALKENE"<<endl;
                        parent_chain=last_num;
                        if (nums.empty()) nums.push_back(1);
                        double_bonds=nums;
                        nums.resize(0);
                        ++i;
                        if (i<name.size()-1 && curr_dict.CBI[1].size()>2 && name[i+1]==curr_dict.CBI[1][2]) ++i;
                    }
                    else if (name.substr(i,2)==curr_dict.CBI[3].substr(0,2))
                    {
                        //cerr<<"ALKYNE"<<endl;
                        parent_chain=last_num;
                        if (nums.empty()) nums.push_back(1);
                        triple_bonds=nums;
                        nums.resize(0);
                        ++i;
                        if (i<name.size()-1 && curr_dict.CBI[1].size()>2 && name[i+1]==curr_dict.CBI[1][2]) ++i;
                    }
                    else if (name.substr(i,2)==curr_dict.SS[1].substr(0,2))
                    {
                        //cerr<<"ALKYL(IDEN)"<<endl;
                        if (nums.empty()) nums.push_back(0);
                        ++i;
                        if (i<name.size()-5 && name.substr(i+1,4)==curr_dict.SS[2].substr(2,4))
                        {
                            //cerr<<"ALKYLIDENE"<<endl;
                            sub=make_tuple(nums,2,last_num);
                            i+=4;
                            if (i<name.size()-1 && curr_dict.SS[2].size()>6 && name[i+1]==curr_dict.CBI[1][6]) ++i;
                        }
                        else
                        {
                            //cerr<<"ALKYL"<<endl;
                            sub=make_tuple(nums,1,last_num);
                        }
                        subs.push_back(sub);
                        nums.resize(0);
                    }
                    else if (was_last_num && nums.size()>1)
                    {
                        //cerr<<"Last was number"<<endl;
                        i+=curr_dict.getNP(nums.size()).size()-1;
                    }
                    else if (name[i]=='a' || name[i]=='а')
                    {
                        //cerr<<"TVA E A"<<endl;
                        continue;
                    }
                    else
                    {

                        //cerr<<"Tva e korena"<<endl;
                        curr="";
                        curr_num=0;
                        i2=i;
                        while (i2<name.size() && name[i2]!='-' && name[i2]!=',')
                        {
                            curr+=name[i2];
                            for (int j=1;j<halogen_symbol.size();++j)
                            {
                                p=curr_dict.HP[halogen_symbol[j]];
                                if (curr==p)
                                {
                                    curr_num=j;
                                    break;
                                }
                            }
                            if (curr_num!=0) break;
                            ++i2;
                        }
                        if (curr_num!=0)
                        {
                            i=i2;
                            if (nums.empty()) nums.push_back(1);
                            hal_subs.push_back(make_pair(nums,curr_num));
                            nums.resize(0);
                            continue;
                        }

                        curr="";
                        curr_num=0;
                        while (i<name.size() && name[i]!='-' && name[i]!=',')
                        {
                            curr+=name[i];
                            for (int j=1;j<curr_dict.NP.size()+1;++j)
                            {
                                p=curr_dict.getCNP(j);
                                if (curr==p.substr(0,p.size()-1))
                                {
                                    curr_num=j;
                                    break;
                                }
                            }
                            if (curr_num!=0) break;
                            ++i;
                        }
                        //cerr<<curr_num<<endl;
                        if (curr_num!=0) last_num=curr_num;
                    }
                }
                was_last_num=0;
            }
        }
        int pr;
        for (int i=0;i<subs.size();++i)
        {
            sub=subs[i];
            nums=get<0>(sub);
            pr=get<1>(sub);
            curr_num=get<2>(sub);
            for (int j=0;j<nums.size();++j)
            {
                //cerr<<nums[j]<<" "<<pr<<" "<<curr_num<<endl;
                if (!nums[j]) nums[j]=min(parent_chain,2);
                subs2.push_back(make_tuple(nums[j],pr,curr_num));
            }
        }
        sort(subs2.begin(),subs2.end());
        for (int i=0;i<hal_subs.size();++i)
        {
            hal_sub=hal_subs[i];
            nums=hal_sub.first;
            curr_num=hal_sub.second;
            for (int j=0;j<nums.size();++j)
            {
                hal_subs2.push_back(make_pair(nums[j],curr_num));
            }
        }
        sort(hal_subs2.begin(),hal_subs2.end());
        //cerr<<"PC: "<<parent_chain<<endl;
        vector<int> PC;
        vector<int> been_up;
        pr=0;
        for (int i=0;i<parent_chain;++i)
        {
            curr_num=addAtom(element_symbol[1],element_valence[1],-(parent_chain/2)*distx+i*distx+x,y);
            if (i) connectAtoms(pr,curr_num);
            pr=curr_num;
            PC.push_back(curr_num);
            been_up.push_back(0);
        }

        for (int i=0;i<double_bonds.size();++i)
        {
            connectAtoms(PC[double_bonds[i]-1],PC[double_bonds[i]]);
        }

        for (int i=0;i<triple_bonds.size();++i)
        {
            connectAtoms(PC[triple_bonds[i]-1],PC[triple_bonds[i]]);
            connectAtoms(PC[triple_bonds[i]-1],PC[triple_bonds[i]]);
        }

        double x2,y2;

        for (int i=0;i<subs2.size();++i)
        {
            pr=PC[get<0>(subs2[i])-1];
            x2=-(parent_chain/2)*distx+(get<0>(subs2[i])-1)*distx+x;
            y2=y;
            if (been_up[get<0>(subs2[i])-1])
            {
                disty=-disty;
                been_up[get<0>(subs2[i])-1]=2;
            }
            for (int j=0;j<get<2>(subs2[i]);++j)
            {
                y2+=disty;
                curr_num=addAtom(element_symbol[1],element_valence[1],x2,y2);
                connectAtoms(curr_num,pr);
                if (!j && get<1>(subs2[i])==2) connectAtoms(curr_num,pr);
                pr=curr_num;
            }
            if (been_up[get<0>(subs2[i])-1]) disty=-disty;
            else been_up[get<0>(subs2[i])-1]=1;
        }

        for (int i=0;i<hal_subs2.size();++i)
        {
            //cerr<<"Halogen substituent: "<<hal_subs2[i].first<<" "<<halogen_symbol[hal_subs2[i].second]<<endl;
            //cerr<<been_up[hal_subs2[i].first-1]<<endl;
            pr=PC[hal_subs2[i].first-1];
            x2=-(parent_chain/2)*distx+(hal_subs2[i].first-1)*distx+x;
            y2=y;
            if (been_up[hal_subs2[i].first-1]<2)
            {
                if (been_up[hal_subs2[i].first-1])
                {
                    disty=-disty;
                    been_up[hal_subs2[i].first-1]=2;
                }
                y2+=disty;
                curr_num=addAtom(halogen_symbol[hal_subs2[i].second],1,x2,y2);
                connectAtoms(curr_num,pr);
                if (been_up[get<0>(hal_subs2[i])-1]) disty=-disty;
                else been_up[get<0>(hal_subs2[i])-1]=1;
            }
            else
            {
                if ((hal_subs2[i].first-1!=0 && hal_subs2[i].first-1!=parent_chain-1) || (parent_chain>1 && been_up[hal_subs2[i].first-1]==3)) y2-=disty;
                if (hal_subs2[i].first-1==0) distx=-distx;
                if (been_up[hal_subs2[i].first-1]>2)
                {
                    distx=-distx;
                    been_up[hal_subs2[i].first-1]=4;
                }
                x2+=distx;
                curr_num=addAtom(halogen_symbol[hal_subs2[i].second],1,x2,y2);
                connectAtoms(curr_num,pr);
                if (hal_subs2[i].first-1==0) distx=-distx;
                if (been_up[hal_subs2[i].first-1]>2) distx=-distx;
                else been_up[hal_subs2[i].first-1]=3;
            }

        }
        //getch();
    }
};

bool selected_element[16];

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
}
double drawSymbol(string symbol, double x, double y, bool centered)
{
    double scale=1.0*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT;
    //cerr<<symbol<<endl;
    double nextpos=x;
    for (int i=0;i<symbol.size();++i)
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
    for (int i=digits.size()-1;i>=0;--i)
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
            nextpos=drawSymbol(element_symbol[0],nextpos,a.y*scale,0);
            drawIndex(a.free_bonds.size(),nextpos-0.02*scale,a.y*scale-0.08*scale);
        }
        else
        {
            nextpos=drawSymbol(element_symbol[0],a.x*scale,a.y*scale,1);
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
    for (int i=0;i<num;++i)
    {
        offset=((num-1)*0.5-i)*0.02;
        glVertex2f((x1+offset*x)*scale,(y1+offset*y)*scale);
        glVertex2f((x2+offset*x)*scale,(y2+offset*y)*scale);
    }
    glEnd();
}
void drawCompound(GLFWwindow* w, compound& c)
{
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

    for (int i=0;i<c.atoms.size();++i)
    {
        a=c.atoms[i];
        if (a.symbol=="") continue;
        for (int i=0;i<a.bonds.size();++i)
        {
            if (a.bonds[i].to!=-1)
            {
                a2=c.atoms[a.bonds[i].to];
                drawBond(a.x,a.y,a2.x,a2.y,a.bonds[i].spots_taken.size());
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
    //cerr<<pressed<<" "<<mxpos<<" "<<mypos<<'\n';
}
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
void help()
{
    system("cls");
    cout<<curr_dict.help<<endl;
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
    deque<compound> history={c};
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
        for (int i=0;i<7;++i)
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
void setDictionaries()
{
    English.CNP={"alka","metha","etha","propa","buta"};
    English.CBI={English.error,"ane","ene","yne"};
    English.NP={English.error,"mono","di","tri","tetra","penta","hexa","hepta","octa","nona",
    "deca","undeca","dodeca","trideca","tetradeca","pentadeca","hexadeca","heptadeca","octadeca","nonadeca",
    "icosa","henicosa","docosa","tricosa","tetracosa","pentacosa","hexacosa","heptacosa","octacosa","nonacosa",
    "triaconta","hentriaconta","hentriaconta","tritriaconta"};
    English.SS={English.error,"yl","ylidene","ylidyne"};
    English.CP="cyclo";

    English.HP["F"]="fluoro";
    English.HP["Cl"]="chloro";
    English.HP["Br"]="bromo";
    English.HP["I"]="iodo";

    English.FGS["alcohol"]="ol";
    English.FGP["alcohol"]="hydroxy";
    English.FGS["ketone"]="one";
    English.FGP["ketone"]="oxo";
    English.FGS["aldehyde"]="al";
    English.FGP["aldehyde"]="oxo";
    English.FGS["carboxylic acid"]="oic acid";
    English.FGP["carboxylic acid"]="carboxy";

    English.benzene="benzene";
    English.phen="phen";

    English.NC="Not Connected";
    English.help="Use the Middle Mouse Button to start or continue chains.\nUse the Left Mouse Button to end chains and to move atoms.";
    English.help+="\nUse the Right mouse button to cancel the current action and remove atoms.\nHold down O, F, C, B or I in order to place aн Oxygen, Fluorine, Chlorine,\nBromine or Iodine atom respectively.";
    English.help+="\nPress Backspace to undo.\nPress Shift to toggle snapping on and off.\nPress R to reset.\nPress L to change the language.";
    English.help+="\nPress N to enter a name of a compound.\nPress H for help.\nPress Escape to end the program.";
    English.PACTC="Press any key to continue.";
    dictionaries.push_back(English);

    Bulgarian.CNP={"алка","мета","ета","пропа","бута"};
    Bulgarian.CBI={Bulgarian.error,"ан","ен","ин"};
    Bulgarian.NP={Bulgarian.error,"моно","ди","три","тетра","пента","хекса","хепта","окта","нона",
    "дека","ундека","додека","тридека","тетрадека","пентадека","хексадека","хептадека","октадека","нонадека"};
    Bulgarian.SS={Bulgarian.error,"ил","илиден","илиден"};
    Bulgarian.CP="цикло";

    Bulgarian.HP["F"]="флуоро";
    Bulgarian.HP["Cl"]="хлоро";
    Bulgarian.HP["Br"]="бромо";
    Bulgarian.HP["I"]="йодо";

    Bulgarian.FGS["alcohol"]="ол";
    Bulgarian.FGP["alcohol"]="хидрокси";
    Bulgarian.FGS["ketone"]="он";
    Bulgarian.FGP["ketone"]="оксо";
    Bulgarian.FGS["aldehyde"]="ал";
    Bulgarian.FGP["aldehyde"]="оксо";
    Bulgarian.FGS["carboxylic acid"]="ова киселина";
    Bulgarian.FGP["carboxylic acid"]="карбокси";

    Bulgarian.benzene="бензен";
    Bulgarian.phen="фен";

    Bulgarian.NC="Не са свързани";
    Bulgarian.help="Използвайте средния бутон на мишката, за да започнете или продължите вериги.\nИзползвайте левия бутон на мишката, за да завършите вериги\nили да местите атоми.";
    Bulgarian.help+="\nИзплозвайте десния бутон на мишката, за да откажете текущото действие\nили да махнете атом.\nЗадръжте O, F, C, B или I, за да сложите кислороден, флуорен, хлорен,\nбромен или йоден атом съответно.";
    Bulgarian.help+="\nНатиснете Backspace, за да върнете последното действие.\nНатиснете Shift, за да включите или изключите\nавтоматичното наместване на атомите.\nНатиснете R, за да рестартирате.";
    Bulgarian.help+="\nНатиснете L, за да смените езика.\nНатиснете H за инструкции.\nНатиснете N, за да въведете името на съединение.";
    Bulgarian.help+="\nНатиснете Escape, за да изключите програмата.";
    Bulgarian.PACTC="Натиснете който и да е клавиш, за да продължите.";
    dictionaries.push_back(Bulgarian);

    curr_dict_N=0;

    curr_dict=dictionaries[curr_dict_N];
}
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

    stopGraphics();
    return 0;
}
