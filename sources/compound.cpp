#include<iostream>
#include<tuple>
#include<algorithm>
#include "../headers/compound.h"
using namespace std;

const int halogen_N=4;
const string halogen_symbols[halogen_N]= {"F","Cl","Br","I"};
unordered_map<string, int> substituent_priorities= {{"alkyl",-1},{"alkyl halide",-1},{"alcohol",1},{"ketone",2},{"aldehyde",3},{"carboxylic acid",4},{"attachment",100},{"carbon dioxide",1000},{"carbonic acid",1000}};

compound::compound()
{
    changed=1;
}
compound::compound(dictionary new_dict)
{
    changed=1;
    curr_dict=new_dict;
}
compound::compound(dictionary new_dict, string new_symbol, int new_valance, double new_x, double new_y)
{
    changed=1;
    curr_dict=new_dict;
    addAtom(new_symbol,new_valance,new_x,new_y);
}

void compound::setDictionary(dictionary new_dict)
{
    changed=1;
    curr_dict=new_dict;
}

int compound::addAtom(atom& a)
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
int compound::addAtom(string new_symbol, int new_valance, double new_x, double new_y)
{
    changed=1;
    if (new_valance<0) return -1;
    atom a(new_symbol,new_valance,new_x,new_y);
    return addAtom(a);
}
int compound::addAtom(string new_symbol, int new_valance, double new_x, double new_y, int new_bond)
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
int compound::connectAtoms(int a1, int a2)
{
    int a=findAtomInCycle(-1,-1);
    vector<int> prev;
    changed=1;
    if (a1<atoms.size() && a2<atoms.size() && a1!=a2 && atoms[a1].symbol!="" && atoms[a2].symbol!="" && atoms[a1].canConnect(a2) && atoms[a2].canConnect(a1))
    {
        /*if (a!=-1)
        {
            prev=findPathFrom(a1,-1,0);
            if (prev[a2]!=-2 && atoms[a1].isConnected(a2)==-1)
            {
                return -1;
            }
        }*/
        atoms[a1].connect(a2);
        atoms[a2].connect(a1);
        return 1;
    }
    return -1;
}
int compound::removeAtom(int a)
{

    changed=1;
    atom a1;
    if (a<atoms.size())
    {
        a1=atoms[a];
        if (a1.symbol!="")// && a1.bonds.size()-a1.free_bonds.size()<=1)
        {
            for (int i=0; i<a1.bonds.size(); ++i)
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
int compound::moveAtom(int a, double x, double y)
{
    if (a<atoms.size() && atoms[a].symbol!="")
    {
        atoms[a].changeXY(x,y);
        return 1;
    }
    return -1;
}
int compound::findAtom(double x, double y)
{
    atom a;
    int minDistInd=-1;
    double minDist=-1;
    double dist;
    for (int i=0; i<atoms.size(); ++i)
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

string compound::addSuffix(string name, string suffix)
{
    bool first_letter_vowel=0;
    for (int i=0; i<suffix.size(); ++i)
    {
        if ((suffix[i]>='a' && suffix[i]<='z') || (suffix[i]>='à' && suffix[i]<='ÿ'))
        {
            if (suffix[i]=='a' || suffix[i]=='e' || suffix[i]=='i' || suffix[i]=='o'
                    || suffix[i]=='u' || suffix[i]=='y' || suffix[i]=='à' || suffix[i]=='ú'
                    || suffix[i]=='î' || suffix[i]=='ó' || suffix[i]=='å' || suffix[i]=='è'
                    || suffix[i]=='é' || suffix[i]=='þ' || suffix[i]=='ÿ' || suffix[i]=='ü')
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
                || name[ns]=='u' || name[ns]=='y' || name[ns]=='à' || name[ns]=='ú'
                || name[ns]=='î' || name[ns]=='ó' || name[ns]=='å' || name[ns]=='è'
                || name[ns]=='é' || name[ns]=='þ' || name[ns]=='ÿ' || name[ns]=='ü')
        {
            name.resize(ns);
        }
    }
    if (suffix[0]>='0' && suffix[0]<='9' && name[name.size()-1]!='-') name+='-';
    name+=suffix;
    return name;
}
string compound::onlyLetters(string a)
{
    string a2;
    for (int i=0; i<a.size(); ++i)
    {
        if ((a[i]>='a' && a[i]<='z') || (a[i]>='à' && a[i]<='ÿ'))
        {
            a2+=a[i];
        }
    }
    return a2;
}
bool compound::cmpBySubName(tuple<int, string, bool> a, tuple<int, string, bool> b)
{
    string a2,b2;
    a2=onlyLetters(get<1>(a));
    b2=onlyLetters(get<1>(b));
    return a2<b2;
}
string compound::intToString(int num)
{
    string ans2="",ans;
    do
    {
        ans2+=num%10+'0';
        num/=10;
    }
    while (num!=0);
    ans=ans2;
    for (int i=0; i<ans.size(); ++i)
    {
        ans[i]=ans2[ans2.size()-1-i];
    }
    return ans;
}
int compound::cmpVectors(vector<int> a, vector<int> b)
{
    if (a.size()>b.size()) return 1;
    if (a.size()<b.size()) return -1;
    for (int i=0; i<a.size(); ++i)
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
vector<int> compound::convertVector(vector<tuple<int, string, bool> > a)
{
    vector<int> ans(a.size());
    for (int i=0; i<a.size(); ++i)
    {
        ans[i]=get<0>(a[i]);
    }
    return ans;
}
vector<int> compound::convertVector(vector<pair<int,int> > a)
{
    vector<int> ans(a.size());
    for (int i=0; i<a.size(); ++i)
    {
        ans[i]=a[i].first;
    }
    return ans;
}

bool compound::isHalogen(string s)
{
    for (int i=0; i<halogen_N; ++i)
    {
        if (s==halogen_symbols[i]) return 1;
    }
    return 0;
}
string compound::findGN(vector<int> pos, string name, bool carbon, int parent_chain_length, bool cyclic, bool benzene, bool most_important) //find group's name
{
    string prefix="";
    bool show_locants=1;
    if (pos.size()==1 && pos[0]==1)
    {
        for (int i=1; i<=3; ++i)
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
    if (most_important && (name==curr_dict.FGS["carboxylic acid"] || name==curr_dict.FGS["aldehyde"]))
    {
        show_locants=0;
    }
    if ((!cyclic && pos.size()==parent_chain_length-1) || pos.size()==parent_chain_length || (most_important && pos.size()==1 && parent_chain_length==3))
    {
        for (int i=2; i<=3; ++i)
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
    if (benzene && pos.size()==parent_chain_length)
    {
        show_locants=0;
    }
    if (show_locants)
    {
        prefix+='-';
        for (int i=0; i<pos.size(); ++i)
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
string compound::findGNs(vector<tuple<int, string, bool> > subs, int parent_chain_length, bool cyclic, bool benzene, bool most_important)//find groups' names
{
    string prefixes="";
    vector<int> curr_pos;
    if (subs.size()==0) return prefixes;
    sort(subs.begin(),subs.end(),compound::cmpBySubName);
    string curr=get<1>(subs[0]);
    bool curr_carbon=get<2>(subs[0]);
    for (int i=0; i<subs.size(); ++i)
    {
        if (get<1>(subs[i])!=curr)
        {
            most_important=0;
            prefixes+=findGN(curr_pos,curr,curr_carbon,parent_chain_length,cyclic,benzene,most_important);
            curr_pos.resize(0);
            curr_carbon=get<2>(subs[i]);
            curr=get<1>(subs[i]);
        }
        curr_pos.push_back(get<0>(subs[i]));
    }
    prefixes+=findGN(curr_pos,curr,curr_carbon,parent_chain_length,cyclic,benzene,most_important);
    return prefixes;
}

bool compound::isConnected()
{
    atom a;
    stack<int> st;
    vector<bool> vis;
    int s=-1;
    vis.resize(atoms.size());
    for (int i=0; i<atoms.size(); ++i)
    {
        vis[i]=0;
        if (s==-1 && atoms[i].symbol!="")
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
        for (int i=0; i<a.bonds.size(); ++i)
        {
            s=a.bonds[i].to;
            if (s!=-1 && !vis[s])
            {
                if (atoms[s].symbol!="" ) st.push(s);
                vis[s]=1;
            }
        }
    }
    for (int i=0; i<atoms.size(); ++i)
    {
        if (atoms[i].symbol!="" && !vis[i]) return 0;
    }
    return 1;
}
int compound::findAtomInCycle(int in, int out)
{
    int vc;
    atom a;
    stack<int> st;
    vector<bool> vis;
    int s=-1,s2;
    vis.resize(atoms.size());
    if (in==-1)
    {
        for (int i=0; i<atoms.size(); ++i)
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
        for (int i=0; i<a.bonds.size(); ++i)
        {
            s2=a.bonds[i].to;
            if (s2!=-1 && atoms[s2].symbol=="C" && s2!=out && atoms_unions[s2]==atoms_unions[s])
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
vector<string> compound::findFunctionalGroups(atom a, int out)
{
    vector<string> ans;
    if (out!=-1 && a.isConnected(out)!=-1) ans.push_back("attachment");
    int carbonyl=0;
    int hydroxyl=0;
    int carbon=0;
    for (int i=0; i<a.bonds.size(); ++i)
    {
        if (a.bonds[i].to==-1 || a.bonds[i].to==out) continue;
        if (atoms[a.bonds[i].to].symbol=="O")
        {
            if (a.bonds[i].spots_taken.size()>1)
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

    if (carbonyl==2 && hydroxyl==0 && out==-1)
    {
        ans.push_back("carbon dioxide");
        return ans;
    }
    if (carbonyl==1 && hydroxyl==2 && out==-1)
    {
        ans.push_back("carbonic acid");
        return ans;
    }

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
pair<vector<int>, vector<int> > compound::findFarthest(int s, int out)
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
    for (int i=0; i<vis.size(); ++i) vis[i]=0;
    curr.first=s;
    curr.second.resize(5);
    max_dist.resize(5);
    for (int i=0; i<5; ++i)
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
            for (int i=1; i<groups.size(); ++i)
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
        for (int i=0; i<a.bonds.size(); ++i)
        {
            s=a.bonds[i].to;
            if (s!=-1 && atoms[s].symbol=="C" && !vis[s] && s!=out && atoms_unions[s]==atoms_unions[curr.first])
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
vector<int> compound::findPathFrom(int s, int out, bool cycle)
{
    atom a;
    int s2;
    stack<pair<int, int> > st;
    pair<int, int> curr;
    vector<int> prev;
    prev.resize(atoms.size());
    for (int i=0; i<prev.size(); ++i) prev[i]=-2;
    st.push(make_pair(s,-1));
    while (!st.empty())
    {
        curr=st.top();
        st.pop();
        s=curr.first;
        if (prev[s]!=-2) continue;
        prev[s]=curr.second;
        a=atoms[s];
        for (int i=0; i<a.bonds.size(); ++i)
        {
            s2=a.bonds[i].to;
            if (s2!=-1 && atoms[s2].symbol=="C" && s2!=out && atoms_unions[s2]==atoms_unions[s])
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
vector<pair<int, int> > compound::findComplexBonds(vector<int> parent_chain)
{
    atom a;
    vector<pair<int, int> > comp_bonds;
    pair<int, int> comp_bond;
    //cerr<<"PC & CB: ";
    for (int i=0; i<parent_chain.size(); ++i)
    {
        //cerr<<" "<<parent_chain[i];
        a=atoms[parent_chain[i]];
        for (int j=0; j<a.bonds.size(); ++j)
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
int compound::findHighestPriority(vector<int> parent_chain, int out)
{
    //cerr<<"TO CALC MAX P"<<endl;
    atom a;
    int max_priority=0;
    int curr_p;
    int PCS=parent_chain.size();
    vector<string> groups;
    for (int i=0; i<PCS; ++i)
    {
        a=atoms[parent_chain[i]];
        groups=findFunctionalGroups(a,out);
        for (int j=0; j<groups.size(); ++j)
        {
            curr_p=substituent_priorities[groups[j]];
            if (curr_p>max_priority) max_priority=curr_p;
        }
    }
    return max_priority;
}
vector<tuple<int, string, bool> > compound::findSubstituents(vector<int> parent_chain, int out, bool prefix, bool want_names)
{
    atom a;
    vector<tuple<int, string, bool> > subs;
    tuple<int, string, bool> sub;
    vector<string> groups;
    int PCS=parent_chain.size();
    int max_priority=findHighestPriority(parent_chain,out);
    //cerr<<"MAX P: "<<max_priority<<endl;
    int curr_p;
    for (int i=0; i<PCS; ++i)
    {
        a=atoms[parent_chain[i]];
        for (int j=0; j<a.bonds.size(); ++j)
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
        for (int j=0; j<groups.size(); ++j)
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
bool compound::isParentChainCyclic(vector<int> parent_chain)
{
    if (parent_chain.size()<3) return 0;
    atom a;
    a=atoms[parent_chain[parent_chain.size()-1]];
    if (a.isConnected(parent_chain[0])==-1) return 0;
    return 1;
}
vector<int> compound::directAcyclicParentChain(vector<int> parent_chain, int out)
{
    vector<int> parent_chain2;
    vector<pair<int, int> > complex_bonds1;
    vector<pair<int, int> > complex_bonds2;
    vector<tuple<int, string, bool> > subs1;
    vector<tuple<int, string, bool> > subs2;
    vector<int> a1,a2;
    int f;
    parent_chain2.resize(parent_chain.size());
    for (int i=0; i<parent_chain.size(); ++i)
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

    for (int i=0; i<complex_bonds1.size(); ++i)
    {
        if (complex_bonds1[i].second==2) a1.push_back(complex_bonds1[i].first);
    }
    for (int i=0; i<complex_bonds2.size(); ++i)
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
    sort(subs1.begin(),subs1.end(),compound::cmpBySubName);
    sort(subs2.begin(),subs2.end(),compound::cmpBySubName);
    a1=convertVector(subs1);
    a2=convertVector(subs2);
    f=cmpVectors(a1,a2);
    if (f==1) return parent_chain;
    if (f==-1) return parent_chain2;
    return parent_chain;
}
vector<int> compound::directCyclicParentChain(vector<int> parent_chain, int out)
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
    for (int i=0; i<PCS; ++i) parent_chains[1][i]=parent_chain[PCS-1-i];
    for (int i=1; i<PCS; ++i)
    {
        parent_chains[i*2].resize(PCS);
        parent_chains[i*2+1].resize(PCS);
        parent_chains[i*2][0]=parent_chains[i*2-2][PCS-1];
        parent_chains[i*2+1][0]=parent_chains[i*2-1][PCS-1];
        for (int j=1; j<PCS; ++j)
        {
            parent_chains[i*2][j]=parent_chains[i*2-2][j-1];
            parent_chains[i*2+1][j]=parent_chains[i*2-1][j-1];
        }
    }

    parent_chains2.resize(0);
    maxx.resize(0);
    for (int i=0; i<parent_chains.size(); ++i)
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
    for (int i=0; i<parent_chains2.size(); ++i)
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
    for (int i=0; i<parent_chains.size(); ++i)
    {
        parent_chain=parent_chains[i];
        complex_bonds=findComplexBonds(parent_chain);
        a1.resize(0);
        for (int i=0; i<complex_bonds.size(); ++i)
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
    for (int i=0; i<parent_chains2.size(); ++i)
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
    for (int i=0; i<parent_chains.size(); ++i)
    {
        parent_chain=parent_chains[i];
        subs=findSubstituents(parent_chain,out,1,0);
        sort(subs.begin(),subs.end(),compound::cmpBySubName);
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
vector<int> compound::directParentChain(vector<int> parent_chain, int out)
{
    if (isParentChainCyclic(parent_chain))
    {
        parent_chain=directCyclicParentChain(parent_chain,out);
    }
    else parent_chain=directAcyclicParentChain(parent_chain,out);
    return parent_chain;
}
void compound::findCyclicUnions(int in)
{
    if (atoms[in].symbol!="C") return;
    if (atoms_unions[in]!=-1) return;

    atom a;
    int s=in;
    int s2;
    stack<pair<int, int> > st;
    pair<int, int> curr;
    vector<int> prev;
    vector<int> in_union;

    int cycle=-1;

    prev.resize(atoms.size());
    for (int i=0; i<prev.size(); ++i) prev[i]=-2;
    st.push(make_pair(s,-1));


    while (!st.empty())
    {
        curr=st.top();
        st.pop();
        s=curr.first;
        if (prev[s]!=-2) continue;
        prev[s]=curr.second;
        a=atoms[s];
        for (int i=0; i<a.bonds.size(); ++i)
        {
            s2=a.bonds[i].to;
            if (s2!=-1 && atoms[s2].symbol=="C" && atoms_unions[s2]==-1)
            {
                if (prev[s2]!=-2 && prev[s]!=s2)
                {
                    prev[s2]=s;
                    cycle=s2;
                    break;
                }
                else if (prev[s2]==-2)
                {
                    st.push(make_pair(s2,s));
                }
            }
        }
        if (cycle!=-1) break;
    }
    if (cycle==-1) return;
    s2=atom_in_union.size();
    atom_in_union.push_back(cycle);
    s=cycle;
    do
    {
        in_union.push_back(s);
        atoms_unions[s]=s2;
        s=prev[s];
    }
    while (s>=0 && s!=cycle);
    ///TODO Deal with multicylce unions

}
void compound::findAcyclicUnions(int in)
{
    //cerr<<endl<<" IN: "<<in<<endl;

    if (atoms[in].symbol!="C") return;
    if (atoms_unions[in]!=-1) return;

    atom a;
    int s=in;
    int s2;
    stack<int> st;
    vector<bool> vis;
    vector<int> touched_unions(atom_in_union.size());
    vector<int> in_union;

    vis.resize(atoms.size());
    for (int i=0; i<vis.size(); ++i) vis[i]=0;
    for (int i=0; i<touched_unions.size(); ++i) touched_unions[i]=0;
    vis[s]=1;
    st.push(s);
    in_union.push_back(s);


    while (!st.empty())
    {
        s=st.top();
        st.pop();
        a=atoms[s];
        for (int i=0; i<a.bonds.size(); ++i)
        {
            s2=a.bonds[i].to;
            if (s2!=-1 && atoms[s2].symbol=="C" && !vis[s2])
            {
                if (atoms_unions[s2]==-1)
                {
                    vis[s2]=1;
                    st.push(s2);
                    in_union.push_back(s2);
                }
                else
                {
                    ++touched_unions[atoms_unions[s2]];
                }
            }
        }
    }
    s2=-1;
    for (int i=0; i<touched_unions.size(); ++i)
    {
        if (touched_unions[i]>=2)
        {
            s2=i;
            break;
        }
    }
    if (s2==-1)
    {
        s2=atom_in_union.size();
        atom_in_union.push_back(in_union[0]);
    }
    for (int i=0; i<in_union.size(); ++i)
    {
        atoms_unions[in_union[i]]=s2;
    }
}
void compound::findUnions()
{
    atoms_unions.resize(atoms.size());
    atom_in_union.resize(0);
    for (int i=0; i<atoms_unions.size(); ++i)
    {
        atoms_unions[i]=-1;
    }
    for (int i=0; i<atoms.size(); ++i)
    {
        findCyclicUnions(i);
    }
    for (int i=0; i<atoms.size(); ++i)
    {
        findAcyclicUnions(i);
    }
}
vector<int> compound::findParentChain(int in, int out)
{
    pair<vector<int>, vector<int>> candidates,candidates2;
    vector<pair<int, int> > finalCandidates;
    pair<int, int> currCandidate;
    vector<int> parent_chain= {0};
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
    for (int i=0; i<5; ++i)
    {
        maxDist[i]=0;
    }
    int starting_atom;


    if (in==-1)
    {
        findUnions();
    }
    else
    {
        //stuff
    }

    starting_atom=findAtomInCycle(in,out);
    if (starting_atom==-1)
    {
        if (in==-1)
        {
            for (int i=0; i<atoms.size(); ++i)
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
        for (int i=0; i<candidates.first.size(); ++i)
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
                for (int j=0; j<candidates2.first.size(); ++j)
                {
                    //cerr<<"FC1: "<<candidates.first[i]<<" "<<candidates2.first[j]<<endl;
                    currCandidate.second=candidates2.first[j];
                    finalCandidates.push_back(currCandidate);
                }
            }
        }
        for (int i=0; i<finalCandidates.size(); ++i)
        {
            if (finalCandidates[i].first>finalCandidates[i].second) swap(finalCandidates[i].first,finalCandidates[i].second);
        }
        sort(finalCandidates.begin(),finalCandidates.end());
        candidateParentChains2.resize(0);
        for (int i=0; i<finalCandidates.size(); ++i)
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
        for (int i=0; i<candidateParentChains2.size(); ++i)
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
        for (int i=0; i<candidateParentChains.size(); ++i)
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
        for (int i=0; i<candidateParentChains2.size(); ++i)
        {
            parent_chain=candidateParentChains2[i];
            complex_bonds=findComplexBonds(parent_chain);
            a1.resize(0);
            for (int i=0; i<complex_bonds.size(); ++i)
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
        for (int i=0; i<candidateParentChains.size(); ++i)
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
        for (int i=0; i<candidateParentChains2.size(); ++i)
        {
            parent_chain=candidateParentChains2[i];
            subs=findSubstituents(parent_chain,out,1,1);
            sort(subs.begin(),subs.end(),compound::cmpBySubName);
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
        /*if (in!=-1)
        {
            parent_chain.push_back(in);
            return parent_chain;
        }*/
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
string compound::generateName(int in, int out)
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
    for (int i=0; i<complex_bonds.size(); ++i)
    {
        //cerr<<"CB: "<<complex_bonds[i].first<<" "<<complex_bonds[i].second<<endl;
        if (complex_bonds[i].second==2) double_bonds.push_back(complex_bonds[i].first);
        else if (complex_bonds[i].second==3) triple_bonds.push_back(complex_bonds[i].first);
    }

    if (out==-1 && parent_chain.size()==1 && suffix_subs.size()==1 && prefix_subs.empty())
    {
        if (get<1>(suffix_subs[0])==curr_dict.FGS["carbonic acid"])
        {
            name=curr_dict.FGS["carbonic acid"];
            goto isInorganic;
        }
        if (get<1>(suffix_subs[0])==curr_dict.FGS["carbon dioxide"])
        {
            name=curr_dict.FGS["carbon dioxide"];
            goto isInorganic;
        }
    }

    benzene=0;
    if (parent_chain.size()==6 && double_bonds.size()==3 && triple_bonds.size()==0 && cyclic)
    {
        benzene=1;
        for (int i=1; i<3; ++i)
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
    prefix=findGNs(prefix_subs,parent_chain.size(),cyclic,benzene,most_important);
    if (!prefix.empty() && prefix[0]=='-') prefix=prefix.substr(1,prefix.size()-1);

    name=prefix+name;

    if (benzene)
    {
        if (suffix_subs.empty()) name+=curr_dict.benz+curr_dict.CBI[2];
        else name+=curr_dict.pheno;
        goto isBenzene;
    }

    if (!double_bonds.empty())
    {
        most_important=0;
        if (suffix_subs.empty() && triple_bonds.empty()) most_important=1;
        suffix=findGN(double_bonds,curr_dict.CBI[2],0,parent_chain.size(),cyclic,benzene,most_important);
        suffixes.push_back(suffix);
    }

    if (!triple_bonds.empty())
    {
        most_important=0;
        if (suffix_subs.empty() && double_bonds.empty()) most_important=1;
        suffix=findGN(triple_bonds,curr_dict.CBI[3],0,parent_chain.size(),cyclic,benzene,most_important);
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

    suffix=findGNs(suffix_subs,parent_chain.size(),cyclic,benzene,1);
    suffixes.push_back(suffix);

    for (int i=0; i<suffixes.size(); ++i)
    {
        name=addSuffix(name,suffixes[i]);
    }

    if (out!=-1)
    {
        if (!prefix_subs.empty()) name="("+name+")";
        else
        {
            for (int i=0; i<name.size(); ++i)
            {
                if (name[i]=='-')
                {
                    name="("+name+")";
                    break;
                }
            }
        }
    }

    isInorganic:

    //cerr<<"FOR IN/OUT: "<<in<<" "<<out<<" name: "<<name<<endl;
    //cerr<<"FOR IN/OUT: "<<in<<" "<<out<<" name: "<<name<<endl;
    return name;
}
string compound::getName()
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
        if (name[i]>='à' && name[i]<='ÿ')
        {
            name[i]+='À'-'à';
            break;
        }
    }*/
    return name;
}
void compound::setName(string name, double x, double y, double distx, double disty)
{
    for (int i=0; i<atoms.size(); ++i) removeAtom(i);
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
    for (int i=0; i<name.size(); ++i)
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
                else if (name[i]=='a' || name[i]=='à')
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
                        for (int j=1; j<halogen_N; ++j)
                        {
                            p=curr_dict.HP[halogen_symbols[j]];
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
                        for (int j=1; j<curr_dict.NP.size()+1; ++j)
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
    for (int i=0; i<subs.size(); ++i)
    {
        sub=subs[i];
        nums=get<0>(sub);
        pr=get<1>(sub);
        curr_num=get<2>(sub);
        for (int j=0; j<nums.size(); ++j)
        {
            //cerr<<nums[j]<<" "<<pr<<" "<<curr_num<<endl;
            if (!nums[j]) nums[j]=min(parent_chain,2);
            subs2.push_back(make_tuple(nums[j],pr,curr_num));
        }
    }
    sort(subs2.begin(),subs2.end());
    for (int i=0; i<hal_subs.size(); ++i)
    {
        hal_sub=hal_subs[i];
        nums=hal_sub.first;
        curr_num=hal_sub.second;
        for (int j=0; j<nums.size(); ++j)
        {
            hal_subs2.push_back(make_pair(nums[j],curr_num));
        }
    }
    sort(hal_subs2.begin(),hal_subs2.end());
    //cerr<<"PC: "<<parent_chain<<endl;
    vector<int> PC;
    vector<int> been_up;
    pr=0;
    for (int i=0; i<parent_chain; ++i)
    {
        curr_num=addAtom("C",4,-(parent_chain/2)*distx+i*distx+x,y);
        if (i) connectAtoms(pr,curr_num);
        pr=curr_num;
        PC.push_back(curr_num);
        been_up.push_back(0);
    }

    for (int i=0; i<double_bonds.size(); ++i)
    {
        connectAtoms(PC[double_bonds[i]-1],PC[double_bonds[i]]);
    }

    for (int i=0; i<triple_bonds.size(); ++i)
    {
        connectAtoms(PC[triple_bonds[i]-1],PC[triple_bonds[i]]);
        connectAtoms(PC[triple_bonds[i]-1],PC[triple_bonds[i]]);
    }

    double x2,y2;

    for (int i=0; i<subs2.size(); ++i)
    {
        pr=PC[get<0>(subs2[i])-1];
        x2=-(parent_chain/2)*distx+(get<0>(subs2[i])-1)*distx+x;
        y2=y;
        if (been_up[get<0>(subs2[i])-1])
        {
            disty=-disty;
            been_up[get<0>(subs2[i])-1]=2;
        }
        for (int j=0; j<get<2>(subs2[i]); ++j)
        {
            y2+=disty;
            curr_num=addAtom("C",4,x2,y2);
            connectAtoms(curr_num,pr);
            if (!j && get<1>(subs2[i])==2) connectAtoms(curr_num,pr);
            pr=curr_num;
        }
        if (been_up[get<0>(subs2[i])-1]) disty=-disty;
        else been_up[get<0>(subs2[i])-1]=1;
    }

    for (int i=0; i<hal_subs2.size(); ++i)
    {
        //cerr<<"Halogen substituent: "<<hal_subs2[i].first<<" "<<halogen_symbols[hal_subs2[i].second]<<endl;
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
            curr_num=addAtom(halogen_symbols[hal_subs2[i].second],1,x2,y2);
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
            curr_num=addAtom(halogen_symbols[hal_subs2[i].second],1,x2,y2);
            connectAtoms(curr_num,pr);
            if (hal_subs2[i].first-1==0) distx=-distx;
            if (been_up[hal_subs2[i].first-1]>2) distx=-distx;
            else been_up[hal_subs2[i].first-1]=3;
        }

    }
    //getch();
}
