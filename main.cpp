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
const double DEG2RAD=3.14159/180.0;
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
    vector<string> CNP; //carbon_number_prefixes
    vector<string> SNP; //substituent_number_prefixes
    vector<string> FGTS; //functional_group_type_suffixes
    vector<string> SS; //substituent_suffix
    string CP; //cycle_prefix
    string NC; //not_connected

    string getCNP(int CN)
    {
        if (CN>=CNP.size())
        {
            if (CN>=SNP.size()) CN=0;
            else return SNP[CN];
        }
        return CNP[CN];
    }
    string getSNP(int ST)
    {
        if (ST>=SNP.size()) ST=0;
        return SNP[ST];
    }
};

dictionary Bulgarian,English,curr_dict;
vector<dictionary> dictionaries;
int curr_dict_N;

string onlyLetters(string a)
{
    string a2;
    for (int i=0;i<a.size();++i)
    {
        if ((a[i]>='a' && a[i]<='z') || (a[i]>='�' && a[i]<='�'))
        {
            a2+=a[i];
        }
    }
    cout<<"OL: "<<a<<" "<<a2<<endl;
    return a2;
}
bool cmpBySubName(pair<int, string> a, pair<int, string> b)
{
    string a2,b2;
    a2=onlyLetters(a.second);
    b2=onlyLetters(b.second);
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
            if (nc>=0 && symbol=="C" && bonds[nc].spots_taken.size()==bonds.size()-1) return 0;
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
    int addAtom(string new_symbol, int new_valance, double new_x, double new_y, int new_bond)
    {
        int ind=-1;
        if (new_valance<=0) return -1;
        if (new_bond>=atoms.size() || atoms[new_bond].symbol=="") return -1;
        atom a(new_symbol,new_valance,new_x,new_y,new_bond);
        if (!atoms[new_bond].free_bonds.empty())
        {
            ind=addAtom(a);
            if (ind!=-1) atoms[new_bond].connect(ind);
        }
        return ind;
    }
    int connectAtoms(int a1, int a2)
    {
        if (a1<atoms.size() && a2<atoms.size() && a1!=a2 && atoms[a1].symbol!="" && atoms[a2].symbol!="" && atoms[a1].canConnect(a2) && atoms[a2].canConnect(a1))
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
            if (x<a.x-0.05 || (x>a.x+0.05 && a.free_bonds.empty()) || x>a.x+0.15 || y<a.y-0.08 || y>a.y+0.08) minDistInd=-1;
        }
        return minDistInd;
    }
    compound() {}
    compound(string new_symbol, int new_valance, double new_x, double new_y)
    {
        addAtom(new_symbol,new_valance,new_x,new_y);
    }


    string addSuffix(string name, string suffix)
    {
        bool first_letter_vowel=0;
        for (int i=0;i<suffix.size();++i)
        {
            if ((suffix[i]>='a' && suffix[i]<='z') || suffix[i]>='�' && suffix[i]<='�')
            {
                if (suffix[i]=='a' || suffix[i]=='e' || suffix[i]=='i' || suffix[i]=='o'
                    || suffix[i]=='u' || suffix[i]=='y' || suffix[i]=='�' || suffix[i]=='�'
                    || suffix[i]=='�' || suffix[i]=='�' || suffix[i]=='�' || suffix[i]=='�'
                    || suffix[i]=='�' || suffix[i]=='�' || suffix[i]=='�' || suffix[i]=='�')
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
                || name[ns]=='u' || name[ns]=='y' || name[ns]=='�' || name[ns]=='�'
                || name[ns]=='�' || name[ns]=='�' || name[ns]=='�' || name[ns]=='�'
                || name[ns]=='�' || name[ns]=='�' || name[ns]=='�' || name[ns]=='�')
            {
                name.resize(ns);
            }
        }
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
                if (s!=-1 && !vis[s] && atoms[s].symbol=="C")
                {
                    st.push(s);
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
    int findMaxCon()
    {
        atom a;
        int maxcon=1;
        for (int i=0;i<atoms.size();++i)
        {
            a=atoms[i];
            if (a.symbol=="C") for (int j=0;j<a.bonds.size();++j)
            {
                if (a.bonds[j].to!=-1 && atoms[a.bonds[j].to].symbol=="C" && a.bonds[j].spots_taken.size()>maxcon) maxcon=a.bonds[j].spots_taken.size();
            }
        }
        return maxcon;
    }
    pair<vector<int>, int> findFarthest(int s, int out)
    {
        vector<int> maxDistAtoms;
        pair<vector<int>, int> ans;
        int maxDist=-1;
        atom a;
        stack<pair<int, int> > st;
        vector<bool> vis;
        pair<int, int> curr;
        int cons;
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
            a=atoms[curr.first];
            curr.second+=atoms.size()*atoms.size();
            cons=0;
            for (int i=0;i<a.bonds.size();++i)
            {
                s=a.bonds[i].to;
                if (s!=-1 && atoms[s].symbol=="C")
                {
                    ++cons;
                    if (s==out) curr.second+=atoms.size()*atoms.size()*atoms.size();
                }
            }
            cons=max(0,cons-2);
            curr.second+=cons;
            if (curr.second>maxDist)
            {
                maxDist=curr.second;
                maxDistAtoms.resize(0);
            }
            if (curr.second==maxDist) maxDistAtoms.push_back(curr.first);
            for (int i=0;i<a.bonds.size();++i)
            {
                s=a.bonds[i].to;
                if (s!=-1 && atoms[s].symbol=="C" && !vis[s] && s!=out)
                {
                    curr.first=s;
                    vis[s]=1;
                    if (a.bonds[i].spots_taken.size()>1)
                    {
                        curr.second+=atoms.size();
                        st.push(curr);
                        curr.second-=atoms.size();
                    }
                    else
                    {
                        st.push(curr);
                    }
                }
            }
        }
        ans.first=maxDistAtoms;
        ans.second=maxDist;
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
    pair<int, int> findAttachment(vector<int> parent_chain, int out)
    {
        atom a;
        pair<int, int> attachment;
        attachment.first=0;
        attachment.second=0;
        for (int i=0;i<parent_chain.size();++i)
        {
            a=atoms[parent_chain[i]];
            for (int j=0;j<a.bonds.size();++j)
            {
                if (a.bonds[j].to==out)
                {
                    attachment.first=i+1;
                    attachment.second=a.bonds[j].spots_taken.size();
                    return attachment;
                }
            }
        }
        return attachment;
    }
    vector<pair<int, int> > findComplexBonds(vector<int> parent_chain)
    {
        atom a;
        vector<pair<int, int> > comp_bonds;
        pair<int, int> comp_bond;
        for (int i=0;i<parent_chain.size()-1;++i)
        {
            a=atoms[parent_chain[i]];
            for (int j=0;j<a.bonds.size();++j)
            {
                if (a.bonds[j].to==parent_chain[i+1] && a.bonds[j].spots_taken.size()>1)
                {
                    comp_bond.first=i+1;
                    comp_bond.second=a.bonds[j].spots_taken.size();
                    comp_bonds.push_back(comp_bond);
                }
            }
        }
        return comp_bonds;
    }
    vector<pair<int, string> > findSubstituents(vector<int> parent_chain, int out, bool wantNames)
    {
        atom a;
        vector<pair<int, string> > subs;
        pair<int, string> sub;
        int PCS=parent_chain.size();
        for (int i=0;i<PCS;++i)
        {
            a=atoms[parent_chain[i]];
            for (int j=0;j<a.bonds.size();++j)
            {
                if (a.bonds[j].to!=-1 && a.bonds[j].to!=parent_chain[(i-1+PCS)%PCS] && a.bonds[j].to!=parent_chain[(i+1)%PCS] && a.bonds[j].to!=out)
                {
                    sub.first=i+1;
                    if (wantNames) sub.second=generateName(a.bonds[j].to,parent_chain[i]);
                    else sub.second="";
                    subs.push_back(sub);
                }
            }
        }
        return subs;
    }
    bool isParentChainCyclic(vector<int> parent_chain)
    {
        if (parent_chain.size()<3) return 0;
        atom a;
        a=atoms[parent_chain[parent_chain.size()-1]];
        for (int i=0;i<a.bonds.size();++i)
        {
            if (a.bonds[i].to==parent_chain[0]) return 1;
        }
        return 0;
    }
    string findSP(vector<int> pos, string name) //find substituents_prefix
    {
        string prefix="";
        for (int i=0;i<pos.size();++i)
        {
            if (i>0) prefix+=',';
            prefix+=intToString(pos[i]);
        }
        prefix+='-';
        if (pos.size()>1) prefix+=curr_dict.getSNP(pos.size());
        prefix+=name;
        return prefix;
    }
    string findSPs(vector<pair<int, string> > subs)//find substituents_prefixes
    {
        string prefixes="";
        vector<int> currPos;
        if (subs.size()==0) return prefixes;

        cerr<<"Before sort: "<<endl;
        for (int i=0;i<subs.size();++i)
        {
            cerr<<subs[i].first<<" "<<subs[i].second<<endl;
        }

        sort(subs.begin(),subs.end(),cmpBySubName);

        cerr<<"After sort: "<<endl;
        for (int i=0;i<subs.size();++i)
        {
            cerr<<subs[i].first<<" "<<subs[i].second<<endl;
        }

        string curr=subs[0].second;
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
    vector<int> directAcyclicParentChain(vector<int> parent_chain, int out)
    {
        atom a;
        vector<int> parent_chain2;
        vector<pair<int, string> > subs1;
        vector<pair<int, string> > subs2;
        vector<pair<int, int> > complex_bonds1;
        vector<pair<int, int> > complex_bonds2;
        vector<int> double_bonds1;
        vector<int> double_bonds2;
        pair<int, int> attachment1;
        pair<int, int> attachment2;
        parent_chain2.resize(parent_chain.size());
        for (int i=0;i<parent_chain.size();++i)
        {
            parent_chain2[parent_chain.size()-1-i]=parent_chain[i];
        }
        attachment1=findAttachment(parent_chain,out);
        attachment2=findAttachment(parent_chain2,out);
        if (attachment1.first<attachment2.first)
        {
            return parent_chain;
        }
        if (attachment1.first>attachment2.first)
        {
            return parent_chain2;
        }
        complex_bonds1=findComplexBonds(parent_chain);
        complex_bonds2=findComplexBonds(parent_chain2);
        for (int i=0;i<complex_bonds1.size();++i)
        {
            if (complex_bonds1[i].first<complex_bonds2[i].first)
            {
                return parent_chain;
            }
            if (complex_bonds1[i].first>complex_bonds2[i].first)
            {
                return parent_chain2;
            }
        }
        for (int i=0;i<complex_bonds1.size();++i)
        {
            if (complex_bonds1[i].second==2) double_bonds1.push_back(complex_bonds1[i].first);
        }
        for (int i=0;i<complex_bonds2.size();++i)
        {
            if (complex_bonds2[i].second==2) double_bonds2.push_back(complex_bonds2[i].first);
        }
        for (int i=0;i<double_bonds1.size();++i)
        {
            if (double_bonds1[i]<double_bonds2[i])
            {
                return parent_chain;
            }
            if (double_bonds1[i]>double_bonds2[i])
            {
                return parent_chain2;
            }
        }
        subs1=findSubstituents(parent_chain,out,0);
        subs2=findSubstituents(parent_chain2,out,0);
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
    vector<int> directCyclicParentChain(vector<int> parent_chain, int out)
    {
        atom a;
        vector<vector<int> > parent_chains;
        vector<vector<int> > parent_chains2;
        vector<pair<int, string> > subs;
        vector<pair<int, string> > max_subs;
        vector<pair<int, int> > complex_bonds;
        vector<pair<int, int> > max_complex_bonds;
        vector<int> double_bonds;
        vector<int> max_double_bonds;
        pair<int, int> attachment;
        pair<int, int> max_attachment;

        max_subs.resize(0);
        max_complex_bonds.resize(0);
        max_double_bonds.resize(0);
        max_attachment=make_pair(-1,-1);


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

        for (int i=0;i<parent_chains.size();++i)
        {
            parent_chain=parent_chains[i];
            attachment=findAttachment(parent_chain,out);
            if ((max_attachment.first==-1 && max_attachment.second==-1) || attachment.first<max_attachment.first)
            {
                parent_chains2.resize(0);
                max_attachment=attachment;
            }
            if (attachment.first==max_attachment.first)
            {
                parent_chains2.push_back(parent_chain);
            }
        }

        int f;

        for (int i=0;i<parent_chains2.size();++i)
        {
            parent_chain=parent_chains2[i];
            complex_bonds=findComplexBonds(parent_chain);
            f=0;
            if (max_complex_bonds.empty() && !complex_bonds.empty()) f=1;
            else
            {
                for (int i=0;i<complex_bonds.size();++i)
                {
                    if (complex_bonds[i].first<max_complex_bonds[i].first)
                    {
                        f=1;
                        break;
                    }
                    if (complex_bonds[i].first>max_complex_bonds[i].first)
                    {
                        f=-1;
                        break;
                    }
                }
            }
            if (f==1)
            {
                f=0;
                parent_chains.resize(0);
                max_complex_bonds=complex_bonds;
            }
            if (f==0)
            {
                parent_chains.push_back(parent_chain);
            }
        }

        for (int i=0;i<parent_chains.size();++i)
        {
            parent_chain=parent_chains[i];
            complex_bonds=findComplexBonds(parent_chain);
            for (int i=0;i<complex_bonds.size();++i)
            {
                if (complex_bonds[i].second==2) double_bonds.push_back(complex_bonds[i].first);
            }
            f=0;
            if (max_double_bonds.empty() && !double_bonds.empty()) f=1;
            else
            {
                for (int i=0;i<double_bonds.size();++i)
                {
                    if (double_bonds[i]<max_double_bonds[i])
                    {
                        f=1;
                        break;
                    }
                    if (double_bonds[i]>max_double_bonds[i])
                    {
                        f=-1;
                        break;
                    }
                }
            }
            if (f==1)
            {
                f=0;
                parent_chains2.resize(0);
                max_double_bonds=double_bonds;
            }
            if (f==0)
            {
                parent_chains2.push_back(parent_chain);
            }
        }

        for (int i=0;i<parent_chains2.size();++i)
        {
            parent_chain=parent_chains2[i];
            subs=findSubstituents(parent_chain,out,0);
            f=0;
            if (max_subs.empty() && !subs.empty()) f=1;
            else
            {
                for (int i=0;i<subs.size();++i)
                {
                    if (subs[i].first<max_subs[i].first)
                    {
                        f=1;
                        break;
                    }
                    if (subs[i].first>max_subs[i].first)
                    {
                        f=-1;
                        break;
                    }
                }
            }
            if (f==1)
            {
                f=0;
                parent_chains.resize(0);
                max_subs=subs;
            }
            if (f==0)
            {
                parent_chains.push_back(parent_chain);
            }
        }

        parent_chain=parent_chains[0];
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
                    candidateParrentChains.push_back(parent_chain);
                }
            }
            parent_chain=candidateParrentChains[0];
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
        string name="";
        vector<string> suffixes;
        string suffix;
        vector<int> parent_chain;
        vector<pair<int, string> > subs;
        vector<pair<int, int> > complex_bonds;
        vector<int> double_bonds;
        vector<int> triple_bonds;
        pair<int, int> attachment;
        parent_chain=findParentChain(in,out);
        parent_chain=directParentChain(parent_chain,out);
        subs=findSubstituents(parent_chain,out,1);
        name+=findSPs(subs);
        complex_bonds=findComplexBonds(parent_chain);
        for (int i=0;i<complex_bonds.size();++i)
        {
            if (complex_bonds[i].second==2) double_bonds.push_back(complex_bonds[i].first);
            else if (complex_bonds[i].second==3) triple_bonds.push_back(complex_bonds[i].first);
        }
        if (!double_bonds.empty())
        {
            suffix="";
            if (parent_chain.size()>2)
            {
                suffix+='-';
                for(int i=0;i<double_bonds.size();++i)
                {
                    if (i) suffix+=',';
                    suffix+=intToString(double_bonds[i]);
                }
                suffix+='-';
            }
            if (double_bonds.size()>1) suffix+=curr_dict.getSNP(double_bonds.size());
            suffix+=curr_dict.FGTS[2];
            suffixes.push_back(suffix);
        }
        if (!triple_bonds.empty())
        {
            suffix="";
            if (parent_chain.size()>2)
            {
                suffix+='-';
                for(int i=0;i<triple_bonds.size();++i)
                {
                    if (i) suffix+=',';
                    suffix+=intToString(triple_bonds[i]);
                }
                suffix+='-';
            }
            if (triple_bonds.size()>1) suffix+=curr_dict.getSNP(triple_bonds.size());
            suffix+=curr_dict.FGTS[3];
            suffixes.push_back(suffix);
        }
        if (complex_bonds.empty())
        {
            suffix=curr_dict.FGTS[1];
            suffixes.push_back(suffix);
        }
        if (in!=-1)
        {
            suffix="";
            attachment=findAttachment(parent_chain,out);
            if (attachment.first!=1)
            {
                suffix+='-';
                suffix+=intToString(attachment.first);
                suffix+='-';
            }
            suffix+=curr_dict.SS[attachment.second];
            if (suffixes[suffixes.size()-1]==curr_dict.FGTS[1] && attachment.first==1)
            {
                suffixes[suffixes.size()-1]=suffix;
            }
            else
            {
                suffixes.push_back(suffix);
            }
        }
        if (isParentChainCyclic(parent_chain)) name+=curr_dict.CP;
        name+=curr_dict.getCNP(parent_chain.size());
        for (int i=0;i<suffixes.size();++i)
        {
            name=addSuffix(name,suffixes[i]);
        }
        if (in!=-1)
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
        return name;
    }
    string getName()
    {
        string name;
        int CS=isConnected(); //bond_status
        if (!CS) return curr_dict.NC;
        name=generateName(-1,-1);

        for (int i=0;i<name.size();++i)
        {
            if (name[i]>='a' && name[i]<='z')
            {
                name[i]+='A'-'a';
                break;
            }
            if (name[i]>='�' && name[i]<='�')
            {
                name[i]+='�'-'�';
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
    if (key==GLFW_KEY_L && action==GLFW_PRESS) pressed=-6;
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
    double nextpos;
    if (a.free_bonds.empty())
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
        drawIndex(a.free_bonds.size(),nextpos-0.02*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT,a.y*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT-0.08*ORIGINAL_WINDOWS_HEIGHT/WINDOWS_HEIGHT);
    }
}
void drawBond(double x1, double y1, double x2, double y2, int num)
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
void run(GLFWwindow* w)
{
    string carbonSymbol="C";
    int carbonValance=4;
    double sx,sy;
    sx=0;
    sy=0;
    if (snappingEnabled) snap(sx,sy);
    compound c(carbonSymbol,carbonValance,sx,sy);
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
        if (pressed==-6)
        {
            ++curr_dict_N;
            curr_dict_N%=dictionaries.size();
            curr_dict=dictionaries[curr_dict_N];
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
                        if (pressed==GLFW_MOUSE_BUTTON_LEFT || c.atoms[last].free_bonds.empty()) last2=-1;
                        else last2=last;
                    }
                    last=-1;

                }
                else
                {
                    result=c.addAtom(carbonSymbol,carbonValance,mxpos,mypos,last2);
                    if (result!=-1)
                    {
                        if (pressed==GLFW_MOUSE_BUTTON_LEFT || c.atoms[result].free_bonds.empty()) last2=-1;
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
                    result=c.addAtom(carbonSymbol,carbonValance,mxpos,mypos);
                    if (result!=-1)
                    {
                        if (pressed==GLFW_MOUSE_BUTTON_LEFT || c.atoms[result].free_bonds.empty()) last2=-1;
                        else last2=result;
                    }
                }
                else
                {
                    if (pressed==GLFW_MOUSE_BUTTON_LEFT) toMove=1;
                    else toMove=0;
                    if (!c.atoms[last].free_bonds.empty())
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
    English.CNP={"alka","metha","etha","propa","buta"};
    English.FGTS={"error","ane","ene","yne"};
    English.SNP={"error","mono","di","tri","tetra","penta","hexa","hepta","octa","nona",
    "deca","undeca","dodeca","trideca","tetradeca","pentadeca","hexadeca","heptadeca","octadeca","nonadeca",
    "icosa","henicosa","docosa","tricosa","tetracosa","pentacosa","hexacosa","heptacosa","octacosa","nonacosa",
    "triaconta","hentriaconta","hentriaconta","tritriaconta"};
    English.SS={"error","yl","ylidene"};
    English.CP="cyclo";
    English.NC="Not Connected";
    dictionaries.push_back(English);

    Bulgarian.CNP={"����","����","���","�����","����"};
    Bulgarian.FGTS={"������","��","��","��"};
    Bulgarian.SNP={"������","����","��","���","�����","�����","�����","�����","����","����",
    "����","������","������","�������","���������","���������","���������","���������","��������","��������"};
    Bulgarian.SS={"������","��","������"};
    Bulgarian.CP="�����";
    Bulgarian.NC="�� �� ��������";
    dictionaries.push_back(Bulgarian);

    curr_dict_N=0;

    curr_dict=dictionaries[curr_dict_N];
}
int main()
{
    system("chcp 1251");
    system("cls");

    setDictionaries();

    cout<<"Use the Middle Mouse Button to start or continue chains.\nUse the Left Mouse Button to end chains and to move atoms,"
    <<"\nUse the Right mouse button to cancel chains and remove atoms.\nPress Backspace to undo.\nPress Shift to toggle snapping on and off."
    <<"\nPress R to reset.\nPress Escape to end the program.\nPress L to change the language.\nPress any key to continue."<<endl;
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
