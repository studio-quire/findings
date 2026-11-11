// five_critical.cpp: the extension of one five-form core S5 as Crapo-Rota's critical problem (apn/README.md, "Toward a proof").
// Written by Claude Fable 5.2 in a consult for Rowen, 2026-11-09 (runs/consult/fable1/core.cpp); copied unchanged below this header.
// Usage: five_critical F1 F2 F3 F4 F5 [seed] [minimize] [a0] [hmode]   (forms as 64-bit hex, byte i = row i)
//   hmode 0: is there pi: span(G5) -> F_2^3 nonzero on every S5-isotropic plane (SAT = an APN S8 over S5)? If not, shrink to a core.
//   hmode 1: the same question on the planes inside each of the 255 hyperplanes of V; last line 'sat X unsat Y cap Z'.
// Critical-problem core finder. Input: 5 hex forms (S5). G5 = S5-isotropic planes as Pluecker vectors.
// Decide: is there pi: span(G) -> F_2^3 nonzero on every point of G? (yes <=> an APN S8 over S5.)
// If not, shrink G to a small unsatisfiable subset (core) and print it.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <random>
using namespace std;
typedef uint32_t u32; typedef uint64_t u64;
static int PIDX[8][8];
static u32 wedge(int a,int b){u32 w=0;for(int i=0;i<8;i++)for(int j=i+1;j<8;j++){int t=((a>>i&1)&(b>>j&1))^((a>>j&1)&(b>>i&1));if(t)w|=1u<<PIDX[i][j];}return w;}
struct Pt{u32 w;int a,b;};
static vector<Pt> ALL;
struct Ech{ u32 lead[28]; u32 combo[28]; int n; Ech(){for(int i=0;i<28;i++)lead[i]=0;n=0;}
  // reduce w; returns residual, sets c to combo over basis indices
  u32 red(u32 w,u32&c)const{c=0;for(int b=27;b>=0;b--)if((w>>b&1)&&lead[b]){w^=lead[b];c^=combo[b];}return w;}
  void add(u32 res,u32 c){int b=31-__lzcnt(res);lead[b]=res;combo[b]=c^(1u<<n);n++;
    // keep fully reduced: clear bit b from others
    for(int k=0;k<28;k++)if(k!=b&&lead[k]&&(lead[k]>>b&1)){lead[k]^=lead[b];combo[k]^=combo[b];}}
};
#include <intrin.h>
static u64 NODECAP=2000000000ULL;
// solve on subset S (indices into ALL); forced: first three points of S are a 3-dim W_a basis (values 1,2,4).
// prio: points to check first. returns 1 sat, 0 unsat, -1 cap. used[] marks failing points.
static u64 nodes; 
static int solve(const vector<int>&S,int nforced,const vector<char>&prio,vector<char>&used,vector<int>*solOut=nullptr){
  int m=S.size(); Ech E; vector<char> done(m,0); vector<u32> mask(m,0); vector<int> basisPt;
  vector<vector<int>> chk; // per level list of subset positions
  auto attach=[&](){ for(int i=0;i<m;i++)if(!done[i]){u32 c;u32 r=E.red(ALL[S[i]].w,c);if(!r){done[i]=1;mask[i]=c;int lv=31-__lzcnt(c);chk[lv].push_back(i);}} };
  for(int k=0;k<nforced;k++){u32 c;u32 r=E.red(ALL[S[k]].w,c);if(!r){fprintf(stderr,"forced dependent\n");exit(1);}E.add(r,c);done[k]=1;mask[k]=1u<<k;chk.push_back(vector<int>());chk[k].push_back(k);basisPt.push_back(k);}
  attach();
  while(true){ // greedy next basis point
    vector<pair<u32,int>> rs; for(int i=0;i<m;i++)if(!done[i]){u32 c;rs.push_back({E.red(ALL[S[i]].w,c),i});}
    if(rs.empty())break; sort(rs.begin(),rs.end());
    int best=-1,bc=0; for(size_t i=0;i<rs.size();){size_t j=i;while(j<rs.size()&&rs[j].first==rs[i].first)j++; if((int)(j-i)>bc){bc=j-i;best=rs[i].second;} i=j;}
    // refine: among candidates, prefer max count of points landing in span after adding (already = same residual class)
    u32 c;u32 r=E.red(ALL[S[best]].w,c);int lv=E.n;E.add(r,c);done[best]=1;mask[best]=1u<<lv;chk.push_back(vector<int>());chk[lv].push_back(best);basisPt.push_back(best);
    attach();
  }
  int nb=E.n;
  for(int lv=0;lv<nb;lv++) stable_sort(chk[lv].begin(),chk[lv].end(),[&](int x,int y){return prio[S[x]]>prio[S[y]];});
  vector<int> val(nb,0); 
  // iterative DFS
  int lv=0; vector<int> cur(nb,0);
  auto ok=[&](int L)->bool{ for(int i:chk[L]){u32 mk=mask[i];int v=0;while(mk){int b=_tzcnt_u32(mk);mk&=mk-1;v^=val[b];} if(!v){used[S[i]]=1;return false;}} return true; };
  for(int k=0;k<nforced;k++){val[k]=1<<k; if(!ok(k)){return 0;}}
  lv=nforced; if(lv==nb){return 1;} cur[lv]=0;
  while(true){
    if(cur[lv]>7){ lv--; if(lv<nforced)return 0; continue; }
    val[lv]=cur[lv]++; nodes++; if(nodes>NODECAP)return -1;
    if(ok(lv)){ if(lv==nb-1){ if(solOut){solOut->assign(val.begin(),val.end());} return 1;} lv++; cur[lv]=0; }
  }
}
int main(int argc,char**argv){
  int n=0;for(int i=0;i<8;i++)for(int j=i+1;j<8;j++)PIDX[i][j]=PIDX[j][i]=n++;
  u64 f[5]; for(int k=0;k<5;k++)f[k]=strtoull(argv[1+k],0,16);
  int seed=argc>6?atoi(argv[6]):1; int minimize=argc>7?atoi(argv[7]):1; int a0req=argc>8?atoi(argv[8]):-1;
  int hmode=argc>9?atoi(argv[9]):0;
  u32 bv[5]; for(int k=0;k<5;k++){bv[k]=0;for(int i=0;i<8;i++)for(int j=i+1;j<8;j++)if(f[k]>>(8*i+j)&1)bv[k]|=1u<<PIDX[i][j];}
  for(int a=1;a<256;a++)for(int b=a+1;b<256;b++){int c=a^b;if(c<b)continue;u32 w=wedge(a,b);bool z=true;for(int k=0;k<5;k++)if(__popcnt(bv[k]&w)&1){z=false;break;}if(z)ALL.push_back({w,a,b});}
  vector<Pt> FULL=ALL;
  int N=ALL.size(); printf("G5 %d\n",N);
  if(hmode){ // hmode=1: all hyperplanes; 2: all codim-2 subspaces
    int nsat=0,nun=0,ncap=0; NODECAP=300000000ULL;
    for(int h1=1;h1<256;h1++)for(int h2=(hmode==2?h1+1:0);h2<(hmode==2?256:1);h2++){
      if(hmode==2&&((h1^h2)<h2))continue;
      ALL.clear(); for(auto&p:FULL){ if((__popcnt(p.a&h1)&1)||(__popcnt(p.b&h1)&1))continue; if(hmode==2&&((__popcnt(p.a&h2)&1)||(__popcnt(p.b&h2)&1)))continue; ALL.push_back(p);}
      int n2=ALL.size(); vector<vector<int>> th(256); for(int i=0;i<n2;i++){th[ALL[i].a].push_back(i);th[ALL[i].b].push_back(i);th[ALL[i].a^ALL[i].b].push_back(i);}
      int a0=-1;size_t bs=0;for(int a=1;a<256;a++)if(th[a].size()>bs&&(th[a].size()==7||th[a].size()==3||th[a].size()==1)){bs=th[a].size();a0=a;}
      vector<int> forced; {Ech E;for(int i:th[a0]){u32 c;u32 r=E.red(ALL[i].w,c);if(r){E.add(r,c);forced.push_back(i);}}}
      vector<int> S=forced; for(int i:th[a0])if(find(S.begin(),S.end(),i)==S.end())S.push_back(i);
      for(int i=0;i<n2;i++)if(find(S.begin(),S.end(),i)==S.end())S.push_back(i);
      vector<char> prio(n2,0),used(n2,0); nodes=0; int r=solve(S,forced.size(),prio,used);
      Ech E; for(auto&p:ALL){u32 c;u32 rr=E.red(p.w,c);if(rr)E.add(rr,c);}
      if(r==1)nsat++;else if(r==0)nun++;else ncap++;
      if(hmode==1||r!=1)printf("H %d %d planes %d rank %d forced %d -> %s nodes %llu\n",h1,h2,n2,E.n,(int)forced.size(),r==1?"SAT":r==0?"UNSAT":"cap",nodes);
    }
    printf("sat %d unsat %d cap %d\n",nsat,nun,ncap); return 0; }
  // per point a: planes through a
  vector<vector<int>> thru(256); for(int i=0;i<N;i++){thru[ALL[i].a].push_back(i);thru[ALL[i].b].push_back(i);thru[ALL[i].a^ALL[i].b].push_back(i);}
  int a0=-1; for(int a=1;a<256;a++)if(thru[a].size()==7&&(a0req<0||a==a0req)){a0=a;break;}
  if(a0<0){printf("no 3-dim W_a\n");return 1;}
  // basis of W_a0: three independent among the 7
  vector<int> wa=thru[a0]; vector<int> forced; {Ech E;for(int i:wa){u32 c;u32 r=E.red(ALL[i].w,c);if(r){E.add(r,c);forced.push_back(i);}}}
  vector<int> mand=forced; for(int i:wa)if(find(mand.begin(),mand.end(),i)==mand.end())mand.push_back(i);
  vector<char> inS(N,1),prio(N,0),used(N,0);
  auto build=[&](){vector<int>S=mand;for(int i=0;i<N;i++)if(inS[i]&&find(mand.begin(),mand.end(),i)==mand.end())S.push_back(i);return S;};
  mt19937 rng(seed);
  int size=N;
  for(int it=0;it<50;it++){
    vector<int>S=build(); fill(used.begin(),used.end(),0); nodes=0; vector<int> sol;
    int r=solve(S,3,prio,used,&sol);
    if(r==1){printf("SAT! size %d nodes %llu\n",(int)S.size(),nodes);return 0;}
    if(r<0){printf("cap\n");return 2;}
    for(int i:mand)used[i]=1; int cnt=0;for(int i=0;i<N;i++){inS[i]=used[i];prio[i]=used[i];cnt+=used[i];}
    printf("iter %d subset %d -> core %d nodes %llu\n",it,(int)S.size(),cnt,nodes); fflush(stdout);
    if(cnt==size)break; size=cnt;
  }
  if(minimize){
    vector<int> order;for(int i=0;i<N;i++)if(inS[i]&&find(mand.begin(),mand.end(),i)==mand.end())order.push_back(i);
    shuffle(order.begin(),order.end(),rng); NODECAP=(u64)minimize*1000000ULL;
    for(int p:order){ if(!inS[p])continue; inS[p]=0; vector<int>S=build(); fill(used.begin(),used.end(),0); nodes=0;
      int r=solve(S,3,prio,used);
      if(r==0){for(int i:mand)used[i]=1;for(int i=0;i<N;i++)inS[i]=used[i];}
      else inS[p]=1; }
    int cnt=0;for(int i=0;i<N;i++)cnt+=inS[i]; printf("minimized core %d\n",cnt);
  }
  {vector<int>S=build(); Ech E; for(int i:S){u32 c;u32 r=E.red(ALL[i].w,c);if(r)E.add(r,c);} printf("core rank %d a0 %d\n",E.n,a0);
   for(int i:S)printf("P %d %d\n",ALL[i].a,ALL[i].b);}
  return 0;
}
