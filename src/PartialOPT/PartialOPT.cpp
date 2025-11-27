#include<bits/stdc++.h>
#pragma GCC optimize(2) 
#define maxm 30000005
#define maxn 15000005
using namespace std;
const int para_h = 3;
const double inf = 1e18;
int n,m,k;
struct Group
{
	int sz;
	vector<int> a;
}grp[maxn];
bool operator < (Group A,Group B){return A.sz<B.sz;}
struct edge
{
	int to,id;
	double w;
	edge(int to=0,double w=0,int id=0):to(to),w(w),id(id){}
};
vector<int> U,V;
vector<double> W; 
vector<vector<edge>> g;
int idcnt;
std::string graph_file;
void readgraph()
{
	freopen(("data/" + graph_file + "/graph.txt").c_str(), "r", stdin);
	scanf("%d%d",&n,&m); 
	int u,v;
	double w;
	idcnt=0;
	U.clear();V.clear();W.clear();g.clear();
	U.resize(m+5);V.resize(m+5);W.resize(m+5);
	g.resize(n+5);
	while(~scanf("%d%d%lf",&u,&v,&w))
	{
		++idcnt;
		g[u].push_back(edge(v,w,idcnt));
		g[v].push_back(edge(u,w,idcnt));
		U[idcnt]=u;V[idcnt]=v;W[idcnt]=w;
	}
	fclose(stdin);
}
bool inq[maxn];
vector<int> used;
vector<double> dis2;
vector<int> preu2,preid2;
void sssp(vector<int> s,vector<double> &dis,vector<int> &preu,vector<int> &preid)
{
	used.clear();used.resize(maxn);
	for(int i=1;i<=n;++i)dis[i]=inf,used[i]=0,preu[i]=0,preid[i]=0;
	priority_queue < pair<double,int> > pq;
	for(int x:s)pq.push(make_pair(0,x)),dis[x]=0;
	while(!pq.empty())
	{
		int u=pq.top().second;pq.pop();
		if(used[u])continue;
		used[u]=1;
		for(auto p:g[u])
		{
			int v=p.to;
			double w=p.w;
			int id=p.id;
			if(dis[v]>dis[u]+w)
			{
				dis[v]=dis[u]+w;
				preu[v]=u;
				preid[v]=id;
				pq.push(make_pair(-dis[v],v));
			}
		}
	}
}

void addpath2(vector<int> &e,int u)
{
	int x=u;
	while(preu2[x])
	{
		e.push_back(preid2[x]);
		x=preu2[x]; 
	}
}

bool cmp (const int &x,const int &y)
{
	return W[x]<W[y];
}
vector<int> par;
int find(int x)
{
	if(par[x]==x)return x;
	return par[x]=find(par[x]);
}
double exlhlerA(int r,vector<Group> grp,vector<int> &anse)
{
	anse.clear();
	vector<int> s;
	s.push_back(r);
	sssp(s,dis2,preu2,preid2);
	vector<int> path;
	path.clear();
	bool ok=1;
	for(int i=0;i<grp.size();++i)
	{
		int mnid=0;
		double mn=inf;
		for(int x:grp[i].a)
		{
			if(mn>dis2[x])
			{
				mn=dis2[x];
				mnid=x;
			}
		}
		if(mn>1e12)
		{
			ok=0;
			break;
		}
		addpath2(path,mnid);
	}
	if(!ok)return 1e12;
	sort(path.begin(),path.end());
	path.erase(unique(path.begin(),path.end()),path.end());
	double ans=0;
	vector<int> edg;
	edg.clear();
	for(int e:path)ans+=W[e],edg.push_back(e);
	int m=edg.size();
	sort(edg.begin(),edg.end(),cmp);
	par.clear();
	par.resize(n+5); 
	for(int i=1;i<=n;++i)par[i]=i;
	double res=0;
	for(int i=0;i<m;++i)
	{
		int u=U[edg[i]],v=V[edg[i]];
		if(find(u)!=find(v))
		{
			par[find(u)]=find(v);
			res+=W[edg[i]];
			anse.push_back(edg[i]);
		}
	}
	return res;
}

#define maxS 17
vector< vector<int> > vis;
vector< vector<double> > dp;
vector< vector< vector<int> > > F;

double DPBF(int r,vector<Group> grp,vector<int> &anse,int h)
{

 	for(int i=0;i<maxn;++i)
    	for(int j=0;j<maxS;++j)
    	{
    		vis[i][j]=0;
    		dp[i][j]=0;
    		F[i][j].clear();
		}
	anse.clear();
	priority_queue< pair<double,pair<int,int> > > q;
	for(int i=0;i<h-1;++i)
		for(int x:grp[i].a)
		{
			vis[x][1<<i]=1;
			dp[x][1<<i]=0;
			F[x][1<<i].clear();
			q.push(make_pair(0,make_pair(x,1<<i)));
		}
	vis[r][1<<(h-1)]=1;
	dp[r][1<<(h-1)]=0;
	F[r][1<<(h-1)].clear();
	q.push(make_pair(0,make_pair(r,1<<(h-1))));
	while(!q.empty())
	{
		int u=q.top().second.first,S=q.top().second.second;
		q.pop();
		if(S==(1<<h)-1)
		{
			for(int e:F[u][S])anse.push_back(e);
			return dp[u][S];
		}
		for(auto p:g[u])
		{
			int v=p.to,id=p.id;
			double w=p.w;
			if(!vis[v][S]||dp[v][S]>dp[u][S]+w)
			{
				vis[v][S]=1;
				dp[v][S]=dp[u][S]+w;
				F[v][S].clear();
				for(int o:F[u][S])F[v][S].push_back(o);
				F[v][S].push_back(id);
				q.push(make_pair(-dp[v][S],make_pair(v,S)));
			}
		}
		int T=(1<<h)-1;
		T^=S;
		for(int p=T;p;p=(p-1)&T)if(vis[u][p])
		{
			if(!vis[u][S|p]||dp[u][S|p]>dp[u][S]+dp[u][p])
			{
				vis[u][S|p]=1;
				dp[u][S|p]=dp[u][S]+dp[u][p];
				F[u][S|p].clear();
				for(int o:F[u][S])F[u][S|p].push_back(o);
				for(int o:F[u][p])F[u][S|p].push_back(o);
				q.push(make_pair(-dp[u][S|p],make_pair(u,S|p)));
			}
		}
	}
	return inf;
}

int d[maxn],g_cnt[maxn];
vector< pair<int,int> > G[maxn];
vector<int> has[maxn];
double reduce_tree(vector<int> Ans2,vector<int> &Ans3)
{
	for(int i=1;i<=n;++i)d[i]=-1,has[i].clear(),G[i].clear();
	for(int e:Ans2)
	{
		int u=U[e],v=V[e];
		G[u].push_back(make_pair(v,e));
		G[v].push_back(make_pair(u,e));
		d[u]++;d[v]++;
	}
	for(int i=1;i<=k;++i)g_cnt[i]=0;
	for(int i=1;i<=k;++i)
		for(int x:grp[i].a)if(d[x]>=0) 
		{
			g_cnt[i]++;
			has[x].push_back(i);
		}
	priority_queue < pair<double,pair<int,int> > > q;
	for(int i=1;i<=n;++i)if(d[i]==0)
	{
		int nx=0;
		for(int j:has[i])if(g_cnt[j]==1)++nx;
		if(nx)continue;
		for(auto p:G[i])
		{
			int x=p.first,e=p.second;
			if(d[x]<=0)continue;
			q.push(make_pair(W[e],make_pair(i,e)));
		}
	}
	unordered_map<int,int> del;
	while(!q.empty())
	{
		int e=q.top().second.second,v=q.top().second.first,u;
		q.pop();
		if(U[e]==v)u=V[e];else u=U[e];
		int nx=0;
		for(int j:has[v])if(g_cnt[j]==1)++nx;
		if(d[v]!=0)continue;
		if(nx)continue;
		for(int j:has[v])g_cnt[j]--;
		--d[u];
		del[e]=1;
		if(d[u]==0)
		{
			int ny=0;
			for(int j:has[u])if(g_cnt[j]==1)++ny;
			if(ny)continue;
			for(auto p:G[u])
			{
				int x=p.first,e=p.second;
				if(d[x]<=0)continue;
				q.push(make_pair(W[e],make_pair(u,e)));
			}
		}
	}
	for(int e:Ans2)if(!del.count(e))Ans3.push_back(e);
	double rr=0;
	for(int e:Ans3)rr+=W[e];
	return rr;
}

double PartialOPT(int h)
{
	time_t st=time(nullptr);
	sort(grp+1,grp+k+1);
	double ans=inf;
	vector<int> Ans;
	for(int r:grp[1].a)
	{
		if(difftime(time(nullptr),st)>500)return -1;
		vector<Group> grp1,grp2;
		for(int i=2;i<=k;++i)
		{
			if(grp1.size()<h-1)grp1.push_back(grp[i]);
			else grp2.push_back(grp[i]);
		}
		vector<int> e1,e2,ed;
		int H=min(h,(int)grp1.size()+1);
		double w1=DPBF(r,grp1,e1,H);
		double w2=exlhlerA(r,grp2,e2);
		if(w1>1e12||w2>1e12)continue;
		for(int e:e1)ed.push_back(e);
		for(int e:e2)ed.push_back(e);
		sort(ed.begin(),ed.end(),cmp);
		vector<int> tote,ne;
		par.clear();
		par.resize(n+5);
		for(int i=1;i<=n;++i)par[i]=i;
		int m=ed.size();
		for(int i=0;i<m;++i)
		{
			int u=U[ed[i]],v=V[ed[i]];
			if(find(u)!=find(v))
			{
				par[find(u)]=find(v);
				tote.push_back(ed[i]);
			}
		}
		double res=reduce_tree(tote,ne);
		if(ans>res)
		{
			ans=res;
      		Ans.clear();
			for(int e:ne)Ans.push_back(e);
		}
	}
	if(ans>1e12)return -1;
	return ans;
}

int numg[1005];
double Ans[1005],Tim[1005];

void work()
{
	vis.clear(),dp.clear(),F.clear();
	vis.resize(maxn);dp.resize(maxn);F.resize(maxn);
	preu2.resize(maxn+5),preid2.resize(maxn+5),dis2.resize(maxn+5);
  	for(int i=0;i<maxn;++i)vis[i].resize(maxS),dp[i].resize(maxS),F[i].resize(maxS);
	freopen(("data/" + graph_file + "/query.txt").c_str(), "r", stdin);
	int Q;
	scanf("%d",&Q);
	for(int cas=1;cas<=Q;++cas)
	{
		cerr<<cas<<" start"<<endl;
		scanf("%d",&k);
		numg[cas]=k;
		cerr<<cas<<" g="<<k<<endl;
		for(int i=1;i<=k;++i)
		{
			scanf("%d",&grp[i].sz);
			for(int j=1;j<=grp[i].sz;++j)
			{
				int x;
				scanf("%d",&x);
				grp[i].a.push_back(x); 
			}
      		sort(grp[i].a.begin(),grp[i].a.end());
      		grp[i].a.erase(unique(grp[i].a.begin(),grp[i].a.end()),grp[i].a.end());
		}
		if(k==1)
		{
			for(int i=1;i<=k;++i)grp[i].a.clear(),grp[i].sz=0;
			Ans[cas]=0;Tim[cas]=0;
			cerr<<"time = "<<Tim[cas]<<" "<<Ans[cas]<<endl;
			continue;
		}
    	auto st = chrono::high_resolution_clock::now();
    	Ans[cas]=PartialOPT(para_h);
    	Tim[cas]=chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - st).count()/1000000.0;
    	cerr<<"time = "<<Tim[cas]<<" "<<Ans[cas]<<endl;
		for(int i=1;i<=k;++i)grp[i].a.clear(),grp[i].sz=0;
	}
	fclose(stdin);
	freopen(("results/" + graph_file + "_PartialOPT_result.txt").c_str(), "w", stdout);
	for(int cas=1;cas<=Q;++cas)printf("%.10f %.10f\n",Tim[cas],Ans[cas]);
}
int main(int argc, char* argv[])
{
	if(argc != 2)
	{
		cout << "Usage: " << argv[0] << " <graphname>" << endl;
		return 1;
	}
	graph_file = argv[1];
	readgraph();
	work();
}
