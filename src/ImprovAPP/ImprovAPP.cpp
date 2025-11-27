#include<bits/stdc++.h>
#pragma GCC optimize(2) 
#define maxm 30000005
#define maxn 15000005
using namespace std;
const double inf = 1e18;
string graph_file;
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
	for(int i=0;i<n+5;++i)g[i].clear();
	while(~scanf("%d%d%lf",&u,&v,&w))
	{
		++idcnt;
		g[u].push_back(edge(v,w,idcnt));
		g[v].push_back(edge(u,w,idcnt));
		U[idcnt]=u;V[idcnt]=v;W[idcnt]=w;
	}
	fclose(stdin);
}
vector<int> used;
vector<double> dis2[12];
vector<int> preu2[12],preid2[12];
void sssp(vector<int> s,vector<double> &dis,vector<int> &preu,vector<int> &preid)
{
	used.clear();
	used.resize(n+5);
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
void init()
{
	sort(grp+1,grp+k+1);
	for(int i=2;i<=k;++i)
	{
		vector<int> st;
		st.clear();
		for(int x:grp[i].a)st.push_back(x);
		sssp(st,dis2[i],preu2[i],preid2[i]);
	}
}

void addpath2(vector<int> &e,int r,int u)
{
	int x=u;
	while(preu2[r][x])
	{
		e.push_back(preid2[r][x]);
		x=preu2[r][x]; 
	}
}

int cov[12],bel[12];
int eid[maxm];
bool cmp(const int &x,const int &y)
{
	return W[x]<W[y];
}
int par[maxn];
int find(int x)
{
	if(par[x]==x)return x;
	return par[x]=find(par[x]);
}
int d[maxn],g_cnt[maxn];
vector< pair<int,int> > G[maxn];
vector<int> has[maxn];
double ImprovAPP()
{
	init();
	double ans=inf;
	vector<int> Ans,Ans2,Ans3;
	double mn=inf;
	for(int r:grp[1].a)
	{
		bool ok=1;
		for(int i=2;i<=k;++i)if(dis2[i][r]>1e12)ok=0;
		if(!ok)continue;
		for(int i=1;i<=k;++i)cov[i]=0;
		priority_queue< pair<double,int> > pq;
		for(int i=2;i<=k;++i)pq.push(make_pair(-dis2[i][bel[i]=r],i));
		vector<int> path;
		while(!pq.empty())
		{
			int o=pq.top().second;
			pq.pop();
			if(cov[o])continue;
			addpath2(path,o,bel[o]);
			cov[o]=1;
			int x=bel[o];
			while(x)
			{
				for(int j=2;j<=k;++j)if(!cov[j])
				{
					if(dis2[j][x]<dis2[j][bel[j]])
					{
						pq.push(make_pair(-dis2[j][x],j));
						bel[j]=x;
					}
				}
				x=preu2[o][x];
			}
		}
		sort(path.begin(),path.end());
		path.erase(unique(path.begin(),path.end()),path.end());
		double totw=0;
		for(int e:path)totw+=W[e];
		if(totw<ans)
		{
			ans=totw;
			Ans.clear();
			for(int e:path)Ans.push_back(e);
		}
	}
	if(ans>1e12)return -1;
	int m=Ans.size();
	for(int i=0;i<m;++i)eid[i]=Ans[i];
	for(int i=1;i<=n;++i)par[i]=i;
	sort(eid,eid+m,cmp);
	for(int i=0;i<m;++i)
	{
		int u=U[eid[i]],v=V[eid[i]];
		if(find(u)!=find(v))
		{
			par[find(u)]=find(v);
			Ans2.push_back(eid[i]);
		}
	}
	for(int i=1;i<=n;++i)d[i]=-1,G[i].clear(),has[i].clear();
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

int numg[1005]; 
double Ans[1005],Tim[1005];
void work()
{
	for(int i=0;i<=10;++i)preu2[i].resize(maxn+5),preid2[i].resize(maxn+5),dis2[i].resize(maxn+5);
	freopen(("data/" + graph_file + "/query.txt").c_str(), "r", stdin);
	int Q;
	scanf("%d",&Q);
	for(int cas=1;cas<=Q;++cas)
	{
		cerr<<cas<<" start"<<endl;
		scanf("%d",&k);
		cerr<<cas<<" g="<<k<<endl;
		numg[cas]=k;
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
		Ans[cas]=ImprovAPP();
		Tim[cas]=chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - st).count()/1000000.0;
		cerr<<"time = "<<Tim[cas]<<" "<<Ans[cas]<<endl;
		for(int i=1;i<=k;++i)grp[i].a.clear(),grp[i].sz=0;
	}
	fclose(stdin);
	freopen(("results/" + graph_file + "_ImprovAPP_result.txt").c_str(), "w", stdout);
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
