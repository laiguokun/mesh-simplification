#include "modelsimplify.h"
#include <map>
#include <cmath>
#include <iostream>
#include <set>
using namespace std;
using namespace SimpleOBJ;
mysimple::mysimple(const char* fn,const char* fout,double ratio)
{
	LoadFromObj(fn);
	work(ratio);
	SaveToObj(fout);
}
bool tri[500000];
bool vertex[500000];
double Q[500000][4][4];
T data[500000];
double cnt=0;
set<T> heap;
void mysimple::print(Vec3f vec)
{
	cout<<vec.x<<" "<<vec.y<<" "<<vec.z<<endl;
}
Vec3f mysimple::cj(Vec3f v1,Vec3f v2)
{
	return (Vec3f(v1.y*v2.z-v1.z*v2.y,v1.z*v2.x-v1.x*v2.z,v1.x*v2.y-v1.y*v2.x));
}
double *q=new double[4];
void mysimple::add(int x,double *q)
{
	for (int i=0;i<4;i++)
		for (int j=0;j<4;j++)
			Q[x][i][j]+=q[i]*q[j];
}
vector< vector<int> > p_tri;
void mysimple::pre()
{
	for (int i=0;i<m_nTriangles;i++)
		for (int j=0;j<4;j++)
			for (int k=0;k<4;k++)
				Q[i][j][k]=0;
	for (int i=0;i<m_nTriangles;i++)
	{
		if (tri[i]!=true)
			continue;
		Vec3f t1=m_pVertexList[m_pTriangleList[i][0]];
		Vec3f t2=m_pVertexList[m_pTriangleList[i][1]];
		Vec3f t3=m_pVertexList[m_pTriangleList[i][2]];
		Vec3f vcj=cj(t2-t1,t3-t2);
		vcj.Normalize();
		q[0]=vcj.x;q[1]=vcj.y;q[2]=vcj.z;
		q[3]=-(t1.x*vcj.x+t1.y*vcj.y+t1.z*vcj.z);
		add(m_pTriangleList[i][0],q);
		add(m_pTriangleList[i][1],q);
		add(m_pTriangleList[i][2],q);
	}
	Vec3f vec;
	int x=0;int y=0;
	for (int i=0;i<m_nVertices;i++)
	{
//		if (vertex[i]==false) continue;
		double tmp=calc(i,x,y,vec);
		T t(tmp,x,y,vec);
		heap.insert(t);
		data[i]=t;
	}
}
double Mat[5][5];
double _p[4];
double mysimple::opt(int v1,int v2,Vec3f& vec)
{
	for (int i=0;i<4;i++)
	{
		for (int j=0;j<4;j++)
		{
			Mat[i][j]=Q[v1][i][j]+Q[v2][i][j];
		}
	}
	for (int i=0;i<3;i++) 
	{
		Mat[3][i]=0;
		Mat[i][4]=0;
	}
	Mat[3][3]=1;Mat[3][4]=1;
	for (int i=0;i<3;i++)
	{
		int p=i;
		for (int j=i+1;j<4;j++)
			if (fabs(Mat[p][i])<(fabs(Mat[j][i])))
				p=j;
		if (p!=i)
			for (int j=0;j<=4;j++)
			{
				double t=Mat[i][j];
				Mat[i][j]=Mat[p][j];
				Mat[p][j]=t;
			}
		for (int j=i+1;j<4;j++)
		{
			double t=Mat[j][i]/Mat[i][i];
			for (int k=i;k<=4;k++)
				Mat[j][k]=Mat[j][k]-t*Mat[i][k];
		}
	}
	for (int i=0;i<4;i++) _p[i]=0;
	_p[3]=Mat[3][4]/Mat[3][3];
	for (int i=2;i>=0;i--)
	{
		double sum=Mat[i][4];
		for (int j=3;j>i;j--)
			sum-=_p[j]*Mat[i][j];
		_p[i]=sum/Mat[i][i];
	}
	vec.x=_p[0];vec.y=_p[1];vec.z=_p[2];
	double res=0;
	for (int i=0;i<4;i++)
		for (int j=0;j<4;j++)
			res+=_p[i]*_p[j]*(Q[v1][i][j]+Q[v2][i][j]);
	return res;
}
map<int,int> Map;
map<int,int>::iterator it;
map<int,int> update_map;
map<int,int> m_;
double mysimple::calc(int x,int& t1,int& t2,Vec3f& opt_vec)
{
	Map.clear();
	int s=p_tri[x].size();
	double min=1e+18;
	Vec3f vec;
	for (int i=0;i<s;i++)
	{
		int t=p_tri[x][i];
		if (tri[t]==false) continue;
		for (int j=0;j<3;j++)
		{
			int now=m_pTriangleList[t][j];
			if (now==x) continue;
			if (vertex[now]==false) continue;
			if (x==90102 && now==89068)
			{
//				cout<<"------------"<<endl;
//				for (int k=0;k<p_tri[89068].size();k++)
//					cout<<p_tri[89068][k]<<endl;
//				cout<<t<<endl;
			}
			if (Map.find(now)!=Map.end()) continue;
			Map[now]=1;
			double tmp=opt(x,now,vec);
			if (tmp<min)
			{
				t1=x;
				t2=now;
				min=tmp;
				opt_vec=vec;
			}
		}
	}
	return min;
}

void mysimple::getmin(int& v1,int& v2,Vec3f& opt_vec)
{
	set<T>::iterator it;
	it=heap.begin();
	v1=it->v1;v2=it->v2;opt_vec=it->vec;
//	cout<<it->f<<endl;
}
void mysimple::update(int v)
{
	for (int j=0;j<4;j++)
		for (int k=0;k<4;k++)
			Q[v][j][k]=0;
	for (int i=0;i<p_tri[v].size();i++)
	{
		int now=p_tri[v][i];
		if (tri[now]==false) continue;
		Vec3f t1=m_pVertexList[m_pTriangleList[now][0]];
		Vec3f t2=m_pVertexList[m_pTriangleList[now][1]];
		Vec3f t3=m_pVertexList[m_pTriangleList[now][2]];
		Vec3f vcj=cj(t2-t1,t3-t2);
		vcj*=10000;
		vcj.Normalize();
		q[0]=vcj.x;q[1]=vcj.y;q[2]=vcj.z;
		q[3]=-(t1.x*vcj.x+t1.y*vcj.y+t1.z*vcj.z);
		add(v,q);
	}
}	
void mysimple::del(int v1,int v2,Vec3f vec)
{
	int new_v=v1;
//	if (v2==90102)
//		cout<<v1<<" "<<v2<<endl;
/*	if (v1==91992 && v2==89068)
	{
		cout<<cnt<<endl;
		cout<<data[90102].v2<<endl;
	}*/
	vertex[v2]=false;
	heap.erase(data[v2]);
	int t=0;
	m_.clear();
	for (int i=0;i<p_tri[v2].size();i++)
	{
		int now=p_tri[v2][i];
		if (tri[now]==false) continue;
		t=0;
		for (int j=0;j<3;j++)
		{
			if (m_pTriangleList[now][j]==v1 or m_pTriangleList[now][j]==v2)
			{
				t++;
				m_pTriangleList[now][j]=new_v;
			}
		}
		if (t==2) 
		{
			tri[now]=false;
			cnt--;
		}
		else
			p_tri[v1].push_back(now);
/*		if (t==0)
		{
			for (int k=0;k<3;k++)
				cout<<m_pTriangleList[54895][k]<<endl;
			cout<<cnt<<endl;
			cout<<now<<endl;
			system("pause");
		}*/
	}
	m_pVertexList[v1]=vec;
	update_map.clear();
	for (int i=0;i<p_tri[v1].size();i++)
	{
		int now=p_tri[v1][i];
		for (int j=0;j<3;j++)
		{
			int t=m_pTriangleList[now][j];
			if (vertex[t]==false) continue;
			if (update_map.find(t)!=update_map.end())
				continue;
			update_map[t]=1;
			update(t);
		}
	}
	for (it=update_map.begin();it!=update_map.end();it++)
	{
		int v=it->first;
		for (int i=0;i<p_tri[v].size();i++)
		{
			int now=p_tri[v][i];
			for (int j=0;j<3;j++)
			{
				int t=m_pTriangleList[now][j];
				if (vertex[t]==false) continue;
				if (m_.find(t)!=m_.end()) continue;
				m_[t]=1;
//				if (v1==91992 && v2==89068 && t==90102)
//					cout<<"fuck"<<data[t].v2<<endl;
				heap.erase(data[t]);
				int x,y;Vec3f vec;
				double tmp=calc(t,x,y,vec);
				T tt(tmp,x,y,vec);
				heap.insert(tt);
				data[t]=tt;
			}
		}
	}			
}
void mysimple::work(double ratio)
{
	cnt=m_nTriangles;
	for (int i=0;i<m_nVertices;i++)
	{
		vector<int> tmp;
		p_tri.push_back(tmp);
		vertex[i]=true;
	}
	for (int i=0;i<m_nTriangles;i++) 
	{
		tri[i]=true;
		p_tri[m_pTriangleList[i][0]].push_back(i);
		p_tri[m_pTriangleList[i][1]].push_back(i);
		p_tri[m_pTriangleList[i][2]].push_back(i);
	}
//	for (int i=0;i<3;i++)
//		cout<<m_pTriangleList[54895][i]<<endl;
//	cout<<"yes"<<endl;
	pre();
	while (cnt/m_nTriangles>ratio)
	{
		int v1,v2;
		Vec3f vec;
		getmin(v1,v2,vec);
//		if (cnt==204501)
//			cout<<v1<<" "<<v2<<endl;
		del(v1,v2,vec);
//		cout<<cnt<<endl;
	}
//	cout<<cnt<<endl;
	int now=0;
	for (int i=0;i<m_nTriangles;i++)
	{
		if (tri[i]==false) continue;
		m_pTriangleList[now]=m_pTriangleList[i];
		now++;
	}
	m_nTriangles=now;	
}
