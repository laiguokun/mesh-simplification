#include "SimpleObject.h"
using namespace std;
using namespace SimpleOBJ;
class mysimple : public CSimpleObject
{
public:
	mysimple(const char* fn,const char* fout,double ratio);
	~mysimple(){};
private:
	void work(double);
	Vec3f cj(Vec3f ,Vec3f);
	void add(int ,double*);
	void pre();
	double opt(int,int,Vec3f&);
	double calc(int,int&,int&,Vec3f&);
	void getmin(int&,int&,Vec3f& vec);
	void del(int,int,Vec3f);
	void print(Vec3f vec);
	void update(int);
};
struct T
{
	double f;
	int v1,v2;
	Vec3f vec;
	T(double f_,int v1_,int v2_,Vec3f vec_)
	{
		f=f_;v1=v1_;v2=v2_;vec=vec_;
	}
	T(){}
	bool operator <(const T r) const
	{
		if (f<r.f) return true;
		if (f==r.f && v1<r.v1) return true;
		return false;
	}
};
