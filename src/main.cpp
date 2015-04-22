#include "modelsimplify.h"
#include <ctime>
#include <iostream>
using namespace std;
using namespace SimpleOBJ;
int main()
{
	double start=clock();
//	mysimple A("kitten.50k.obj","output.obj",0.5);
	mysimple A("fixed.perfect.dragon.100K.0.07.obj","output.obj",0.01);
	double end=clock();
	cout<<start-end<<endl;
//	mysimple A("cube.obj","output.obj",0.4);
}
