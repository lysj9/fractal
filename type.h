#ifndef _TYPE_H_
#define _TYPE_H_

struct node{
	double m;
	double x,y,z,vx,vy,vz;
	struct node *next;
};

struct vector_s{
	double m,x,y,z,vx,vy,vz;
};

#endif
