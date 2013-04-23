#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NN 10

void quick_sort( double *a, int n )
{
	int i,j,k,l=0,r=n-1;
	double temp,a0;
	int nstack=2*8*sizeof(n);
//	int nstack=2*log2(n);
//	if ( nstack<50 ) nstack=50;
	int istack[nstack],jstack=0;
//	int NN=10;

	while ( 1 ){
		if ( r-l<NN ){
			for ( j=l+1;j<=r;++j ){
				a0 = a[j];
				for ( i=j-1;i>=l;--i ){
					if ( a[i]<a0 ) break;
					a[i+1] = a[i];
				}
				a[i+1] = a0;
			}
			if ( jstack <= 0 ) break;
			jstack -= 2;
			l = istack[jstack];
			r = istack[jstack+1];
			continue;
		}
		k = ( l+r )/2;
		temp=a[k];
		a[k]=a[l+1];
		a[l+1]=temp;
		if ( a[l]>a[r] ){
			temp=a[l];
			a[l]=a[r];
			a[r]=temp;
		}
		if ( a[l+1]>a[r] ){
			temp=a[l+1];
			a[l+1]=a[r];
			a[r]=temp;
		}
		if ( a[l]>a[l+1] ){
			temp=a[l];
			a[l]=a[l+1];
			a[l+1]=temp;
		}
		i = l+1;
		j = r;
		a0 = a[l+1];
		while ( 1 ){
			for ( ; (a[++i]<a0) ; );
			for ( ; (a[--j]>a0) ; );
			if ( j<i ) break;
			temp=a[i];
			a[i]=a[j];
			a[j]=temp;
		}
		a[l+1] = a[j];
		a[j] = a0;
		if ( r-i+1 >= j-l ){
			istack[jstack]=i;
			istack[jstack+1]=r;
			r=j-1;
		} else {
			istack[jstack]=l;
			istack[jstack+1]=j-1;
			l=i;
		}
		jstack += 2;
	}

	return;
}
