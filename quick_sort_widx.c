#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NN 10

void quick_sort_widx( double *a, int *idx, int n )
{
	int i,j,k,l=0,r=n-1;
	double temp,a0;
	int idt,id0;
	int nstack=2*8*sizeof(nstack);
//	int nstack=2*log2(n);
//	if ( nstack<50 ) nstack=50;
	int istack[nstack],jstack=0;
	
	while ( 1 ){
		if ( r-l<NN ){
			for ( j=l+1;j<=r;++j ){
				a0 = a[j];
				id0 = idx[j];
				for ( i=j-1;i>=l;--i ){
					if ( a[i]<a0 ) break;
					a[i+1] = a[i];
					idx[i+1] = idx[i];
				}
				a[i+1] = a0;
				idx[i+1] = id0;
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
			idt=idx[k];
			idx[k]=idx[l+1];
			idx[l+1]=idt;
		if ( a[l]>a[r] ){
			temp=a[l];
			a[l]=a[r];
			a[r]=temp;
			idt=idx[l];
			idx[l]=idx[r];
			idx[r]=idt;
		}
		if ( a[l+1]>a[r] ){
			temp=a[l+1];
			a[l+1]=a[r];
			a[r]=temp;
			idt=idx[l+1];
			idx[l+1]=idx[r];
			idx[r]=idt;
		}
		if ( a[l]>a[l+1] ){
			temp=a[l];
			a[l]=a[l+1];
			a[l+1]=temp;
			idt=idx[l];
			idx[l]=idx[l+1];
			idx[l+1]=idt;
		}
		i = l+1;
		j = r;
		a0 = a[l+1];
		id0 = idx[l+1];
		while ( 1 ){
			for ( ; (a[++i]<a0) ; );
			for ( ; (a[--j]>a0) ; );
			if ( j<i ) break;
			temp=a[i];
			a[i]=a[j];
			a[j]=temp;
			idt=idx[i];
			idx[i]=idx[j];
			idx[j]=idt;
		}
		a[l+1] = a[j];
		a[j] = a0;
		idx[l+1] = idx[j];
		idx[j] = id0;
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
