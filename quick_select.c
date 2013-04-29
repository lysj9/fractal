#define SWAP(a,b,type) do\
{\
	type temp;\
	temp=a;\
	a=b;\
	b=temp;\
} while(0)

double quick_select(double *a, int k, int n)
{
	double a0;
//	double ftemp;
	double ak;
	int i,j,l,r;

	l=0;
	r=n-1;
	while (1){
		if (r-l<=1){
			if (1==r-l){
				if (a[l]>a[r]){
					SWAP(a[l],a[r],double);
				}
			}
			ak = a[k];
			break;
		} else {
			i = (l+r)/2;
			SWAP(a[i],a[l+1],double);
			if (a[l]>a[r]){
				SWAP(a[l],a[r],double);
			}
			if (a[l+1]>a[r]){
				SWAP(a[l+1],a[r],double);
			}
			if (a[l]>a[l+1]){
				SWAP(a[l],a[l+1],double);
			}
			i=l+1;
			j=r;
			a0=a[l+1];
			while (1){
				while (a[++i]<=a0);
				while (a[--j]>=a0);
				if (j<i) break;
				SWAP(a[i],a[j],double);
			}
			a[l+1]=a[j];
			a[j]=a0;
			if (j>=k) r=j-1;
			if (j<=k) l=i;
		}
	}
	return ak;
}
