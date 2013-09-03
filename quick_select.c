#define SWAP(a,b,type) \
({\
	type temp;\
	temp=a;\
	a=b;\
	b=temp;\
})

double quick_select(double *a, int k, int n)
{
	int i,j;
	int l=0,r=n-1;
	double a0;
	double ak;

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
				while (a[++i]<a0);
				while (a[--j]>a0);
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

double quick_select_widx(double *a, int *idx, int k, int n, int *ik)
{
	int i,j;
	int l=0,r=n-1;
	double a0;
//	double ftemp;
	int id0;
//	int itemp;
	double ak;

	while (1){
		if (r-l<=1){
			if (1==r-l){
				if (a[l]>a[r]){
					SWAP(a[l],a[r],double);
					SWAP(idx[l],idx[r],int);
				}
			}
			ak = a[k];
			*ik = idx[k];
			break;
		} else {
			i = (l+r)/2;
			SWAP(a[i],a[l+1],double);
			SWAP(idx[i],idx[l+1],int);
			if (a[l]>a[r]){
				SWAP(a[l],a[r],double);
				SWAP(idx[l],idx[r],int);
			}
			if (a[l+1]>a[r]){
				SWAP(a[l+1],a[r],double);
				SWAP(idx[l+1],idx[r],int);
			}
			if (a[l]>a[l+1]){
				SWAP(a[l],a[l+1],double);
				SWAP(idx[l],idx[l+1],int);
			}
			i=l+1;
			j=r;
			a0=a[l+1];
			id0=idx[l+1];
			while (1){
				while (a[++i]<a0);
				while (a[--j]>a0);
				if (j<i) break;
				SWAP(a[i],a[j],double);
				SWAP(idx[i],idx[j],int);
			}
			a[l+1]=a[j];
			a[j]=a0;
			idx[l+1]=idx[j];
			idx[j]=id0;
			if (j>=k) r=j-1;
			if (j<=k) l=i;
		}
	}
	return ak;
}
