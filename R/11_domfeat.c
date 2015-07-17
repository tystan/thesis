#include <R.h> 

void domfeat(int *n, double *obja, double *objb, int *domvec)
{
	int i, j, nonDomI;
	for (i=0; i<*n; i++)
	{
		j=0;
		nonDomI=1;
		while(nonDomI && j<*n)
		{
			if((obja[i]<obja[j]) && (objb[i]<=objb[j])) nonDomI=0;
			else if((objb[i]<objb[j]) && (obja[i]<=obja[j])) nonDomI=0;
			j+=1;
		}
		domvec[i]=-nonDomI+1;
	}
}


