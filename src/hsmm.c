#include "stdio.h"
#include <R.h>
#include "stdlib.h"

#define int8 unsigned char

void checkmem(void *x) {
  if(x==NULL) error("Out of memory.");
}

int min(int a,int b) {
  if(a<b) return(a);
  else return(b);	
}


void sim_mc(double *start, double *a,int *nstates,int *state,int *T,int *nseq) {
  int t,i,n,*s=NULL;	
  int K = *nstates;
  int N = *nseq;
  double tmp;
	
  GetRNGstate();

  for(n=0;n<N;n++) {
    if(n==0) s = state;
    else s = s+T[n-1];
    tmp=unif_rand();
    i=0;
    while(tmp>start[i])	i++;
    s[0]=i+1;	

    for(t=1;t<T[n];t++) {
      tmp=unif_rand();
      i=0;
      while(tmp>a[i*K+s[t-1]-1])	i++;
      s[t]=i+1;	
    }	
  }

  PutRNGstate();
}

void **alloc_matrix(int nrow,int ncol,int size) {
  int i;
  void **x = malloc(sizeof(void *)*nrow);
  checkmem(x);
  for(i=0;i<nrow;i++)	{
    x[i]=malloc(size*ncol);
    checkmem(x[i]);
  }
  return(x);
}

void free_matrix(int nrow,int ncol,void **x) {
  int i;
  for(i=0;i<nrow;i++)	free(x[i]);
  free(x);
}

void print_matrix(int nrow,int ncol,double *x) {
  int i,j;
  for(i=0;i<nrow;i++) {
    for(j=0;j<ncol;j++)		
      Rprintf("%.3g\t",x[i*ncol+j]);
    //			Rprintf("%.3f\t",x[j*nrow+i]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

void print_matrix2(int nrow,int ncol,double **x) {
  int i,j;
  for(i=0;i<nrow;i++) {
    for(j=0;j<ncol;j++)		
      Rprintf("%.3g\t",x[i][j]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

void print_imatrix2(int nrow,int ncol,int **x) {
  int i,j;
  for(i=0;i<nrow;i++) {
    for(j=0;j<ncol;j++)		
      Rprintf("%d\t\t",x[i][j]);
    Rprintf("\n");
  }
  Rprintf("\n");
}


void forward(double *a,double *start,double *p0,double *d,double *D,int *timelength,int *nstates,int *M,
	     double **F0,double *N0,double **si0,int *nsequences,int *totallength) {
  double obs;
  int u,i,j,t,n;
  int T;
  int J = *nstates;	
  double **p,**F,**si,*N;
  int offset;
  int ln = *totallength;
  int nseq = *nsequences;

  p = (double **)malloc(sizeof(double *)*J);
  F = (double **)malloc(sizeof(double *)*J);
  si = (double **)malloc(sizeof(double *)*J);

  /*
    printf("J = %d     T = %d\n",J,T);
    print_matrix(J,J,a);
    print_matrix(1,J,start);	
	
    for(j=0;j<J;j++)  printf("M[%d] = %d\t",j,M[j]);
    printf("\nd = \n");
    print_matrix(J,100,d);
    printf("p = \n");
    print_matrix(J,T,p);	   
    printf("D = \n");
    print_matrix(J,100,D);	
  */
  for(j=0;j<J;j++) {
    p[j]  = p0 + j*ln;
    F[j]  = F0[j];
    si[j] = si0[j];
  }
  N = N0;

  for(n=0;n<nseq;n++) {	
    T = timelength[n];
    if(n>0){
      offset = timelength[n-1];
      for(j=0;j<J;j++) {
	p[j]  += offset;
	F[j]  += offset;
	si[j] += offset;
      }
      N += offset;
    }
    for(t=0;t<T;t++) {
      //		printf("\nt = %d\n",t);
      N[t]=0;
      for(j=0;j<J;j++) {
	F[j][t]=0;
	obs=p[j][t];
	if(t<T-1) {
	  //				printf("j = %d   ",j);
	  for(u=1;u<=min(t+1,M[j]);u++) {
	    if(u<t+1) {
	      F[j][t] += obs*d[j*M[j]+u-1]*si[j][t-u+1];
	      N[t] += obs*D[j*M[j]+u-1]*si[j][t-u+1];
	      obs*= p[j][t-u]/N[t-u];
	      //						printf("%.3g (%.3g) \t",obs,N[t-u]);
	    }
	    else {
	      F[j][t] += obs*d[j*M[j]+t]*start[j];
	      N[t] += obs*D[j*M[j]+t]*start[j];						
	    }
	  }
	  //				printf("\n");
	}
	else{
	  for(u=1;u<=min(t+1,M[j]);u++) {
	    if(u<T) {
	      F[j][T-1] += obs*D[j*M[j]+u-1]*si[j][T-u];
	      obs *= p[j][T-1-u]/N[T-1-u];
	    }
	    else F[j][T-1] += obs*D[j*M[j]+T-1]*start[j];
	  }
	  N[T-1] += F[j][T-1];
	}
      }

      for(j=0;j<J;j++){
	F[j][t] /= N[t];
	F[j][t]+=1e-300;
      }

      if(t<T-1) {
	for(j=0;j<J;j++){
	  si[j][t+1]=0;
	  for(i=0;i<J;i++) 
	    si[j][t+1]+=a[j*J+i]*F[i][t];
	}
      }
    }
  }
  free(si);
  free(p);
  free(F);
  /*
    printf("\nF=\n");
    print_matrix2(J,T,F);
    print_matrix(1,T,N);
  */
}

void backward(double *a,double *start,double *p0,double *d,double *D,int *timelength,int *nstates,int *M,double *L10,
	      double *N0,double *eta,double *F1,double *statein,double *gamma,int *nsequences,int *totallength,double *Gret) {
	
  double **L0,**G0,**F0,**si0,**L,**G,obs,**si,**F,**num,*den,*N,**p,**L1;
  int u,i,j,k,t,T,n;	
  int J = *nstates;
  int offset;
  int ln = *totallength;
  int nseq = *nsequences;

  F0 =   (double **)alloc_matrix(J,ln,sizeof(double));
  si0  = (double **)alloc_matrix(J,ln,sizeof(double));   
  den  = (double *)malloc(sizeof(double)*J);
  num  = (double **)alloc_matrix(J,J,sizeof(double));   

  forward(a,start,p0,d,D,timelength,nstates,M,F0,N0,si0,nsequences,totallength);
  /*
    for(j=0;j<J;j++) memcpy(F1+j*T,F[j],T*sizeof(double));
    for(j=0;j<J;j++) memcpy(statein+j*T,si[j],T*sizeof(double));
    printf("J = %d     T = %d\n",J,T);	
    for(j=0;j<J;j++)  printf("M[%d] = %d\t",j,M[j]);
    printf("\n");
    print_matrix(1,T,N);
    print_matrix2(J,50,F);
    printf("\n");
    print_matrix2(J,50,si);
  */

  //	L1 =  (double **)alloc_matrix(J,T,sizeof(double));
  L0  =  (double **)alloc_matrix(J,ln,sizeof(double));
  G0  =  (double **)alloc_matrix(J,ln,sizeof(double));	

  p = (double **)malloc(sizeof(double *)*J);
  F = (double **)malloc(sizeof(double *)*J);
  G = (double **)malloc(sizeof(double *)*J);
  L = (double **)malloc(sizeof(double *)*J);
  si = (double **)malloc(sizeof(double *)*J);
  L1 = (double **)malloc(sizeof(double *)*J);

  /*
    printf("J = %d     T = %d\n",J,T);
    print_matrix(J,J,a);
    print_matrix(1,J,start);	
	
    for(j=0;j<J;j++)  printf("M[%d] = %d\t",j,M[j]);
    printf("\nd = \n");
    print_matrix(J,100,d);
    printf("p = \n");
    print_matrix(J,T,p);	   
    printf("D = \n");
    print_matrix(J,100,D);	
  */
  for(j=0;j<J;j++) {
    p[j]  = p0 + j*ln;
    F[j]  = F0[j];
    si[j] = si0[j];
    G[j] = G0[j];
    L[j] = L0[j];
    L1[j] = L10+j*ln;
  }
  N = N0;
	
  for(n=0;n<nseq;n++) {	
    T = timelength[n];
    if(n>0) {
      offset = timelength[n-1];
      for(j=0;j<J;j++) {
	p[j]  += offset;
	F[j]  += offset;
	si[j] += offset;
	G[j] += offset;
	L[j] += offset;
	L1[j] += offset;
      }
      N += offset;
    }

    for(j=0;j<J;j++) L[j][T-1] = F[j][T-1];

    if(n==0) {
      for(j=0;j<J;j++) 
	for(u=1;u<M[j];u++)
	  eta[j*M[j]+u-1] = 0;
    }

    for(t=T-2;t>=0;t--) {
      //		printf("\nt = %d",t);
      for(j=0;j<J;j++) {
	//			printf("\nj = %d:   ",j);
	G[j][t+1] = 0;
	obs = 1;
	for(u=1;u<=min(T-1-t,M[j]);u++) {				
	  obs *= p[j][t+u]/N[t+u];
	  if(u<T-1-t) {
	    G[j][t+1] += L1[j][t+u]/F[j][t+u]*obs*d[j*M[j]+u-1];
	    eta[j*M[j]+u-1] += L1[j][t+u]/F[j][t+u]*obs*d[j*M[j]+u-1]*si[j][t+1];
	  }
	  else {
	    G[j][t+1] += obs*D[j*M[j]+T-1-t-1];
	    eta[j*M[j]+u-1] += obs*d[j*M[j]+u-1]*si[j][t+1];
	  }
	  if(t==0){
	    if(u>T-1) eta[j*M[j]+u-1] += L1[j][t+u]/F[j][t+u]*obs*d[j*M[j]+u-1]*start[j];
	    else eta[j*M[j]+u-1] += obs*d[j*M[j]+u-1]*start[j];
	  }				
	  //				printf(" %.3g/%.3g = %.3g\t",L1[j][t+u],F[j][t+u],L1[j][t+u]/F[j][t+u]);
	}
      }		
      for(j=0;j<J;j++){
	L1[j][t] = 0;
	for(k=0;k<J;k++) L1[j][t] += G[k][t+1]*a[k*J+j];
	L1[j][t] *= F[j][t];
	L[j][t] = L1[j][t] + L[j][t+1] - G[j][t+1]*si[j][t+1];
      }
    }
  }

  //reset pi
  for(i=0;i<J;i++) {
    start[i]=0;
    den[i] = 0;
    for(j=0;j<J;j++) num[i][j] = 0;
  }
  //new estimates for a and pi

  for(i=0;i<J;i++){
    for(n=0;n<nseq;n++) {	
      T = timelength[n];
      if(n==0) {
	for(j=0;j<J;j++) {
	  F[j]  = F0[j];
	  G[j] = G0[j];
	  L[j] = L0[j];
	  L1[j] = L10+j*ln;
	}
      }
      else {
	offset = timelength[n-1];
	for(j=0;j<J;j++) {
	  F[j]  += offset;
	  G[j] += offset;
	  L[j] += offset;
	  L1[j] += offset;
	}
      }
      start[i]+=L[i][0];			
      for(t=0;t<T-2;t++) den[i] += L1[i][t];
      for(j=0;j<J;j++){
	for(t=0;t<T-2;t++) num[i][j] += G[j][t+1]*a[j*J+i]*F[i][t];				
      }		
    }
  }
	
  for(i=0;i<J;i++) {
    start[i]/=nseq;
    for(j=0;j<J;j++) a[j*J+i]=num[i][j]/den[i];
  }

  /*
    printf("\n\nG=\n");
    print_matrix2(J,T,G);	
    printf("\nL=\n");
    print_matrix2(J,T,L);
    printf("\nL1=\n");
    print_matrix(J,T,L1);
    printf("\n");
  */

  for(j=0;j<J;j++) {
    memcpy(gamma+j*ln,L0[j],ln*sizeof(double));
    memcpy(F1+j*ln,F0[j],ln*sizeof(double));
    memcpy(Gret+j*ln,G0[j],ln*sizeof(double));
  }

  //	free_matrix(J,T,(void **)L1);
  free_matrix(J,ln,(void **)F0);
  free_matrix(J,ln,(void **)G0);
  free_matrix(J,ln,(void **)si0);
  free_matrix(J,ln,(void **)L0);
  free_matrix(J,J,(void **)num);
  free(den);
  free(p);
  free(F);
  free(G);
  free(L);
  free(si);
  free(L1);
}

/*
  void viterbi(double *a,double *start,double *p0,double *d,double *D,int *timelength,int *nstates,int *M,
  double *alpha0,int *nsequences,int *totallength,int *statehat) {
  double obs;
  int u,i,j,t,n,T;
  int J = *nstates;	
  double **p,**alpha,**si,**si0,tmp_max,tmp1;
  int offset,*q;
  int ln = *totallength;
  int nseq = *nsequences;
  int **psi_time;
  int **psi_state;
  int *psi_state0, *psi_time0;


  Rprintf("J = %d     Total = %d\n",J,ln);	
  for(j=0;j<J;j++)  Rprintf("M[%d] = %d\t",j,M[j]);
  Rprintf("\n");	
  for(j=0;j<nseq;j++)  Rprintf("T[%d] = %d\t",j,timelength[j]);
  Rprintf("\n");	


  psi_time0 = (int *)malloc(ln*J*sizeof(int));
  psi_state0 = (int *)malloc(ln*J*sizeof(int));
  si0 = (double **)alloc_matrix(J,ln,sizeof(double));
  psi_time = (int **)malloc(J*sizeof(int));
  psi_state = (int **)malloc(J*sizeof(int));
  p = (double **)malloc(sizeof(double *)*J);
  alpha = (double **)malloc(sizeof(double *)*J);
  si = (double **)malloc(sizeof(double *)*J);

  print_matrix(J,ln,p0);
  print_matrix(J,M[0],d);

  for(j=0;j<J;j++) {
  p[j]  = p0 + j*ln;
  alpha[j]  = alpha0 + j*ln;
  si[j] = si0[j];
  psi_time[j] = psi_time0 + j*ln;
  psi_state[j] = psi_state0 + j*ln;
  }

  for(n=0;n<nseq;n++) {	
  T = timelength[n];
  if(n>0){
  offset = timelength[n-1];
  for(j=0;j<J;j++) {
  p[j]  += offset;
  alpha[j] += offset;
  psi_state[j]  += offset;
  psi_time[j]  += offset;
  si[j] += offset;
  }
  }
  for(t=0;t<T;t++) {
  Rprintf("t = %d\n",t);
  for(j=0;j<J;j++) {				
  obs=0;
  Rprintf("alpha[%d][%d]",j,t);
  if(t<T-1) {
  for(u=1;u<=min(t+1,M[j]);u++) {
  if(u<t+1) {
  tmp1 = obs+d[j*M[j]+u-1]+si[psi_state[j][t-u+1]][t-u+1];
  if(u==1 || tmp_max < tmp1) {
  tmp_max = tmp1;
  psi_time[j][t]=u;
  }
  obs += p[j][t-u];
  }
  else {
  tmp1 = obs+d[j*M[j]+t]+start[j];
  Rprintf(" = max(%.3g,%.3g)",tmp_max,tmp1);
  if(u==1 || tmp_max < tmp1) {
  tmp_max = tmp1;
  psi_time[j][t]=u;
  }
  }
  }
  alpha[j][t] = tmp_max+p[j][t];
  tmp_max=0;
  Rprintf(" + %.3g = %.3g\npsi_time[%d][%d]  = %d\n",p[j][t],alpha[j][t],j,t,psi_time[j][t]);
  }
  else{
  for(u=1;u<=min(t+1,M[j]);u++) {
  if(u<T) {
  //							Rprintf(" = max(%.3g,%.3g)",j,t,tmp_max,tmp1);
  tmp1 = obs+D[j*M[j]+u-1]+si[psi_state[j][t-u+1]][t-u+1];
  if(u==1 || tmp_max < tmp1) {
  tmp_max = tmp1;
  psi_time[j][t]=u;
  }
  obs += p[j][T-1-u];
  }
  else {
  tmp1 = obs+D[j*M[j]+T-1]+start[j];							
  if(u==1 || tmp_max < tmp1) {
  tmp_max = tmp1;
  psi_time[j][t]=u;
  }
  }						
  }
  alpha[j][t] = tmp_max+p[j][t];
  Rprintf(" = %.3g + %.3g = %.3g\npsi_time[%d][%d]  = %d\n",p[j][t],alpha[j][t],alpha[j][t]-p[j][t],j,t,psi_time[j][t]);
  }
  }

  if(t<T-1) {
  for(j=0;j<J;j++){
  i=0;
  si[j][t+1]=a[j*J+i]+alpha[i][t];
  psi_state[j][t+1]=0;
  for(i=1;i<J;i++) 
  if(i!=j) {
  tmp1 = a[j*J+i]+alpha[i][t];
  if(si[j][t+1] < tmp1) {
  si[j][t+1] = tmp1;
  psi_state[j][t+1]=i;
  }
  }
  Rprintf("psi_state[%d][%d] = %d\tsi[%d][%d] = %.3g\n",j,t+1,psi_state[j][t+1],j,t+1,si[j][t+1]);
  }
  }
  }
  }

  //and now we backtrack!

  for(j=0;j<J;j++) {
  p[j]  = p0 + j*ln;
  alpha[j]  = alpha0 + j*ln;
  si[j] = si0[j];
  psi_time[j] = psi_time0 + j*ln;
  psi_state[j] = psi_state0 + j*ln;
  }

  q=statehat;

  for(n=0;n<nseq;n++) {	
  T = timelength[n];
  if(n>0){
  offset = timelength[n-1];
  q+=offset;
  for(j=0;j<J;j++) {
  p[j]  += offset;
  alpha[j] += offset;
  psi_state[j]  += offset;
  psi_time[j]  += offset;
  si[j] += offset;
  }
  }
  q[T-1] = 0;
  for(j=1;j<J;j++) if(alpha[q[T-1]][T-1] < alpha[j][T-1]) q[T-1] = j;
  u=1;
  for(t=T-2;t>=0;t--) {
  if(u < psi_time[q[t+u]][t+u]) {
  q[t] = q[t+u];
  u++;
  }
  else {
  q[t] = psi_state[q[t+u]][t+u];
  u=1;
  }
  }
  }

  free(si);
  free_matrix(J,ln,(void **)si0);
  free(p);
  free(alpha);
  free(psi_time);
  free(psi_state);
  free(psi_state0);
  free(psi_time0);
  }
*/


void viterbi(double *a,double *start,double *p0,double *d0,double *D0,int *timelength,int *nstates,int *M,
	     double *alpha0,int *statehat,int *psi_state0,int *psi_time0) {
  double obs;
  int u,i,j,t,T;
  int J = *nstates;	
  double **p,**d,**D,**alpha,**si,**si0,tmp_max=-1e300,tmp1;
  int *q;
  int **psi_time;
  int **psi_state;
  //	int *psi_state0, *psi_time0;
  T = *timelength;


  //	Rprintf("J = %d\n",J);	
  //	for(j=0;j<J;j++)  Rprintf("M[%d] = %d\t",j,M[j]);
  //	Rprintf("\n");	
	

  //	psi_time0 = (int *)malloc(T*J*sizeof(int));
  //	psi_state0 = (int *)malloc(T*J*sizeof(int));
  si0 = (double **)alloc_matrix(J,T,sizeof(double));
  psi_time = (int **)malloc(J*sizeof(int *));
  psi_state = (int **)malloc(J*sizeof(int *));
  p = (double **)malloc(sizeof(double *)*J);
  d = (double **)malloc(sizeof(double *)*J);
  D = (double **)malloc(sizeof(double *)*J);
  alpha = (double **)malloc(sizeof(double *)*J);
  si = (double **)malloc(sizeof(double *)*J);

  checkmem(si0);
  checkmem(psi_time);
  checkmem(psi_state);
  checkmem(p);
  checkmem(d);
  checkmem(D);
  checkmem(alpha);
  checkmem(si);


  //	print_matrix(J,T,p0);
  //	print_matrix(J,M[0],d);

  for(j=0;j<J;j++) {
    d[j] = d0 + j*M[j];
    D[j] = D0 + j*M[j];		
    p[j]  = p0 + j*T;
    alpha[j]  = alpha0 + j*T;
    si[j] = si0[j];
    psi_time[j] = psi_time0 + j*T;
    psi_state[j] = psi_state0 + j*T;
    //		for(t=0;t<10;t++) Rprintf("%g\t",p[j][t]);
    //		 Rprintf("\n");
  }

  for(t=0;t<T;t++) {
    //			Rprintf("t = %d\n",t);
    for(j=0;j<J;j++) {								
      //				Rprintf("alpha[%d][%d]",j,t);
      obs=0;
      if(t<T-1) {					
	for(u=1;u<=min(t+1,M[j]);u++) {
	  if(u<t+1) {
	    //							tmp1 = obs+d[j][u-1]+si[psi_state[j][t-u+1]][t-u+1];
	    tmp1 = obs+d[j][u-1]+si[j][t-u+1];
	    //							Rprintf("S%d %d %g %g %g %g\n",j,u,obs,d[j][u-1],si[j][t-u+1],tmp1);
	    if(u==1 || tmp_max < tmp1) {
	      tmp_max = tmp1;
	      psi_time[j][t]=u;
	    }
	    obs += p[j][t-u];
	    //							if(obs > 0) 								Rprintf("WOOT t = %d   u = %d   obs = %.3g\n",t,u,obs);

	  }
	  else {
	    tmp1 = obs+d[j][t]+start[j];
	    //							Rprintf("S%d %d %g %g %g %g\n",j,u,obs,d[j][t],start[j],tmp1);
	    if(u==1 || tmp_max < tmp1) {
	      tmp_max = tmp1;
	      psi_time[j][t]=u;
	    }
	  }
	}
	alpha[j][t] = tmp_max+p[j][t];					
	//					Rprintf("psi_time[%d][%d]  = %d\n",j,t,psi_time[j][t]);
      }
      else{
	for(u=1;u<=min(t+1,M[j]);u++) {
	  if(u<T) {
	    //							Rprintf(" = max(%.3g,%.3g)",j,t,tmp_max,tmp1);
	    //							tmp1 = obs+D[j][u-1]+si[psi_state[j][t-u+1]][t-u+1];
	    tmp1 = obs+D[j][u-1]+si[j][t-u+1];
	    if(u < 2000) 						
	      //								Rprintf("S%d %d %g %g %g %g\n",j,u,obs,D[j][u-1],si[j][t-u+1],tmp1);
	      if(u==1 || tmp_max < tmp1) {
		tmp_max = tmp1;
		psi_time[j][t]=u;
	      }
	    obs += p[j][T-1-u];
	  }
	  else {
	    tmp1 = obs+D[j][T-1]+start[j];	
	    //				Rprintf(" = max(%.3g,%.3g +%.3g + %.3g = %.3g)",tmp_max,obs,D[j][t],start[j],tmp1);
	    if(u==1 || tmp_max < tmp1) {
	      tmp_max = tmp1;
	      psi_time[j][t]=u;
	    }
	  }						
	}
	alpha[j][t] = tmp_max+p[j][t];
	//					Rprintf(" = %.3g\npsi_time[%d][%d]  = %d\n",alpha[j][t]-p[j][t],j,t,psi_time[j][t]);
      }
    }

    if(t<T-1) {
      for(j=0;j<J;j++){
	i=0;
	si[j][t+1]=a[j*J+i]+alpha[i][t];
	psi_state[j][t+1]=0;
	for(i=1;i<J;i++) 
	  if(i!=j) {
	    tmp1 = a[j*J+i]+alpha[i][t];
	    if(si[j][t+1] <= tmp1) {
	      si[j][t+1] = tmp1;
	      psi_state[j][t+1]=i;
	    }
	  }
	//					Rprintf("psi_state[%d][%d] = %d\tsi[%d][%d] = %.3g\n",j,t+1,psi_state[j][t+1],j,t+1,si[j][t+1]);
      }
    }
  }
	

  //and now we backtrack!
  q=statehat;

  q[T-1] = 0;
  for(j=1;j<J;j++) if(alpha[q[T-1]][T-1] < alpha[j][T-1]) q[T-1] = j;
  u=1;
  for(t=T-2;t>=0;t--) {
    if(u < psi_time[q[t+u]][t+u]) {
      q[t] = q[t+u];
      u++;
    }
    else {
      q[t] = psi_state[q[t+u]][t+u];
      u=1;
    }
  }

  free(si);
  free_matrix(J,T,(void **)si0);
  free(p);
  free(alpha);
  free(psi_time);
  free(psi_state);
}



void viterbi_online(double *a,double *start,double *p0,double *d0,double *D0,int *timelength,int *nstates,int *M,
		    double *alpha0) {
  double obs;
  int u,i,j,t,T;
  int J = *nstates;	
  double **p,**d,**D,**alpha,**si,**si0,tmp_max=-10000,tmp1;
  int **psi_time;
  int **psi_state;
  int *psi_state0, *psi_time0;
  T = *timelength;

  //	Rprintf("J = %d\n",J);	
  //	for(j=0;j<J;j++)  Rprintf("M[%d] = %d\t",j,M[j]);
  //	Rprintf("\n");	
	

  psi_time0 = (int *)malloc(T*J*sizeof(int));
  psi_state0 = (int *)malloc(T*J*sizeof(int));
  si0 = (double **)alloc_matrix(J,T,sizeof(double));
  psi_time = (int **)malloc(J*sizeof(int));
  psi_state = (int **)malloc(J*sizeof(int));
  p = (double **)malloc(sizeof(double *)*J);
  d = (double **)malloc(sizeof(double *)*J);
  D = (double **)malloc(sizeof(double *)*J);
  alpha = (double **)malloc(sizeof(double *)*J);
  si = (double **)malloc(sizeof(double *)*J);

  //	print_matrix(J,T,p0);
  //	print_matrix(J,M[0],d);

  for(j=0;j<J;j++) {
    d[j] = d0 + j*M[j];
    D[j] = D0 + j*M[j];		
    p[j]  = p0 + j*T;
    alpha[j]  = alpha0 + j*T;
    si[j] = si0[j];
    psi_time[j] = psi_time0 + j*T;
    psi_state[j] = psi_state0 + j*T;
  }

  for(t=0;t<T;t++) {
    for(j=0;j<J;j++) {				
      obs=0;
      if(t<T-1) {
	for(u=1;u<=min(t+1,M[j]);u++) {
	  if(u<t+1) {
	    //							tmp1 = obs+d[j][u-1]+si[psi_state[j][t-u+1]][t-u+1];
	    tmp1 = obs+d[j][u-1]+si[j][t-u+1];
	    if(u==1 || tmp_max < tmp1) {
	      tmp_max = tmp1;
	      psi_time[j][t]=u;
	    }
	    obs += p[j][t-u];
	  }
	  else {
	    tmp1 = obs+d[j][t]+start[j];
	    if(u==1 || tmp_max < tmp1) {
	      tmp_max = tmp1;
	      psi_time[j][t]=u;
	    }
	  }
	}
	alpha[j][t] = tmp_max+p[j][t];
	tmp_max=0;
      }
      obs=0;
      //this is calculating the online predictions, as well as alpha [T-1]
      for(u=1;u<=min(t+1,M[j]);u++) {
	if(u<T) {
	  tmp1 = obs+D[j][u-1]+si[j][t-u+1];
	  if(u==1 || tmp_max < tmp1) {
	    tmp_max = tmp1;
	    if(t==T-1) psi_time[j][t]=u;
	  }
	  obs += p[j][T-1-u];
	}
	else {
	  tmp1 = obs+D[j][T-1]+start[j];	
	  if(u==1 || tmp_max < tmp1) {
	    tmp_max = tmp1;
	    if(t==T-1) psi_time[j][t]=u;
	  }
	}						
      }
      if(t==T-1) alpha[j][t] = tmp_max+p[j][t];
    }

    if(t<T-1) {
      for(j=0;j<J;j++){
	i=0;
	si[j][t+1]=a[j*J+i]+alpha[i][t];
	psi_state[j][t+1]=0;
	for(i=1;i<J;i++) 
	  if(i!=j) {
	    tmp1 = a[j*J+i]+alpha[i][t];
	    if(si[j][t+1] <= tmp1) {
	      si[j][t+1] = tmp1;
	      psi_state[j][t+1]=i;
	    }
	  }
	//					Rprintf("psi_state[%d][%d] = %d\tsi[%d][%d] = %.3g\n",j,t+1,psi_state[j][t+1],j,t+1,si[j][t+1]);
      }
    }
  }
	
  free(si);
  free_matrix(J,T,(void **)si0);
  free(p);
  free(alpha);
  free(psi_time);
  free(psi_state);
  free(psi_state0);
  free(psi_time0);
}

void forward_online(double *a,double *start,double *p0,double *d,double *D,int *timelength,int *nstates,int *M,
		    double **F0,double *N0,double **si0,int *nsequences,int *totallength) {
  double obs;
  int u,i,j,t,n;
  int T;
  int J = *nstates;	
  double **p,**F,**si,*N;
  int offset;
  int ln = *totallength;
  int nseq = *nsequences;

  p = (double **)malloc(sizeof(double *)*J);
  F = (double **)malloc(sizeof(double *)*J);
  si = (double **)malloc(sizeof(double *)*J);

  for(j=0;j<J;j++) {
    p[j]  = p0 + j*ln;
    F[j]  = F0[j];
    si[j] = si0[j];
  }
  N = N0;

  for(n=0;n<nseq;n++) {	
    T = timelength[n];
    if(n>0){
      offset = timelength[n-1];
      for(j=0;j<J;j++) {
	p[j]  += offset;
	F[j]  += offset;
	si[j] += offset;
      }
      N += offset;
    }
    for(t=0;t<T;t++) {
      //		printf("\nt = %d\n",t);
      N[t]=0;
      for(j=0;j<J;j++) {
	F[j][t]=0;
	obs=p[j][t];
	if(t<T-1) {
	  //				printf("j = %d   ",j);
	  for(u=1;u<=min(t+1,M[j]);u++) {
	    if(u<t+1) {
	      F[j][t] += obs*d[j*M[j]+u-1]*si[j][t-u+1];
	      N[t] += obs*D[j*M[j]+u-1]*si[j][t-u+1];
	      obs*= p[j][t-u]/N[t-u];
	      //						printf("%.3g (%.3g) \t",obs,N[t-u]);
	    }
	    else {
	      F[j][t] += obs*d[j*M[j]+t]*start[j];
	      N[t] += obs*D[j*M[j]+t]*start[j];						
	    }
	  }
	  //				printf("\n");
	}
	else{
	  for(u=1;u<=min(t+1,M[j]);u++) {
	    if(u<T) {
	      F[j][T-1] += obs*D[j*M[j]+u-1]*si[j][T-u];
	      obs *= p[j][T-1-u]/N[T-1-u];
	    }
	    else F[j][T-1] += obs*D[j*M[j]+T-1]*start[j];
	  }
	  N[T-1] += F[j][T-1];
	}
      }

      for(j=0;j<J;j++){
	F[j][t] /= N[t];
	F[j][t]+=1e-300;
      }

      if(t<T-1) {
	for(j=0;j<J;j++){
	  si[j][t+1]=0;
	  for(i=0;i<J;i++) 
	    si[j][t+1]+=a[j*J+i]*F[i][t];
	}
      }
    }
  }
  free(si);
  free(p);
  free(F);
  /*
    printf("\nF=\n");
    print_matrix2(J,T,F);
    print_matrix(1,T,N);
  */
}



void forward_hmm(double *a,double *start,double *p,int offset,int *timelength,int *nstates,double ***output) {
  int K = *nstates;
  int T = *timelength;
  int i,k,t;

  double **alpha = *output;
	
  t=0;
  for(i=0;i<K;i++) alpha[i][t]=start[i]*p[offset*K+i];
	
  alpha[K][t] = 1; //no scaling at initialisation
	
  for(t=1;t<T;t++) {
    for(k=0;k<K;k++) {
      alpha[k][t]=0;
      for(i=0;i<K;i++) alpha[k][t]+=alpha[i][t-1]*a[i*K+k];
      alpha[k][t]*=p[offset*K+t*K+k];
    }
    alpha[K][t]=0;
    for(k=0;k<K;k++) alpha[K][t]+=alpha[k][t];
    alpha[K][t] = 1/alpha[K][t];
    for(k=0;k<K;k++) alpha[k][t]*=alpha[K][t];
  }
}

void backward_hmm(double *a,double *start,double *p,int offset,int *timelength,int *nstates,double *c,double ***output) {
  int K = *nstates;
  int T = *timelength;
  int i,j,t;

  double **beta = *output;

  for(i=0;i<K;i++) beta[i][T-1] = 1.0*c[T-1];

  for(t=T-2;t>=0;t--) {
    for(i=0;i<K;i++){
      beta[i][t]=0;
      for(j=0;j<K;j++) beta[i][t]+= a[i*K+j]*p[offset*K+(t+1)*K+j]*beta[j][t+1]*c[t];
    }
  }
}


//multiple observation estep
void mo_estep_hmm(double *a,double *start,double *p,int *T,int *nsequences,int *nstates,double *forward_p,double *backward_p,double *gam,double *loglik) {
  int K = *nstates;
  int N = *nsequences;
  double ***alpha,***beta,*xi,den,num,ll=0;
  int i,j,t,n;
  int *ii = (int *)malloc((N+1)*sizeof(int));
  checkmem(ii);
	
  ii[0]=0;
  for(n=1;n<N+1;n++) ii[n]= T[n-1]+ii[n-1];

  xi = (double *)malloc(K*K*ii[N]*sizeof(double));
  checkmem(xi);

  alpha = (double ***)malloc(N*sizeof(double **));
  //	checkmem(alpha);
  beta = (double ***)malloc(N*sizeof(double **));
  //	checkmem(beta);

  for(n=0;n<N;n++){
    //		Rprintf("Sequence %d\n",n);
    //		xi[n] = (double ***)malloc(T[n]*sizeof(double **));
    //		checkmem(xi[n]);
    /*
      alpha[n] = (double **)alloc_matrix(K+1,T[n],sizeof(double));
      beta[n] = (double **)alloc_matrix(K,T[n],sizeof(double));
    */
    alpha[n] = (double **)malloc((K+1)*sizeof(double *));
    beta[n] =  (double **)malloc((K)*sizeof(double *));
    for(i=0;i<K;i++) {
      alpha[n][i]=forward_p+i*ii[N]+ii[n];
      beta[n][i]=backward_p+i*ii[N]+ii[n];
    }
    i=K;
    alpha[n][i]=forward_p+i*ii[N]+ii[n];

    forward_hmm(a,start,p,ii[n],T+n,nstates,&alpha[n]);
    backward_hmm(a,start,p,ii[n],T+n,nstates,alpha[n][K],&beta[n]);
    //		print_matrix2(K,T[n],alpha[n]);
    for(t=0;t<T[n]-1;t++) {
      //			Rprintf("t = %d\n",t);
      //			xi[n][t] = (double **)alloc_matrix(K,K,sizeof(double));
      den = 0;
      for(i=0;i<K;i++)
	for(j=0;j<K;j++) {
	  //					xi[n][t][i][j]=alpha[n][i][t]*a[i*K+j]*p[ii[n]*K+(t+1)*K+j]*beta[n][j][t+1];
	  xi[ii[n]*K*K + t*K*K + i*K + j]=alpha[n][i][t]*a[i*K+j]*p[ii[n]*K+(t+1)*K+j]*beta[n][j][t+1];
	  //					den+=xi[n][t][i][j];
	  den+=xi[ii[n]*K*K + t*K*K + i*K + j];
	}

      for(i=0;i<K;i++)
	for(j=0;j<K;j++)
	  //					xi[n][t][i][j]/=den;
	  xi[ii[n]*K*K + t*K*K + i*K + j]/=den;
      //			print_matrix2(K,K,xi[n][t]);
      for(i=0;i<K;i++) {
	gam[i*ii[N]+ii[n]+t]=0;
	for(j=0;j<K;j++)
	  gam[i*ii[N]+ii[n]+t]+=xi[ii[n]*K*K + t*K*K + i*K + j];
      }
    }

    den = 0;
    for(i=0;i<K;i++) {
      gam[i*ii[N]+ii[n]+T[n]-1]=alpha[n][i][t]*beta[n][i][t];
      den+=gam[i*ii[N]+ii[n]+T[n]-1];
    }
    for(i=0;i<K;i++) gam[i*ii[N]+ii[n]+T[n]-1]/=den;
  }

  for(i=0;i<K;i++) start[i] = 0;
  for(n=0;n<N;n++)
    for(i=0;i<K;i++) 
      start[i]+=gam[i*ii[N]+ii[n]+0]/N;
		
  for(i=0;i<K;i++) {
    den = 0;		
    for(n=0;n<N;n++) for(t=0;t<T[n]-1;t++) 	den += gam[i*ii[N]+ii[n]+t];

    for(j=0;j<K;j++) {
      num = 0;
      for(n=0;n<N;n++) for(t=0;t<T[n]-1;t++) num+=xi[ii[n]*K*K + t*K*K + i*K + j];
      a[i*K+j]=num/den;
    }		
  }

  //Calculate the loglikelihood
  for(n=0;n<N;n++) for(t=0;t<T[n];t++) {
      //		Rprintf("alpha[%d][K][%d] = %.3g\n",n,t,alpha[n][K][t]);
      if(alpha[n][K][t]<0) {
	error("Negative likelihood\n alpha[%d][K][%d] = %.3g\n",n,t,alpha[n][K][t]);
      }
      else ll+=log(alpha[n][K][t]);
    }
  ll*=-1;
	
  *loglik=ll;
  //	Rprintf("%.5g\n",ll);
  /*
  for(n=0;n<N;n++) free_matrix(K,T[n],(void **)beta[n]);	
  for(n=0;n<N;n++) free_matrix(K+1,T[n],(void **)alpha[n]);
  */
  for(n=0;n<N;n++){
    free(beta[n]);
    free(alpha[n]);
  }
  //	for(n=0;n<N;n++) for(t=0;t<T[n]-1;t++) free_matrix(K,K,(void **)xi[n][t]);
  free(alpha);
  free(beta);
  free(xi);
  free(ii);
}


void viterbi_hmm(double *a,double *start,double *p,int *T,int *nsequences,int *nstates,int *q,double *loglik) {
  int K = *nstates;
  int N = *nsequences;
  int i,j,t,n;
  double P;
  int **psi;
  double **delta;
  double *tmp1;
  int maxind;
  int *ii = (int *)malloc((N+1)*sizeof(int));

  //	Rprintf("%d %d %d\n",K,N,*T);
  checkmem(ii);

  ii[0]=0;
  for(n=1;n<N+1;n++) ii[n]= T[n-1]+ii[n-1];

  if(K<2) error("Invalid number of states (K = %d)\n",K);

  psi = (int **)alloc_matrix(K,ii[N],sizeof(int));
  delta = (double **)alloc_matrix(K,ii[N],sizeof(double));
  tmp1 = (double *)malloc(K*sizeof(double));
  checkmem(tmp1);
  checkmem(psi);
  checkmem(delta);
  /*
    Rprintf("%g %g\n",start[0],start[1]);
    print_matrix(1,K,start);
    print_matrix(K,K,a);
    print_matrix(K,T,p);	
    print_matrix2(K,T,psi);	
    print_matrix2(K,T,delta);	
  */
  for(n=0;n<N;n++) {
    t=ii[n];
    for(i=0;i<K;i++) {
      delta[i][t]=start[i] + p[i];
      psi[i][t]=0;
    }
    for(t=ii[n]+1;t<ii[n+1];t++) {
      for(j=0;j<K;j++) {
	i=0;
	maxind = i;
	tmp1[i] = delta[i][t-1] + a[i*K+j];
	for(i=1;i<K;i++) {			
	  tmp1[i] = delta[i][t-1] + a[i*K+j];
	  if(tmp1[i]>tmp1[maxind]) maxind = i;
	}			
	delta[j][t]=tmp1[maxind] + p[t*K+j];
	psi[j][t]=maxind;
      }
    }				
  }

  P=1;
  *loglik=0.0;
  for(n=1;n<=N;n++) {
    maxind = 0;
    for(i=1;i<K;i++) 
      if(delta[i][ii[n]-1]>delta[maxind][ii[n]-1]) maxind=i;	
    *loglik+=delta[maxind][ii[n]-1];
    q[ii[n]-1]=maxind;
    P += delta[maxind][ii[n]-1];
  }

  for(n=0;n<N;n++)
    for(t=ii[n+1]-2;t>=ii[n];t--) {
      if(q[t+1]<0) {
	error("Invalid state at n = %d and t = %d\n",n,t+1);
	free_matrix(K,ii[N],(void **)psi);	
	free_matrix(K,ii[N],(void **)delta);	
      }	
      else q[t]=psi[q[t+1]][t+1];
    }

  //	print_imatrix2(K,T,psi);	
  //	print_matrix2(K,3050,delta);	
  free_matrix(K,ii[N],(void **)psi);	
  free_matrix(K,ii[N],(void **)delta);	
  free(tmp1);
  free(ii);
}
