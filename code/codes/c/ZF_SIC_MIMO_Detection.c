#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#define MAX 4
#define MAXM 16
#define DETECTOR "ZF-SIC"
#define ML_LIMIT 100000

double grandom(void) {
    double u1=((double)rand()+1.0)/((double)RAND_MAX+2.0);
    double u2=((double)rand()+1.0)/((double)RAND_MAX+2.0);
    return sqrt(-2.0*log(u1))*cos(2.0*acos(-1.0)*u2);
}
void make_constellation(int mode,double complex S[MAXM],int *M) {
    if(mode==1) {
        *M=4; S[0]=(1+I)/sqrt(2.0); S[1]=(-1+I)/sqrt(2.0);
        S[2]=(1-I)/sqrt(2.0); S[3]=(-1-I)/sqrt(2.0);
    } else {
        int lev[4]={-3,-1,1,3},k=0; *M=16;
        for(int q=0;q<4;q++) for(int i=0;i<4;i++) S[k++]=(lev[i]+I*lev[q])/sqrt(10.0);
    }
}
void make_channel(int nr,int nt,int mode,double complex H[MAX][MAX]) {
    for(int r=0;r<nr;r++) for(int c=0;c<nt;c++) H[r][c]=0;
    if(mode==2) { for(int i=0;i<nr&&i<nt;i++) H[i][i]=1; }
    else if(mode==1) {
        for(int r=0;r<nr;r++) for(int c=0;c<nt;c++)
            H[r][c]=(grandom()+I*grandom())/sqrt(2.0);
    } else {
        double K=5.0,los=sqrt(K/(K+1.0))/sqrt(2.0),sc=sqrt(1.0/(K+1.0))/sqrt(2.0);
        for(int r=0;r<nr;r++) for(int c=0;c<nt;c++) H[r][c]=los+sc*(grandom()+I*grandom());
    }
}
int solve_linear(int n,double complex A[MAX][MAX],double complex b[MAX],double complex x[MAX]) {
    double complex a[MAX][MAX+1];
    for(int r=0;r<n;r++){for(int c=0;c<n;c++)a[r][c]=A[r][c];a[r][n]=b[r];}
    for(int c=0;c<n;c++){
        int p=c;for(int r=c+1;r<n;r++)if(cabs(a[r][c])>cabs(a[p][c]))p=r;
        if(cabs(a[p][c])<1e-13)return 0;
        for(int j=c;j<=n;j++){double complex t=a[c][j];a[c][j]=a[p][j];a[p][j]=t;}
        double complex d=a[c][c];for(int j=c;j<=n;j++)a[c][j]/=d;
        for(int r=0;r<n;r++)if(r!=c){double complex f=a[r][c];for(int j=c;j<=n;j++)a[r][j]-=f*a[c][j];}
    }
    for(int i=0;i<n;i++)x[i]=a[i][n];return 1;
}
int weights(int nr,int nt,double complex H[MAX][MAX],double nv,double Es,int mmse,double complex W[MAX][MAX]) {
    double complex G[MAX][MAX]={{0}},rhs[MAX],sol[MAX];
    for(int i=0;i<nt;i++)for(int j=0;j<nt;j++)for(int r=0;r<nr;r++)G[i][j]+=conj(H[r][i])*H[r][j];
    if(mmse)for(int i=0;i<nt;i++)G[i][i]+=nv/Es;
    for(int c=0;c<nr;c++){
        for(int i=0;i<nt;i++)rhs[i]=conj(H[c][i]);
        if(!solve_linear(nt,G,rhs,sol))return 0;
        for(int i=0;i<nt;i++)W[i][c]=sol[i];
    }
    return 1;
}
int nearest(double complex z,double complex S[MAXM],int M) {
    int best=0;double d=pow(cabs(z-S[0]),2);
    for(int i=1;i<M;i++){double v=pow(cabs(z-S[i]),2);if(v<d){d=v;best=i;}}
    return best;
}
void linear_detect(int nr,int nt,double complex y[MAX],double complex H[MAX][MAX],
                   double complex S[MAXM],int M,double nv,double Es,int mmse,
                   double complex xhat[MAX]) {
    double complex W[MAX][MAX]={{0}};
    if(!weights(nr,nt,H,nv,Es,mmse,W)){for(int i=0;i<nt;i++)xhat[i]=S[0];return;}
    for(int i=0;i<nt;i++){double complex z=0;for(int r=0;r<nr;r++)z+=W[i][r]*y[r];xhat[i]=S[nearest(z,S,M)];}
}
int sic_detect(int nr,int nt,double complex y[MAX],double complex H[MAX][MAX],
               double complex S[MAXM],int M,double nv,double Es,int mmse,
               double complex xhat[MAX],int order[MAX]) {
    double complex residual[MAX],Hr[MAX][MAX],W[MAX][MAX],z[MAX];
    int rem[MAX],count=nt;
    for(int r=0;r<nr;r++)residual[r]=y[r];
    for(int i=0;i<nt;i++)rem[i]=i;
    for(int stage=0;stage<nt;stage++){
        for(int r=0;r<nr;r++)for(int c=0;c<count;c++)Hr[r][c]=H[r][rem[c]];
        if(!weights(nr,count,Hr,nv,Es,mmse,W))return 0;
        for(int i=0;i<count;i++){z[i]=0;for(int r=0;r<nr;r++)z[i]+=W[i][r]*residual[r];}
        int pick=0;double best=mmse?1e300:-1.0;
        for(int i=0;i<count;i++){
            double metric=0;
            if(mmse){for(int r=0;r<nr;r++)metric+=pow(cabs(W[i][r]),2.0);if(metric<best){best=metric;pick=i;}}
            else {for(int r=0;r<nr;r++)metric+=pow(cabs(Hr[r][i]),2.0);if(metric>best){best=metric;pick=i;}}
        }
        int stream=rem[pick];xhat[stream]=S[nearest(z[pick],S,M)];
        for(int r=0;r<nr;r++)residual[r]-=H[r][stream]*xhat[stream];
        order[stage]=stream+1;
        for(int i=pick;i<count-1;i++)rem[i]=rem[i+1];count--;
    }
    return 1;
}
int ml_detect(int nr,int nt,double complex y[MAX],double complex H[MAX][MAX],
              double complex S[MAXM],int M,double complex xhat[MAX]) {
    long total=1;for(int i=0;i<nt;i++){if(total>ML_LIMIT/M)return 0;total*=M;}
    double best=1e300;int idx[MAX],bestidx[MAX];
    for(long q=0;q<total;q++){
        long v=q;for(int i=0;i<nt;i++){idx[i]=(int)(v%M);v/=M;}
        double metric=0;
        for(int r=0;r<nr;r++){double complex e=y[r];for(int c=0;c<nt;c++)e-=H[r][c]*S[idx[c]];metric+=pow(cabs(e),2.0);}
        if(metric<best){best=metric;for(int i=0;i<nt;i++)bestidx[i]=idx[i];}
    }
    for(int i=0;i<nt;i++)xhat[i]=S[bestidx[i]];return 1;
}

/* Complex Gram-Schmidt and LLL */
void gs(int m,int n,double complex B[MAX][MAX],double complex star[MAX][MAX],
        double complex mu[MAX][MAX],double norms[MAX]) {
    for(int r=0;r<m;r++)for(int c=0;c<n;c++)star[r][c]=0;
    for(int i=0;i<n;i++){
        norms[i]=0;for(int j=0;j<n;j++)mu[i][j]=0;
        for(int r=0;r<m;r++)star[r][i]=B[r][i];
        for(int j=0;j<i;j++)if(norms[j]>1e-14){
            double complex v=0;for(int r=0;r<m;r++)v+=conj(B[r][i])*star[r][j];
            mu[i][j]=v/norms[j];for(int r=0;r<m;r++)star[r][i]-=mu[i][j]*star[r][j];
        }
        for(int r=0;r<m;r++)norms[i]+=pow(cabs(star[r][i]),2.0);
    }
}
double complex ground(double complex z){return round(creal(z))+I*round(cimag(z));}
int lll(int m,int n,double complex H[MAX][MAX],double delta,
        double complex B[MAX][MAX],double complex T[MAX][MAX]) {
    double complex star[MAX][MAX],mu[MAX][MAX];double norms[MAX];
    for(int r=0;r<m;r++)for(int c=0;c<n;c++)B[r][c]=H[r][c];
    for(int r=0;r<n;r++)for(int c=0;c<n;c++)T[r][c]=(r==c)?1:0;
    if(n<2)return 0;
    int k=1,it=0;
    while(k<n&&it<5000){
        it++;gs(m,n,B,star,mu,norms);
        for(int j=k-1;j>=0;j--){
            double complex q=ground(mu[k][j]);
            if(cabs(q)>0){
                for(int r=0;r<m;r++)B[r][k]-=q*B[r][j];
                for(int r=0;r<n;r++)T[r][k]-=q*T[r][j];
                gs(m,n,B,star,mu,norms);
            }
        }
        if(norms[k-1]<1e-14)return -1;
        if(norms[k]>=(delta-pow(cabs(mu[k][k-1]),2.0))*norms[k-1])k++;
        else{
            for(int r=0;r<m;r++){double complex t=B[r][k];B[r][k]=B[r][k-1];B[r][k-1]=t;}
            for(int r=0;r<n;r++){double complex t=T[r][k];T[r][k]=T[r][k-1];T[r][k-1]=t;}
            k=(k>1)?k-1:1;
        }
    }
    return it;
}
int lr_detect(int nr,int nt,double complex H[MAX][MAX],double complex y[MAX],
              double complex S[MAXM],int M,double nv,double Es,double delta,
              double complex xhat[MAX],double complex Hred[MAX][MAX],
              double complex T[MAX][MAX],int *iterations) {
    *iterations=lll(nr,nt,H,delta,Hred,T);if(*iterations<0)return 0;
    double complex W[MAX][MAX]={{0}},z[MAX],xraw[MAX];
    if(!weights(nr,nt,Hred,nv,Es,1,W))return 0;
    for(int i=0;i<nt;i++){z[i]=0;for(int r=0;r<nr;r++)z[i]+=W[i][r]*y[r];}
    for(int i=0;i<nt;i++){xraw[i]=0;for(int j=0;j<nt;j++)xraw[i]+=T[i][j]*z[j];xhat[i]=S[nearest(xraw[i],S,M)];}
    return 1;
}
int bit_errors(int tx,int rx,int bps){int v=tx^rx,e=0;for(int i=0;i<bps;i++)e+=(v>>i)&1;return e;}
void pv(const char *name,double complex x[MAX],int n){printf("%s: ",name);for(int i=0;i<n;i++)printf("%+.4f%+.4fj ",creal(x[i]),cimag(x[i]));printf("\n");}
void pm(const char *name,double complex A[MAX][MAX],int m,int n){printf("%s:\n",name);for(int r=0;r<m;r++){for(int c=0;c<n;c++)printf("%+.3f%+.3fj ",creal(A[r][c]),cimag(A[r][c]));printf("\n");}}

int main(void) {
    int nt,nr,mc,cc,frames,ns,bps,M;double snr,delta=0.75,Es=0;
    double complex S[MAXM],H[MAX][MAX],Hr[MAX][MAX],T[MAX][MAX],x[MAX],xh[MAX],y[MAX];
    double complex lastx[MAX],lastxh[MAX],lasty[MAX],lastH[MAX][MAX],lastHr[MAX][MAX],lastT[MAX][MAX];
    double snrs[8],ber[8];int nS=0,lastOrder[MAX]={0},lastIt=0;
    srand(41);
    printf("============================================\n%s MIMO DETECTION SIMULATION\n============================================\n",DETECTOR);
    printf("Transmit antennas Nt (1-4): ");if(scanf("%d",&nt)!=1)return 1;
    printf("Receive antennas Nr (1-4): ");if(scanf("%d",&nr)!=1)return 1;
    if(nt<1||nt>MAX||nr<1||nr>MAX){printf("Nt and Nr must be 1..4.\n");return 1;}
    if((strcmp(DETECTOR,"ZF")==0||strcmp(DETECTOR,"ZF-SIC")==0||strcmp(DETECTOR,"LR")==0)&&nr<nt){printf("%s requires Nr >= Nt in this implementation.\n",DETECTOR);return 1;}
    printf("Modulation: 1=QPSK, 2=16-QAM. Choose: ");if(scanf("%d",&mc)!=1||(mc!=1&&mc!=2))return 1;
    make_constellation(mc,S,&M);bps=(mc==1)?2:4;for(int i=0;i<M;i++)Es+=pow(cabs(S[i]),2.0);Es/=M;
    printf("Channel: 1=Rayleigh, 2=AWGN, 3=Rician. Choose: ");if(scanf("%d",&cc)!=1||cc<1||cc>3)return 1;
    printf("Selected SNR dB (-5..40): ");if(scanf("%lf",&snr)!=1)return 1;
    printf("Frames (1..2000): ");if(scanf("%d",&frames)!=1)return 1;
    printf("Symbols/frame (10..1000): ");if(scanf("%d",&ns)!=1)return 1;
    if(strcmp(DETECTOR,"LR")==0){printf("LLL delta (0.51..0.99): ");if(scanf("%lf",&delta)!=1)return 1;}
    if(snr< -5||snr>40||frames<1||frames>2000||ns<10||ns>1000||delta<0.51||delta>0.99){printf("Invalid input range.\n");return 1;}
    for(int s=0;s<=30;s+=5)snrs[nS++]=s;if(snr<0||snr>30)snrs[nS++]=snr;
    printf("\nRunning simulation...\n");
    for(int si=0;si<nS;si++){
        long errors=0,total=0;double nv=pow(10.0,-snrs[si]/10.0);
        for(int f=0;f<frames;f++)for(int t=0;t<ns;t++){
            int txi[MAX],order[MAX]={0},it=0;
            for(int i=0;i<nt;i++){txi[i]=rand()%M;x[i]=S[txi[i]];}
            make_channel(nr,nt,cc,H);
            for(int r=0;r<nr;r++){y[r]=0;for(int c=0;c<nt;c++)y[r]+=H[r][c]*x[c];y[r]+=sqrt(nv/2.0)*(grandom()+I*grandom());}
            int ok=1;
            if(strcmp(DETECTOR,"ML")==0)ok=ml_detect(nr,nt,y,H,S,M,xh);
            else if(strcmp(DETECTOR,"ZF")==0)linear_detect(nr,nt,y,H,S,M,nv,Es,0,xh);
            else if(strcmp(DETECTOR,"MMSE")==0)linear_detect(nr,nt,y,H,S,M,nv,Es,1,xh);
            else if(strcmp(DETECTOR,"ZF-SIC")==0)ok=sic_detect(nr,nt,y,H,S,M,nv,Es,0,xh,order);
            else if(strcmp(DETECTOR,"MMSE-SIC")==0)ok=sic_detect(nr,nt,y,H,S,M,nv,Es,1,xh,order);
            else ok=lr_detect(nr,nt,H,y,S,M,nv,Es,delta,xh,Hr,T,&it);
            if(!ok){printf("Detection failed (singular channel or ML limit).\n");return 1;}
            for(int i=0;i<nt;i++)errors+=bit_errors(txi[i],nearest(xh[i],S,M),bps),total+=bps;
            if(si==nS-1&&f==0&&t==0){
                for(int i=0;i<nt;i++){lastx[i]=x[i];lastxh[i]=xh[i];lastOrder[i]=order[i];}
                for(int r=0;r<nr;r++){lasty[r]=y[r];for(int c=0;c<nt;c++){lastH[r][c]=H[r][c];if(strcmp(DETECTOR,"LR")==0)lastHr[r][c]=Hr[r][c];}}
                if(strcmp(DETECTOR,"LR")==0){for(int r=0;r<nt;r++)for(int c=0;c<nt;c++)lastT[r][c]=T[r][c];lastIt=it;}
            }
        }
        ber[si]=(double)errors/(double)total;
    }
    printf("\nBER versus SNR\n%-12s %-14s\n","SNR(dB)","BER");
    for(int i=0;i<nS;i++)printf("%-12.1f %.8e\n",snrs[i],ber[i]);
    pv("\nTransmitted x",lastx,nt);pv("Received y",lasty,nr);pv("Detected xhat",lastxh,nt);
    if(strcmp(DETECTOR,"ZF-SIC")==0||strcmp(DETECTOR,"MMSE-SIC")==0){printf("Detection order: ");for(int i=0;i<nt;i++)printf("%d ",lastOrder[i]);printf("\n");}
    if(strcmp(DETECTOR,"LR")==0){pm("Original H",lastH,nr,nt);pm("Reduced H*T",lastHr,nr,nt);pm("Transformation T",lastT,nt,nt);printf("LLL iterations: %d\n",lastIt);}
    printf("\nDetector: %s | Modulation: %s | Channel: %s\n",DETECTOR,mc==1?"QPSK":"16-QAM",cc==1?"Rayleigh":cc==2?"AWGN":"Rician");
    printf("Antenna configuration: %d x %d | Selected SNR: %.1f dB\n",nt,nr,snr);
    printf("Simulation completed.\n");return 0;
}
