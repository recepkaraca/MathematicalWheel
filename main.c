#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


#ifdef NAN
/* NAN is supported */
#endif
#ifdef INFINITY
/* INFINITY is supported */
#endif


typedef struct
{
    int m, n;
    double ** v;
} mat_t, *mat;

void bekle();
int cark();
void matriscarp(double [100][100],double [100][100],double,double [100][100]);
mat matrix_new(int,int);
void matrix_delete(mat);
void matrix_transpose(mat);
mat matrix_mul(mat,mat);
mat matrix_minor(mat,int);
double *vmadd(double [],double [],double,double [],int);
mat vmul(double [],int);
double vnorm(double [],int);
double* vdiv(double [],double,double [],int);
double* mcol(mat , double *, int);
void matrix_show(mat,short);
void householder(mat,mat *,mat *);
void ozdegerbul(double [100][100],double,double *,int *,short);
void schur(double [100][100],int,double*,int);
void nilpotent(int [100][100],int);

FILE  *DosyaOku;

void bekle()
{
    clock_t hedefzaman=1000+clock();
    clock_t suankizaman=clock();
    while(hedefzaman>=suankizaman)
    {
        suankizaman=clock();
    };
}
int cark()
{
    int uretileceksayi,i,rastgele;
    int dur=0;

    int temp, durum;

	printf("Uretilecek sayi miktarini giriniz:");
	durum = scanf("%d",&uretileceksayi);
	while(durum!=1){
		while((temp=getchar()) != EOF && temp != '\n');
		printf("Gecersiz giris.Lutfen tekrar giris yapiniz: ");
		durum = scanf("%d",&uretileceksayi);
	}


    for(i=0; i<uretileceksayi; i++)
    {
        rastgele=rand()%241;
        printf("%d.tur : %d\n",i+1,rastgele);
        bekle();
    }
    return rastgele%4;
}

void matriscarp(double matris1[100][100],double matris2[100][100],double boyut,double geridonenmatris[100][100]){
    double carpimmatris[100][100];
    int i,j,k;
    double toplam=0,n=boyut;
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            for(k=0; k<n; k++){
                toplam += matris1[i][k] * matris2[k][j];
            }
            carpimmatris[i][j] = toplam;
            toplam = 0;
        }
    }
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            geridonenmatris[i][j]=carpimmatris[i][j];
        }
    }
}



mat matrix_new(int m, int n)
{
    int i;
    mat x = malloc(sizeof(mat_t));
    x->v = malloc(sizeof(double) * m);
    x->v[0] = calloc(sizeof(double), m * n);
    for (i = 0; i < m; i++)
        x->v[i] = x->v[0] + n * i;
    x->m = m;
    x->n = n;
    return x;
}

void matrix_delete(mat m)
{
    free(m->v[0]);
    free(m->v);
    free(m);
}

void matrix_transpose(mat m)
{
    int i,j;
    for (i = 0; i < m->m; i++)
    {
        for (j = 0; j < i; j++)
        {
            double t = m->v[i][j];
            m->v[i][j] = m->v[j][i];
            m->v[j][i] = t;
        }
    }
}

mat matrix_copy(int n, double a[][n], int m)
{
    mat x = matrix_new(m, n);
    int i,j;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            x->v[i][j] = a[i][j];
    return x;
}

mat matrix_mul(mat x, mat y)
{
    int i,j,k;
    if (x->n != y->m) return 0;
    mat r = matrix_new(x->m, y->n);
    for (i = 0; i < x->m; i++)
        for (j = 0; j < y->n; j++)
            for (k = 0; k < x->n; k++)
                r->v[i][j] += x->v[i][k] * y->v[k][j];
    return r;
}

mat matrix_minor(mat x, int d)
{
    mat m = matrix_new(x->m, x->n);
    int i,j;
    for (i = 0; i < d; i++)
        m->v[i][i] = 1;
    for (i = d; i < x->m; i++)
        for (j = d; j < x->n; j++)
            m->v[i][j] = x->v[i][j];
    return m;
}

/* c = a + b * s */
double *vmadd(double a[], double b[], double s, double c[], int n)
{
    int i;
    for (i = 0; i < n; i++)
        c[i] = a[i] + s * b[i];
    return c;
}

/* m = I - v v^T */
mat vmul(double v[], int n)
{
    int i,j;
    mat x = matrix_new(n, n);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            x->v[i][j] = -2 *  v[i] * v[j];
    for (i = 0; i < n; i++)
        x->v[i][i] += 1;

    return x;
}

/* ||x|| */
double vnorm(double x[], int n)
{
    int i;
    double sum = 0;
    for (i = 0; i < n; i++) sum += x[i] * x[i];
    return sqrt(sum);
}

/* y = x / d */
double* vdiv(double x[], double d, double y[], int n)
{
    int i;
    for (i = 0; i < n; i++) y[i] = x[i] / d;
    return y;
}

/* take c-th column of m, put in v */
double* mcol(mat m, double *v, int c)
{
    int i;
    for (i = 0; i < m->m; i++)
        v[i] = m->v[i][c];
    return v;
}

void matrix_show(mat m,short dosyayayaz)
{
    int i,j;
    for(i = 0; i < m->m; i++)
    {
        for (j = 0; j < m->n; j++)
        {
            printf(" %-8.3f", m->v[i][j]);
            if(dosyayayaz==1)fprintf(DosyaOku, "%-8.3f  ", m->v[i][j]);
        }
        if(dosyayayaz==1)fputs("\n",DosyaOku );
        printf("\n");
    }
    printf("\n");
}

void householder(mat m, mat *R, mat *Q)
{
    int k,i;
    mat q[m->m];
    mat z = m, z1;
    for (k = 0; k < m->n && k < m->m - 1; k++)
    {
        double e[m->m], x[m->m], a;
        z1 = matrix_minor(z, k);
        if (z != m) matrix_delete(z);
        z = z1;

        mcol(z, x, k);
        a = vnorm(x, m->m);
        if (m->v[k][k] > 0) a = -a;

        for (i = 0; i < m->m; i++)
            e[i] = (i == k) ? 1 : 0;

        vmadd(x, e, a, e, m->m);
        vdiv(e, vnorm(e, m->m), e, m->m);
        q[k] = vmul(e, m->m);
        z1 = matrix_mul(q[k], z);
        if (z != m) matrix_delete(z);
        z = z1;
    }
    matrix_delete(z);
    *Q = q[0];
    *R = matrix_mul(q[0], m);
    for (i = 1; i < m->n && i < m->m - 1; i++)
    {
        z1 = matrix_mul(q[i], *Q);
        if (i > 1) matrix_delete(*Q);
        *Q = z1;
        matrix_delete(q[i]);
    }
    matrix_delete(q[0]);
    z = matrix_mul(*Q, m);
    matrix_delete(*R);
    *R = z;
    matrix_transpose(*Q);
}

double in[100][100];


void ozdegerbul(double matris[][100],double boyut,double *ozdegermatris,int *ozdegerboyut,short dosyayayaz)
{
    double size=boyut,d;
    double ozdegerler[100];
    int i,j,k;

    if(size>=3)
    {

        for(i=0; i<size; i++)
        {
            for(j=0; j<size; j++)
            {
                in[i][j]=0;
            }
        }

        mat R, Q;

        for(i=0; i<size; i++)
        {
            for(j=0; j<size; j++)
            {
                in[i][j]=matris[i][j];
            }
        }

        mat x = matrix_copy(size, in, size);
        mat x_kopya =matrix_copy(size,in,size);
        double toplam=0,kosegenharictoplam=0;
        int dur=0,sayac=0,sayac2=0,sayac3=0,infsayac=0;

        for(i=0; i<size; i++)
        {
            for(j=0; j<size; j++)
            {
                x->v[i][j]=in[i][j];
                x_kopya->v[i][j]=in[i][j];
            }
        }

        while(!dur)
        {
            puts("X");
            if(dosyayayaz==1)fputs("X Matrisi = \n",DosyaOku);
            matrix_show(x,dosyayayaz);
            if(dosyayayaz==1)fputs("\n\n",DosyaOku);
            for(i=0; i<size; i++)
            {
                for(j=0; j<size; j++)
                {

                    if(i!=j)
                    {
                        if(x->v[i][j]==0.0)
                        {
                            x->v[i][j]=0.0001;
                        }
                        if(x->v[i][j]<=0.0005 && x->v[i][j]>=-0.0005)
                        {
                            sayac++;
                        }
                    }
                    x_kopya->v[i][j]=x->v[i][j];
                }
            }

            householder(x, &R, &Q);

            x=matrix_mul(R,Q);

            printf("\n");

            if(dosyayayaz==1)fputs("Q Matrisi = \n",DosyaOku);
            puts("Q");
            matrix_show(Q,dosyayayaz);
            if(dosyayayaz==1)fputs("\n\n",DosyaOku);
            if(dosyayayaz==1)fputs("R Matrisi = \n",DosyaOku);
            puts("R");
            matrix_show(R,dosyayayaz);
            if(dosyayayaz==1)fputs("\n\n",DosyaOku);
            if(dosyayayaz==1)fputs("R * Q = X Matrisi = \n",DosyaOku);
            puts("R * Q = X");
            matrix_show(x,dosyayayaz);
            if(dosyayayaz==1)fputs("\n\n",DosyaOku);
            printf("\n\n");

            for(i=0; i<size; i++)
            {
                for(j=0; j<size; j++)
                {
                    if(fabs(x->v[i][j]-x_kopya->v[i][j])<=0.0005)
                    {
                        sayac3++;
                    }
                    if(x->v[i][j]==INFINITY || x->v[i][j]==NAN || x->v[i][j]==(-1)*INFINITY || x->v[i][j]==(-1)*NAN || x->v[i][j]!=x->v[i][j]){infsayac++;}
                }
            }
            if(sayac==(size*size-size) || sayac3==(size*size) || infsayac==(size*size))
            {
                dur=1;
            }
            sayac=0;
            sayac2++;
            //printf("%d",sayac3);
            sayac3=0;
            printf("%d ",infsayac);
            infsayac=0;
            if(sayac2>=2000)
            {
                break;
            }
        }

        int koksayac=0,kokler_eleman=0;
        for(i=0; i<size; i++)
        {
            for(j=0; j<size; j++)
            {
                if(i==j)
                {
                    for(k=i; k<size; k++)
                    {
                        if(x->v[k][j]<=0.1 && x->v[k][j]>=-0.1)
                        {
                            koksayac++;
                        }
                    }
                }
            }
            if(i!=size-1)
            {
                if(koksayac==size-i-1 && (x->v[i][i-1]<=0.1 && x->v[i][i-1]>=-0.1))
                {
                    printf("%d. eleman koktur.",i+1);
                    ozdegerler[kokler_eleman]=x->v[i][i];
                    kokler_eleman++;
                }
            }
            else
            {
                if(x->v[i][i-1]<=0.0001 && x->v[i][i-1]>=-0.0001)
                {
                    printf("%d. eleman koktur.",i+1);
                    ozdegerler[kokler_eleman]=x->v[i][i];
                    kokler_eleman++;
                }
            }
            koksayac=0;
        }
        *ozdegerboyut=kokler_eleman;
        printf("%d adim surdu.",sayac2);

        for(i=0; i<kokler_eleman; i++)
        {
            ozdegermatris[i]=ozdegerler[i];
        }

        matrix_delete(x);
        matrix_delete(R);
        matrix_delete(Q);
    }else if(size==2){
        double a,b,c,delta,kok1,kok2;
        a=1;
        b=-1*(matris[0][0]+matris[1][1]);
        c=(matris[0][0]*matris[1][1])-(matris[0][1]*matris[1][0]);

        delta=(b*b)-(4*a*c);
        if(dosyayayaz==1)fprintf(DosyaOku, "Delta = %4.3f \n", delta);
        if(delta<0){
                if(dosyayayaz==1)fputs("Girilen matrisin koku delta<0 oldugundan yoktur.",DosyaOku);
                printf("Kok yok");
                *ozdegerboyut=0;
        }
        else if(delta>0){
            kok1=((-1*b)-sqrt(delta))/2*a;
            kok2=((-1*b)+sqrt(delta))/2*a;
            ozdegermatris[0]=kok1;
            ozdegermatris[1]=kok2;
            *ozdegerboyut=2;
        }
        else{
            kok1=((-1*b)-sqrt(delta))/2*a;
            if(dosyayayaz==1)fputs("Delta=0 oldugundan kokler cakisiktir.\n",DosyaOku);
            ozdegermatris[0]=kok1;
            *ozdegerboyut=1;
        }
    }
    else if(size==1){
        printf("%4.2f",matris[0][0]);
        ozdegermatris[0]=matris[0][0];
        *ozdegerboyut=1;
    }
}

void schur(double girilenmatris[100][100],int girilenmatrisboyut,double *ozdegerler,int ozdegerboyut){
    int i,j;
    double esitsizlik_sol=0;
    double esitsizlik_sag=0;

    for(i=0;i<ozdegerboyut;i++){
        esitsizlik_sol+=pow(fabs(ozdegerler[i]),2);
    }

    fputs("Girilen matris: \n",DosyaOku);
    for(i=0;i<girilenmatrisboyut;i++){
        for(j=0;j<girilenmatrisboyut;j++){
            fprintf(DosyaOku, "%-8.2f", girilenmatris[i][j]);
            esitsizlik_sag+=pow(fabs(girilenmatris[i][j]),2);
        }
        fputs("\n",DosyaOku);
    }
    printf("\n%4.2f < %4.2f",esitsizlik_sol,esitsizlik_sag);
    fputs("\n",DosyaOku);
    for(i=0;i<ozdegerboyut;i++){
            fprintf(DosyaOku,"%d. ozdeger: %f \n",i+1,ozdegerler[i]);
        }

    fprintf(DosyaOku, "Sinir Deger = %-8.2f < %-8.2f", esitsizlik_sol,esitsizlik_sag);
}

void nilpotent(int rastgelematris[100][100],int boyut){
    double matris[100][100],carpim[100][100],carpan=1.0;

    int ilkbasamaksayisi[100][100],sonbasamaksayisi[100][100];
    int i,j,n=boyut,dur=0,sayac=0,sayac2=0,hepsibuyuk=0;
    int basamaksayac=0;
    for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                matris[i][j]=rastgelematris[i][j];
            }
        }
    matriscarp(matris,matris,n,carpim);
    fputs("Random matris: \n",DosyaOku);
    for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                fprintf(DosyaOku, "%4.0f    ",matris[i][j]);
            }
            fputs("\n",DosyaOku);
        }
        fputs("\n",DosyaOku);
    while(!dur){
        fputs("Asagidaki carpimin basamak sayilari: \n",DosyaOku);
        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                while(carpim[i][j]/carpan>=1.0 || carpim[i][j]/carpan<=-1.0){
                    basamaksayac++;
                    carpan=carpan*10;
                }
                ilkbasamaksayisi[i][j]=basamaksayac;
                carpan=1.0;

                basamaksayac=0;
                fprintf(DosyaOku, "%4d    ", ilkbasamaksayisi[i][j]);
            }
            fputs("\n",DosyaOku);
        }
        fputs("Carpim: \n",DosyaOku);
        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                fprintf(DosyaOku, "%4.0f    ",carpim[i][j]);
            }
            fputs("\n",DosyaOku);
        }

        fputs("\n",DosyaOku);
        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                if(carpim[i][j]==0.0){sayac++;}
            }
        }
        matriscarp(carpim,matris,n,carpim);
        if(sayac==n*n || sayac2>20){dur=1;}
        else{sayac=0;}
        if(hepsibuyuk==n*n){dur=1;}
        sayac2++;

        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                while(carpim[i][j]/carpan>=1.0 || carpim[i][j]/carpan<=-1.0){
                    basamaksayac++;
                    carpan=carpan*10;
                }
                sonbasamaksayisi[i][j]=basamaksayac;
                carpan=1.0;
                basamaksayac=0;
            }
        }
        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                if(sonbasamaksayisi[i][j]>ilkbasamaksayisi[i][j]){hepsibuyuk++;}
            }
        }
        if(hepsibuyuk==n*n){dur=1;}
        hepsibuyuk=0;
    }
    printf("\nRandom matris: \n");
    for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                printf("%-8.2f",matris[i][j]);
            }
            printf("\n");
        }
    printf("\n");
    printf("\nCarpim matris: \n");
    for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                printf("%-8.2f  ",carpim[i][j]);
            }
            printf("\n");
        }
    if(sayac==n*n){
            printf("\nBu matris nilpotent matristir.\n");
            fputs("\nBu matris nilpotent matristir.\n",DosyaOku);
    }
    else{
            printf("\nBu matris nilpotent matris degildir.\n");
            fputs("\nBu matris nilpotent matris degildir.\n",DosyaOku);
    }


}

int main()
{
    srand(time(NULL));
    int gelenislem,i,j;
    gelenislem=cark();
    printf("%d",gelenislem);
    if(gelenislem==0)
    {
        int n,ozdegerboyut=100;
        printf("Gelen islem: Oz Deger Bul\n");
        printf("nxn boyutundaki kare matrisin n sayisini giriniz: ");

        int temp,durum,pozitifmi=pozitifmi=1;
        durum = scanf("%d",&n);
        if(n<=0){pozitifmi=0;}
        else if(n>0){pozitifmi=1;}
        while(durum!=1 || pozitifmi==0){
            while(((temp=getchar()) != EOF && temp != '\n') );
            printf("Gecersiz giris.Lutfen tekrar giris yapiniz: ");
            durum = scanf("%d",&n);
            if(n<=0){pozitifmi=0;}
            else if(n>0){pozitifmi=1;}
        }



        double matris[100][100];
        double ozdegerler[100];
        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                printf("[%d][%d]. eleman: ",i,j);
                durum = scanf("%lf",&matris[i][j]);
                while(durum!=1){
                    while(((temp=getchar()) != EOF && temp != '\n') );
                    printf("Gecersiz giris.Lutfen tekrar giris yapiniz: ");
                    durum = scanf("%lf",&matris[i][j]);
                }

            }
        }


        DosyaOku = fopen("ozdeger.txt","w");
        remove("ozdeger.txt");
        DosyaOku = fopen("ozdeger.txt","w");

        fputs("Girilen matris:\n",DosyaOku);
        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                fprintf(DosyaOku,"%5.3f    ",matris[i][j]);

            }
            fputs("\n",DosyaOku);
        }
        fputs("\n",DosyaOku);
        ozdegerbul(matris,n,ozdegerler,&ozdegerboyut,1);
        if(ozdegerboyut>=1){printf("\nOzdegerler:\n\n|");}
        for(i=0;i<ozdegerboyut;i++){
            printf("%4.2f | ",ozdegerler[i]);
            fprintf(DosyaOku, "%d. Ozdeger = %5.3f\n", i+1,ozdegerler[i]);
        }

        fclose(DosyaOku);

    }else if(gelenislem==1)
    {
        int n,ozdegerboyut=100;
        printf("Gelen islem: Schur teoremi\n");
        printf("nxn boyutundaki kare matrisin n sayisini giriniz: ");
        int temp,durum,pozitifmi=pozitifmi=1;
        durum = scanf("%d",&n);
        if(n<=0){pozitifmi=0;}
        else if(n>0){pozitifmi=1;}
        while(durum!=1 || pozitifmi==0){
            while(((temp=getchar()) != EOF && temp != '\n') );
            printf("Gecersiz giris.Lutfen tekrar giris yapiniz: ");
            durum = scanf("%d",&n);
            if(n<=0){pozitifmi=0;}
            else if(n>0){pozitifmi=1;}
        }

        double matris[100][100];
        double ozdegerler[100];
        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                printf("[%d][%d]. eleman: ",i,j);
                durum = scanf("%lf",&matris[i][j]);
                while(durum!=1){
                    while(((temp=getchar()) != EOF && temp != '\n') );
                    printf("Gecersiz giris.Lutfen tekrar giris yapiniz: ");
                    durum = scanf("%lf",&matris[i][j]);
                }
            }
        }


        DosyaOku = fopen("schur.txt","w");
        remove("schur.txt");
        DosyaOku = fopen("schur.txt","w");

        ozdegerbul(matris,n,ozdegerler,&ozdegerboyut,0);
        schur(matris,n,ozdegerler,ozdegerboyut);
        fclose(DosyaOku);

    }else if(gelenislem==3){
        int n,rastgele,i,j;
        printf("Gelen islem: Nilpotent Matris \n");
        printf("nxn boyutundaki kare matrisin n sayisini giriniz: ");
        int temp,durum,pozitifmi=pozitifmi=1;
        durum = scanf("%d",&n);
        if(n<=0){pozitifmi=0;}
        else if(n>0){pozitifmi=1;}
        while(durum!=1 || pozitifmi==0){
            while(((temp=getchar()) != EOF && temp != '\n') );
            printf("Gecersiz giris.Lutfen tekrar giris yapiniz: ");
            durum = scanf("%d",&n);
            if(n<=0){pozitifmi=0;}
            else if(n>0){pozitifmi=1;}
        }


        DosyaOku = fopen("nilpotent.txt","w");
        remove("nilpotent.txt");
        DosyaOku = fopen("nilpotent.txt","w");

        int rastgelematris[100][100];
        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
            {
                rastgelematris[i][j]=(rand()%30)-15;
            }
        }
        nilpotent(rastgelematris,n);
    }

    return 0;
}
