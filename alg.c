#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <ctype.h>
#include <math.h>
#define NEXT 6
/*
// функция ввода матрицы, возвращает ссылку на нее и её размеры: str stb
double **getm(int *str,int *stb){
    int k;
    printf ("\nчисло строк(>0): ");
    k=scanf("%d",str);
    if (!k || k==EOF)
        return NULL;
    printf ("\nчисло столбцов(>0): ");
    k=scanf("%d",stb);
    if(!k||k==EOF)
        return NULL;
    if (*str<1 || *stb<1)
        return NULL;
    printf("OUTPUT \n");
    double **matr=(double**) malloc((*str)*sizeof(double*));
    if (!matr)
        return NULL;
    for (int i=0;i<(*str);i++){
        matr[i]=(double*) malloc((*stb)*sizeof(double));
        if (!matr[i]){
            for (int j=0;j<i;j++)
                free(matr[j]);
            free(matr);
            return NULL;
        }
    }
    for (int i=0;i<(*str);i++){
        for(int j=0;j<(*stb);j++){
            printf("\na(%d,%d)=",i,j);
            k=scanf("%lf",&matr[i][j]);
            if (!k || k==EOF){
                for (int j=0;j<(*str);j++)
                    free(matr[j]);
                free(matr);
                return NULL;
            }
        }
    }
        return matr;
}
*/
//матрица 1
void getm1(double ***matr, int *m){
    int n,k;
    double alf=0.0;
    printf("Ввод расширенной матрицы: \n");
    printf ("Введите n(параметр размера матрицы,натуральное число): ");
    k=scanf("%d",m);
    if (!k || k==EOF)
        exit(1);
    n=*m;

    *matr=(double**) malloc((n-1)*sizeof(double*));
    if (!(*matr))
        exit(1);
    for (int i=0;i<n-1;i++){
        (*matr)[i]=(double*) malloc(n*sizeof(double));
        if (!(*matr)[i]){
            for (int j=0;j<i;j++)
                free((*matr)[j]);
            free(*matr);
            exit(1);
        }
    }

    for(int i=0;i<n-1;i++){
        (*matr)[i][0]=1.0;
        alf=(double) 1/(i+2);

        (*matr)[i][n-1]=alf;
        for(int p=0;p<n-3;p++)
            (*matr)[i][n-1]*=alf;
        for(int j=1;j<n-1;j++){
            (*matr)[i][j]=alf;
            for(int p=0;p<j-1;p++)
                (*matr)[i][j]*=alf;
        }
    }

}

//матрица 2
void getm2(double ***matr, int *m){
    int n,k;
    printf("Ввод расширенной матрицы: \n");
    printf ("Введите n(параметр размера матрицы,натуральное число): ");
    k=scanf("%d",m);
    if (!k || k==EOF)
        exit(1);
    n=*m;

    *matr=(double**) malloc((n+1)*sizeof(double*));
    if (!(*matr))
        exit(1);
    for (int i=0;i<n+1;i++){
        (*matr)[i]=(double*) malloc((n+2)*sizeof(double));
        if (!(*matr)[i]){
            for (int j=0;j<i;j++)
                free((*matr)[j]);
            free(*matr);
            exit(1);
        }
    }

    for(int i=0;i<n+1;i++){

        for(int j=0;j<i;j++){
            (*matr)[i][j]=(double) -(i-j);
        }
        (*matr)[i][i]=1.0;
        (*matr)[i][n+1]=1.5;
        (*matr)[i][n]=1.0;
        for(int j=i+1;j<n;j++){
            (*matr)[i][j]=0.0;
        }
    }

}
//функция вывода матрицы matr размера: m n
void printm(double **matr,int m,int n){
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++)
            printf("%.3lf ",matr[i][j]);
        printf("\n");
    }
    printf("\n");
}
//метод гаусса с выбором главного элемента
void Gauss(double **mas1, int m, double **solve,long long *counter){
    int k=0,nmax=0,*nsol=NULL,p=0,p1=0;
    double alf=0.0,max=0.0,lf=0.0,**mas=NULL,p2=0.0;
    mas=(double **) malloc(m*sizeof(double *));
    for (int i=0;i<m;i++){
        mas[i]=(double*) malloc((m+1)*sizeof(double));
        if (!mas[i]){
            for (int j=0;j<i;j++)
                free(mas[j]);
            free(mas);
            exit(1);
        }
    }
    for(int i=0;i<m;i++)
        for(int j=0;j<m+1;j++)
            mas[i][j]=mas1[i][j];
    *counter=0;
    //прямой ход
    nsol=malloc(m*sizeof(int));
    for(int i=0;i<m;i++)
        nsol[i]=i;
    for(int i=0;i<m;i++){
        k=i;
        nmax=i;
        max=mas[i][i];
        max=max>0.0?max:max*(-1);
        while(k<m){
            lf=mas[i][k]>0.0?mas[i][k]:mas[i][k]*(-1);
            if(max<lf){
                max=lf;
                nmax=k;
            }
            k++;
        }
        if (nmax!=i){
            p=nsol[i];
            nsol[i]=nsol[nmax];
            nsol[nmax]=p;
            for(int j=0;j<m;j++){
                p2=mas[j][i];
                mas[j][i]=mas[j][nmax];
                mas[j][nmax]=p2;
            }
        }
        alf=mas[i][i];
        for(int j=i;j<m+1;j++)
            mas[i][j]/=alf;
        *counter+=(long long) (m-i+1);
        for(int p=i+1;p<m;p++)
            for(int s=m;s>i-1;s--)
                mas[p][s]-=mas[i][s]*mas[p][i];
        *counter+=(long long) (m-i+1)*(m-i-1);
    }
    //обратный ход
    *solve = malloc(m*sizeof(double));
    for(int i=m-1;i>-1;i--){
        p=nsol[i];
        (*solve)[p]=mas[i][m];
        for(int j=i;j<m-1;j++){
            p1=nsol[j+1];
            (*solve)[p]-=(*solve)[p1]*mas[i][j+1];
        }
        *counter+=(long long)(m-1-i);
    }
    for (int i=0;i<m;i++)
        free(mas[i]);
    free(mas);
    free(nsol);
}

//метод Гаусса
void Gaussl(double **mas1, int m, double **solve,long long *counter){
    int k=0,*nsol=NULL,p=0,p1=0;
    double alf=0.0,**mas=NULL,p2=0.0;
    mas=(double **) malloc(m*sizeof(double *));
    for (int i=0;i<m;i++){
        mas[i]=(double*) malloc((m+1)*sizeof(double));
        if (!mas[i]){
            for (int j=0;j<i;j++)
                free(mas[j]);
            free(mas);
            exit(1);
        }
    }
    for(int i=0;i<m;i++)
        for(int j=0;j<m+1;j++)
            mas[i][j]=mas1[i][j];
    *counter=0;
    //прямой ход
    nsol=malloc(m*sizeof(int));
    for(int i=0;i<m;i++)
        nsol[i]=i;
    for(int i=0;i<m;i++){
        k=i;
        while(k<m && mas[i][k]==0)
            k++;
        if (k!=i){
            p=nsol[i];
            nsol[i]=nsol[k];
            nsol[k]=p;
            for(int j=0;j<m;j++){
                p2=mas[j][i];
                mas[j][i]=mas[j][k];
                mas[j][k]=p2;
            }
        }
        alf=mas[i][i];
        for(int j=i;j<m+1;j++)
            mas[i][j]/=alf;
        *counter+=(long long) (m-i+1);
        for(int p=i+1;p<m;p++)
            for(int s=m;s>i-1;s--)
                mas[p][s]-=mas[i][s]*mas[p][i];

        *counter+=(long long) (m-i+1)*(m-i-1);
    }
    //обратный ход
    *solve = malloc(m*sizeof(double));
    for(int i=m-1;i>-1;i--){
        p=nsol[i];
        (*solve)[p]=mas[i][m];
        for(int j=i;j<m-1;j++){
            p1=nsol[j+1];
            (*solve)[p]-=(*solve)[p1]*mas[i][j+1];
        }
        *counter+=(long long)(m-1-i);
    }
    free(nsol);
    for (int i=0;i<m;i++)
        free(mas[i]);
    free(mas);
}



//метод отражений Хаусхолдера
void refl(double **mas1, int m, double **solve,long long *counter){
    int p2=0,p3=0,fl=0,sign=0;
    double alf=0.0,norm=0.0;
    double **u=NULL,*w=NULL,*w1=NULL,*scal=NULL,**mas=NULL;
    u=(double**) malloc(m*sizeof(double*));
    w=(double*)malloc(m*sizeof(double));
    w1=(double*)malloc(m*sizeof(double));
    scal=(double*)malloc((m+1)*sizeof(double));
    mas=(double **) malloc(m*sizeof(double *));
    for (int i=0;i<m;i++){
        mas[i]=(double*) malloc((m+1)*sizeof(double));
        if (!mas[i]){
            for (int j=0;j<i;j++)
                free(mas[j]);
            free(mas);
            exit(1);
        }
    }
    for(int i=0;i<m;i++)
        for(int j=0;j<m+1;j++)
            mas[i][j]=mas1[i][j];
    *counter=0;
    if (!(u))
        exit(1);
    for (int i=0;i<m;i++){
        u[i]=(double*) malloc((m+1)*sizeof(double));
        if (!u[i]){
            for (int j=0;j<i;j++)
                free(u[j]);
            free(u);
            exit(1);
        }
    }
    //построение верхнетреугольной матрицы
    for(int i=0;i<m-1;i++){
        fl=0;
        alf=0.0;
        for(int  j=i+1;j<m;j++)
            if(mas[j][i]!=0){
                fl=1;
                break;
            }
        if(fl){
        //построение матрицы Хаусхолдера(отражения)
            //построение w
            if (mas[i][i]>=0)
                sign=1;
            else
                sign=-1;
            for(int j=i;j<m;j++)
                alf+=mas[j][i]*mas[j][i];
            *counter+=(long long)(m-i);
            norm=sqrt(alf);
            p3=sign*norm;
            p2=mas[i][i]+p3;
            alf=alf+p2*p2-mas[i][i]*mas[i][i];
            alf=sqrt(alf);
            w[i]=(double) (mas[i][i]+p3)/alf;
            *counter+=4;
            for(int j=i+1;j<m;j++)
                w[j]=mas[j][i]/alf;
            *counter+=(long long)(m-i-1);
            for(int j=i;j<m;j++)
                w1[j]=2*w[j];
            *counter+=(long long)(m-i);
            for (int j=i;j<m+1;j++)
                scal[j]=0;
            for(int p=i;p<m+1;p++)
                for(int k=i;k<m;k++)
                    scal[p]+=mas[k][p]*w[k];
            *counter+=(long long)(m-i)*(m+1-i);
            for(int p=i;p<m;p++)
                for(int k=i;k<m+1;k++)
                    u[p][k]=w1[p]*scal[k];
            *counter+=(long long)(m+1-i)*(m-i);
            //получение новой расширенной матрицы системы на i-й итерации
            for(int p=i;p<m;p++)
                for(int k=i;k<m+1;k++)
                    mas[p][k]=mas[p][k]-u[p][k];
        }
    }
    //обратный ход
    *solve = malloc(m*sizeof(double));
    for(int i=m-1;i>-1;i--){
        (*solve)[i]=mas[i][m];
        for(int j=i;j<m-1;j++){
            (*solve)[i]-=(*solve)[j+1]*mas[i][j+1];
        }
        (*solve)[i]/=mas[i][i];
        *counter+=(long long)(m-i);
    }
    free(scal);
    free(w);
    free(w1);
    for (int i=0;i<m;i++)
        free(u[i]);
    free(u);
    for (int i=0;i<m;i++)
        free(mas[i]);
    free(mas);
}
//поиск нормы невязки
double dscr(double **mas, int m, double *solve){
    double *res=NULL,n=0;
    res=(double *) malloc(m*sizeof(double));
    for(int i=0;i<m;i++)
        res[i]=0;
    for(int i=0;i<m;i++)
        for(int j=0;j<m;j++)
            res[i]+=mas[i][j]*solve[j];
    for(int i=0;i<m;i++)
        res[i]-=mas[i][m];
    for(int i=0;i<m;i++)
        n+=res[i]*res[i];
    n=sqrt(n);
    free(res);
    return n;
}

//Метод Зейделя
void Zeidel(double **mas, int m, double **solve){
    int k=0;
    double *solvepr=NULL,alf1=0.0,alf2=0.0;


    *solve = malloc(m*sizeof(double));
    solvepr = malloc(m*sizeof(double));
    //начальное приближение
    for(int i=0;i<m;i++)
        solvepr[i]=0.0;
    while(k<6){
        for(int i=0;i<m;i++){
            alf1=alf2=0.0;
            for(int j=0;j<i;j++)
                alf1+=mas[i][j]*(*solve)[j];
            for(int j=i+1;j<m;j++)
                alf2+=mas[i][j]*solvepr[j];
            (*solve)[i]=(mas[i][m]-alf1-alf2)/mas[i][i];

        }
        for(int i=0;i<m;i++)
            solvepr[i]=(*solve)[i];
         k++;

        printf("Номер итерации: %d\n",k);
        printf("Норма невязки: %.20lf\n",dscr(mas,m,(*solve)));


    }
    free(solvepr);
}



int main(void){
    int m1=0,m2=0,a,b;
    double **a1=NULL,*x=NULL,**a2=NULL;
    long long k=0;
    while(1){
        printf("         Меню \n");
        printf("1 Первая матрица \n" );
        printf("2 Вторая матрица \n" );
        printf("3 Метод Гаусса  \n" );
        printf("4 Метод Гаусса с выбором главного элемента\n");
        printf("5 Метод отражений \n" );
        printf("6 Метод Зейделя\n");
        printf("Любые другие два числа: Выход из программы \n");
        printf("Введите матрицу и метод (два числа):");
        while(scanf("%d%d",&a,&b)!=2){
            while(getchar()!='\n');
            printf("Некорректный ввод, повторите попытку.\n");
            printf("Введите матрицу и метод (два числа):");
        }
        if ((a!=1 && a!=2)||(b!=3 && b!=4 && b!=5 && b!=6)){
            printf("Завершение программы.\n");
            exit(0);
        }
        if (a==1){
            getm1(&a1,&m1);
            printf("Полученная матрица\n");
            printm(a1,m1-1,m1);
            if(b==3)
                Gaussl(a1,m1-1,&x,&k);
            else if(b==4)
                Gauss(a1,m1-1,&x,&k);
            else if(b==5)
                refl(a1,m1-1,&x,&k);
            else
                Zeidel(a1,m1-1,&x);
            if(b!=6){
                printf("Норма невязки: %.20lf\n",dscr(a1,m1-1,x));
                printf("Количество умножений/делений: %lld\n",k);
            }
            for(int i=0;i<m1-1;i++)
                free(a1[i]);
            free(a1);
        }
        else{
            getm2(&a2,&m2);
            printf("Полученная матрица\n");
            printm(a2,m2+1,m2+2);
            if(b==3)
                Gaussl(a2,m2+1,&x,&k);
            else if(b==4)
                Gauss(a2,m2+1,&x,&k);
            else if(b==5)
                refl(a2,m2+1,&x,&k);
            else
                Zeidel(a2,m2+1,&x);
            if(b!=6){
                printf("Норма невязки: %.20lf\n",dscr(a2,m2+1,x));
                printf("Количество умножений/делений: %lld\n",k);
            }
            for(int i=0;i<m2+1;i++)
                free(a2[i]);
            free(a2);
        }
        free(x);
        printf("Нажмите любую клавишу для продолжения \n");
        while(getchar()!='\n');
        getchar();



    }
}







/*
    getm1(&a1,&m1);
    printm(a1,m1-1,m1);
    printf(" \n");
    Gaussl(a1,m1-1,&x,&k);
   // refl(a1,m1-1,&x,&k);
    printf("Норма невязки: %.20lf\n",dscr(a1,m1-1,x));
 //   printm(a1,m1-1,m1);
  //  for (int i=0;i<m1-1;i++)
    //        printf("%.50lf ",x[i]);
//    printf(" \n");
*/
/*
    getm2(&a2,&m2);
    printm(a2,m2+1,m2+2);
    printf(" \n");
    //Gaussl(a2,m2+1,&x,&k);
    refl(a2,m2+1,&x,&k);
*/
/*
    printf("Solve: \n");

    for (int i=0;i<m2+1;i++)
            printf("%.20lf ",x[i]);
    printf(" \n");
*/
/*
    printf("Норма невязки: %.20lf\n",dscr(a2,m2+1,x));
    printf("Количество умножений/делений: %lld\n",k);
*/
    /*
    printf("Ввод расширенной матрицы: \n");
    getm2(&a2,&m2);
    printm(a2,m2+1,m2+2);
    printf(" \n");
    Gaussl(a2,m2+1,&x1);
    printm(a2,m2+1,m2+2);
    printf("Solve: \n");
    for (int i=0;i<m2+1;i++)
            printf("%.5lf ",x1[i]);
    printf(" \n");
    */

