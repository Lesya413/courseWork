#ifndef NA_ME
#define NA_ME 1
#include <STDLIB.H>

FILE *dat,*tab,*name;

void namdatab(char *s);
double readdat(FILE *f);
long read_l(FILE *f);
int intread(FILE *f);

void namdatab(char *s)
{
    char *ss=new char[200];
    name=fopen(s,"r");
    if(name==NULL)
    {
        printf("Error opening file %s",s);
        exit(0);
    }
    fgets(ss,30,name);
    dat=fopen("warmsloj.dat","r");
    if(dat==NULL)
    {
        printf("Error opening file %s",ss);
        exit(0);
    }
    fgets(ss,30,name);
    tab=fopen("warmsloj.tab","w+");
    if(tab==NULL)
    {
        printf("Error opening file %s",ss);
        exit(0);
    }
    fclose(name);
}

double readdat(FILE *f)
{
    float tmp;
    int rez;
    rez=fscanf(f,"%f",&tmp);
    if(!rez)
    {
        printf("Error of reading data from file!");
        exit(0);
    }
    return tmp;
}

long read_l(FILE *f)
{
    long tmp;
    int rez;
    rez=fscanf(f,"%d",&tmp);
    if(!rez)
    {
        printf("Error of reading data from file!");
        exit(0);
    }
    return tmp;
}

int intread(FILE *f)
{
    int tmp;
    int rez;
    rez=fscanf(f,"%d",&tmp);
    if(!rez)
    {
        printf("Error of reading data from file!");
        exit(0);
    }
    return tmp;
}
#endif

