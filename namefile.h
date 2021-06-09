#ifndef NA_ME
#define NA_ME 1
#include <STDLIB.H>
FILE *readname(char *,FILE *);  // Возвращает указатель name на файл,
FILE *readname(char *,char *);  // Возвращает указатель name на файл,
int readname(char *);           // обработанное имя файла,
double readdat(FILE *);         // Читаем dat, возвращаем прочитанное число
int intread(FILE *);
int simbol(char ch);
char symb(FILE *, char *str = 0, char ch = '#');
char symb(FILE *dat, char ch);
int probel(char *xx, char c = ' '); // Возвращаем номер элемента массива,
                                    // содержащего символ 'c'
double readnamedat(FILE *);     // Читаем name, возвращаем прочитанное число
int countdat(FILE *,int, char *);  // Читаем dat;
long read_l(FILE *);            // Читаем dat, возвращаем прочитанное число
int searchdat(FILE *, char *b = 0, char *e = 0, char *ig = 0, char *s = 0);
                     // Ищем в файле массив, ограниченный двумя заданными
FILE **datab(char *, int narg = 1, char **argc = 0);
FILE **datab(char *, char *, int narg = 1, char **argc = 0);
void tipe(FILE *, FILE *, int r = 0, char c = '!', char s = '*');
void twoline(FILE *tab, char *s, char a = '¤', char b = '2', int p = 0);
char *krap(double *, int *, char *fr = 0); // Цетр. относ. десятичной точки
char *krap(double, int *, char *fr = 0); // Цетр. относ. десятичной точки
char *krap(double, int, int, char fr = 0); // Цетр. относ. десятичной точки
char *krp(double ss, int kkr0 = 0, int kkr1 = -1, char fr = 0);
char *krpk(double w, int n, int m, char fr = ' ');
char *sigkr(double, int *, char *fr = 0);
char *sgn(double, int k = -1);  // в нужном случае отодвигаем знак от числа
char *osgn(double, int k = -1);  // в нужном случае отодвигаем знак от числа
char *sgn(int, int k = -1);
char *spr(double x);  // отодвигаем знак от числа
char *skb(double x);  // отрицательное число заключаем в скобки
char *tbl(int n, char c, char *s = 0, int m = 92);
void namdatab(char *, int narg = 1, char **argc = 0);
void namdatab(char *, char *, int narg = 1, char **argc = 0);
int updown(char *, FILE *, char *up = "", char dn=' '); //строки вверх-вниз
char **drobchr(char *);
int rsdwig(char *s, int k, char ch = ' ');
FILE  *name = 0, *dat = 0, *tab = 0, **ndt;
void namdatab(char *s, int narg, char **argc)
{
  int  c;
  c = (narg < 0) ? -narg : narg;
  ndt = datab(s, c, argc);
  name = ndt[0];
  dat = ndt[1];
  tab = ndt[2];
  if(!dat)
    exit(0);
  if(narg < 0 && name)
    fclose(name);
}
void namdatab(char *ss, char *s, int narg, char **argc)
{
  int  c;
  c = (narg < 0) ? -narg : narg;
  ndt = datab(ss, s, c, argc);
  name = ndt[0];
  dat = ndt[1];
  tab = ndt[2];
  if(!dat)
    exit(0);
  if(narg < 0 && name)
    fclose(name);
}
#endif

