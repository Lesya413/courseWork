#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <STDLIB.H>
#include<graphics.h>
#include "fileop.h"

double zr(double x)
{
  double eps = 1e-20;
  if(x < eps && x > -eps)
    x = 0.;
  return(x);
}

long ordnung(long m, long *a)
{
    long i,izm,stmp;
    for(;;)
    {
        izm=0;
        for(i=1;i<m;i++)
        {
            if(a[i]<a[i-1])
            {
                stmp=a[i];
                a[i]=a[i-1];
                a[i-1]=stmp;
                izm++;
            }
        }
        if(izm==0)break;
    }
    return a[0];
}

int main(int narg, char **argc)
{
  double L, F, c, ll, Tl, Ql, lr, Tr, Qr, dt, lmbd = 0.8, eps = 1e-6,
         err = 1., err0, h, kf,  ex, *u0, *u1, *u2, *LL, *cc, *qq,
	       Dx, Kxy, mx, my, max, *u, *v, *w, *x, f, a, b, sy, maxv;
  long N, M = 100,  iter = 0, *Mit, itr, slength;
  int n, bc, pc,  yes = 0, i, j, k, m, correct = 1;
  char *format=new char[20];
  int x1,y1,y2,wx,hy,kks1,kks2,kks3,dx;

  int GD=DETECT,GM;
  kks3=0;
  initgraph(&GD, &GM, "C:\\TC\\BGI");
  wx=getmaxx();
  hy=getmaxy();

  x1=50;
  y1=hy/4;
  y2=3*hy/4;

  /*if(!(narg>0))
  {
      printf("Few many parameters...");
      return 0;
  }
  if(narg>2)
  {
      printf("Too many parameters...");
      return 0;
  }*/
  namdatab("warmsloj.0"); //argc[0]
  F = readdat(dat);  // Площадь сечения             0.0001  (м¤)
  m = intread(dat);  // Количество стержней            3     (штук)
  if(m < 0.)
  {
    m = -m;
    correct = 0;    // Корректор не участвует
  }
  LL = new double[m];
  cc = new double[m];
  qq = new double[m];
  slength=0;
  for(k = 0; k < m; k++)
  {
    LL[k] = readdat(dat);  // Длины стержней          0.1     (м)
    slength+=LL[k];
  }
  for(k = 0; k < m; k++)
    cc[k] = readdat(dat);  // Удельные теплоемкости 3610000 (дж/м.куб?град)
  for(k = 0; k < m; k++)
    qq[k] = readdat(dat);  // Теплопроводности          46   (дж/м?сек?град)
  ll = readdat(dat); // Теплоотдача в воздух слева     100  (дж/кв.м?сек?град)
  Tl = readdat(dat); // Температура воздуха  слева      20   (град)
  Ql = readdat(dat); // Мощность теплоисточника слева  12000 дж/сек
  lr = readdat(dat); // Теплоотдача в воздух справа    100  (дж/кв.м?сек?град)
  Tr = readdat(dat); // Температура воздуха  справа     20   (град)
  Qr = readdat(dat); // Мощность теплоисточника справа  12000 дж/сек
  eps = readdat(dat); // Точность установившегося процесса 1e-4
  dt = readdat(dat); //  Шаг по времени                 0.001 сек
  lmbd = readdat(dat); //Коэффициент усреднения
  n = intread(dat);  //  Число элементарных участков      10
  N = read_l(dat);   //  Число итераций                    2
  M = read_l(dat);   //  Число записываемых итераций      10
  Mit = new long[M];
  for(itr = 0; itr < M; itr++)
    Mit[itr] =  read_l(dat);
  itr = ordnung(M, Mit);
  bc = intread(dat);  // Число цифр до запятой             5
  pc = intread(dat);  // Число цифр после запятой          4
  fclose(dat);

  sprintf(format," %%%d.%df",bc+pc+1,pc);

  if(n < 0)
  {
    n = -n;
    yes += 1;  // пишем безразмерные величины
  }
  if(dt < 0.)
  {
    dt = -dt;
    yes += 2;  // пишем теплофизические параметры
  }
  if(yes > 1)
  {
    fprintf(tab," Длины стержней:  ");
    for(k = 0; k < m; k++)
      fprintf(tab," %lg", LL[k]);
    fprintf(tab,"\n");
    fprintf(tab," Теплоемкости:    ");
    for(k = 0; k < m; k++)
      fprintf(tab," %lg", cc[k]);
    fprintf(tab,"\n");
    fprintf(tab," Теплопроводности: ");
    for(k = 0; k < m; k++)
      fprintf(tab," %lg", qq[k]);
    fprintf(tab,"\n");
    fprintf(tab," Длины новые:     ");
    fflush(tab);
  }
    // Единообразивание
  for(k = 0, L = 0.; k < m; k++)
  {
    if(qq[k] < eps || cc[k] < eps)
    {
      fprintf(tab," Отрицательный теплофизический параметр.\n");
      goto ABC;
    }
    cc[k] = qq[k]/cc[k];   // Температуропроводности
    err0 = LL[k]*sqrt(cc[0]/cc[k]);
    if(yes > 1)
      fprintf(tab," %lg ", err0);
    L += err0;
  }
  if(yes > 1)
    fprintf(tab,"\n");
     // Обезразмеривание
  h = 1./(double)n;

  Ql *= h*L/F/qq[0];
  Qr *= h*L/F/qq[m-1]*sqrt(cc[m-1]/cc[0]);
  ll *= Tl*h*L/qq[0];
  lr *= h*Tr*L/qq[m-1]*sqrt(cc[m-1]/cc[0]);
  Ql += ll;
  Qr += lr;

  c = L*L/cc[0];    // Размерность: секунда; теперь длина стержня равна 1.
  dt /= c;
  kf = 2./h/h;
  ex = exp(-kf*dt);

  if(yes%2 == 1)
  {
    fprintf(tab,"   Эталон времени %lg сек\n", c);
    fprintf(tab,"     Безразмерные величины\n");
    fprintf(tab,"   Шаг по времени                  %lg\n", dt);
    fprintf(tab,"   Левый коэффициент границы       %lg\n", 1.+ll);
    fprintf(tab,"   Температура воздуха слева       %lg град\n",Tl);
    fprintf(tab,"   Теплоисточник слева             %lg град\n", Ql);
    fprintf(tab,"   Правый коэффициент границы      %lg\n", 1.+lr);
    fprintf(tab,"   Температура воздуха справа      %lg град\n",Tr);
    fprintf(tab,"   Теплоисточник справа            %lg град\n", Qr);
    fprintf(tab,"   Показатель                      %lg\n", kf);
    fprintf(tab,"   Экспонента                      %lg\n", ex);
    fflush(tab);
  }
  u0 = new double[n+1];
  u1 = new double[n+1];
  u2 = new double[n+1];
  x = new double[n+1];
  u = u2;
  w = u1;
  v = u0;
  for(i = 0, x[0] = 0.; i <= n; i++)
  {
    u[i] = v[i] = w[i] = 0.;
    if(i)
      x[i] = x[i-1] + h;
  }
  if(!correct)
    fprintf(tab,"   Работает только предиктор.\n");
  for(iter = 0, itr = 0; iter < N && err > eps; iter++)
  {
    if(kbhit() && getch() == 27)
      break;        //    goto ABC;
    w = u;
    u = v;
    v = w;
    w = u1;
    err0 = err;
    if(correct)
    {
      for(i = 1; i < n; i++)    // Предиктор
      {
        f = (u[i-1] + u[i+1])/2.;
        w[i] = u[i]*ex + f*(1. - ex);
      }
      w[0] = (w[1] + Ql)/(1. + ll);
      w[n] = (w[n-1] + Qr)/(1. + lr);

      for(i = 1; i < n; i++)    // Корректор
      {
        f = (u[i-1] + u[i+1] + w[i-1] + w[i+1])/4.;
        v[i] = u[i]*ex + f*(1. - ex);
      }
    }
    else
      for(i = 1; i < n; i++)    // Предиктор
      {
        f = (u[i-1] + u[i+1])/2.;
        v[i] = u[i]*ex + f*(1. - ex);
      }
    v[0] = (v[1] + Ql)/(1. + ll);
    v[n] = (v[n-1] + Qr)/(1. + lr);
    if(Mit[itr] == iter)
    {
      fprintf(tab,"  Итерация %ld.\n", iter);
      for(i = 0, mx = 0., my = 0.; i <= n; i++)
      {
        if(i && !(i%5))
          fprintf(tab,"\n");
        fprintf(tab,format, v[i]);
      }
      itr++;
      fprintf(tab,"\n");
    if(kks3==0)
    {
        kks3++;
        setcolor(13);
        line(x1,0,x1,hy);
        line(0,y1,wx,y1);
        line(0,y2,wx,y2);
        outtextxy(x1+2,10,"Y");
        outtextxy(wx-20,y1+2,"X");
        outtextxy(wx-20,y2+2,"X");
        kks2=0;
        dx=(wx-50)/n;
        maxv=abs(v[0]);
        sy=0;

        for(kks1=1;kks1<n;kks1++)
        {
            if(abs(v[kks1])>maxv) maxv=abs(v[kks1]);
        }

        if(maxv>0)
        {
            if(maxv>(h/4)) sy=((double)hy/4.0)/((double)maxv*1.1);
            else sy=((double)maxv*0.9)/((double)hy/4.0);
        }
        else sy=1;

        for(kks1=1;kks1<n;kks1++)
        {
            line(x1+kks2,y1-v[kks1-1]*sy,x1+kks2+dx,y1-v[kks1]*sy);
            kks2+=dx;
        }
    }
    }
    for(i = 0, err = 0.; i <= n; i++)
      err += fabs(v[i] - u[i]);
    err = lmbd*err0 + (1.- lmbd)*err;
  }
  fprintf(tab,"  Итерация %ld заключительная. Точность %13.6le.\n",
                                                   iter, err);
  for(i = 0, mx = 0., my = 0.; i <= n; i++)
  {
    if(i && !(i%5))
      fprintf(tab,"\n");
    fprintf(tab,format, v[i]);
    mx += x[i];
    my += v[i];
  }
  mx /= (double)(n+1);
  my /= (double)(n+1);
  for(i = 0, Dx = 0., Kxy = 0.; i <= n; i++)
  {
    a = x[i] - mx;
    Dx += a*a;
    Kxy += a*(v[i] - my);
  }
  fprintf(tab,"\n");
  Kxy /= (double)(n+1);
  Dx /= (double)(n+1);
  a = Kxy/Dx;
  if(a > 0.)
    fprintf(tab," Линейная регрессия: y = %lg + %lg?(x - %lg).\n", my,
		                                 a, mx);
  else
    fprintf(tab," Линейная регрессия: y = %lg - %lg?(x - %lg).\n", my,
		                                 -a, mx);
  for(Dx = 0., max = 0., i = 0; i <= n; i++)
  {
    b = fabs(my + a*(x[i] - mx) - v[i]);
    Dx += b;
    if(b > max)
      max = b;
  }
  Dx /= (double)(n+1);
  fprintf(tab," Максимальное абсолютное отклонение %lg,\n", max);
  fprintf(tab," Среднее   абсолютное    отклонение %lg.\n", Dx);
  f = 2.*h;
  fprintf(tab,"        Центральные отнесенные разности.\n  ");
  for(j = 0; j < bc+pc; j++)
    fprintf(tab," ");
  for(i = 1, mx = 0., my = 0.; i < n; i++)
  {
    if(i && !(i%5))
      fprintf(tab,"\n");
    fprintf(tab,format, (v[i+1] - v[i-1])/f);
  }
  fprintf(tab,"\n");
ABC:
  fflush(tab);
  fclose(tab);

        setcolor(13);
        kks2=0;
        maxv=abs(v[0]);
        sy=0;

        for(kks1=1;kks1<n;kks1++)
        {
            if(abs(v[kks1])>maxv) maxv=abs(v[kks1]);
        }

        if(maxv>0)
        {
            if(maxv>(h/4)) sy=((double)hy/4.0)/((double)maxv*1.1);
            else sy=((double)maxv*0.9)/((double)hy/4.0);
        }
        else sy=1;

        for(kks1=1;kks1<n;kks1++)
        {
            line(x1+kks2,y2-v[kks1-1]*sy,x1+kks2+dx,y2-v[kks1]*sy);
            kks2+=dx;
        }


  if(u0) delete [] u0;
  if(u1) delete [] u1;
  if(u2) delete [] u2;
  if(x) delete [] x;
  if(LL) delete [] LL;
  if(qq) delete [] qq;
  if(cc) delete [] cc;
  if(Mit) delete [] Mit;

  system("pause");
  closegraph();

  return 0;
}

