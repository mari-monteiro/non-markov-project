#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int i, j, i90, ym1, N, n;
double te1[6000000], pr1[6000000];
double r1, r2, r3, r4, r5, r6;
double m1, m2, dm, alpha;
double pr2[6000000];
double pr3[6000000], aux;
double RR1, RR2, RR3, RR4, RR5, L100, L90, DT;

int main()
{

  FILE *fil0, *fil1, *fil2, *fil3;
  char filename0[650], filename1[650];
  char filename2[650], filename3[650];

  printf("alpha:\n ");
  scanf("%lf", &alpha);

  printf("arquivo de entrada:\n ");
  scanf("%s", filename0);
  fil0 = fopen(filename0, "r");
  //   printf("O arquivo escolhido foi: %s\n",filename0);
  N = 100002;

  // arquivos de saida
  sprintf(filename0, "LambdaxalphaN%ialpha%1.2lf.dat", N, alpha);

  fil2 = fopen(filename0, "w");
  sprintf(filename0, "funcaoNxalphaN%ialpha%1.2lf.dat", N, alpha);

  fil3 = fopen(filename0, "w");

  ym1 = 0;

  while (fscanf(fil0, "%lf %lf %lf %lf ", &r1, &r2, &r3, &r4) != EOF)
  {
    ym1 = ym1 + 1;
    te1[ym1] = r1;
    pr1[ym1] = r2;
  }

  aux = (double)(ym1);
  aux = 0.95 * aux;
  ym1 = (int)ym1;

  DT = te1[2] - te1[1];

  aux = (double)(ym1);
  aux = aux / 2.5;
  i90 = (int)aux;

  for (n = 1; n <= ym1; n = n + 1)
  {
    pr2[n] = log(pr1[n]);
  }
  RR1 = RR2 = RR3 = RR4 = RR5 = 0.;
  for (n = i90; n <= ym1; n = n + 1)
  {
    RR1 = RR1 + te1[n] * pr2[n];
    RR2 = RR2 + te1[n];
    RR3 = RR3 + pr2[n];
    RR4 = RR4 + te1[n] * te1[n];
  }

  L90 = (double)(ym1 - i90);

  L100 = L90 * RR1 - RR2 * RR3;
  L100 = L100 / (L90 * RR4 - (RR2) * (RR2));

  fprintf(fil2, "%15.6lf %17.7g   \n", alpha, fabs(L100));

  // calculo de N

  for (n = 1; n <= ym1; n = n + 1)
  {
    pr2[n] = pr1[n] * pr1[n];
  }

  // calculo da derivada numérica

  for (n = 2; n <= ym1 - 1; n = n + 1)
  {
    pr3[n] = (pr2[n + 1] - pr2[n - 1]) / (2. * DT);
  }

  RR1 = RR2 = 0.;
  for (n = 2; n <= ym1 - 1; n = n + 1)
  {

    if (pr3[n] < 0.)
    {

      RR1 = RR1 + (pr3[n] + pr3[n + 1]) * 0.5 * DT;
    }

    if (pr3[n] > 0.)
    {
      RR2 = RR2 + (pr3[n] + pr3[n + 1]) * 0.5 * DT;
    }
  }

  fprintf(fil3, "%15.6lf %17.7g   \n", alpha, RR2 / fabs(RR1));

  return 0;
}
