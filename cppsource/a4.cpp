#include "a4.h"

vector load_data(const char *fp)
{
  FILE *file = fopen(fp, "r");
  int r, i, N = 0;
  while (r != EOF)
  {
    r = fscanf(file, "%*i %i\n", &i);
    if (r == 1)
      N++;
  }
  rewind(file);
  printf("N = %i\n", N);

  vector v = null_vector(N);
  for (i = 0; i < N; i++)
  {
    int val;
    fscanf(file, "%*i %i\n", val);
    VectorSET(v, i, (double) val);
  }
  return v;
}

int main ()
{
  vector dat = load_data("../dat/Agdecay.dat");
  return 0;
}
