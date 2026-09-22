/* Exact check of the cubic Borwein conjecture (class m = 2 mod 3):
 *   c(m) = [q^m] P_n(q)^3,  P_n = prod_{k<=3n, 3 not| k} (1-q^k).
 *   Verifies: c(m) = 0 for m = 2 mod 3, m < 3n;  c(m) <= 0 (in fact < 0) for 3n <= m <= floor(9n^2/2), m = 2 mod 3.
 * Only coefficients up to D = floor(9n^2/2) are needed (palindromic symmetry handles the rest).
 * Each factor (1-q^k)^3 = 1 - 3q^k + 3q^{2k} - q^{3k} is applied in one in-place pass (high to low index).
 * Usage: ./cubic_check n1 n2      (prints one line per n; flush after each)
 * Build: gcc -O2 -o cubic_check cubic_check.c -lgmp
 */
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
int main(int argc, char **argv){
  if(argc<3){fprintf(stderr,"usage: %s n1 n2\n",argv[0]);return 1;}
  long n1=atol(argv[1]), n2=atol(argv[2]);
  for(long n=n1;n<=n2;n++){
    clock_t t0=clock();
    long D=(9*n*n)/2;
    mpz_t *c=malloc((D+1)*sizeof(mpz_t));
    for(long i=0;i<=D;i++) mpz_init(c[i]);
    mpz_set_ui(c[0],1);
    mpz_t tmp; mpz_init(tmp);
    for(long k=1;k<=3*n;k++){
      if(k%3==0) continue;
      for(long i=D;i>=k;i--){
        /* c[i] = c[i] - 3c[i-k] + 3c[i-2k] - c[i-3k]  using OLD values: process high->low so lower entries are old */
        mpz_mul_ui(tmp,c[i-k],3); mpz_sub(c[i],c[i],tmp);
        if(i>=2*k){ mpz_mul_ui(tmp,c[i-2*k],3); mpz_add(c[i],c[i],tmp); }
        if(i>=3*k){ mpz_sub(c[i],c[i],c[i-3*k]); }
      }
    }
    long viol=0, nonzero_below=0, firstneg=-1, zeros_above=0;
    for(long m=2;m<=D;m+=3){
      int s=mpz_sgn(c[m]);
      if(m<3*n){ if(s!=0) nonzero_below++; }
      else { if(s>0) viol++; if(s==0 && 2*m<=D) zeros_above++; if(s<0 && firstneg<0) firstneg=m; }
    }
    printf("n=%ld D=%ld violations=%ld zeros_in_3n_to_D/2=%ld nonzero_below_3n=%ld first_negative=%ld time=%.1fs\n",
           n,D,viol,zeros_above,nonzero_below,firstneg,(double)(clock()-t0)/CLOCKS_PER_SEC);
    fflush(stdout);
    for(long i=0;i<=D;i++) mpz_clear(c[i]); free(c); mpz_clear(tmp);
  }
  return 0;
}
