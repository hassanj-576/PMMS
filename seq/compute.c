#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include "compute.h"

/// Replace/erase the following line:
#include "ref1.c"

void do_compute(const struct parameters* p, struct results *r)
{
/// Replace/erase the following line:
#include "ref2.c"
	printf("TESTING\n%zd\n%zd\n%zd\n%zd\n%zd\n%d\n%d\n",p->N,p->M,p->maxiter,p->period,p->printreports,p->tinit,p->conductivity);

}
