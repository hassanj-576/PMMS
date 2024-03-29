/* Syed Hassan Jalil - 11161604
 * Jedda Boyle - 11398221
 */

#include "../src/input.h"
#include "../src/output.h"

#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "math.h"
#include "float.h"

// Vector instructions.
#include <immintrin.h>

// The constansts used to weight the current cells neighbours.
#define WEIGHT_ADJACENT_NEIGHBOURS 0.5857864376269 // 0495119831127579030192
#define WEIGHT_DIAGONAL_NEIGHBOURS 0.4142135623730 // 9504880168872420969807

typedef int bool;
#define TRUE 1
#define FALSE 0

// Use this to store the three conductivity weights for a single cell.
typedef struct conductivity conductivity;
struct conductivity {
    // x := same cell conductivity, y := adjacent cells conductivity, z:= diagonal cells conductivity.
    double x, y, z;
};

typedef struct heat_cylinder heat_cylinder;
struct heat_cylinder {
    size_t N;
    size_t M;
    double **heat;
    double **temp;
    conductivity **cond;
};

heat_cylinder* init_heat_cylinder (const struct parameters *p) {

    heat_cylinder *hc = malloc(sizeof(heat_cylinder));

    size_t N;
    size_t M;
    double **heat;
    double **temp;
    conductivity **cond;

    // Rotate the grid, which is why N and M are swapped.
    // Add 2 to M for padding
    N = p->M;
    M = p->N + 2;

    // Allocate the memory required for heat, temp and conductivity arrays.
    heat = malloc(sizeof(double*) * N);
    temp = malloc(sizeof(double*) * N);
    cond = malloc(sizeof(conductivity*) * N);
    for (size_t i = 0; i < N; i++) {
        heat[i] = malloc(sizeof(double) * M);
        temp[i] = malloc(sizeof(double) * M);
        cond[i] = malloc(sizeof(conductivity) * M);
    }

    // Calculate each cells three conductivity co-efficients.
    double x;
    for (size_t i = 0; i < p->N; i++) {
        for (size_t j = 0; j < p->M; j++) {
            heat[j][i+1] = p->tinit[i * p->M + j];
            x = p->conductivity[i * p->M + j];
            cond[j][i+1] = (conductivity) {
                .x = x,
                .y = WEIGHT_ADJACENT_NEIGHBOURS * (1 - x) * .25,
                .z = WEIGHT_DIAGONAL_NEIGHBOURS * (1 - x) * .25
            };
        }
    }

    // Pad the edge that doesn't wrap around.
    for (size_t i = 0; i < N; i++) {
        heat[i][0] = heat[i][1];
        temp[i][0] = heat[i][1];
        heat[i][M-1] = heat[i][M-2];
        temp[i][M-1] = heat[i][M-2];
    }

    hc->N = N;
    hc->M = M;
    hc->heat = heat;
    hc->temp = temp;
    hc->cond = cond;

    return hc;
}

void free_heat_cyclinder(heat_cylinder *hc) {
    for (size_t i = 0; i < hc->N; i++) {
            free(hc->heat[i]);
            free(hc->temp[i]);
            free(hc->cond[i]);
    }
    free(hc->heat);
    free(hc->temp);
    free(hc->cond);
    free(hc);
}

int mod(int a, int b) {
    int r = a % b;
    return r < 0 ? r + b : r;
}

#define AVERAGES(x1, x2, x3, y1, y2, y3) {\
    adjacent_sum = heat[x2][y1] + heat[x1][y2] + heat[x2][y3] + heat[x3][y2]; \
    diagonal_sum = heat[x1][y1] + heat[x1][y3] + heat[x3][y1] + heat[x3][y3]; \
    temp[x2][y2] = cond[x2][y2].x * heat[x2][y2] + \
                   cond[x2][y2].y * adjacent_sum + \
                   cond[x2][y2].z * diagonal_sum ; \
}

// #define AVERAGES(x1, x2, x3, y1, y2, y3) {\
//     adjacent_sum = heat[x2][y1] + heat[x1][y2] + heat[x2][y3] + heat[x3][y2]; \
//     diagonal_sum = heat[x1][y1] + heat[x1][y3] + heat[x3][y1] + heat[x3][y3]; \
//     temp[x2][y2] = (cond[x2][y2].x * heat[x2][y2]) + \
//                    ((( (sqrt(2.)) / (sqrt(2) + 1) )   )     * (1 - cond[x2][y2].x) * adjacent_sum) + \
//                    (   (( 1. / (sqrt(2) + 1) )   ) * (1 - cond[x2][y2].x) * diagonal_sum); \
// }

void __dissipate (size_t N,
                  size_t M,
                  double **heat,
                  double **temp,
                  conductivity **cond) {

    double adjacent_sum;
    double diagonal_sum;

    // Update the central area, not considering the edge cases.
    for (size_t i = 1; i < N - 1; i++) {
        for (size_t j = 1; j < M - 1; j++) {
            AVERAGES(i-1, i, i+1, j-1, j, j+1);
        }
    }
    for (size_t j = 2; j < M - 2; j++) {
        // Top edge.
        AVERAGES(N-1, 0, 1, j-1, j, j+1);

        // Bottom edge.
        AVERAGES(N-2, N-1, 0, j-1, j, j+1);
    }

    // Top left corner.
    AVERAGES(N-1, 0, 1, 0, 1, 2);

    // Top right corner
    AVERAGES(N-1, 0, 1, M-3, M-2, M-1);

    // Bottom left corner.
    AVERAGES(N-2, N-1, 0, 0, 1, 2);

    // Bottom right corner.
    AVERAGES(N-2, N-1, 0, M-3, M-2, M-1);
}

void dissipate(heat_cylinder *hc) {
    __dissipate(hc->N, hc->M, hc->heat, hc->temp, hc->cond);
    double **swap = hc->heat;
    hc->heat = hc->temp;
    hc->temp = swap;
}

void __vec_dissipate (size_t N,
                  size_t M,
                  double **heat,
                  double **temp,
                  conductivity **cond) {

    double adjacent_sum;
    double diagonal_sum;

    // Update the central area, not considering the edge cases.
    int j;
    __m256d adj_values, diag_values, reg, values;
    for (size_t i = 1; i < N - 1; i++) {
        for (j = 1; j <= M - 5 ; j+=4) {
            /****************/
            adj_values  = _mm256_setr_pd(heat[i][j-1], heat[i][j], heat[i][j+1], heat[i][j+2]); //left

            reg = _mm256_setr_pd(heat[i-1][j], heat[i-1][j+1], heat[i-1][j+2], heat[i-1][j+3]); //above
            adj_values = _mm256_add_pd(reg, adj_values);

            reg = _mm256_setr_pd(heat[i][j+1], heat[i][j+2], heat[i][j+3], heat[i][j+4]);//right
            adj_values = _mm256_add_pd(reg, adj_values);

            reg = _mm256_setr_pd(heat[i+1][j], heat[i+1][j+1], heat[i+1][j+2], heat[i+1][j+3]);//below
            adj_values = _mm256_add_pd(reg, adj_values);

            reg = _mm256_setr_pd(cond[i][j].y, cond[i][j+1].y , cond[i][j+2].y , cond[i][j+3].y);
            adj_values = _mm256_mul_pd(reg, adj_values);

            /****************/
            diag_values  = _mm256_setr_pd(heat[i-1][j-1], heat[i-1][j], heat[i-1][j+1], heat[i-1][j+2]); //left

            reg = _mm256_setr_pd(heat[i-1][j+1], heat[i-1][j+2], heat[i-1][j+3], heat[i-1][j+4]); //above
            diag_values = _mm256_add_pd(reg, diag_values);

            reg = _mm256_setr_pd(heat[i+1][j-1], heat[i+1][j], heat[i+1][j+1], heat[i+1][j+2]);//right
            diag_values = _mm256_add_pd(reg, diag_values);
            //
            reg = _mm256_setr_pd(heat[i+1][j+1], heat[i+1][j+2], heat[i+1][j+3], heat[i+1][j+4]);//below
            diag_values = _mm256_add_pd(reg, diag_values);
            //
            reg = _mm256_setr_pd(cond[i][j].z, cond[i][j+1].z, cond[i][j+2].z, cond[i][j+3].z);
            diag_values = _mm256_mul_pd(reg, diag_values);

            /****************/
            values = _mm256_setr_pd(heat[i][j], heat[i][j+1], heat[i][j+2], heat[i][j+3]);
            reg = _mm256_setr_pd(cond[i][j].x, cond[i][j+1].x, cond[i][j+2].x, cond[i][j+3].x);
            values = _mm256_mul_pd(reg, values);

            values = _mm256_add_pd(values, diag_values);
            values = _mm256_add_pd(values, adj_values);

            _mm256_store_pd (&temp[i][j], values);


        }
        // Iterate over the remaining values that don't fit into a 256 bit register.
        for ( ; j < M - 1; j++) {

            AVERAGES(i-1, i, i+1, j-1, j, j+1);
        }

    }
    for (size_t j = 2; j < M - 2; j++) {
        // Top edge.
        AVERAGES(N-1, 0, 1, j-1, j, j+1);

        // Bottom edge.
        AVERAGES(N-2, N-1, 0, j-1, j, j+1);
    }

    // Top left corner.
    AVERAGES(N-1, 0, 1, 0, 1, 2);

    // Top right corner
    AVERAGES(N-1, 0, 1, M-3, M-2, M-1);

    // Bottom left corner.
    AVERAGES(N-2, N-1, 0, 0, 1, 2);

    // Bottom right corner.
    AVERAGES(N-2, N-1, 0, M-3, M-2, M-1);

}

void vec_dissipate(heat_cylinder *hc) {
    __vec_dissipate(hc->N, hc->M, hc->heat, hc->temp, hc->cond);
    double **swap = hc->heat;
    hc->heat = hc->temp;
    hc->temp = swap;
}

void __dissipate2 (size_t N,
                   size_t M,
                   double **heat,
                   double **temp,
                   conductivity **cond) {

    double adjacent_sum;
    double diagonal_sum;

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 1; j < M - 1; j++) {
            AVERAGES(mod(i-1, N), i, mod(i+1, N), mod(j-1, N), j, mod(j+1, N));
        }
    }

}

void dissipate2(heat_cylinder *hc) {
    __dissipate2(hc->N, hc->M, hc->heat, hc->temp, hc->cond);
    double **swap = hc->heat;
    hc->heat = hc->temp;
    hc->temp = swap;
}

double calc_time(struct timespec start, struct timespec end) {
	struct timespec temp;
	if ((end.tv_nsec - start.tv_nsec) < 0) {
		temp.tv_sec = end.tv_sec - start.tv_sec - 1;
		temp.tv_nsec = 1e+9 + end.tv_nsec - start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec - start.tv_sec;
		temp.tv_nsec = end.tv_nsec - start.tv_nsec;
	}
	return temp.tv_sec + (temp.tv_nsec * 1e-9);
}

void __compute_results(struct results *r,
                       size_t N,
                       size_t M,
                       double **heat,
                       double **temp,
                       struct timespec *start_timer,
                       struct timespec *end_timer,
                       size_t iter_num) {

    double total_heat = 0.f;
    double min_heat = DBL_MAX;
    double max_heat = DBL_MIN;
    double max_diff = DBL_MIN;
    double diff;
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 1; j < M-1; j++) {
            diff = fabs(temp[i][j] - heat[i][j]);
            total_heat += heat[i][j];

            max_diff =  diff <= max_diff ? max_diff : diff;

            max_heat =  heat[i][j] <= max_heat ? max_heat : heat[i][j];

            min_heat = heat[i][j] < min_heat ? heat[i][j] : min_heat;
        }
    }

    r->niter = iter_num;
    r->tmin = min_heat;
    r->tmax = max_heat;
    r->maxdiff = max_diff;
    r->tavg = total_heat / (N * (M - 2)); /* average temperature */
    r->time = calc_time(*start_timer, *end_timer);

}

void compute_results(struct results *r,
                     heat_cylinder *hc,
                     struct timespec *start_timer,
                     struct timespec *end_timer,
                     size_t iter_num) {
    __compute_results(r, hc->N, hc->M, hc->heat, hc->temp, start_timer, end_timer, iter_num);
}

int main(int argc, char **argv) {
    // Read input and init datastructures.
    struct parameters p;
    struct results r;
    read_parameters(&p, argc, argv);
    heat_cylinder *hc = init_heat_cylinder(&p);

    // Perform the main loop.
    struct timespec start_timer, end_timer;
    clock_gettime(CLOCK_MONOTONIC, &start_timer);
    for (size_t i = 0; i < p.maxiter; i++) {

        vec_dissipate(hc);

        // Check to see if the results should be printed.
        if ((i % p.period) == 0 && i != 0) {
            clock_gettime(CLOCK_MONOTONIC, &end_timer);
            compute_results(&r, hc, &start_timer, &end_timer, i);
            report_results(&p, &r);
            // Check if the threshold exit condition has been met.
            if (r.maxdiff < p.threshold)
                goto EXIT;

        }
    }

    clock_gettime(CLOCK_MONOTONIC, &end_timer);
    compute_results(&r, hc, &start_timer, &end_timer, p.maxiter);
    report_results(&p, &r);

    EXIT:free_heat_cyclinder(hc);
	return 0;
}
