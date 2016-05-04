#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "util.h"

bool load_data(char *fname, double **matrix, size_t nrow,
        size_t ncol) {
	FILE *ifile = fopen(fname, "r");
	if(!ifile) {
		return false;
	}
	size_t i;
	size_t j;
	for(i = 0; i < nrow; ++i) {
		for(j = 0; j < ncol; ++j) {
			if(fscanf(ifile, "%lf", &matrix[i][j]) == EOF) {
				fclose(ifile);
				return false;
			}
		}
	}
	fclose(ifile);
	return true;
}

void print_mtx_d(double **matrix, size_t nrow, size_t ncol) {
	size_t i;
	size_t j;
	for(i = 0; i < nrow; ++i) {
		for(j = 0; j < ncol - 1; ++j) {
			printf("%lf ", matrix[i][j]);
		}
		printf("%lf\n", matrix[i][j]);
	}
}

void fprint_mtx_d(FILE *file, double **matrix, size_t nrow,
					size_t ncol) {
	size_t i;
	size_t j;
	for(i = 0; i < nrow; ++i) {
		for(j = 0; j < ncol - 1; ++j) {
			fprintf(file, "%.4lf ", matrix[i][j]);
		}
		fprintf(file, "%.4lf\n", matrix[i][j]);
	}
}

void print_mtx_size_t(size_t **matrix, size_t nrow, size_t ncol) {
	size_t i;
	size_t j;
	for(i = 0; i < nrow; ++i) {
		for(j = 0; j < ncol - 1; ++j) {
			printf("%u ", matrix[i][j]);
		}
		printf("%u\n", matrix[i][j]);
	}
}

bool deq(double a, double b) {
    return (a < (b + FPOINT_OFFSET) && a > (b - FPOINT_OFFSET));
}

bool dgt(double a, double b) {
    return a > (b + FPOINT_OFFSET);
}

bool dlt(double a, double b) {
    return a < (b - FPOINT_OFFSET);
}

void mtxcpy_d(double **destination, double **source, size_t nrow,
        size_t ncol) {
    size_t i;
    for(i = 0; i < nrow; ++i) {
        memcpy(destination[i], source[i], sizeof(double) * ncol);
    }
}

void mtxcpy_size_t(size_t **destination, size_t **source, size_t nrow,
        size_t ncol) {
    size_t i;
    for(i = 0; i < nrow; ++i) {
        memcpy(destination[i], source[i], sizeof(size_t) * ncol);
    }
}
