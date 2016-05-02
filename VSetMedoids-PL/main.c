#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define FPOINT_OFFSET 1e-10

bool debug;
size_t max_iter;
int objc;
size_t clustc;
double mfuz;
double mfuzval;
double ***dmatrix;
int dmatrixc;
size_t ***medoids;
size_t medoids_card;
double **weights;
double **memb;
double epsilon;
double theta;
double *parc_obj_adeq;
double *parc_cluster_adeq;
double prev_adeq;

bool load_data(char *fname, double **matrix) {
	FILE *ifile = fopen(fname, "r");
	if(!ifile) {
		return false;
	}
	size_t i;
	size_t j;
	for(i = 0; i < objc; ++i) {
		for(j = 0; j < objc; ++j) {
			if(fscanf(ifile, "%lf", &matrix[i][j]) == EOF) {
				fclose(ifile);
				return false;
			}
		}
	}
	fclose(ifile);
	return true;
}

void print_mtx_d(double **matrix, size_t nrow, size_t ncol, bool lnum) {
	size_t i;
	size_t j;
	size_t last = ncol - 1;
	for(i = 0; i < nrow; ++i) {
	    if(lnum) {
	        printf("%d: ", i + 1);
	    }
		for(j = 0; j < last; ++j) {
			printf("%lf ", matrix[i][j]);
		}
		printf("%lf\n", matrix[i][j]);
	}
}

void fprint_mtx_d(FILE *file, double **matrix, size_t nrow,
					size_t ncol) {
	size_t i;
	size_t j;
	size_t last = ncol - 1;
	for(i = 0; i < nrow; ++i) {
		for(j = 0; j < last; ++j) {
			fprintf(file, "%.4lf ", matrix[i][j]);
		}
		fprintf(file, "%.4lf\n", matrix[i][j]);
	}
}

void print_mtx_size_t(size_t **matrix, size_t nrow, size_t ncol, bool lnum) {
	size_t i;
	size_t j;
	size_t last = ncol - 1;
	for(i = 0; i < nrow; ++i) {
	    if(lnum) {
	        printf("%d: ", i + 1);
	    }
		for(j = 0; j < last; ++j) {
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

void mtxcpy(void **destination, const void **source, size_t nrow,
        size_t num) {
    size_t i;
    for(i = 0; i < nrow; ++i) {
        memcpy(destination[i], source[i], num);
    }
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

void print_weights() {
	printf("Weights:\n");
	size_t j;
	size_t k;
	double prod;
	for(k = 0; k < clustc; ++k) {
		prod = 1.0;
		for(j = 0; j < dmatrixc; ++j) {
			if(weights[k][j] < 0.0) {
				printf("!");
			}
			printf("%lf ", weights[k][j]);
			prod *= weights[k][j];
		}
		printf("[%lf]", prod);
		if(!deq(prod, 1.0)) {
			printf(" =/= 1.0?\n");
		} else {
			printf("\n");
		}
	}
}

void init_medoids() {
	size_t i;
    size_t j;
	size_t e;
	size_t k;
	int obj;
	bool chosen[objc];
	for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            for(i = 0; i < objc; ++i) {
                chosen[i] = false;
            }
            for(e = 0; e < medoids_card; ++e) {
                do {
                    obj = rand() % objc;
                } while(chosen[obj]);
                medoids[k][j][e] = obj;
                chosen[obj] = true;
            }
        }
	}
}

void print_medoids() {
    printf("Medoids:\n");
    size_t e;
    size_t j;
    size_t k;
    for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            printf("{ ");
            for(e = 0; e < medoids_card; ++e) {
                printf("%d ", medoids[k][j][e]);
            }
            printf("} ");
        }
        printf("\n");
    }
}

void print_memb() {
	printf("Membership:\n");
	size_t i;
	size_t k;
	double sum;
	for(i = 0; i < objc; ++i) {
		printf("%u: ", i);
		sum = 0.0;
		for(k = 0; k < clustc; ++k) {
			printf("%lf ", memb[i][k]);
			sum += memb[i][k];
		}
		printf("[%lf]", sum);
		if(!deq(sum, 1.0)) {
			printf("*\n");
		} else {
			printf("\n");
		}
	}
}

void update_memb() {
	size_t e;
	size_t h;
	size_t i;
	size_t j;
	size_t k;
	size_t a_card;
	double sums[clustc];
	double sumd;
	double val;
	for(i = 0; i < objc; ++i) {
		a_card = 0;
		for(k = 0; k < clustc; ++k) {
			sums[k] = 0.0;
			for(j = 0; j < dmatrixc; ++j) {
				sumd = 0.0;
				for(e = 0; e < medoids_card; ++e) {
					sumd += dmatrix[j][i][medoids[k][j][e]];
				}
				sums[k] += weights[k][j] * sumd;
			}
			if(deq(sums[k], 0.0)) {
				++a_card;
			}
		}
		if(a_card) {
			printf("Object %u has zero val.\n", i);
			val = 1.0 / ((double) a_card);
			for(k = 0; k < clustc; ++k) {
				if(deq(sums[k], 0.0)) {
					memb[i][k] = val;
				} else {
					memb[i][k] = 0.0;
				}
			}
		} else {
			for(k = 0; k < clustc; ++k) {
				memb[i][k] = 0.0;
				for(h = 0; h < clustc; ++h) {
					memb[i][k] += pow(sums[k] / sums[h], mfuzval);
				}
				memb[i][k] = 1.0 / memb[i][k];
			}
		}
	}
}

double adequacy_cluster(bool check) {
    size_t e;
    size_t i;
    size_t j;
    size_t k;
    double sum;
    double sumd;
	double sumw;
    double adeq = 0.0;
	for(k = 0; k < clustc; ++k) {
		sum = 0.0;
		for(i = 0; i < objc; ++i) {
			sumw = 0.0;
            for(j = 0; j < dmatrixc; ++j) {
                sumd = 0.0;
                for(e = 0; e < medoids_card; ++e) {
                    sumd += dmatrix[j][i][medoids[k][j][e]];
                }
                sumw += weights[k][j] * sumd;
            }
            sum += pow(memb[i][k], mfuz) * sumw;
        }
		if(check) {
			if(dlt(parc_cluster_adeq[k], sum)) {
				printf("Msg: adequacy for cluster %d is greater "
                        "than previous (%.15lf).\n", k,
						sum - parc_cluster_adeq[k]);
			}
		}
		parc_cluster_adeq[k] = sum;
		adeq += sum;
    }
    if(check) {
        if(dlt(prev_adeq, adeq)) {
            printf("Msg: overall adequacy is greater than "
                    "previous (%.15lf).\n", adeq - prev_adeq);
        }
    }
    prev_adeq = adeq;
    return adeq;
}

double adequacy_obj(bool check) {
    size_t e;
    size_t i;
    size_t j;
    size_t k;
    double sum;
    double sumd;
	double sumw;
    double adeq = 0.0;
	for(i = 0; i < objc; ++i) {
		sum = 0.0;
		for(k = 0; k < clustc; ++k) {
			sumw = 0.0;
            for(j = 0; j < dmatrixc; ++j) {
                sumd = 0.0;
                for(e = 0; e < medoids_card; ++e) {
                    sumd += dmatrix[j][i][medoids[k][j][e]];
                }
                sumw += weights[k][j] * sumd;
            }
            sum += pow(memb[i][k], mfuz) * sumw;
        }
		if(check) {
			if(dlt(parc_obj_adeq[i], sum)) {
				printf("Msg: adequacy for object %d is greater "
                        "than previous (%.15lf).\n", i,
						sum - parc_obj_adeq[i]);
			}
		}
		parc_obj_adeq[i] = sum;
		adeq += sum;
    }
    if(check) {
        if(dlt(prev_adeq, adeq)) {
            printf("Msg: overall adequacy is greater than "
                    "previous (%.15lf).\n", adeq - prev_adeq);
        }
    }
    prev_adeq = adeq;
    return adeq;
}

typedef struct objnval {
	size_t obj;
	double val;
} objnval;

static int objnval_cmp(const void *p1, const void *p2) {
	const objnval *a = (const objnval *) p1;
	const objnval *b = (const objnval *) p2;
	return (a->val > b->val) - (a->val < b->val);
}

void update_medoids() {
	size_t h;
	size_t i;
	size_t j;
	size_t k;
	objnval candidates[objc];
	for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            for(h = 0; h < objc; ++h) {
                candidates[h].obj = h;
                candidates[h].val = 0.0;
                for(i = 0; i < objc; ++i) {
                    candidates[h].val += pow(memb[i][k], mfuz) * dmatrix[j][i][h];
                }
            }
            qsort(candidates, objc, sizeof(objnval), objnval_cmp);
            for(h = 0; h < medoids_card; ++h) {
                medoids[k][j][h] = candidates[h].obj;
            }
        }
	}
}

void update_weights() {
	size_t e;
	size_t i;
	size_t j;
	size_t k;
	size_t v_card;
	double sums[dmatrixc];
	double sumd;
	double chi;
	double prod_num;
	double exp;
	for(k = 0; k < clustc; ++k) {
		v_card = 0;
		for(j = 0; j < dmatrixc; ++j) {
			sums[j] = 0.0;
			for(i = 0; i < objc; ++i) {
				sumd = 0.0;
				for(e = 0; e < medoids_card; ++e) {
					sumd += dmatrix[j][i][medoids[k][j][e]];
				}
				sums[j] += pow(memb[i][k], mfuz) * sumd;
			}
			if(sums[j] <= theta) {
				++v_card;
			}
		}
		exp = 1.0 / ((double)(dmatrixc - v_card));
		chi = 1.0;
		prod_num = 1.0;
		for(j = 0; j < dmatrixc; ++j) {
			if(sums[j] <= theta) {
				chi *= pow(weights[k][j], exp);
			} else {
				prod_num *= pow(sums[j], exp);
			}
		}
		prod_num = (1.0 / chi) * prod_num;
		for(j = 0; j < dmatrixc; ++j) {
			if(sums[j] > theta) {
				weights[k][j] = prod_num / sums[j];
			}
		}
	}
}

double run() {
	size_t i;
	size_t j;
	size_t k;
	printf("Initialization.\n");
	init_medoids();
    print_medoids();
	for(k = 0; k < clustc; ++k) {
		for(j = 0; j < dmatrixc; ++j) {
			weights[k][j] = 1.0;
		}
	}
	print_weights();
	update_memb();
    //memb_adequacy(false);
	print_memb();
	double prev_adeq = 0.0;
	double adeq = adequacy_obj(false);
	printf("Adequacy: %.20lf\n", adeq);
    double diff = fabs(adeq - prev_adeq);
	for(i = 1; i <= max_iter && diff > epsilon; ++i) {
        printf("Iteration %d.\n", i);
        prev_adeq = adeq;
		adequacy_cluster(false);
        update_medoids();
		adeq = adequacy_cluster(true);
        print_medoids();
        printf("Adequacy1: %.20lf\n", adeq);
		adequacy_cluster(false);
        update_weights();
		adeq = adequacy_cluster(true);
        print_weights();
        printf("Adequacy2: %.20lf\n", adeq);
		adequacy_obj(false);
        update_memb();
		adeq = adequacy_obj(true);
        print_memb();
        printf("Adequacy: %.20lf\n", adeq);
        if(dgt(adeq, prev_adeq)) {
            printf("Warn: current adequacy is greater than "
                    "previous iteration (%.20lf)\n",
                    adeq - prev_adeq);
        }
        diff = fabs(adeq - prev_adeq);
	}
    printf("Adequacy difference threshold reached (%.20lf).\n",
            diff);
    return adeq;
}

// lenergy - J value, calculated using the medoids used by the
// algorithm
// genergy - T value, calculated using the global medoids computed
// in this function
void global_energy() {
    size_t e;
    size_t h;
    size_t i;
    size_t j;
    size_t k;
    // Defining global medoids
    objnval candidates[objc];    
    size_t global_medoids[dmatrixc][medoids_card];
    for(j = 0; j < dmatrixc; ++j) {
        for(h = 0; h < objc; ++h) {
            candidates[h].obj = h;
            candidates[h].val = 0.0;
            for(k = 0; k < clustc; ++k) {
                for(i = 0; i < objc; ++i) {
                    candidates[h].val += pow(memb[i][k], mfuz) *
                        weights[k][j] * dmatrix[j][i][h];
                }
            }
        }
        qsort(candidates, objc, sizeof(objnval), objnval_cmp);
        for(h = 0; h < medoids_card; ++h) {
            global_medoids[j][h] = candidates[h].obj;
        }
    }
    printf("Global medoids:\n");
    for(j = 0; j < dmatrixc; ++j) {
        printf("%d:", j + 1);
        for(e = 0; e < medoids_card; ++e) {
            printf(" %d", global_medoids[j][e]);
        }
        printf("\n");
    }
    printf("\n");
    // Global energy per cluster and matrix
    double genergy[clustc][dmatrixc];
    double lenergy[clustc][dmatrixc];
    double sumd_global;
    double sumd_local;
    double val;
    double mtx_genergy[dmatrixc];
    double mtx_lenergy[dmatrixc];
    for(j = 0; j < dmatrixc; ++j) {
        mtx_genergy[j] = 0.0;
        mtx_lenergy[j] = 0.0;
    }
    double clust_genergy[clustc];
    double clust_lenergy[clustc];
    for(k = 0; k < clustc; ++k) {
        clust_genergy[k] = 0.0;
        clust_lenergy[k] = 0.0;
    }
    double total_genergy = 0.0;
    double total_lenergy = 0.0;
    for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            genergy[k][j] = 0.0;
            lenergy[k][j] = 0.0;
            for(i = 0; i < objc; ++i) {
                sumd_global = 0.0;
                sumd_local = 0.0;
                for(e = 0; e < medoids_card; ++e) {
                    sumd_global +=
                        dmatrix[j][i][global_medoids[j][e]];
                    sumd_local += dmatrix[j][i][medoids[k][j][e]];
                }
                val = pow(memb[i][k], mfuz) * weights[k][j];
                genergy[k][j] += val * sumd_global;
                lenergy[k][j] += val * sumd_local;
            }
            clust_genergy[k] += genergy[k][j];
            clust_lenergy[k] += lenergy[k][j];
            mtx_genergy[j] += genergy[k][j];
            mtx_lenergy[j] += lenergy[k][j];
            total_genergy += genergy[k][j];
            total_lenergy += lenergy[k][j];
        }
    }
    printf("Global energy matrix:\n");
    for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            printf("%lf ", genergy[k][j]);
        }
        printf("\n");
    }
    printf("Local energy matrix:\n");
    for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            printf("%lf ", lenergy[k][j]);
        }
        printf("\n");
    }
    for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            if(dgt(lenergy[k][j], genergy[k][j])) {
                printf("Msg: lenergy > genergy for cluster %d and matrix %d.\n", k, j);
            }
        }
    }
    printf("\n");
    double sum = 0.0;
    printf("Cluster global energy:\n");
    for(k = 0; k < clustc; ++k) {
        sum += clust_genergy[k];
        printf("%lf ", clust_genergy[k]);
    }
    printf("[%lf]\n", sum);
    sum = 0.0;
    printf("Cluster local energy:\n");
    for(k = 0; k < clustc; ++k) {
        sum += clust_lenergy[k];
        printf("%lf ", clust_lenergy[k]);
    }
    printf("[%lf]\n", sum);
    for(k = 0; k < clustc; ++k) {
        if(dgt(clust_lenergy[k], clust_genergy[k])) {
            printf("Warn: clust_lenergy > clust_genergy for cluster %d.\n", k);
        }
    }
    printf("\n");
    sum = 0.0;
    printf("Matrix global energy:\n");
    for(j = 0; j < dmatrixc; ++j) {
        sum += mtx_genergy[j];
        printf("%lf ", mtx_genergy[j]);
    }
    printf("[%lf]\n", sum);
    sum = 0.0;
    printf("Matrix local energy:\n");
    for(j = 0; j < dmatrixc; ++j) {
        sum += mtx_lenergy[j];
        printf("%lf ", mtx_lenergy[j]);
    }
    printf("[%lf]\n", sum);
    for(j = 0; j < dmatrixc; ++j) {
        if(dgt(mtx_lenergy[j], mtx_genergy[j])) {
            printf("Warn: mtx_lenergy > mtx_genergy for matrix %d.\n", j);
        }
    }
    printf("\n");
    printf("Global energy: %lf\n", total_genergy);
    printf("Local energy: %lf\n", total_lenergy);
    if(dgt(total_lenergy, total_genergy)) {
        printf("Warn: total_lenergy > total_genergy.\n");
    }
    printf("\n");
    printf("Matrix global heterogeneity index:\n");
    for(j = 0; j < dmatrixc; ++j) {
        printf("%lf ", (1.0 - (mtx_lenergy[j] / mtx_genergy[j])));
    }
    printf("\n");
    printf("\n");
    printf("Global heterogeneity index: %lf\n",
            (1.0 - (total_lenergy / total_genergy)));
    printf("\n");
    printf("Cluster heterogeinety indexes:\n");
    printf("T:\n");
    for(k = 0; k < clustc; ++k) {
        printf("%lf ", clust_genergy[k] / total_genergy);
    }
    printf("\n");
    printf("J:\n");
    for(k = 0; k < clustc; ++k) {
        printf("%lf ", clust_lenergy[k] / total_lenergy);
    }
    printf("\n");
    printf("B:\n");
    for(k = 0; k < clustc; ++k) {
        printf("%lf ", (clust_genergy[k] - clust_lenergy[k]) /
                    (total_genergy - total_lenergy));
    }
    printf("\n");
    printf("Q:\n");
    for(k = 0; k < clustc; ++k) {
        printf("%lf ", (1.0 - (clust_lenergy[k] / clust_genergy[k])));
    }
    printf("\n");
    printf("\n");
    printf("Cluster heterogeinety indexes for matrices:\n");
    for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            printf("%lf ", (1.0 - (lenergy[k][j] / genergy[k][j])));
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {
	debug = true;
	size_t argpos = 1;
    int val;
	clustc = 3;
	medoids_card = 3;
	max_iter = 100;
	epsilon = 0.000001;
	theta = 0.02;
	mfuz = 2;
	int insts = 10;
    for(; argpos < argc; ++argpos) {
        if(!strcmp(argv[argpos], "-k")) {
            val = atoi(argv[++argpos]);
            if(val <= 0) {
                printf("Error: k <= 0.\n");
                return 1;
            }
            clustc = val;
        } else if(!strcmp(argv[argpos], "-q")) {
            val = atoi(argv[++argpos]);
            if(val <= 0) {
                printf("Error: q <= 0.\n");
                return 1;
            }
            medoids_card = val;
        } else if(!strcmp(argv[argpos], "-T")) {
            val = atoi(argv[++argpos]);
            if(val <= 0) {
                printf("Error: T <= 0.\n");
                return 1;
            }
            max_iter = val;
        } else if(!strcmp(argv[argpos], "-e")) {
            epsilon = atof(argv[++argpos]);
            if(epsilon < 0) {
                printf("Error: e <= 0.\n");
                return 1;
            }
        } else if(!strcmp(argv[argpos], "-t")) {
            theta = atof(argv[++argpos]);
            if(theta < 0) {
                printf("Error: t <= 0.\n");
                return 1;
            }
        } else if(!strcmp(argv[argpos], "-m")) {
            val = atoi(argv[++argpos]);
            if(val <= 0) {
                printf("Error: m <= 0.\n");
                return 1;
            }
            mfuz = val;
        } else if(!strcmp(argv[argpos], "-i")) {
            val = atoi(argv[++argpos]);
            if(val <= 0) {
                printf("Error: i <= 0.\n");
                return 1;
            }
            insts = val;
        } else {
            break;
        }
    }
    if((argc - argpos) < 2) {
        printf("Error: not enough args.\n");
        return 1;
    }
	mfuzval = 1.0 / (mfuz - 1.0);
	objc = atoi(argv[argpos++]);
	if(objc <= 0) {
		printf("Error: objc <= 0.\n");
		return 1;
	}
	dmatrixc = atoi(argv[argpos++]);
	if(dmatrixc <= 0) {
		printf("Error: dmatrixc <= 0.\n");
		return 1;
	}
    if((argc - argpos) < dmatrixc) {
        printf("Error: not enough args.\n");
        return 1;
    }
    printf("######Config summary:######\n");
    printf("Number of clusters: %d.\n", clustc);
    printf("Medoids cardinality: %d.\n", medoids_card);
    printf("Number of iterations: %d.\n", max_iter);
    printf("Epsilon: %.15lf.\n", epsilon);
    printf("Theta: %.15lf.\n", theta);
    printf("Parameter m: %.15lf.\n", mfuz);
    printf("Number of instances: %d.\n", insts);
    printf("###########################\n");
	size_t i;
	size_t j;
	size_t k;
	// Allocating memory start
	parc_cluster_adeq = malloc(sizeof(double) * clustc);
    parc_obj_adeq = malloc(sizeof(double) * objc);
	dmatrix = malloc(sizeof(double **) * dmatrixc);
	for(j = 0; j < dmatrixc; ++j) {
		dmatrix[j] = malloc(sizeof(double *) * objc);
		for(i = 0; i < objc; ++i) {
			dmatrix[j][i] = malloc(sizeof(double) * objc);
		}
	}
	medoids = malloc(sizeof(size_t **) * clustc);
	size_t ***best_medoids = malloc(sizeof(size_t **) * clustc);
	for(k = 0; k < clustc; ++k) {
		medoids[k] = malloc(sizeof(size_t *) * dmatrixc);
		best_medoids[k] = malloc(sizeof(size_t *) * dmatrixc);
        for(j = 0; j < dmatrixc; ++j) {
            medoids[k][j] = malloc(sizeof(size_t) * medoids_card);
            best_medoids[k][j] = malloc(sizeof(size_t) * medoids_card);
        }
	}
	weights = malloc(sizeof(double *) * clustc);
	double **best_weights = malloc(sizeof(double *) * clustc);
	for(k = 0; k < clustc; ++k) {
		weights[k] = malloc(sizeof(double) * dmatrixc);
		best_weights[k] = malloc(sizeof(double) * dmatrixc);
	}
	memb = malloc(sizeof(double *) * objc);
	double **best_memb = malloc(sizeof(double *) * objc);
	for(i = 0; i < objc; ++i) {
		memb[i] = malloc(sizeof(double) * clustc);
		best_memb[i] = malloc(sizeof(double) * clustc);
	}
	// Allocating memory end
	for(j = 0; j < dmatrixc; ++j, ++argpos) {
		if(!load_data(argv[argpos], dmatrix[j])) {
			printf("Error: could not load %s.\n", argv[argpos]);
			goto END;
		}
	}
    size_t best_inst;
    double best_inst_adeq;
    double cur_inst_adeq;
	srand(time(NULL));
	for(i = 1; i <= insts; ++i) {
		printf("Instance %u:\n", i);
		cur_inst_adeq = run();
		printf("\n");
		global_energy();
		printf("\n");
        if(i == 1 || cur_inst_adeq < best_inst_adeq) {
            mtxcpy_d(best_memb, memb, objc, clustc);
            mtxcpy_d(best_weights, weights, clustc, dmatrixc);
            for(k = 0; k < clustc; ++k) {
                mtxcpy_size_t(best_medoids[k], medoids[k], dmatrixc,
                        medoids_card);
            }
            best_inst_adeq = cur_inst_adeq;
            best_inst = i;
        }
	}
	printf("\n");
    printf("Best adequacy %.15lf on instance %d.\n",
            best_inst_adeq, best_inst);
    printf("\n");
    size_t ***swp2 = medoids;
    medoids = best_medoids;
    best_medoids = swp2;
    print_medoids();
    printf("\n");
    double **swp = memb;
	memb = best_memb;
	best_memb = swp;
	print_memb();
	printf("\n");
	swp = weights;
	weights = best_weights;
	best_weights = swp;
	print_weights();
	printf("\n");
    global_energy();
END:
	for(i = 0; i < dmatrixc; ++i) {
		for(j = 0; j < objc; ++j) {
			free(dmatrix[i][j]);
		}
		free(dmatrix[i]);
	}
	free(dmatrix);
	for(k = 0; k < clustc; ++k) {
        for(j = 0; j < dmatrixc; ++j) {
            free(medoids[k][j]);
            free(best_medoids[k][j]);
        }
		free(medoids[k]);
		free(best_medoids[k]);
		free(weights[k]);
		free(best_weights[k]);
	}
	free(medoids);
	free(best_medoids);
	free(weights);
	free(best_weights);
	for(i = 0; i < objc; ++i) {
		free(memb[i]);
		free(best_memb[i]);
	}
	free(memb);
	free(best_memb);
	free(parc_cluster_adeq);
    free(parc_obj_adeq);
	return 0;
}
