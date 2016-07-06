#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "util.h"

#define BUFF_SIZE 1024

bool debug;
size_t max_iter;
int objc;
size_t clustc;
double mfuz;
double mfuzval;
double ***dmatrix;
int dmatrixc;
size_t **medoids;
size_t medoids_card;
double **weights;
double **memb;
double epsilon;
double theta;
double *parc_obj_adeq;
double *parc_cluster_adeq;
double prev_adeq;

void print_weights(double **weights) {
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
	size_t e;
	size_t k;
	int obj;
	bool chosen[objc];
	for(k = 0; k < clustc; ++k) {
		for(i = 0; i < objc; ++i) {
			chosen[i] = false;
		}
		for(e = 0; e < medoids_card; ++e) {
			do {
				obj = rand() % objc;
			} while(chosen[obj]);
			medoids[k][e] = obj;
			chosen[obj] = true;
		}
	}
}

void print_medoids(size_t **medoids) {
    printf("Medoids:\n");
    size_t e;
    size_t k;
    for(k = 0; k < clustc; ++k) {
        for(e = 0; e < medoids_card; ++e) {
            printf("%d ", medoids[k][e]);
        }
        printf("\n");
    }
}

void print_memb(double **memb) {
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
					sumd += dmatrix[j][i][medoids[k][e]];
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
                    sumd += dmatrix[j][i][medoids[k][e]];
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
                    sumd += dmatrix[j][i][medoids[k][e]];
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
	double sumweights;
	for(k = 0; k < clustc; ++k) {
		for(h = 0; h < objc; ++h) {
			candidates[h].obj = h;
			candidates[h].val = 0.0;
			for(i = 0; i < objc; ++i) {
				sumweights = 0.0;
				for(j = 0; j < dmatrixc; ++j) {
					sumweights += weights[k][j] * dmatrix[j][i][h];
				}
				candidates[h].val += pow(memb[i][k], mfuz) * sumweights;
			}
		}
		qsort(candidates, objc, sizeof(objnval), objnval_cmp);
        for(h = 0; h < medoids_card; ++h) {
            medoids[k][h] = candidates[h].obj;
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
					sumd += dmatrix[j][i][medoids[k][e]];
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
	print_medoids(medoids);
	for(k = 0; k < clustc; ++k) {
		for(j = 0; j < dmatrixc; ++j) {
			weights[k][j] = 1.0;
		}
	}
	print_weights(weights);
	update_memb();
    //memb_adequacy(false);
	print_memb(memb);
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
        print_medoids(medoids);
        printf("Adequacy1: %.20lf\n", adeq);
		adequacy_cluster(false);
        update_weights();
		adeq = adequacy_cluster(true);
        print_weights(weights);
        printf("Adequacy2: %.20lf\n", adeq);
		adequacy_obj(false);
        update_memb();
		adeq = adequacy_obj(true);
        print_memb(memb);
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

int main(int argc, char **argv) {
	debug = true;
	int insts;
    FILE *cfgfile = fopen(argv[1], "r");
    if(!cfgfile) {
        printf("Error: could not open config file.\n");
        return 1;
    }
    fscanf(cfgfile, "%d", &objc);
    if(objc <= 0) {
        printf("Error: objc <= 0.\n");
        return 2;
    }
    // ignore labels
    fscanf(cfgfile, "%*d");
	size_t i;
    for(i = 0; i < objc; ++i) {
        fscanf(cfgfile, "%*d");
    }
    // ignore labels end
    fscanf(cfgfile, "%d", &dmatrixc);
    if(dmatrixc <= 0) {
        printf("Error: dmatrixc <= 0.\n");
        return 2;
    }
    char dmtx_file_name[dmatrixc][BUFF_SIZE];
	size_t j;
    for(j = 0; j < dmatrixc; ++j) {
        fscanf(cfgfile, "%s", dmtx_file_name[j]);
    }
    char out_file_name[BUFF_SIZE];
    fscanf(cfgfile, "%s", out_file_name);
    fscanf(cfgfile, "%d", &clustc);
    if(clustc <= 0) {
        printf("Error: clustc <= 0.\n");
        return 2;
    }
    fscanf(cfgfile, "%d", &medoids_card);
    if(medoids_card <= 0) {
        printf("Error: medoids_card <= 0.\n");
        return 2;
    }
    fscanf(cfgfile, "%d", &insts);
    if(insts <= 0) {
        printf("Error: insts <= 0.\n");
        return 2;
    }
    fscanf(cfgfile, "%lf", &theta);
    if(dlt(theta, 0.0)) {
        printf("Error: theta < 0.\n");
        return 2;
    }
    fscanf(cfgfile, "%d", &max_iter);
    fscanf(cfgfile, "%lf", &epsilon);
    if(dlt(epsilon, 0.0)) {
        printf("Error: epsilon < 0.\n");
        return 2;
    }
    fscanf(cfgfile, "%lf", &mfuz);
    if(!dgt(mfuz, 0.0)) {
        printf("Error: mfuz <= 0.\n");
        return 2;
    }
    fclose(cfgfile);
    freopen(out_file_name, "w", stdout);
	mfuzval = 1.0 / (mfuz - 1.0);
    printf("######Config summary:######\n");
    printf("Number of objects: %d.\n", objc);
    printf("Number of clusters: %d.\n", clustc);
    printf("Medoids cardinality: %d.\n", medoids_card);
    printf("Number of iterations: %d.\n", max_iter);
    printf("Epsilon: %.15lf.\n", epsilon);
    printf("Theta: %.15lf.\n", theta);
    printf("Parameter m: %.15lf.\n", mfuz);
    printf("Number of instances: %d.\n", insts);
    printf("###########################\n");
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
	medoids = malloc(sizeof(size_t *) * clustc);
	size_t **best_medoids = malloc(sizeof(size_t *) * clustc);
	for(k = 0; k < clustc; ++k) {
		medoids[k] = malloc(sizeof(size_t) * medoids_card);
		best_medoids[k] = malloc(sizeof(size_t) * medoids_card);
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
	for(j = 0; j < dmatrixc; ++j) {
		if(!load_data(dmtx_file_name[j], dmatrix[j], objc, objc)) {
			printf("Error: could not load %s.\n", dmtx_file_name[j]);
			goto END;
		}
	}
	srand(time(NULL));
    size_t best_inst;
    double best_inst_adeq;
    double cur_inst_adeq;
	for(i = 1; i <= insts; ++i) {
		printf("Instance %u:\n", i);
		cur_inst_adeq = run();
        if(i == 1 || cur_inst_adeq < best_inst_adeq) {
            mtxcpy_d(best_memb, memb, objc, clustc);
            mtxcpy_d(best_weights, weights, clustc, dmatrixc);
            mtxcpy_size_t(best_medoids, medoids, clustc,
                    medoids_card);
            best_inst_adeq = cur_inst_adeq;
            best_inst = i;
        }
	}
	printf("\n");
    printf("Best adequacy %.15lf on instance %d.\n",
            best_inst_adeq, best_inst);
    printf("\n");
    print_medoids(best_medoids);
    printf("\n");
	print_memb(best_memb);
	printf("\n");
	print_weights(best_weights);
END:
    fclose(stdout);
	for(i = 0; i < dmatrixc; ++i) {
		for(j = 0; j < objc; ++j) {
			free(dmatrix[i][j]);
		}
		free(dmatrix[i]);
	}
	free(dmatrix);
	for(k = 0; k < clustc; ++k) {
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
