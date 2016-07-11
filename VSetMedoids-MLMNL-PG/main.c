#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "util.h"
#include "matrix.h"
#include "stex.h"

#define BUFF_SIZE 1024
#define HEADER_SIZE 51

typedef struct int_vec {
	int *get;
	size_t size;
} int_vec;

void int_vec_init(int_vec *v, size_t size) {
	v->get = malloc(sizeof(int) * size);
	v->size = 0;
}

void int_vec_free(int_vec *v) {
	free(v->get);
}

void int_vec_push(int_vec *v, int val) {
	v->get[v->size++] = val;
}

typedef struct constraint {
	int_vec *ml;
	int_vec *mnl;
} constraint;

void constraint_init(constraint *c, size_t mlsize, size_t mnlsize) {
//	c = malloc(sizeof(constraint));
    c->ml = malloc(sizeof(int_vec));
	int_vec_init(c->ml, mlsize);
    c->mnl = malloc(sizeof(int_vec));
	int_vec_init(c->mnl, mnlsize);
}

void constraint_free(constraint *c) {
	int_vec_free(c->ml);
    free(c->ml);
	int_vec_free(c->mnl);
    free(c->mnl);
	free(c);
}

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
double *weights;
double **memb;
double epsilon;
double theta;
double alpha;
double *parc_obj_adeq;
double *parc_cluster_adeq;
double prev_adeq;
double prev_rest_adeq;
double prev_free_adeq;
int_vec *class;
int classc;
int_vec *sample;
constraint **constraints;
size_t constsc;

void print_weights(double *weights) {
	printf("Weights:\n");
	size_t j;
	double prod = 1.0;
    for(j = 0; j < dmatrixc; ++j) {
        if(weights[j] < 0.0) {
            printf("!");
        }
        printf("%lf ", weights[j]);
        prod *= weights[j];
    }
    printf("[%lf]", prod);
    if(!deq(prod, 1.0)) {
        printf(" =/= 1.0?\n");
    } else {
        printf("\n");
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

void print_medoids(size_t ***medoids) {
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

void constrained_update_memb() {
	size_t e;
	size_t h;
	size_t i;
	size_t j;
	size_t k;
	size_t a_card[objc];
	double sumd;
	double val;
    double mtx_a[objc][clustc];
    double mtx_b[objc][clustc];
    bool mtx_v[clustc][clustc];
	for(i = 0; i < objc; ++i) {
		a_card[i] = 0;
		for(k = 0; k < clustc; ++k) {
			mtx_a[i][k] = 0.0;
			for(j = 0; j < dmatrixc; ++j) {
				sumd = 0.0;
				for(e = 0; e < medoids_card; ++e) {
					sumd += dmatrix[j][i][medoids[k][j][e]];
				}
				mtx_a[i][k] += weights[j] * sumd;
			}
			if(deq(mtx_a[i][k], 0.0)) {
				a_card[i] += 1;
			}
		}
	}
    size_t m;
    size_t obj;
    for(i = 0; i < objc; ++i) {
        if(!a_card[i]) {
            for(k = 0; k < clustc; ++k) {
                mtx_v[i][k] = true;
                mtx_a[i][k] *= 2.0;
                mtx_b[i][k] = 0.0;
                if(constraints[i]) {
                    for(m = 0; m < constraints[i]->ml->size; ++m) {
                        obj = constraints[i]->ml->get[m];
                        for(h = 0; h < clustc; ++h) {
                            if(h != k) {
                                mtx_b[i][k] += memb[obj][h];
                            }
                        }
                    }
                    for(m = 0; m < constraints[i]->mnl->size; ++m) {
//                        obj = constraints[i]->mnl->get[m];
                        mtx_b[i][k] += memb[constraints[i]->mnl->get[m]][k];
                    }
                    mtx_b[i][k] = alpha * (2.0 * mtx_b[i][k]);
                }
            }
        }
    }
    double sum_num;
    double sum_den;
    double gamma;
    bool test;
    do {
        test = false;
        for(i = 0; i < objc; ++i) {
            if(!a_card[i]) {
                sum_num = 0.0;
                sum_den = 0.0;
                for(k = 0; k < clustc; ++k) {
                    if(mtx_v[i][k]) {
                        sum_num += mtx_b[i][k] / mtx_a[i][k];
                        sum_den += 1.0 / mtx_a[i][k];
                    }
                }
                gamma = (1.0 + sum_num) / sum_den;
                for(k = 0; k < clustc; ++k) {
                    if(mtx_v[i][k]) {
                        memb[i][k] = (gamma - mtx_b[i][k]) / mtx_a[i][k];
                        if(memb[i][k] <= 0.0) {
                            memb[i][k] = 0.0;
                            mtx_v[i][k] = false;
                            test = true;
                        }
                    }
                }
            }
        }
    } while(test);
    for(i = 0; i < objc; ++i) {
        if(a_card[i]) {
            val = 1.0 / a_card[i];
            for(k = 0; k < clustc; ++k) {
                if(!deq(mtx_a[i][k], 0.0)) {
                    memb[i][k] = 0.0;
                } else {
                    memb[i][k] = val;
                }
            }
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
				sums[k] += weights[j] * sumd;
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
    size_t m;
    size_t r;
    size_t s;
	double sumweights;
	double sumd;
	double free_adeq = 0.0; // overall free adeq
    size_t obj;
    double rest_adeq1 = 0.0; // ML overall rest adeq
    double rest_adeq2 = 0.0; // MNL overall rest adeq
    double cluster_adeq; // overall adeq for a cluster
    for(k = 0; k < clustc; ++k) {
        cluster_adeq = 0.0;
        for(i = 0; i < objc; ++i) {
			sumweights = 0.0;
			for(j = 0; j < dmatrixc; ++j) {
				sumd = 0.0;
				for(e = 0; e < medoids_card; ++e) {
					sumd += dmatrix[j][i][medoids[k][j][e]];
				}
				sumweights += weights[j] * sumd;
			}
            cluster_adeq += pow(memb[i][k], mfuz) * sumweights;
            free_adeq += pow(memb[i][k], mfuz) * sumweights;
            if(constraints[i]) {
                for(m = 0; m < constraints[i]->ml->size; ++m) {
                    obj = constraints[i]->ml->get[m];
                    for(r = 0; r < clustc; ++r) {
                        for(s = 0; s < clustc; ++s) {
                            if(s != r) {
                                cluster_adeq +=
                                            memb[i][r] * memb[obj][s];
                                rest_adeq1 +=
                                            memb[i][r] * memb[obj][s];
                            }
                        }
                    }
                }
                for(m = 0; m < constraints[i]->mnl->size; ++m) {
                    obj = constraints[i]->mnl->get[m];
                    for(r = 0; r < clustc; ++r) {
                        cluster_adeq += memb[i][r] * memb[obj][r];
                        rest_adeq2 += memb[i][r] * memb[obj][r];
                    }
                }
            }
		}
        if(check) {
            if(dlt(parc_cluster_adeq[k], cluster_adeq)) {
                printf("Msg: adequacy for cluster %d is greater "
                        "than previous (%.15lf).\n", k,
						cluster_adeq - parc_cluster_adeq[k]);
            }
        }
        parc_cluster_adeq[k] = cluster_adeq;
	}
    double rest_adeq = alpha * (rest_adeq1 + rest_adeq2);
//    printf("constraint adeq: %.20lf\n", alpha * (sum1 + sum2));
//    printf("adeq: %.20lf\n", adeq);
    double adeq = free_adeq + rest_adeq;
    if(check) {
        if(dlt(prev_free_adeq, free_adeq)) {
            printf("Msg: unconstrained adequacy is greater than "
                    "previous (%.15lf).\n",
                    free_adeq - prev_free_adeq);
        }
        if(dlt(prev_rest_adeq, rest_adeq)) {
            printf("Msg: constrained adequacy is greater than "
                    "previous (%.15lf).\n",
                    rest_adeq - prev_rest_adeq);
        }
        if(dlt(prev_adeq, adeq)) {
            printf("Msg: overall adequacy is greater than "
                    "previous (%.15lf).\n", adeq - prev_adeq);
        }
    }
    prev_rest_adeq = rest_adeq;
    prev_free_adeq = free_adeq;
    prev_adeq = adeq;
    return adeq;
}

double adequacy_obj(bool check) {
	size_t e;
	size_t i;
	size_t j;
	size_t k;
    size_t m;
    size_t r;
    size_t s;
	double sumweights;
	double sumd;
	double free_adeq = 0.0;
    size_t obj;
    double rest_adeq1 = 0.0;
    double rest_adeq2 = 0.0;
    double obj_adeq;
    for(i = 0; i < objc; ++i) {
        obj_adeq = 0.0;
        for(k = 0; k < clustc; ++k) {
			sumweights = 0.0;
			for(j = 0; j < dmatrixc; ++j) {
				sumd = 0.0;
				for(e = 0; e < medoids_card; ++e) {
					sumd += dmatrix[j][i][medoids[k][j][e]];
				}
				sumweights += weights[j] * sumd;
			}
            obj_adeq += pow(memb[i][k], mfuz) * sumweights;
            free_adeq += pow(memb[i][k], mfuz) * sumweights;
		}
//        free_adeq += obj_adeq;
        if(constraints[i]) {
            for(m = 0; m < constraints[i]->ml->size; ++m) {
                obj = constraints[i]->ml->get[m];
                for(r = 0; r < clustc; ++r) {
                    for(s = 0; s < clustc; ++s) {
                        if(s != r) {
                            obj_adeq += memb[i][r] * memb[obj][s];
                            rest_adeq1 += memb[i][r] * memb[obj][s];
                        }
                    }
                }
            }
            for(m = 0; m < constraints[i]->mnl->size; ++m) {
                obj = constraints[i]->mnl->get[m];
                for(r = 0; r < clustc; ++r) {
                    obj_adeq += memb[i][r] * memb[obj][r];
                    rest_adeq2 += memb[i][r] * memb[obj][r];
                }
            }
        }
        if(check) {
            if(dlt(parc_obj_adeq[i], obj_adeq)) {
                if(constraints[i]) {
                    printf("Msg: adequacy for constrained object %d "
                            "is greater than previous (%.15lf).\n", i,
    						obj_adeq - parc_obj_adeq[i]);
                } else {
                    printf("Msg: adequacy for unconstrained object %d"
                            " is greater than previous (%.15lf).\n",
                            i, obj_adeq - parc_obj_adeq[i]);
			    }
            }
        }
        parc_obj_adeq[i] = obj_adeq;
	}
    double rest_adeq = alpha * (rest_adeq1 + rest_adeq2);
//    printf("constraint adeq: %.20lf\n", alpha * (sum1 + sum2));
//    printf("adeq: %.20lf\n", adeq);
    double adeq = free_adeq + rest_adeq;
    if(check) {
        if(dlt(prev_free_adeq, free_adeq)) {
            printf("Msg: unconstrained adequacy is greater than "
                    "previous (%.15lf).\n",
                    free_adeq - prev_free_adeq);
        }
        if(dlt(prev_rest_adeq, rest_adeq)) {
            printf("Msg: constrained adequacy is greater than "
                    "previous (%.15lf).\n",
                    rest_adeq - prev_rest_adeq);
        }
        if(dlt(prev_adeq, adeq)) {
            printf("Msg: overall adequacy is greater than "
                    "previous (%.15lf).\n", adeq - prev_adeq);
        }
    }
    prev_rest_adeq = rest_adeq;
    prev_free_adeq = free_adeq;
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
                    candidates[h].val += pow(memb[i][k], mfuz) *
                                            dmatrix[j][i][h];
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
    v_card = 0;
    for(j = 0; j < dmatrixc; ++j) {
        sums[j] = 0.0;
        for(k = 0; k < clustc; ++k) {
            for(i = 0; i < objc; ++i) {
                sumd = 0.0;
                for(e = 0; e < medoids_card; ++e) {
                    sumd += dmatrix[j][i][medoids[k][j][e]];
                }
                sums[j] += pow(memb[i][k], mfuz) * sumd;
            }
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
            chi *= pow(weights[j], exp);
        } else {
            prod_num *= pow(sums[j], exp);
        }
    }
    prod_num = (1.0 / chi) * prod_num;
    for(j = 0; j < dmatrixc; ++j) {
        if(sums[j] > theta) {
            weights[j] = prod_num / sums[j];
        }
    }
}

double run() {
	size_t i;
	size_t j;
	printf("Initialization.\n");
	init_medoids();
    print_medoids(medoids);
    for(j = 0; j < dmatrixc; ++j) {
        weights[j] = 1.0;
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
        constrained_update_memb();
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

bool load_class_file(char *filename) {
	bool ret = true;
	FILE *ifile = fopen(filename, "r");
	if(ifile) {
		fscanf(ifile, "%d", &classc);
		if(classc) {
			size_t i;
			class = malloc(sizeof(int_vec) * classc);
			for(i = 0; i < classc; ++i) {
				int_vec_init(&class[i], objc);
			}
			int obj_class;
			for(i = 0; i < objc; ++i) {
				fscanf(ifile, "%d", &obj_class);
				if(obj_class < 0 || obj_class >= classc) {
					printf("Error: invalid object class.\n");
					for(i = 0; i < classc; ++i) {
						int_vec_free(&class[i]);
					}
					free(class);
                    ret = false;
					break;
				}
				int_vec_push(&class[obj_class], i);
			}
		} else {
			ret = false;
		}
		fclose(ifile);
	} else {
		ret = false;
	}
	return ret;
}

void print_class() {
	size_t i;
	size_t k;
	for(k = 0; k < classc; ++k) {
		printf("Class %d (%d members):\n", k, class[k].size);
		for(i = 0; i < class[k].size; ++i) {
			printf("%u ", class[k].get[i]);
		}
		printf("\n");
	}
}

void gen_sample(size_t size) {
	printf("sample size: %d\n", size);
	sample = malloc(sizeof(int_vec) * classc);
	size_t per_class = size / classc;
	constsc = per_class * classc;
	int pos;
	size_t i;
	size_t k;
	for(k = 0; k < classc; ++k) {
		int_vec_init(&sample[k], per_class);
		bool chosen[class[k].size];
        for(i = 0; i < class[k].size; ++i) {
            chosen[i] = false;
        }
		for(i = 0; i < per_class; ++i) {
			do {
				pos = rand() % class[k].size;
			} while(chosen[pos]);
			chosen[pos] = true;
			int_vec_push(&sample[k], class[k].get[pos]);
		}
	}
}

void print_sample() {
	size_t i;
	size_t k;
	for(k = 0; k < classc; ++k) {
		printf("Sample %d (%d members):\n", k, sample[k].size);
		for(i = 0; i < sample[k].size; ++i) {
			printf("%u ", sample[k].get[i]);
		}
		printf("\n");
	}
}

void gen_constraints() {
	constraints = calloc(objc, sizeof(constraint *));
	size_t i;
	size_t e;
	size_t h;
	size_t k;
	size_t obj;
	size_t obj2;
	for(k = 0; k < classc; ++k) {
		for(i = 0; i < sample[k].size; ++i) {
			obj = sample[k].get[i];
            constraints[obj] = malloc(sizeof(constraint));
			constraint_init(constraints[obj], sample[k].size,
							constsc - sample[k].size);
			for(h = 0; h < classc; ++h) {
				for(e = 0; e < sample[h].size; ++e) {
					obj2 = sample[h].get[e];
					if(obj != obj2) {
						if(h == k) {
							int_vec_push(constraints[obj]->ml, obj2);
						} else {
							int_vec_push(constraints[obj]->mnl, obj2);
						}
					}
				}
			}
		}
	}
}

void print_constraints() {
    printf("Constraints:\n");
    size_t e;
    size_t i;
    size_t k;
    size_t obj;
	for(k = 0; k < classc; ++k) {
		for(i = 0; i < sample[k].size; ++i) {
			obj = sample[k].get[i];
            printf("Obj %d:\n", obj);
            printf("ML:");
            for(e = 0; e < constraints[obj]->ml->size; ++e) {
                printf(" %d", constraints[obj]->ml->get[e]);
            }
            printf("\n");
            printf("MNL:");
            for(e = 0; e < constraints[obj]->mnl->size; ++e) {
                printf(" %d", constraints[obj]->mnl->get[e]);
            }
            printf("\n");
		}
	}
}

st_matrix* medoid_dist(double *weights, size_t ***medoids) {
    st_matrix *ret = malloc(sizeof(st_matrix));
    init_st_matrix(ret, objc, clustc);
    size_t e;
    size_t i;
    size_t j;
    size_t k;
    double val;
    double sumd;
    for(i = 0; i < objc; ++i) {
        for(k = 0; k < clustc; ++k) {
            val = 0.0;
            for(j = 0; j < dmatrixc; ++j) {
                sumd = 0.0;
                for(e = 0; e < medoids_card; ++e) {
                    sumd += dmatrix[j][i][medoids[k][j][e]];
                }
                val += weights[j] * sumd;
            }
            set(ret, i, k, val);
        }
    }
    return ret;
}

st_matrix* agg_dmatrix(double *weights) {
    st_matrix *ret = malloc(sizeof(st_matrix));
    init_st_matrix(ret, objc, objc);
    size_t e;
    size_t i;
    size_t j;
    double val;
    for(i = 0; i < objc; ++i) {
        for(e = 0; e < objc; ++e) {
            val = 0.0;
            for(j = 0; j < dmatrixc; ++j) {
                val += weights[j] * dmatrix[j][i][e];
            }
            set(ret, i, e, val);
        }
    }
    return ret;
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
    // read labels
    int labels[objc];
    fscanf(cfgfile, "%d", &classc);
	size_t i;
    for(i = 0; i < objc; ++i) {
        fscanf(cfgfile, "%d", &labels[i]);
        if(labels[i] < 0 || labels[i] >= classc) {
            printf("Error: invalid object class.\n");
            return 2;
        }
    }
    // read labels end
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
    double sample_perc;
    fscanf(cfgfile, "%lf", &sample_perc);
    if(dlt(sample_perc, 0.0)) {
        printf("Error: sample_perc < 0.\n");
        return 2;
    }
    fclose(cfgfile);
    freopen(out_file_name, "w", stdout);
	mfuzval = 1.0 / (mfuz - 1.0);
    printf("######Config summary:######\n");
    printf("Number of clusters: %d.\n", clustc);
    printf("Medoids cardinality: %d.\n", medoids_card);
    printf("Number of iterations: %d.\n", max_iter);
    printf("Epsilon: %.15lf.\n", epsilon);
    printf("Theta: %.15lf.\n", theta);
    printf("Parameter m: %.15lf.\n", mfuz);
    printf("Number of instances: %d.\n", insts);
    printf("Sample percentage: %lf.\n", sample_perc);
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
	medoids = malloc(sizeof(size_t **) * clustc);
	size_t ***best_medoids = malloc(sizeof(size_t **) * clustc);
	for(k = 0; k < clustc; ++k) {
		medoids[k] = malloc(sizeof(size_t *) * dmatrixc);
		best_medoids[k] = malloc(sizeof(size_t *) * dmatrixc);
        for(j = 0; j < dmatrixc; ++j) {
            medoids[k][j] = malloc(sizeof(size_t) * medoids_card);
            best_medoids[k][j] = malloc(sizeof(size_t) *
                                        medoids_card);
        }
	}
    weights = malloc(sizeof(double) * dmatrixc);
    double *best_weights = malloc(sizeof(double) * dmatrixc);
	memb = malloc(sizeof(double *) * objc);
	double **best_memb = malloc(sizeof(double *) * objc);
	for(i = 0; i < objc; ++i) {
		memb[i] = malloc(sizeof(double) * clustc);
		best_memb[i] = malloc(sizeof(double) * clustc);
	}
    class = malloc(sizeof(int_vec) * classc);
    for(i = 0; i < classc; ++i) {
        int_vec_init(&class[i], objc);
    }
	// Allocating memory end
    // Loading labels
    for(i = 0; i < objc; ++i) {
        int_vec_push(&class[labels[i]], i);
    }
	print_class();
    // Loading matrices
	for(j = 0; j < dmatrixc; ++j) {
		if(!load_data(dmtx_file_name[j], dmatrix[j], objc, objc)) {
			printf("Error: could not load %s.\n", dmtx_file_name[j]);
			goto END;
		}
	}
    double avg_partcoef;
    double avg_modpcoef;
    double avg_partent;
    double avg_aid;
    silhouet *csil;
    silhouet *fsil;
    silhouet *ssil;
    silhouet *avg_csil;
    silhouet *avg_fsil;
    silhouet *avg_ssil;
    int *pred;
    st_matrix *groups;
    st_matrix *dists;
    st_matrix *agg_dmtx;
    st_matrix *memb_mtx;
	srand(time(NULL));
    size_t best_inst;
    double best_inst_adeq;
    double cur_inst_adeq;
	gen_sample(sample_perc * objc);
	print_sample();
	gen_constraints();
    print_constraints();
	for(i = 1; i <= insts; ++i) {
		printf("Instance %u:\n", i);
		cur_inst_adeq = run();
        memb_mtx = convert_mtx(memb, objc, clustc);
        pred = defuz(memb_mtx);
        groups = asgroups(pred, objc, classc);
        agg_dmtx = agg_dmatrix(weights);
        dists = medoid_dist(weights, medoids);
        csil = crispsil(groups, agg_dmtx);
        fsil = fuzzysil(csil, groups, memb_mtx, mfuz);
        ssil = simplesil(pred, dists);
        if(i == 1) {
            avg_partcoef = partcoef(memb_mtx);
            avg_modpcoef = modpcoef(memb_mtx);
            avg_partent = partent(memb_mtx);
            avg_aid = avg_intra_dist(memb_mtx, dists, mfuz);
            avg_csil = csil;
            avg_fsil = fsil;
            avg_ssil = ssil;
        } else {
            avg_partcoef = (avg_partcoef + partcoef(memb_mtx)) / 2.0;
            avg_modpcoef = (avg_modpcoef + modpcoef(memb_mtx)) / 2.0;
            avg_partent = (avg_partent + partent(memb_mtx)) / 2.0;
            avg_aid = (avg_aid +
                        avg_intra_dist(memb_mtx, dists, mfuz)) / 2.0;
            avg_silhouet(avg_csil, csil);
            avg_silhouet(avg_fsil, fsil);
            avg_silhouet(avg_ssil, ssil);
            free_silhouet(csil);
            free(csil);
            free_silhouet(fsil);
            free(fsil);
            free_silhouet(ssil);
            free(ssil);
        }
        free_st_matrix(memb_mtx);
        free(memb_mtx);
        free(pred);
        free_st_matrix(groups);
        free(groups);
        free_st_matrix(agg_dmtx);
        free(agg_dmtx);
        free_st_matrix(dists);
        free(dists);
        if(i == 1 || cur_inst_adeq < best_inst_adeq) {
            mtxcpy_d(best_memb, memb, objc, clustc);
            memcpy(best_weights, weights, sizeof(double) * dmatrixc);
            for(k = 0; k < clustc; ++k) {
                mtxcpy_size_t(best_medoids[k], medoids[k], dmatrixc,
                        medoids_card);
            }
            best_inst_adeq = cur_inst_adeq;
            best_inst = i;
        }
	}
    printf("Best adequacy %.15lf on instance %d.\n",
            best_inst_adeq, best_inst);
    printf("\n");
    print_medoids(best_medoids);
    printf("\n");
	print_memb(best_memb);
	printf("\n");
	print_weights(best_weights);

    st_matrix *best_memb_mtx = convert_mtx(best_memb, objc, clustc);
    pred = defuz(best_memb_mtx);
    groups = asgroups(pred, objc, classc);
    print_header("Partitions", HEADER_SIZE);
    print_groups(groups);

    dists = medoid_dist(best_weights, best_medoids);
    print_header("Best instance indexes", HEADER_SIZE);
    printf("\nPartition coefficient: %.10lf\n",
            partcoef(best_memb_mtx));
    printf("Modified partition coefficient: %.10lf\n",
            modpcoef(best_memb_mtx));
    printf("Partition entropy: %.10lf (max: %.10lf)\n",
            partent(best_memb_mtx), log(clustc));
    printf("Average intra cluster distance: %.10lf\n",
            avg_intra_dist(best_memb_mtx, dists, mfuz));

    print_header("Average indexes", HEADER_SIZE);
    printf("\nPartition coefficient: %.10lf\n", avg_partcoef);
    printf("Modified partition coefficient: %.10lf\n", avg_modpcoef);
    printf("Partition entropy: %.10lf (max: %.10lf)\n", avg_partent,
            log(clustc));
    printf("Average intra cluster distance: %.10lf\n", avg_aid);

    print_header("Averaged crisp silhouette", HEADER_SIZE);
    print_silhouet(avg_csil);
    print_header("Averaged fuzzy silhouette", HEADER_SIZE);
    print_silhouet(avg_fsil);
    print_header("Averaged simple silhouette", HEADER_SIZE);
    print_silhouet(avg_ssil);

    agg_dmtx = agg_dmatrix(best_weights);
    csil = crispsil(groups, agg_dmtx);
    print_header("Best instance crisp silhouette", HEADER_SIZE);
    print_silhouet(csil);
    fsil = fuzzysil(csil, groups, best_memb_mtx, mfuz);
    print_header("Best instance fuzzy silhouette", HEADER_SIZE);
    print_silhouet(fsil);
    ssil = simplesil(pred, dists);
    print_header("Best instance simple silhouette", HEADER_SIZE);
    print_silhouet(ssil);

    free_st_matrix(best_memb_mtx);
    free(best_memb_mtx);
    free(pred);
    free_st_matrix(groups);
    free(groups);
    free_st_matrix(dists);
    free(dists);
    free_st_matrix(agg_dmtx);
    free(agg_dmtx);
    free_silhouet(avg_csil);
    free(avg_csil);
    free_silhouet(avg_fsil);
    free(avg_fsil);
    free_silhouet(avg_ssil);
    free(avg_ssil);
    free_silhouet(csil);
    free(csil);
    free_silhouet(fsil);
    free(fsil);
    free_silhouet(ssil);
    free(ssil);
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
	}
	free(medoids);
	free(best_medoids);
	free(weights);
	free(best_weights);
	for(i = 0; i < objc; ++i) {
		free(memb[i]);
		free(best_memb[i]);
		if(constraints[i]) {
            constraint_free(constraints[i]);
        }
	}
    free(constraints);
	free(memb);
	free(best_memb);
	for(i = 0; i < classc; ++i) {
		int_vec_free(&class[i]);
		int_vec_free(&sample[i]);
	}
	free(class);
	free(sample);
	free(parc_cluster_adeq);
    free(parc_obj_adeq);
	return 0;
}
