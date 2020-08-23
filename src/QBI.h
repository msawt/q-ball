#ifndef QBI_H_
#define QBI_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"nifti1_io.h"
#include"Interpolation.h"
#include"ConvolutionAndPrediction.h"
#include"DirOfInterest.h"
#include"Matrices.h"


#ifdef _WIN32
const char DIRSEP = '\\';
#else
const char DIRSEP = '/';
#endif

#define DO_TENSOR 0

const double DECO_P        = 2.0;
const double DECO_EVALS[3] = {1.5, 0.4, 0.4};
const double DELTA_SM      = 16.9100;
const double DELTA_LG      = 29.7300;
const double DIFF_RADIUS   = 10.0;
const double PROB_THRESH   = 0.5;

typedef struct ICOS_TESS {
    int num_vertices;
    float **vertices;
    int **connectivity;
} ICOS_TESS;


typedef struct OUTPUT_DATA {
    nifti_image **nim;
    nifti_brick_list **nbl;
    int num_images;
} OUTPUT_DATA;

typedef struct MAXIMA {
    double value;
    int index;
} MAXIMA;

typedef struct SYMMAT33 {
    float data[6];
} SYMMAT33;

typedef struct MAT33 {
    float data[9];
} MAT33;

typedef struct DIFF_DATA {
    nifti_brick_list *nii_brick_list;
    nifti_image *nii_image;
    nifti_image *mask;
    double *single_voxel_storage;
    int n_volumes;
    float *bvecs;
    int *b_high_ind;
    int *b_low_ind;

} DIFF_DATA;

typedef struct QBI_RECON {
    DIFF_DATA *diff;
    ICOS_TESS *reco_tess; //One of these tess is equal to U/V, need to ask which one
    ICOS_TESS *restart_tess;
    double deco_p; //Don't know if this is needed, have to ask
    double deco_evals[3];
    double prob_thresh;
    char *output_directory;
    char *data_filename; //e
    char *mask_filename;
    char *bvec_filename; //Q
    char *datadir;
    int log_bad_voxels;
    int num_output_files;
} QBI_RECON;

void qbi_initialize_opts(QBI_RECON *qbi, int argc, char **argv);
ICOS_TESS *load_tess_from_file(const char *tess_file);
void read_diff_data_from_file(char *filename, DIFF_DATA *diff);
void __tokenize_line(char *line, const char *tokens, float *split, int *nsplit);
int read_acsii_file_to_float_array(char *filename, float *data, int data_size);
int maxima_compare(const void *m1, const void *m2);
int find_local_maxima(ICOS_TESS *domain, double *values, double min_value, ICOS_TESS *restart, MAXIMA *maxima_list);
void qbi_print_usage(void);
OUTPUT_DATA *initialize_output (nifti_image *template, int nfiles);
inline int nii_voxel3_index(nifti_image *nim, const int x, const int y, const int z);
void read_bvecs_from_file(char *filename, DIFF_DATA *diff);
int load_voxel_double_all(DIFF_DATA *diff, int x, int y, int z, double *dest);
void add_maxima_to_output(OUTPUT_DATA *output, int x, int y, int z, double **vertlist, MAXIMA *maxima_list, int n_maxima);
double read_nii_voxel_anytype(void *src, int index, int datatype);
int nii_recast_to_int32 (nifti_image *nim);
int save_output (const char *basedir, OUTPUT_DATA *output);
void QBall(double* e, double ** Q, double ** U, double** V, int k,int n, int m, int p,double* ODF);




int temp_find_max(double* ODF, int ODFLen, MAXIMA* maxima_list);

#endif
