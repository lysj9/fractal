void randomz_seed(int seed);
double randomz();
double gaussrand1( double,double );
double gaussrand2( double,double );
void randqueue( int,int*,int,int* );
void randqueue2(int n, int *idx, struct vector_s *star, int randnum, int *randidx, double eps);
//void fractal(int,double,double*,struct vector_s*);
void fractal(int StarNum, double D, double mlow, double mhigh, struct vector_s *star);
int getopt(int,char*[],char*);
void generate_mass ( int N, double mlow, double mhigh, double *mass );
void generate_binaries( int N_star, int nbin,double *mass, struct vector_s *star );
double make_mass ( double mlow, double mhigh );
double nbody_scale(int N_node, double q, struct vector_s *star, int *nnbmax_out, double *rs0_out);

void mempool_init(size_t);
void* mempool_malloc();
void mempool_free(void*);
void mempool_destroy();

void quick_sort_widx(double *a, int *idx, int n);
void quick_sort(double *a, int n);
double quick_select(double *a, int k, int n);

void sort_radius(int N_cm, struct vector_s *star, double *r2_sort, int *idx);
double get_radius(int N_cm, struct vector_s *star, double truncate, double *r2_sort, int *idx);

//void output_nbody6(char *outname, int N_star, int nbin, double *r2_sort, int seed, double r_virial, double m_mean, double mlow, double mhigh, double q, int *NNBMAX, double *RS0);
void output_nbody6(char *outname, int N_star, int nbin, double *r2_sort, int seed, double r_virial, double m_mean, double mlow, double mhigh, double q, int NNBMAX, double RS0);
