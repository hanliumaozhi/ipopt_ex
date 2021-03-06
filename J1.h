/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

int jac_f1(casadi_real* arg, casadi_real* res, casadi_int* iw, casadi_real* w, int mem);
int jac_f1_alloc_mem(void);
int jac_f1_init_mem(int mem);
void jac_f1_free_mem(int mem);
int jac_f1_checkout(void);
void jac_f1_release(int mem);
void jac_f1_incref(void);
void jac_f1_decref(void);
casadi_int jac_f1_n_out(void);
casadi_int jac_f1_n_in(void);
casadi_real jac_f1_default_in(casadi_int i);
const char* jac_f1_name_in(casadi_int i);
const char* jac_f1_name_out(casadi_int i);
const casadi_int* jac_f1_sparsity_in(casadi_int i);
const casadi_int* jac_f1_sparsity_out(casadi_int i);
int jac_f1_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w);
#ifdef __cplusplus
} /* extern "C" */
#endif
