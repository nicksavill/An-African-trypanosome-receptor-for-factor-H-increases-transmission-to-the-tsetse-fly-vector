#include <model.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>

#define NVSG 3
#define IMMUNOCOMPETENT (trt <= 2)

void Variables(struct variables *V) 
{
  int i;
  Var(V, 0, "day");
  Var(V, 1, "log-concentration");
  Var(V, 2, "proportion slender");
  Var(V, 3, "proportion stumpy");
  Var(V, 4, "SIF");
  Var(V, 5, "I_0");
  Var(V, 6, "I_1");
  for (i = 0; i < NVSG; i++)
    Var(V, i+7, "log-variant concentration");
}

void Parameters(int trt, int ind, double *p, struct parameterset *Q) 
{
  // Par(Q, 0,  Real, "alpha",   TruncNormal, 0.1,   V(0.1));
  Par(Q, 0,  Real, "alpha",   Cst, 0.14);
  Par(Q, 1,  Real, "beta",    TruncNormal, 1.e-10, V(1.e-10));
  Par(Q, 2,  Real, "gamma",   TruncNormal, 0.05,   V(0.05));
  Par(Q, 6,  Real,  "tau_{l}", TruncNormal, 24.,    V(4.));
  Par(Q, 7,  Real,  "tau_{s}", TruncNormal, 48.,    V(24.));
  if (IMMUNOCOMPETENT) {
    Par(Q, 12, Real, "l_{0}", Normal,      3.,    V(0.5));
    Par(Q, 3,  Real, "delta", Uniform,     0.0,   1.0);
    Par(Q, 13, Real, "psi",   TruncNormal, 1.e-10, V(1.e-10));
    Par(Q, 20, Real, "omega0", Normal,      -6.,   2.);
    Par(Q, 21, Real, "omega1", Normal,      -6.,   2.);
  }
  else
    Par(Q, 12, Real, "l_{0}", Normal,      5.5,    V(0.5));

  Par(Q, 22, Real, "mean_s", Der);
  // Par(Q, 41, Real, "b", Der);
  // Par(Q, 43, Real, "c", Der);
  // Par(Q, 44, Int,  "d", Der);
}

void DerivedParameters(int trt, int ind, int nhours, double *p, double **V, TimePoints *TP)
{
  int t;
  double mean_s = 0;
  Model(trt, ind, nhours, p, V, TP);

  for (t = 0; t < nhours; t++)
    mean_s += V[t][1]*V[t][3];
  p[22] = mean_s/(double)nhours;
}

int Model(const uint trt, const uint ind, const uint nhours, double *p, double **V, TimePoints *TP)
{
  int v, t, a, n = NVSG;
  double c, l, l1, s, f, fp, dva;
  double alpha = p[0];
  double beta  = p[1];
  double gamma = p[2];
  double delta = p[3];
  double l0    = exp10(p[12]);
  double psi   = p[13];
  int    tau_l = lrint(p[6]);
  int    tau_s = tau_l+lrint(p[7]);

  double *omega = doublevector(n);  // immune response of variant v
  double *I     = doublevector(n);  // immune response of variant v
  double *lv    = doublevector(n);  // total non-committed slender of variant v
  double *l1v   = doublevector(n);  // total committed slender of variant v
  double *sv    = doublevector(n);  // total stumpy of variant v
  double *lvp   = doublevector(n);  // storage for non-committed slender
  double **d    = doublematrix(n, nhours); // age structure

  if (IMMUNOCOMPETENT) {
    gamma = p[2];
    omega[0] = exp10(p[20]);
    omega[1] = exp10(p[21]);
  }
  // else
  //   gamma = 0.05;

  // initial non-committed slender of variant 0
  lv[0] = l0;
  // SIF = 0
  f = 0;
  // no immune response against variant 0 at t=0

  for (t = 0; t < nhours; t++) {
    l = l1 = s = 0;

    // calculated total non-committed, committed slender and stumpy for each variant v
    for (v = 0; v < n; v++) {
      l1v[v] = sv[v] = 0;
      // total non-committed slender 
      if (lv[v]) 
        l += lv[v];
      
      for (a = 0; a < min(tau_s, t); a++)
        if ((dva = d[v][a])) {
          if (a < tau_l) {
            // total committed slender (and each variant)
            l1v[v] += dva;
            l1     += dva;
          }
          else {
            // total stumpy (and each variant)
            sv[v] += dva;
            s     += dva;
          }
        }
    }
    // total cells per ml
    c = l+l1+s;
    
    // save variables
    V[t][0] = t;
    V[t][1] = c;
    V[t][4] = f;
    V[t][5] = I[0];
    V[t][6] = I[1];
    if (c) {
      V[t][2] = (l+l1)/c;
      V[t][3] = s/c;
    } else {
      V[t][2] = 1;
      V[t][3] = 0;
    }
    // total cells/ml of each variant
    for (v = 0; v < n; v++)
      V[t][v+7] = lv[v]+l1v[v]+sv[v];

    // SIF produced by non-committed and committed slenders, lost at rate gamma
    fp = l+l1+f*exp(-gamma);

    // update committed slenders and stumpies by 1 hour
    for (v = 0; v < n; v++) {
      for (a = min(tau_s-1, t-1); a >= 1; a--)
        if (d[v][a-1] > 0) {
          // committed cells are killed by immune response at rate I[v]
          d[v][a] = d[v][a-1]*exp(-I[v]);
          // committed slenders reproduce at rate alpha
          if (a < tau_l)
            d[v][a] *= exp(alpha);
        }
      
      // slenders become committed due to SIF at rate beta*SIF, and enter at age 0
      if (lv[v])
        d[v][0] = lv[v]*(1.-exp(-beta*f));

      // slenders reproduce at rate alpha, differentiate at rate beta*f and 
      // killed at rate delta*I[v]
      lvp[v] = lv[v]*exp(alpha-delta*I[v]-beta*f);
    }

    // variant switch 0->1 and 1->2, 2 does not switch
    if (IMMUNOCOMPETENT)
      for (v = 0; v < n-1; v++)
        if (lv[v]) {
          lvp[v]   -= omega[v]*lv[v];
          lvp[v+1] += omega[v]*lv[v];
        }

    // increase immune response against variants 0 and 1 up to a max of 1
    if (IMMUNOCOMPETENT)
      for (v = 0; v < n-1; v++)
        I[v] = 1.-(1.-I[v])*exp(-psi*lv[v]);

    // update non-committed slenders from temporary storage
    for (v = 0; v < n; v++)
      if (lvp[v]) 
        lv[v] = fmax(0, lvp[v]);

    // update SIF
    f = fp;
  }
  
  free(omega); free(I);
  free(d[0]); free(d); free(lv); free(lvp); free(l1v); free(sv);
  return SUCCESS;
}

#define ML_PER_FIELD (32./pow(10, 8.1)) // reference values from herbert (1976)

double logLikelihood(const uint trt, const uint ind, const double *p, const double *V, const double *Data, const uint *value)
{
  if (value[1] == NA)
    return 0;

  static unsigned int categories[] = {6, 8, 12, 16, 24, 32, 48, 64, 128, 192, 256, 1024};
  // static double range[] = {          5, 7, 10, 14, 20, 28, 40, 56, 96, 160, 224, 1024};
  int category = lrint(Data[1]), index;
  double l, vc = V[1]*ML_PER_FIELD;  //expected cells per field (cells/ml*ml/field)
  gsl_sf_result g1, g2;
  int n1, n2;

  if      (category == -4) l = poissonPMF(1, 20.*vc);
  else if (category == -3) l = poissonPMF(2, 20.*vc)+poissonPMF(3, 20.*vc);
  else if (category == -2) l = poissonPMF(2, 10.*vc)+poissonPMF(3, 10.*vc);
  else if (category == -1) l = poissonPMF(2,  5.*vc)+poissonPMF(3,  5.*vc);
  else if (category ==  1) l = poissonPMF(4,  5.*vc)+poissonPMF(5,  5.*vc);
  else if (category <=  4) l = poissonPMF(category, vc);
  else {
    // find the index in categories of the recorded cell count
    for (index = 0; index < 12; index++)
      if (category == categories[index]) break;
    // n1 = range[index];
    // n2 = range[index+1];
    n1 = categories[max(0, index-1)];
    n2 = categories[min(11, index+1)];
    // find the probable range of cells per field
    // range[index] is the lower value of the range
    // range[index+1] is the upper value of the range+1
    if (gsl_sf_gamma_inc_Q_e(n1, vc, &g1)) {
      printf("error\n");
      return -INFINITY;
    }
    if (gsl_sf_gamma_inc_Q_e(n2, vc, &g2)) {
      printf("error\n");
      return -INFINITY;
    }
    if (g2.val-g1.val == 0)
      // prevents -infinity from halting the sim
      l = -100;
    else
      l = log(g2.val-g1.val);
  }
  // if (!isfinite(l))
    // printf("%d %d %f %d %f %f %f\n", trt, ind, Data[0], category, V[1], vc, l);
  return l;
}
    
void SaturatedModel(int trt, int ind, double *V, double *p, double *Data, uint *value)
{
  int category = (int)Data[1];

  if      (category == -4) V[1] = 1.0/(20.*ML_PER_FIELD);
  else if (category == -3) V[1] = 2.5/(20.*ML_PER_FIELD);
  else if (category == -2) V[1] = 2.5/(10.*ML_PER_FIELD);
  else if (category == -1) V[1] = 2.5/( 5.*ML_PER_FIELD);
  else if (category ==  1) V[1] = 4.5/( 5.*ML_PER_FIELD);
  else V[1] = (double)category/ML_PER_FIELD;
}

void OutputData(int trt, int ind, double *Output, double *Var, double *Data, uint *value)
{
  int category = (int)Data[1];

  if      (category == -4) Output[1] = 1.0/(20.*ML_PER_FIELD);
  else if (category == -3) Output[1] = 2.5/(20.*ML_PER_FIELD);
  else if (category == -2) Output[1] = 2.5/(10.*ML_PER_FIELD);
  else if (category == -1) Output[1] = 2.5/( 5.*ML_PER_FIELD);
  else if (category ==  1) Output[1] = 4.5/( 5.*ML_PER_FIELD);
  else Output[1] = (double)category/ML_PER_FIELD;
  Var[1] = 0;
  Output[1] = log10(Output[1]);
}

void OutputModel(int trt, int ind, double *Output, double *V)
{
  int i;
  for (i = 1; i < 10; i++)
    if (i == 1)
      Output[i] = log10(V[i]);
    else if (i == 4)
      Output[i] = V[i]/1e9;
    else if (i >= 7)
      Output[i] = V[i]/1e8;
    else
      Output[i] = V[i];
}

double timestep(void) {return 1;}
void PredictData(int trt, int ind, double *Output, double *V, double *p, gsl_rng *stream) {}
void SimulateData(int trt, int ind, double *Output, double *V, double *p, double *Data, uint *value, gsl_rng *stream) {}
void GlobalParameters() {}
double UserPDF(double x) {return 0;}
void HyperParameters(Treatments T, Hyperparameters *H) {}
int function(int trt, int ind, double *x, double *y, double *p, TimePoints *TP) {return 0; }
void Residual(int trt, int ind, double *R, double *V, double *p, double *Data, uint *value) {}
void WAIC(int trt, int ind, double *lnL, double *V, double *p, double *Data, uint *value) {}
