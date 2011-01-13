//R CMD SHLIB raremethods.cpp

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <iostream>
using namespace std;

#include "thmalloc.h"

/** Random number generation **/

// attaching and detaching the random number generator
extern "C"{
  void rndAttach()
  {
    GetRNGstate();
  }
  void rndDetach()
  {
    PutRNGstate();
  }
}


// Uniform random number U[0,1)
double RUnif() { return unif_rand(); }
//double RUnif() { return( (double)rand() / ( (double)(RAND_MAX) + 1 ) ); }
double RUnifExt( double min, double max ){
  return( unif_rand() * (max-min) + min );
}

// N(0,1)
double RNormal() { return norm_rand(); }

// N(mean,sigma^2)
double RNormalExt( double mean, double sigma ){
  return( mean + RNormal()*sigma );
}

// Gets a random integer b/n min & max, including both endpoints
int RandInt( int min, int max ){
  int range = (max - min) + 1;
  return(  (int)(RUnif() * range)  + min );
}

/** Computing the allele frequency **/

extern "C" {
  void gAfreq(const double *g, const int *gcols,
              const int *nindiv, const int *ng,
              double *ret_afreq) {
    // computes the allele frequency, but only at the gcols (others are left as -1)
    int N = *nindiv, G = *ng;
    for(int c=0; c<G; c++){
      if(gcols[c] == 1){
        // then compute this one
        ret_afreq[c] = 0;
        int n = 0;
        for(int r=0; r<N; r++){
          if(!isnan(g[r + c*N])){
            n++;
            ret_afreq[c] += g[r + c*N];
          }
        }
        ret_afreq[c] /= (2*n);
      }else{
        ret_afreq[c] = -1;
      }
    }
  } // completely debugged

  void gAfreqCACO(const double *g, const int *gcols,
                  const double *aff,
                  const int *nindiv, const int *ng,
                  double *ret_afreq_case, double *ret_afreq_cont){
    int N = *nindiv, G = *ng;
    for(int c=0; c<G; c++){
      if(gcols[c] == 1){
        ret_afreq_case[c] = ret_afreq_cont[c] = 0.0;
        int ncase = 0, ncont = 0;
        for(int r=0; r<N; r++){
          if(!isnan(g[r + c*N])){
            if(aff[r] == 1){
              ncase++;
              ret_afreq_case[c] += g[r + c*N];
            }else{
              ncont++;
              ret_afreq_cont[c] += g[r + c*N];
            }
          }
        }
        ret_afreq_case[c] /= (2*ncase);
        ret_afreq_cont[c] /= (2*ncont);
      }else{
        ret_afreq_case[c] = -1;
        ret_afreq_cont[c] = -1;
      }
    }
  } // not debugged -- intended for the 'sign' routine
}

extern "C" {
  void sortRemDupl(double *v, int *n){
    // simple bubble sort, but removing zeroes...
    for(int i=0; i<*n; i++){
      for(int j=i+1; j<*n; j++){
        if(v[i] > v[j]){
          //swap
          double temp = v[i];
          v[i] = v[j];
          v[j] = temp;
        }else if(v[i] == v[j]){
          // swap v[j] with last, and decrement n
          v[j] = v[(*n) - 1];
          (*n)--;
        }
      }
    }
  }

  void afreq_uafreq(const double *g, const int *gcols,
                    const int *nindiv, const int *ng,
                    double *afreq, double *uafreq, int *nuafreq){ // returns, really should be named that way...
    // get all unique allele frequencies (sorted as well for debug) --> afreq, nafreq
    int N = *nindiv, G = *ng;
    for(int i=0; i<G; i++)
      afreq[i] = 0;
    gAfreq(g, gcols, nindiv, ng, afreq);
    for(int i=0; i<G; i++)
      uafreq[i] = afreq[i];

    //` << "Raw afreq ";
    //for(int i=0; i<G; i++)
    //  cout << afreq[i] << " ";
    //cout << endl;

    // remove any duplicates
    *nuafreq = G;
    for(int i=0; i<G; i++)
      uafreq[i] = afreq[i];
    sortRemDupl(uafreq, nuafreq);

    // and remove the '0' element
    if(*nuafreq > 0 && uafreq[0]==0){
      for(int i=0; i<(*nuafreq)-1; i++)
        uafreq[i] = uafreq[i+1];
      (*nuafreq)--;
    }
  } // seems debugged

  void fix_gcols_afreq(int *gcols, const int *ng, const double *afreq){
    int G = *ng;
    for(int i=0; i<G; i++)
      if(afreq[i] == 0)
        gcols[i] = 0;
  }

  void permute_double(const double *in, int in_size, double *out){
    for(int i=0; i<in_size; i++)
      out[i] = in[i];
    for(int i=0; i<in_size; i++){
      // I mixed up out instead of in here, but this actually doesn't screw things up too horribly... still permutes, except for the first one...
      int randn = RandInt(0, in_size-1);
      double temp = out[i];
      out[i] = out[randn];
      out[randn] = temp;
    }
  }
}

// returns true if there is another, false when done
bool next_mask(bool *mask, int n){
  // essentially does binary addition
  int i = 0;
  for(i=0; i<n && mask[i]; i++)
    mask[i] = false;
  if(i < n){
    mask[i] = true;
    return(true);
  }
  return(false);
}

bool next_gcols_mask(int *mask, int *gcols, int n){
  //adjust n
  //while(gcols[n-1]==0 && n>0)
  //  n--;

  // essentially binary addition, skipping where gcols != 1
  int i=0;
  for(i=0; i<n && (mask[i] || !gcols[i]); i++)
    mask[i] = false;

  //printv(mask, n);

  if(i < n){
    while(!gcols[i])
      i++;
  }
  if(i < n){
    mask[i] = 1;
    return(true);
  }
  return(false);
}

extern "C"{
  // actually calculate the squared value of this one...
  void zstat(const double *g, const int *m, const int *ng,
             const int *s, const double *w, const double *aff,
             const int *naff,
             double *res){
    // calculate the means
    int K = *ng, I = *naff;

    double tbar = 0; // could make extrenal? speed?
    for(int i=0; i<I; i++)
      tbar += aff[i];
    tbar /= I;

    ///int KN = 0;
    ///for(int k=0; k<K; k++)
    ///  KN += m[k];
    ///double xbar = 0; // could make external? speed?
    ///for(int k=0; k<K; k++)
    ///  if(m[k] == 1)
    ///    for(int i=0; i<I; i++)
    ///      xbar += g[i + k*I];
    ///xbar /= (I * KN); // need to be careful here -- not all columns are used!!!

    ///double xbar_num = 0.0, xbar_den = 0.0;
    ///for(int k=0; k<K; k++){
    ///  if(m[k] == 1){
    ///    for(int i=0; i<I; i++){
    ///      xbar_num += w[k] * g[i + k*I];
    ///      xbar_den += w[k];
    ///    }
    ///  }
    ///}
    ///double xbar = xbar_num / xbar_den;

    THMALLOC(double, _xbar, xbar, K); //double xbar[K];
    for(int k=0; k<K; k++){
      xbar[k] = 0.0;

      int k_x_I = k*I;
      if(m[k] == 1) {
        for(int i=0; i<I; i++)
          xbar[k] += g[i + k_x_I];
        xbar[k] /= I;
      }
    }

    // and then calculate the test statistic
    double num=0.0, den=0.0;
    for(int i=0; i<I; i++){
      double temp = 0;
      for(int k=0; k<K; k++)
        if(m[k] == 1)
          temp += s[k] * w[k] * (g[i + k*I] - xbar[k]) * (aff[i] - tbar);
//          temp += s[k] * w[k] * (g[i + k*I] - xbar) * (aff[i] - tbar);
      num += temp;
      den += temp * temp;
    }

    *res = (num*num) / den;
    //cout << "zstat " << *res << endl; exit(1);
  }
}

void zstat2(const double *g, const int *m, const int *ng,
            const double *aff, const int *naff,
            const int *use_sign, const int *use_weight,
            double *res){
  int K = *ng, I = *naff, G = *ng, N = I;
  const int *gcols = m;

  // 1) compute the weights
  THMALLOC(double, _wi, wi, K); //double wi[K];

  if(*use_weight == 0){
    for(int k=0; k<K; k++)
      wi[k] = 1.0;
  }else if(*use_weight == 1 || *use_weight == 2){
    // dichotomous
    //double miu[G], niu[G], qi[G], ni[G]; //, wi[G];  // you stupid, fucking idiot!!!!! No wonder weighting wasn't working, they weren't really being set!!!!!!!!!!!! ARGHHHHHHHHHH!!!!!!!!!!!!!!!
    THMALLOC(double, _miu, miu, G);
    THMALLOC(double, _niu, niu, G);
    THMALLOC(double, _qi, qi, G);
    THMALLOC(double, _ni, ni, G);
    for(int c=0; c<G; c++){
      miu[c] = niu[c] = qi[c] = ni[c] = wi[c] = 0.0;

      if(gcols[c] == 1){
        for(int r=0; r<N; r++){
          if(aff[r] == 0) {
            miu[c] += g[r + c*N];
            niu[c] += !isnan(g[r + c*N]);
          }//fi
          ni[c] += !isnan(g[r + c*N]);
        }//rof
        qi[c] = (miu[c] + 1.0) / (2.0*niu[c] + 2.0);
        //wi[c] = 1.0 / sqrt(ni[c] * qi[c] * (1-qi[c])); // the weights... --> major fix, the weights were bloody inverted!!!
        wi[c] = 1.0 / sqrt(qi[c] * (1-qi[c])); // the weights... --> major fix, the weights were bloody inverted!!!
        //cout << wi[c] << ' ';
      }//fi
    }//rof
    //cout << endl;
  }else if(*use_weight == 3 || *use_weight == 4){
    // continuous
    for(int k=0; k<K; k++){
      wi[k] = 0;
      if(gcols[k] == 1){
        double pk = 0.0;
        for(int r=0; r<N; r++)
          pk += g[r + k*N];
        pk = (pk + 1.0) / (2.0*I + 2.0);
        wi[k] = 1.0 / sqrt(pk * (1-pk));
      }
    }
  }else{
    // error
    cout << "ERROR: zstat2, (*use_weight) value is not possible (" << *use_weight << ")" << endl;
    exit(1);
  }

  // 2) compute the signs
  THMALLOC(int, _s, s, K); //int s[K];
  if(*use_sign == 0){
    for(int k=0; k<K; k++)
      s[k] = 1;
  }else if(*use_sign == 1){
    ///double afreq_case[K], afreq_cont[K];
    ///gAfreqCACO(g, m,  aff,  naff, ng,  afreq_case, afreq_cont);
    ///for(int k=0; k<K; k++){
    ///  if(afreq_case[k] < afreq_cont[k]){
    ///    s[k] = -1;
    ///  }else if(afreq_case[k] == afreq_cont[k]){
    ///    s[k] = 0;
    ///  }else{
    ///    s[k] = 1;
    ///  }
    ///}

    // That doesn't work for continuous traits!!!
    THMALLOC(double, _xbar, xbar, K); //double xbar[K];
    for(int k=0; k<K; k++){
      xbar[k] = 0;
      if(m[k] == 1){
        for(int i=0; i<I; i++)
          xbar[k] += g[i + k*I];
        xbar[k] /= I;
      }
    }
    double affbar = 0.0;
    for(int i=0; i<I; i++)
      affbar += aff[i];
    affbar /= I;
    THMALLOC(double, _cov, cov, K); //double cov[K];
    for(int k=0; k<K; k++){
      cov[k] = 0.0;
      s[k] = 0.0;
      if(m[k] == 1){
        for(int i=0; i<I; i++)
          cov[k] += (g[i + k*I] - xbar[k]) * (aff[i] - affbar);
      }
      if(cov[k] < 0){
        s[k] = 1.0;
      }else if(cov[k]>0){
        s[k] = -1.0;
      }else{
        s[k] = 0.0;
      }
      //cout << cov[k] << '(' << s[k] << ')' << ' ';
    }
    //cout << endl;
  }else{
    cout << "ERROR: zstat2, (*use_sign) value is not possible (" << *use_sign << ")" << endl;
    exit(1);
  }

  // 3) and then pass these along and compute the test statistic!
  zstat(g, m, ng,  s, wi, aff,  naff,  res);

  if(*use_weight == 2 || *use_weight == 4){
    // then we really wanted _both_; weighted has been run, so reset weights to 1, and give the best_loc
    for(int k=0; k<K; k++)
      wi[k] = 1.0;
    double res2 = 0.0;
    zstat(g, m, ng,  s, wi, aff,  naff,  &res2);
    if(res2 > *res)
      *res = res2;
  }
}


extern "C" {
  void zstat_perm(const double *g, int *m, const int *ng, // NEW: NO LONGER int *m
                  const double *aff, const int *naff, // warning naff is the length of aff, _not_ the number of affected
                  const double *thresh, // allele frequency threshhold (for hard)
                  const int *gsubsetmatrix, const int *ngsubset,
                  int *use_sign, int *use_weight, int *strategy, // strategy = 1 (HARD THRESH), 2 (ALL AFREQ), 3 (STEP UP), 4 (ALL SUBSETS)
                  int *nperm,
                  double *ret_pvalue) {
    rndAttach();

    int P = *nperm, I = *naff, G = *ng, NGS=*ngsubset;

    // sometimes we will need allele frequencies
    //double afreq[G], uafreq[G];
    THMALLOC(double, _afreq, afreq, G);
    THMALLOC(double, _uafreq, uafreq, G);
    int nuafreq;
    afreq_uafreq(g, m, naff, ng, afreq, uafreq, &nuafreq);

    // gear up to compute the test statistic
    THMALLOC(double, _z, z, P+1); //double z[P + 1]; // z[0] is the observed...
    THMALLOC(double, _aff_perm, aff_perm, I); //double aff_perm[I];
    THMALLOC(int, _mbest, mbest, G); //int mbest[G];

    // for the initial, we want to calculate the observed
    //memcpy(aff, aff_perm, sizeof(double)*I);
    for(int i=0; i<I; i++)
      aff_perm[i] = aff[i];

    // need to alter this, and work in gsubsetmatrix!!!
    for(int p=0; p<=P; p++){
      z[p] = 0;

      for(int ngs=0; ngs<NGS; ngs++){ // over gsubsetmatrix
        THMALLOC(int, _msub, msub, G); //int msub[G];
        for(int k=0; k<G; k++)
          msub[k] = gsubsetmatrix[ngs + k*NGS] && m[k] && (afreq[k]<=(*thresh));

        // compute under each strategy
        double znew = 0.0;

        if(*strategy == 1){
          // HARD THRESHHOLD
          zstat2(g, msub, ng,  aff_perm, naff,  use_sign, use_weight,  &znew);
          if(znew > z[p]) z[p] = znew; // store it if better
        }else if(*strategy == 2){
          // ALL ALLELE FREQUENCIES
          double znew;
          for(int ua=0; ua<nuafreq; ua++){
            THMALLOC(int, _m2, m2, G); //int m2[G];
            for(int k=0; k<G; k++)
              m2[k] = (int)((msub[k] != 0) & (afreq[k] <= uafreq[ua]));
            zstat2(g, m, ng,  aff_perm, naff,  use_sign, use_weight,  &znew);
            if(znew > z[p]) z[p] = znew; // store it if better
          }
        }else if(*strategy == 3 || *strategy == 33){
          // STEP UP
          THMALLOC(int, _mask, mask, G); //int mask[G];
          for(int k=0; k<G; k++)
            mask[k] = 0.0;
          THMALLOC(int, _mask_diff, mask_diff, G); //int mask_diff[G];

          bool keep_going = true;
          while(keep_going){
            // generate the mask difference
            bool any_nonzero = false;
            for(int k=0; k<G; k++){
              mask_diff[k] = msub[k] - mask[k];
              any_nonzero = any_nonzero || (mask_diff[k] == 1);
            }

            // try all the new spots in the mask difference
            int best_loc = -1;
            for(int k=0; k<G; k++){
              if(mask_diff[k] == 1){
                mask[k] = 1; // twiddle it on

                zstat2(g, mask, ng,  aff_perm, naff,  use_sign, use_weight,  &znew);
                if(znew > z[p]) {
                  z[p] = znew; // store it if better
                  best_loc = k;
                }

                mask[k] = 0; // twiddle it back off
              }
            }
            // now, was there a difference?
            if(best_loc == -1){
              // best has already come to pass -- stop!
              keep_going = false;
            }else{
              // we found something better -- keep going!
              mask[best_loc] = 1; // set it on
            }

          }
          if(p == 0 && *strategy==33){
            for(int k=0; k<G; k++)
              m[k] = mask[k];
          }
        }else if(*strategy == 4){
          // ALL POSSIBLE SUBSETS
          THMALLOC(int, _mask, mask, G); //int mask[G];
          for(int k=0; k<G; k++)
            mask[k] = 0;
          while(next_gcols_mask(mask, msub, G)){
            zstat2(g, mask, ng,  aff_perm, naff,  use_sign, use_weight,  &znew);
            if(znew > z[p]) z[p] = znew; // store it if better
          }
        }else{
          // error!!!
          cout << "ERROR: zstat_meta, (*strategy) is not possible (" << *use_sign << ")" << endl;
          exit(1);
        }
      }

      // then permute --> note we do this after the initial, so that z[0] is observed, and z[1] to z[P] is the permuted test statistic...
      permute_double(aff, I, aff_perm);
    }

    rndDetach();

    // basically mask stores the best mask...

    // you need to be _very_ careful with += and /= with pointers -- don't use them, they actually add to the pointer, not the value!!!
    double pvalue = 0.0;
    for(int p=1; p<=P; p++)
      pvalue += (int)(z[p] >= z[0]);
    pvalue /= (double)P;
    *ret_pvalue = pvalue;

    //cout << "z[p] values: ";
    //for(int p=0; p<=P; p++)
    //  cout << z[p] << " ";
    //cout << endl;
    //cout << "pvalue " << *ret_pvalue << endl;
    //exit(1);
  }
}

// nonzero m now indicates the gene... -- should be numeric
void zstat_pathway_stat(const double *g, const int *m, const int *ng,
                  const double *aff, const int *naff, // warning naff is the length of aff, _not_ the number of affected
                  const double *thresh, // allele frequency threshhold (for hard)
                  const int *gsubsetmatrix, const int *ngsubset,
                  int *use_sign, int *use_weight, int *strategy, // strategy = 1 (HARD THRESH), 2 (ALL AFREQ), 3 (STEP UP), 4 (ALL SUBSETS)
                  int *nperm,
                  double *ret_pvalue,
                  bool print) {

  // First we need to know how many unique genes there are...
  int P = *nperm, I = *naff, G = *ng, NGS=*ngsubset;

  //cout << "P " << P << " I " << I << " G " << G << " NGS " << NGS << endl;

  int numGenes = 0;
  THMALLOC(int, _genes, genes, G);
  ///int genes[G];///
  for(int k=0; k<G; k++){
    bool geneFound = false;
    for(int cg=0; cg<numGenes; cg++)
      if(m[k]!=0 && genes[cg]==m[k])
        geneFound = true;
    if(!geneFound){
      genes[numGenes] = m[k];
      numGenes++;
    }
  }
  //cout << "**** There are " << numGenes << " genes found." << endl; // DEBUG ONLY

  ///int s[G], w[G];
  THMALLOC(int, _s, s, G);
  THMALLOC(int, _w, w, G);
  for(int k=0; k<G; k++){
    s[k] = 1;
    w[k] = 1;
  }

  // For each individual gene, make a note of the best genes in that step-up routine...
  //cout << "pathway c++ code -- making note of best genes" << endl;
  THMALLOC2(int, _mask, mask, numGenes, G);
  ///int mask[numGenes][G]; ///
  for(int gene=0; gene<numGenes; gene++){
    //cout << gene << "/" << numGenes << endl;
    // set up the current mask to that gene
    THMALLOC(int, _curmask, curmask, G);
    ///int curmask[G];
    for(int k=0; k<G; k++)
      curmask[k] = (genes[gene] == m[k]);
    // run that gene...
    int local_nperm = 1;
    double local_pvalue = 1;
    int strategy_varselect = 33; // needs to be set so that it variable selects...
    zstat_perm(g, curmask, ng,  aff, naff,  thresh,  gsubsetmatrix, ngsubset,  use_sign, use_weight, &strategy_varselect, &local_nperm, &local_pvalue);
    // copy that mask in...
    for(int k=0; k<G; k++)
      mask[gene][k] = curmask[k];
  }
  //cout << "done with making best genes" << endl;


  /// THE ALTERNATIVE is to allow the pathway to then choose from markers found by stepup in all of the genes...
  //int pathway_mask[G];
  //for(int gene=0; gene<numGenes; gene++){
  //  if(gene==0){
  //    for(int k=0; k<G; k++)
  //      pathway_mask[k] = mask[gene][k];
  //  }else{
  //    for(int k=0; k<G; k++)
  //      pathway_mask[k] = pathway_mask[k] || mask[gene][k];
  //  }
  //}
  //int nnperm = 0;
  //zstat_perm(g, pathway_mask, ng,  aff, naff,  thresh, gsubsetmatrix, ngsubset,  use_sign, use_weight, strategy,  &nnperm,  ret_pvalue);
  //*ret_pvalue = 1 - (*ret_pvalue); // set up for stat, not pvalue, so flip
  //return;
  ///THE_ALTERNATIVE



  // then allow the pathway to choose from markers in any of the genes? how many is this??? or perhaps just a gene at a time <-- YES
  THMALLOC(bool, _geneUsed, geneUsed, numGenes);
  ///bool geneUsed[numGenes];///
  for(int gene=0; gene<numGenes; gene++)
    geneUsed[gene] = false;
  double pvalue = 1.0;
  THMALLOC(int, _finalmask, finalmask, G); /// THIS APPEARS TO BE THE OFFENDING LINE...
  ///int finalmask[G];
  for(int k=0; k<G; k++)
    finalmask[k] = 0;   /// here was the mistake
  bool keepGoing = true;
  while(keepGoing){
    int bestGene = -1;
    double bestGenePvalue = 0; // really a 'stat', _not_ a p-value...
    //cout << "pvalue " << pvalue << endl;

    for(int gene=0; gene<numGenes; gene++){
      if(!geneUsed[gene]){
        //cout << "Gene " << gene << " mask ";
        //for(int k=0; k<G; k++)
        //  cout << mask[gene][k];
        //cout << endl;

        // toggle this gene on (according to it's mask)
        for(int k=0; k<G; k++)
          finalmask[k] = finalmask[k] || mask[gene][k];
        // compute the pvalue
        double newPvalue = 0.0;
        zstat2(g, finalmask, ng,  aff, naff,  use_sign, use_weight,  &newPvalue);
        if(newPvalue > bestGenePvalue){
          //cout << " (current newPvalue (" << newPvalue << ") beat best (" << bestGenePvalue << ")" << endl;
          bestGenePvalue = newPvalue;
          bestGene = gene;
        }
        // toggle this gene off
        for(int k=0; k<G; k++)
          finalmask[k] = finalmask[k] && !mask[gene][k];

        ///cout << "gene " << gene << " pvalue " << newPvalue << endl;
      }
    }

    if(bestGenePvalue > pvalue){
      // set the new pvalue, the gene used, and toggle that gene to be on...
      pvalue = bestGenePvalue;
      geneUsed[bestGene] = 1;
      for(int k=0; k<G; k++)
        finalmask[k] = finalmask[k] || mask[bestGene][k];
    }else{
      keepGoing = false;
    }
  }

  if(print){
    cout << "Genes/masks chosen (gene starts at 0):" << endl;
    for(int gene=0; gene<numGenes; gene++){
      if(geneUsed[gene]){
        cout << "Gene " << gene << ": ";
        for(int k=0; k<G; k++)
          cout << mask[gene][k];
        cout << endl;
      }
    }
    cout << "End of masks chosen." << endl;
  }

  *ret_pvalue = pvalue;
}

extern "C" {
  void zstat_pathway_perm(const double *g, const int *m, const int *ng,
                    const double *aff, const int *naff, // warning naff is the length of aff, _not_ the number of affected
                    const double *thresh, // allele frequency threshhold (for hard)
                    const int *gsubsetmatrix, const int *ngsubset,
                    int *use_sign, int *use_weight, int *strategy, // strategy = 1 (HARD THRESH), 2 (ALL AFREQ), 3 (STEP UP), 4 (ALL SUBSETS)
                    int *nperm,
                    double *ret_pvalue) {
    rndAttach();

    int P = *nperm, I = *naff;

    //cout << "P" << P << " I" << I << endl;
    //cout << "*ng" << *ng << endl;

    THMALLOC(double, _z, z, P+1); //double z[P+1];
    THMALLOC(double, _aff_perm, aff_perm, I); //double aff_perm[I];
    for(int i=0; i<I; i++)
      aff_perm[i] = aff[i];

    for(int p=0; p<=P; p++){
      ///cout << "\r" << (double)p/(double)P*100 << "       ";
      zstat_pathway_stat(g, m, ng,  aff_perm, naff,  thresh,  gsubsetmatrix, ngsubset,  use_sign, use_weight, strategy,  nperm,  &z[p], p==0);
      //cout << "  pathway_stat z[" << p << "] = " << z[p] << endl;
      permute_double(aff, I, aff_perm);
    }

    rndDetach();

    double pvalue = 0.0;
    for(int p=1; p<=P; p++)
      pvalue += (int)(z[p] >= z[0]);
    pvalue /= (double)P;
    *ret_pvalue = pvalue;
    //cout << "zstat_pathway_perm pvalue " << pvalue << endl;
  }
}

// next todo is to envelope this within a routine that also does permutation... then later an adjustable permutation???
