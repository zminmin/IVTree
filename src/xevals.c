/*
 * The cross validation evaluation function
 */
#include<stdio.h>
#include <math.h>
#include "IVTree.h"
#include "IVTreeproto.h"

double CTH_xpred(double *y, double wt, double treatment, double tr_mean,
                 double con_mean, double trs, double cons, double alpha, 
                 double xtrain_to_est_ratio, double propensity) {
   double res;
   double tr_var;
   double con_var;
   double tmp;
   //Rprintf("Entered xeval\n");
   if (treatment == 0) {
       // con
       con_var = wt * (y[0] - con_mean) *  (y[0] - con_mean);
       tmp = con_var / ((1 - propensity) * cons);
   } else {
       // tr
       tr_var = wt * (y[0] - tr_mean) * (y[0] - tr_mean);
       tmp = tr_var / (propensity * trs);
   } 
   double effect = tr_mean - con_mean;
   
   res = 4 * ct.max_y * ct.max_y - alpha *  effect * effect + (1 + xtrain_to_est_ratio / (ct.NumXval - 1)) 
       * (1 - alpha) *  tmp; 
   
   return res;
}

double CTA_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                     double tree_tr_mean, double tree_con_mean, double alpha) {
    double res;
    double effect_tr = tree_tr_mean - tree_con_mean;
    double effect_te = tr_mean - con_mean;
    res = 2 * ct.max_y * ct.max_y + effect_tr * effect_tr  -  2 *  effect_tr * effect_te;

    return res;
}
