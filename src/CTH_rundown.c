/*
 * This rundown function for CTH.
 *
 */
#include "IVTree.h"
#include "node.h"
#include "IVTreeproto.h"

#ifdef NAN
/* NAN is supported */
#endif

void
CTH_rundown(pNode tree, int obs, double *cp, double *xpred, double *xtemp, int k, double alpha, 
            double xtrain_to_est_ratio, double propensity)
{
    int i, obs2 = (obs < 0) ? -(1 + obs) : obs;
    int my_leaf_id;
    pNode otree =  tree;
    pNode otree_tmp = tree;
    pNode tree_tmp = tree;
    
    int opnumber = 0;
    int j, s;
    int tmp_obs, tmp_id;
    double tr_mean, con_mean;
    double tr_sqr_sum, con_sqr_sum;
    double consums, trsums, cons, trs;
    // double tr_var, con_var;
    double xz_sum, xy_sum, x_sum, y_sum, z_sum;
    double yz_sum, xx_sum, yy_sum, zz_sum;
    int n;

    /*
     * Now, repeat the following: for the cp of interest, run down the tree
     *   until I find a node with smaller complexity.  The parent node will
     *   not have collapsed, but this split will have, so this is my
     *   predictor.
     */
    for (i = 0; i < ct.num_unique_cp; i++) {
        cons = 0.;
        trs = 0.;
        consums = 0.;
        trsums = 0.;
        tr_sqr_sum = 0.;
        con_sqr_sum = 0.;

    	n = 0;
    	xz_sum = 0.;
    	xy_sum = 0.;
    	x_sum = 0.;
    	y_sum = 0.;
    	z_sum = 0.;
        yz_sum = 0.;
    	xx_sum = 0.;
    	yy_sum = 0.;
    	zz_sum = 0.;
        
        while (cp[i] < tree->complexity) {
	        tree = branch(tree, obs);
	        if (tree == 0)
		        goto oops;
	        otree = tree;
	    }
	    xpred[i] = tree->response_est[0];
        my_leaf_id = tree->id;
        
        for (s = k; s < ct.n; s++) {
            tree_tmp = otree_tmp;
            j = ct.sorts[0][s];
            tmp_obs = (j < 0) ? -(1 + j) : j;
            while (cp[i] < tree_tmp->complexity) {
                tree_tmp = branch(tree_tmp, tmp_obs);
            }
            tmp_id = tree_tmp->id;
            if (tmp_id == my_leaf_id) {
                if (ct.treatment[tmp_obs] == 0) {
                    cons += ct.wt[tmp_obs];
                    consums += *ct.ydata[tmp_obs] * ct.wt[tmp_obs];
                    con_sqr_sum += (*ct.ydata[tmp_obs]) * (*ct.ydata[tmp_obs]) * ct.wt[tmp_obs];
                } else {
                    trs += ct.wt[tmp_obs];
                    trsums += *ct.ydata[tmp_obs] * ct.wt[tmp_obs];
                    tr_sqr_sum += (*ct.ydata[tmp_obs]) * (*ct.ydata[tmp_obs]) * ct.wt[tmp_obs];
                }
                
                ++n;
                xz_sum += ct.IV[tmp_obs] * *ct.ydata[tmp_obs];
                xy_sum += ct.IV[tmp_obs] * ct.treatment1[tmp_obs];
                x_sum += ct.IV[tmp_obs];
                y_sum += ct.treatment1[tmp_obs];
                z_sum += *ct.ydata[tmp_obs];
                yz_sum += *ct.ydata[tmp_obs] * ct.treatment1[tmp_obs];
                xx_sum += ct.IV[tmp_obs] * ct.IV[tmp_obs];
                yy_sum += ct.treatment1[tmp_obs] * ct.treatment1[tmp_obs];
                zz_sum += *ct.ydata[tmp_obs] * *ct.ydata[tmp_obs];
            }
        }

        if (trs == 0) {
            tr_mean = tree->parent->xtreatMean[0];
        } else {
            tr_mean = trsums / trs;
            tree->xtreatMean[0] = tr_mean;
        }
        
        if (cons == 0) {
            con_mean = tree->parent->xcontrolMean[0];
        } else {
            con_mean = consums / cons;
            tree->xcontrolMean[0] = con_mean;
        }
        
        double alpha_1;
// PARAMETER!
        double b1_hat = (n * xy_sum - x_sum * y_sum) / (n * xx_sum - x_sum * x_sum);
        double y_mean = y_sum / n;
        double b0_hat = y_mean - b1_hat * x_sum / n;
        double MSM_temp = b1_hat * b1_hat * xx_sum + 2.0 * b0_hat * b1_hat * x_sum - 2.0 * n * b0_hat * y_mean 
                          - 2.0 * b1_hat * y_mean * x_sum + n * b0_hat * b0_hat + n * y_mean * y_mean;
        double MSE_temp = n * b0_hat * b0_hat + 2.0 * b0_hat * b1_hat * y_sum - 2.0 * b0_hat * y_sum
                          - 2.0 * b1_hat * yy_sum + b1_hat * b1_hat * yy_sum + yy_sum;
        MSE_temp = MSE_temp / (n-2); 
        if(MSM_temp / MSE_temp < F_test_threshold){
            Rprintf("Entered CTH_rundown (a week IV).\n"); 
            alpha_1 = (n * yz_sum - y_sum * z_sum) / (n * yy_sum - y_sum * y_sum);
            
            double effect = alpha_1;
            double alpha_0 = (z_sum - alpha_1 * y_sum) / n;

            double mu = zz_sum - 2.0 * alpha_0 * z_sum - 2.0 * alpha_1 * yz_sum
                        + 2.0 * alpha_0 * alpha_1 * y_sum + n * alpha_0 * alpha_0 + alpha_1 * alpha_1 * yy_sum;
            mu = mu / (n-2);
            double var_alpha_1_est = mu / (yy_sum - y_sum * y_sum / n);

            xtemp[i] = 4 * ct.max_y * ct.max_y - alpha * effect * effect + (1 + xtrain_to_est_ratio / (ct.NumXval - 1)) * (1 - alpha) * var_alpha_1_est;
        } 
        else{
            alpha_1 = (n * xz_sum - x_sum * z_sum) / (n * xy_sum - x_sum * y_sum);

            double effect = alpha_1;
            double alpha_0 = (z_sum - alpha_1 * y_sum) / n;

            double numerator = (zz_sum + n * alpha_0 * alpha_0 + alpha_1 * alpha_1 * yy_sum - 2 * alpha_0 * z_sum - 2 * alpha_1 * yz_sum + 2 * alpha_0 * alpha_1 * y_sum)/n;
            double denominator =  1 / (xx_sum / n - (x_sum / n) * (x_sum / n)) * (xy_sum / n - x_sum/n * y_sum / n) * (xy_sum / n - x_sum/n * y_sum / n) * n;   
            double tmp;
            if (n > 2 && denominator!=0) {
                    tmp = numerator / denominator / (n - 2);
            } else {
                tmp = 0.;
            }  
            
            xtemp[i] = 4 * ct.max_y * ct.max_y - alpha * effect * effect + (1 + xtrain_to_est_ratio / (ct.NumXval - 1)) * (1 - alpha) * tmp;
        } 
    	// if (fabs(n * xy_sum - x_sum * y_sum) <= 0.1 * n * n){
    	// 	alpha_1 = 0.;
    	// }
    	// else{
    	// 	alpha_1 = (n * xz_sum - x_sum * z_sum) / (n * xy_sum - x_sum * y_sum);
    	// }
        // double effect = alpha_1;
        // double alpha_0 = (z_sum - alpha_1 * y_sum) / n;

        // double numerator = (zz_sum + n * alpha_0 * alpha_0 + alpha_1 * alpha_1 * yy_sum - 2 * alpha_0 * z_sum - 2 * alpha_1 * yz_sum + 2 * alpha_0 * alpha_1 * y_sum)/n;
        // double denominator =  1 / (xx_sum / n - (x_sum / n) * (x_sum / n)) * (xy_sum / n - x_sum/n * y_sum / n) * (xy_sum / n - x_sum/n * y_sum / n) * n;   
        // double tmp;
        // if (n > 2 && denominator!=0) {
        //         tmp = numerator / denominator / (n - 2);
        // } else {
        //     tmp = 0.;
        // }  
        
        // xtemp[i] = 4 * ct.max_y * ct.max_y - alpha * effect * effect + (1 + xtrain_to_est_ratio / (ct.NumXval - 1)) * (1 - alpha) * tmp;
    }
    return;

oops:;
    if (ct.usesurrogate < 2) {  /* must have hit a missing value */
    	Rprintf("Entered CTH_rundown.c. Double check.\n");
    	for (i = 0; i < ct.num_unique_cp; i++)
    	    xpred[i] = otree->response_est[0];

    	xtemp[i] = (*ct_xeval)(ct.ydata[obs2], ct.wt[obs2], ct.treatment[obs2], tr_mean, con_mean);
    	Rprintf("oops number %d.\n", opnumber++);
    return;
    }
    warning("Warning message--see rundown.c");
}
