#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix new_func_sd(NumericMatrix X) {
    int m = X.nrow();
    int d = X.ncol();
    NumericMatrix M(m, m);
    NumericVector vec(2*d);
    NumericVector x1(d);
    NumericVector x2(d);

    M(m-1, m-1) = 1;
    for (int i= 0; i < (m-1); i++) {
        M(i, i) = 1;
        for (int j = (i+1); j < m; j++) {
            NumericVector vec(2*d);
            for (int k = 0; k < d; k++) {
                vec[k] = x1[k] = X(i, k);
                vec[d+ k] = x2[k] = X(j, k);
            }

            std::sort( vec.begin(), vec.end() );
            vec.erase( std::unique( vec.begin(), vec.end() ), vec.end() );

            int nobs = vec.size();
            NumericVector F1(nobs);
            NumericVector F2(nobs);
            std::sort(x1.begin(), x1.end());
            std::sort(x2.begin(), x2.end());
            for (int l = 0; l < nobs; ++l) {
                F1[l] =  (std::upper_bound(x1.begin(), x1.end(), vec[l]) - x1.begin());
                F2[l] =  (std::upper_bound(x2.begin(), x2.end(), vec[l]) - x2.begin());
            }

            F1 = F1 / ((double) nobs);
            F2 = F2 / ((double) nobs);
            int test = 1;
            for (int k = 0; k< (nobs-1); k++) {
                if (F2[k] > F1[k]) {
                    test = 0;
                }
            }
            if (test == 1) {
                M(i, j) = 1;
            }
            if (M(i, j) == 0) {
                int test2 = 1;
                for (int k = 0; k< (nobs-1); k++) {
                    if (F1[k] > F2[k]) {
                        test2 = 0;
                    }
                }

                if (test2 == 1) {
                    M(j, i) = 1;
                }
            }

        }
    }
    return M;

}


// [[Rcpp::export]]
List new_func_eps(NumericMatrix X) {
    int m = X.nrow();
    int d = X.ncol();
    NumericMatrix vec_check(m*m, 2);
    NumericMatrix sums_mat(m*m, 3);
    NumericVector vec(2*d);
    NumericVector x1(d);
    NumericVector x2(d);
    NumericMatrix M(m, m);

    M(m-1, m-1) = 1;
    int sd = 0;
    int count = 0;
    for (int i= 0; i < (m-1); i++) {
        M(i, i) = 1;
        for (int j = (i+1); j < m; j++) {
            NumericVector vec(2*d);
            for (int k = 0; k < d; k++) {
                vec[k] = x1[k] = X(i, k);
                vec[d+ k] = x2[k] = X(j, k);
            }

            std::sort( vec.begin(), vec.end() );
            vec.erase( std::unique( vec.begin(), vec.end() ), vec.end() );

            int nobs = vec.size();
            NumericVector F1(nobs);
            NumericVector F2(nobs);
            std::sort(x1.begin(), x1.end());
            std::sort(x2.begin(), x2.end());
            int d2 = 1;
            int d3 = 1;
            for (int l = 0; l < nobs; ++l) {
                F1[l] =  (std::upper_bound(x1.begin(), x1.end(), vec[l]) - x1.begin());
                F2[l] =  (std::upper_bound(x2.begin(), x2.end(), vec[l]) - x2.begin());
                if (F1[l] < F2[l]){
                    d2 = 0;
                }
                if (F1[l] > F2[l]){
                    d3 = 0;
                }

            }
            sd = sd + d2 + d3;
            M(i, j) = d2;
            M(j, i) = d3;

            if (d2+d3 == 0){
                F1 = F1 / ((double) nobs);
                F2 = F2 / ((double) nobs);
                double sum_s_left = 0;
                double sum_s_right = 0;
                double sum_all = 0;
                for (int k = 0; k< (nobs-1); k++) {
                    sum_all += std::abs(F1[k] - F2[k])*(vec[k+1] - vec[k]);
                    if (F2[k] > F1[k]) {
                        sum_s_left += (F2[k] - F1[k])*(vec[k+1] - vec[k]);
                    }
                    if (F1[k] > F2[k]) {
                        sum_s_right += (F1[k] - F2[k])*(vec[k+1] - vec[k]);
                    }
                }
                vec_check(count, 0) = i;
                vec_check(count, 1) = j;
                sums_mat(count, 0) = sum_all;
                sums_mat(count, 1) = sum_s_left;
                sums_mat(count, 2) = sum_s_right;
                count = count+1;
            }
        }
    }

    List list_d_mat;
    list_d_mat["d"] = sd;
    list_d_mat["k"] = count;
    list_d_mat["mat"] = sums_mat;
    list_d_mat["M"] = M;
    list_d_mat["indx"] = vec_check;
    return list_d_mat;
}


// [[Rcpp::export]]
List normal_comp_eps(NumericMatrix X) {
    int m = X.nrow();

    NumericMatrix vec_check(m*m, 2);
    NumericMatrix sums_mat(m*m, 5);

    NumericMatrix M(m, m);

    M(m-1, m-1) = 1;
    int sd = 0;
    int count = 0;
    for (int i= 0; i < (m-1); i++) {
        M(i, i) = 1;
        for (int j = (i+1); j < m; j++) {
            if (std::abs(X(i, 1) - X(j, 1)) < 1e-9){
                if (X(i, 0) >= X(j, 0)){
                    sd = sd + 1;
                    M(j, i) = 1;
                }
                if (X(j, 0) >= X(i, 0)) {
                    M(i, j) = 1;
                    sd = sd + 1;
                }
            }
            else {
                if (X(i, 1) < X(j, 1)) {
                    static const float inv_sqrt_2pi = 0.3989422804014327;
                    double sigmu1;
                    double sigmu2;
                    sigmu1 = X(i, 0)*X(j, 1);
                    sigmu2 = X(j, 0)*X(i, 1);
                    double v = (sigmu1 - sigmu2) / (X(j, 1) - X(i, 1));

                    double a = (v - X(i,0)) / X(i, 1);
                    double b = (v - X(j,0)) / X(j, 1);
                    double g2 = inv_sqrt_2pi / X(j, 1) * std::exp(-0.5f * b * b);
                    double f1 = inv_sqrt_2pi / X(i, 1) * std::exp(-0.5f * a * a);
                    double G2 = R::pnorm(v,  X(j, 0),  X(j, 1), 1, 0);
                    double F1 = R::pnorm(v,  X(i, 0), X(i, 1), 1, 0);
                    double sum_s_left = (X(j, 1)*X(j, 1)*g2- X(i, 1)*X(i, 1)*f1);

                    vec_check(count, 0) = i;
                    vec_check(count, 1) = j;
                    sums_mat(count, 0) = sum_s_left;
                    sums_mat(count, 1) = (X(i, 0) - X(j, 0));
                    sums_mat(count, 2) = F1;
                    sums_mat(count, 3) = G2;
                    sums_mat(count, 4) = 0;
                    count = count+1;
                }
                else if (X(i, 1) > X(j, 1)) {
                    static const float inv_sqrt_2pi = 0.3989422804014327;
                    double sigmu1;
                    double sigmu2;
                    sigmu1 = X(i, 0)*X(j, 1);
                    sigmu2 = X(j, 0)*X(i, 1);
                    double v = (sigmu1 - sigmu2) / (X(j, 1) - X(i, 1));
                    double a = (v - X(i,0)) / X(i, 1);
                    double b = (v - X(j,0)) / X(j, 1);
                    double g2 = inv_sqrt_2pi / X(j, 1) * std::exp(-0.5f * b * b);
                    double f1 = inv_sqrt_2pi / X(i, 1) * std::exp(-0.5f * a * a);
                    double G2 = R::pnorm(v,  X(j, 0),  X(j, 1), 1, 0);
                    double F1 = R::pnorm(v, X(i, 0), X(i, 1), 1, 0);
                    double sum_s_left = (X(i, 1)*X(i, 1)*f1- X(j, 1)*X(j, 1)*g2);

                    vec_check(count, 0) = i;
                    vec_check(count, 1) = j;
                    sums_mat(count, 0) = sum_s_left;
                    sums_mat(count, 1) = (X(i, 0) - X(j, 0));
                    sums_mat(count, 2) = F1;
                    sums_mat(count, 3) = G2;
                    sums_mat(count, 4) = 1;
                    count = count+1;
                }

            }
        }
    }
    List list_d_mat;
    list_d_mat["d"] = sd;
    list_d_mat["k"] = count;
    list_d_mat["mat"] = sums_mat;
    list_d_mat["M"] = M;
    list_d_mat["indx"] = vec_check;
    return list_d_mat;

}

// [[Rcpp::export]]
NumericMatrix normal_comp_sd(NumericMatrix X) {
    int m = X.nrow();
    NumericMatrix M(m, m);
    M(m-1, m-1) = 1;
    for (int i= 0; i < (m-1); i++) {
        M(i, i) = 1;
        for (int j = (i+1); j < m; j++) {
            if (std::abs(X(i, 1) - X(j, 1)) < 1e-9){
                if (X(i, 0) >= X(j, 0)){
                    M(j, i) = 1;
                }
                if (X(j, 0) >= X(i, 0)) {
                    M(i, j) = 1;
                }
            }
        }
    }
    return M;
}

// [[Rcpp::export]]
NumericMatrix normal_comp_sd_ab(NumericMatrix X, double a, double b) {
    int m = X.nrow();
    NumericMatrix M(m, m);
    M(m-1, m-1) = 1;
    for (int i= 0; i < (m-1); i++) {
        M(i, i) = 1;
        for (int j = (i+1); j < m; j++) {
            if (std::abs(X(i, 1) - X(j, 1)) < 1e-9){
                if (X(i, 0) >= X(j, 0)){
                    M(j, i) = 1;
                }
                if (X(j, 0) >= X(i, 0)) {
                    M(i, j) = 1;
                }
            }
            else {
                double diff = X(j, 1) - X(i, 1);
                double num = X(i, 0)*X(j, 1) - X(j, 0)*X(i, 1);
                double test = num / diff;

                if (test < a){
                    if (X(i, 0) >= X(j, 0)){
                        M(j, i) = 1;
                    }
                    if (X(j, 0) >= X(i, 0)) {
                        M(i, j) = 1;
                    }
                }
                if (test > b){
                    if (X(i, 0) >= X(j, 0)){
                        M(j, i) = 1;
                    }
                    if (X(j, 0) >= X(i, 0)) {
                        M(i, j) = 1;
                    }
                }
            }
        }
    }
    return M;
}



// [[Rcpp::export]]
List indx_norm_sd(NumericMatrix X, NumericMatrix x) {
    int mX = X.nrow();
    int mx = x.nrow();

    IntegerMatrix smaller_indx(mx, mX);
    IntegerMatrix greater_indx(mx, mX);

    List ret;


    for (int i = 0; i < mx; i++) {
        double sigma1;
        double mu1;
        sigma1 = x(i, 1);
        mu1 = x(i, 0);
        for (int j = 0; j < mX; j++) {
            double sigma2;
            double mu2;
            sigma2 = X(j,1);
            mu2 = X(j,0);

            if (std::abs(sigma2 - sigma1) < 1e-9){
                if (mu1 >= mu2){
                    smaller_indx(i, j) = 1;
                }
                if (mu2 >= mu1){
                    greater_indx(i, j) = 1;
                }
            }
        }
    }
    ret["smaller"] = smaller_indx;
    ret["greater"] = greater_indx;
    return ret;
}

// [[Rcpp::export]]
List new_func2_sd(NumericMatrix X, NumericMatrix x) {
    int mX = X.nrow();
    int d = X.ncol();
    int mx = x.nrow();

    IntegerMatrix smaller_indx(mx, mX);
    IntegerMatrix greater_indx(mx, mX);
    NumericVector vec(2*d);
    NumericVector x1(d);
    NumericVector x2(d);
    List ret;


    for (int i = 0; i < mx; i++) {
        for (int j = 0; j < mX; j++) {
            NumericVector vec(2*d);
            for (int k = 0; k < d; k++) {
                vec[k] = x1[k] = x(i, k);
                vec[k + d] = x2[k] = X(j, k);
            }
            std::sort( vec.begin(), vec.end() );
            vec.erase( std::unique( vec.begin(), vec.end() ), vec.end() );

            int nobs = vec.size();
            NumericVector F1(nobs);
            NumericVector F2(nobs);
            std::sort(x1.begin(), x1.end());
            std::sort(x2.begin(), x2.end());
            for (int l = 0; l < nobs; ++l) {
                F1[l] =  (std::upper_bound(x1.begin(), x1.end(), vec[l]) - x1.begin());
                F2[l] =  (std::upper_bound(x2.begin(), x2.end(), vec[l]) - x2.begin());
            }

            F1 = F1 / ((double) nobs);
            F2 = F2 / ((double) nobs);
            int test = 1;
            for (int k = 0; k< (nobs-1); k++) {

                if (F2[k] > F1[k]) {
                    test = 0;
                }
            }
            if (test == 1) {
                greater_indx(i, j) = 1;
            }

            int test2 = 1;
            for (int k = 0; k< (nobs-1); k++) {
                if (F1[k] > F2[k]) {
                    test2 = 0;
                }
            }
            if (test2 == 1) {
                smaller_indx(i, j) = 1;
            }
        }
    }

    ret["smaller"] = smaller_indx;
    ret["greater"] = greater_indx;
    return ret;
}


// [[Rcpp::export]]
List new_func_single_grid_sd(NumericMatrix X, NumericMatrix x, NumericVector vec) {
    int mX = X.nrow();
    int d = X.ncol();
    int mx = x.nrow();

    IntegerMatrix smaller_indx(mx, mX);
    IntegerMatrix greater_indx(mx, mX);

    List ret;

    for (int i = 0; i < mx; i++){
        for (int j = 0; j < mX; j++){
            int test = 1;
            for (int k = 0; k< (d-1); k++) {

                if (X(j, k) > x(i, k)) {
                    test = 0;
                }
            }
            if (test == 1) {
                greater_indx(i, j) = 1;
            }

            int test2 = 1;
            for (int k = 0; k< (d-1); k++) {
                if (x(i, k) > X(j, k)) {
                    test2 = 0;
                }
            }
            if (test2 == 1) {
                smaller_indx(i, j) = 1;
            }

        }
    }

    ret["smaller"] = smaller_indx;
    ret["greater"] = greater_indx;
    return ret;
}

// [[Rcpp::export]]
List new_func_mat_sd(NumericMatrix X, NumericMatrix x, NumericVector gridx, NumericVector gridX) {
    int mX = X.nrow();
    int dX = X.ncol();
    int dx = x.ncol();
    int mx = x.nrow();
    int dXx = dX + dx;

    IntegerMatrix smaller_indx(mx, mX);
    IntegerMatrix greater_indx(mx, mX);

    List ret;
    NumericVector vec(dXx);
    for (int k = 0; k < dx; k++) {
        vec[k] = gridx[k];
    }
    for (int k = dx; k < dXx; k++){
        vec[k] = gridX[k-dx];
    }

    std::sort( vec.begin(), vec.end() );
    vec.erase( std::unique( vec.begin(), vec.end() ), vec.end() );
    int d = vec.size();

    NumericMatrix F1(mx, d);
    NumericMatrix F2(mX, d);

    for(int k = 0; k < d; k++){
        int indx1 = 0;
        int indx2 = 0;
        while (vec[k] >= gridx[indx1] && indx1 < dx){
            indx1 +=1;
        }
        if (indx1 == 0){
            for (int i = 0; i < mx; i++){
                F1(i, k) = 0;
            }
        } else {
            for (int i = 0; i < mx; i++){
                F1(i, k) = x(i,(indx1-1));
            }
        }
        while (vec[k] >= gridX[indx2] && indx2 < dX){
            indx2 +=1;
        }
        if (indx2 == 0){
            for (int i = 0; i < mX; i++){
                F2(i, k) = 0;
            }
        } else {
            for (int i = 0; i < mX; i++){
                F2(i, k) = X(i, (indx2-1));
            }
        }
    }

    for (int i = 0; i < mx; i++){
        for (int j = 0; j < mX; j++){
            int test = 1;
            for (int k = 0; k< (d-1); k++) {
                if (F2(j, k) > F1(i, k)) {
                    test = 0;
                }
            }
            if (test == 1) {
                greater_indx(i, j) = 1;
            }

            int test2 = 1;
            for (int k = 0; k< (d-1); k++) {
                if (F1(i, k) > F2(j, k)) {
                    test2 = 0;
                }
            }
            if (test2 == 1) {
                smaller_indx(i, j) = 1;
            }
        }
    }

    ret["smaller"] = smaller_indx;
    ret["greater"] = greater_indx;
    return ret;
}


// [[Rcpp::export]]
List new_func_mat_list_sd(List X, NumericMatrix x, NumericVector gridx, List gridX) {
    int mX = X.size();
    int mx = x.nrow();
    int t_l1 = gridx.size();

    IntegerMatrix smaller_indx(mx, mX);
    IntegerMatrix greater_indx(mx, mX);

    List ret;

    for (int i = 0; i < mx; i++){
        for (int j = 0; j < mX; j++){
            NumericVector t2 = as<NumericVector>(gridX[j]);
            NumericVector x2 = as<NumericVector>(X[j]);
            int t_l2 = t2.size();
            int dim_thresh = t_l1 + t_l2;
            NumericVector thresholds(dim_thresh);
            for (int l = 0; l < t_l1; l++){
                thresholds[l] = gridx[l];
            }
            for (int l = t_l1; l < dim_thresh; l++){
                thresholds[l] = t2[(l-t_l1)];
            }
            std::sort( thresholds.begin(), thresholds.end() );
            thresholds.erase( std::unique( thresholds.begin(), thresholds.end() ), thresholds.end() );
            int d = thresholds.size();

            NumericVector F1(d);
            NumericVector F2(d);

            for(int k = 0; k < d; k++){
                int indx1 = 0;
                int indx2 = 0;
                while (thresholds[k] >= gridx[indx1] && indx1 < t_l1){
                    indx1 +=1;
                }
                if (indx1 == 0){
                    F1[k] = 0;
                } else {
                    F1[k] = x(i, (indx1-1));
                }
                while (thresholds[k] >= t2[indx2] && indx2 < t_l2){
                    indx2 +=1;
                }
                if (indx2 == 0){
                    F2[k] = 0;
                } else {
                    F2[k] = x2[indx2-1];
                }
            }

            int test = 1;
            for (int k = 0; k< (d-1); k++) {
                if (F2[k] > F1[k]) {
                    test = 0;
                }
            }
            if (test == 1) {
                greater_indx(i, j) = 1;
            }

            int test2 = 1;
            for (int k = 0; k< (d-1); k++) {
                if (F1[k] > F2[k]) {
                    test2 = 0;
                }
            }
            if (test2 == 1) {
                smaller_indx(i, j) = 1;
            }
        }
    }

    ret["smaller"] = smaller_indx;
    ret["greater"] = greater_indx;
    return ret;
}


// [[Rcpp::export]]
List new_func_list_sd_back(List X, List x, List gridx, List gridX) {
    int mX = X.size();
    int mx = x.size();


    IntegerMatrix smaller_indx(mx, mX);
    IntegerMatrix greater_indx(mx, mX);

    List ret;

    for (int i = 0; i < mx; i++){

        NumericVector t1 = as<NumericVector>(gridx[i]);
        NumericVector x1 = as<NumericVector>(x[i]);
        int t_l1 = t1.size();
        for (int j = 0; j < mX; j++){
            NumericVector t2 = as<NumericVector>(gridX[j]);
            NumericVector x2 = as<NumericVector>(X[j]);

            int t_l2 = t2.size();
            int dim_thresh = t_l1 + t_l2;
            NumericVector thresholds(dim_thresh);
            for (int l = 0; l < t_l1; l++){
                thresholds[l] = t1[l];
            }
            for (int l = t_l1; l < dim_thresh; l++){
                thresholds[l] = t2[(l-t_l1)];
            }
            std::sort( thresholds.begin(), thresholds.end() );
            thresholds.erase( std::unique( thresholds.begin(), thresholds.end() ), thresholds.end() );
            int d = thresholds.size();

            NumericVector F1(d);
            NumericVector F2(d);

            for(int k = 0; k < d; k++){
                int indx1 = 0;
                int indx2 = 0;
                while (thresholds[k] >= t1[indx1] && indx1 < t_l1){
                    indx1 +=1;
                }
                if (indx1 == 0){
                    F1[k] = 0;
                } else {
                    F1[k] = x1[indx1-1];
                }
                while (thresholds[k] >= t2[indx2] && indx2 < t_l2){
                    indx2 +=1;
                }
                if (indx2 == 0){
                    F2[k] = 0;
                } else {
                    F2[k] = x2[indx2-1];
                }
            }

            int test = 1;
            for (int k = 0; k< (d-1); k++) {
                if (F2[k] > F1[k]) {
                    test = 0;
                }
            }

            if (test == 1) {
                greater_indx(i, j) = 1;
            }

            int test2 = 1;
            for (int k = 0; k< (d-1); k++) {

                if (F1[k] > F2[k]) {
                    test2 = 0;
                }
            }
            if (test2 == 1) {
                smaller_indx(i, j) = 1;
            }
        }
    }

    ret["smaller"] = smaller_indx;
    ret["greater"] = greater_indx;
    return ret;
}



// [[Rcpp::export]]
List new_func_list_sd(List X, List x, List gridx, List gridX) {
    int mX = X.size();
    int mx = x.size();

    IntegerMatrix smaller_indx(mx, mX);
    IntegerMatrix greater_indx(mx, mX);

    List ret;

    for (int i = 0; i < mx; i++) {
        NumericVector x1 = as<NumericVector>(x[i]);
        NumericVector t1 = as<NumericVector>(gridx[i]);
        int t_l1 = t1.size();

        for (int j = 0; j < mX; j++) {
            NumericVector x2 = as<NumericVector>(X[j]);
            NumericVector t2 = as<NumericVector>(gridX[j]);
            int t_l2 = t2.size();
            int dim_thresh = t_l1 + t_l2;
            NumericVector thresholds(dim_thresh);

            for (int l = 0; l < t_l1; l++) {
                thresholds[l] = t1[l];
            }
            for (int l = t_l1; l < dim_thresh; l++) {
                thresholds[l] = t2[l - t_l1];
            }

            std::sort(thresholds.begin(), thresholds.end());
            thresholds.erase(std::unique(thresholds.begin(), thresholds.end()), thresholds.end());
            int d = thresholds.size();

            NumericVector F1(d);
            NumericVector F2(d);

            for (int k = 0; k < d; k++) {
                int indx1 = 0;
                int indx2 = 0;

                while (indx1 < t_l1 && thresholds[k] >= t1[indx1]) {
                    indx1 += 1;
                }

                if (indx1 == 0) {
                    F1[k] = 0;
                } else {
                    F1[k] = x1[indx1 - 1];
                }

                while (indx2 < t_l2 && thresholds[k] >= t2[indx2]) {
                    indx2 += 1;
                }

                if (indx2 == 0) {
                    F2[k] = 0;
                } else {
                    F2[k] = x2[indx2 - 1];
                }
            }

            int test = 1;
            for (int k = 0; k < (d - 1); k++) {
                if (F2[k] > F1[k]) {
                    test = 0;
                }
            }

            if (test == 1) {
                greater_indx(i, j) = 1;
            }

            int test2 = 1;
            for (int k = 0; k < (d - 1); k++) {
                if (F1[k] > F2[k]) {
                    test2 = 0;
                }
            }

            if (test2 == 1) {
                smaller_indx(i, j) = 1;
            }
        }
    }

    ret["smaller"] = smaller_indx;
    ret["greater"] = greater_indx;
    return ret;
}


// [[Rcpp::export]]
List new_func_list_mat_sd(NumericMatrix X, List x, List gridx, NumericVector gridX) {
    int mX = X.nrow();
    int mx = x.size();

    int t_l2 = gridX.size();
    IntegerMatrix smaller_indx(mx, mX);
    IntegerMatrix greater_indx(mx, mX);

    List ret;

    for (int i = 0; i < mx; i++){

        NumericVector t1 = as<NumericVector>(gridx[i]);
        NumericVector x1 = as<NumericVector>(x[i]);
        int t_l1 = t1.size();
        for (int j = 0; j < mX; j++){

            int dim_thresh = t_l1 + t_l2;
            NumericVector thresholds(dim_thresh);
            for (int l = 0; l < t_l1; l++){
                thresholds[l] = t1[l];
            }
            for (int l = t_l1; l < dim_thresh; l++){
                thresholds[l] = gridX[(l-t_l1)];
            }
            std::sort( thresholds.begin(), thresholds.end() );
            thresholds.erase( std::unique( thresholds.begin(), thresholds.end() ), thresholds.end() );
            int d = thresholds.size();

            NumericVector F1(d);
            NumericVector F2(d);

            for(int k = 0; k < d; k++){
                int indx1 = 0;
                int indx2 = 0;
                while (thresholds[k] >= t1[indx1] && indx1 < t_l1){
                    indx1 +=1;
                }
                if (indx1 == 0){
                    F1[k] = 0;
                } else {
                    F1[k] = x1[indx1-1];
                }
                while (thresholds[k] >= gridX[indx2] && indx2 < t_l2){
                    indx2 +=1;
                }
                if (indx2 == 0){
                    F2[k] = 0;
                } else {
                    F2[k] = X(j, (indx2-1));
                }
            }

            int test = 1;
            for (int k = 0; k< (d-1); k++) {

                if (F2[k] > F1[k]) {
                    test = 0;
                }
            }
            if (test == 1) {
                greater_indx(i, j) = 1;
            }

            int test2 = 1;
            for (int k = 0; k< (d-1); k++) {

                if (F1[k] > F2[k]) {
                    test2 = 0;
                }
            }
            if (test2 == 1) {
                smaller_indx(i, j) = 1;
            }
        }
    }

    ret["smaller"] = smaller_indx;
    ret["greater"] = greater_indx;
    return ret;
}


// [[Rcpp::export]]
List ecdf_list_comp_class_sd(List X, List t) {
    int m = X.size();

    NumericMatrix M(m, m);
    NumericVector classY(m);
    int class_count = 1;
    M(m-1, m-1) = 1;

    for (int i= 0; i < (m-1); i++) {
        M(i, i) = 1;
        bool class_check = false;
        if (classY[i] == 0){
            classY[i] = class_count;
            class_count += 1;
            class_check = true;
        }
        NumericVector t1 = as<NumericVector>(t[i]);
        NumericVector x1 = as<NumericVector>(X[i]);
        int t_l1 = t1.size();
        for (int j = (i+1); j < m; j++) {
            NumericVector t2 = as<NumericVector>(t[j]);
            NumericVector x2 = as<NumericVector>(X[j]);

            int t_l2 = t2.size();
            int dim_thresh = t_l1 + t_l2;
            NumericVector thresholds(dim_thresh);
            for (int l = 0; l < t_l1; l++){
                thresholds[l] = t1[l];
            }
            for (int l = t_l1; l < dim_thresh; l++){
                thresholds[l] = t2[(l-t_l1)];
            }

            std::sort( thresholds.begin(), thresholds.end() );
            thresholds.erase( std::unique( thresholds.begin(), thresholds.end() ), thresholds.end() );
            int d = thresholds.size();

            NumericVector F1(d);
            NumericVector F2(d);

            int d2 = 1;
            int d3 = 1;

            bool check_equal = true;

            for (int k = 0; k < d; k++) {
                int indx1 = 0;
                int indx2 = 0;

                while (indx1 < t_l1 && thresholds[k] >= t1[indx1]) {
                    indx1 += 1;
                }
                if (indx1 == 0) {
                    F1[k] = 0;
                } else {
                    F1[k] = x1[indx1 - 1];
                }

                while (indx2 < t_l2 && thresholds[k] >= t2[indx2]) {
                    indx2 += 1;
                }
                if (indx2 == 0) {
                    F2[k] = 0;
                } else {
                    F2[k] = x2[indx2 - 1];
                }

                if (F1[k] != F2[k]) {
                    check_equal = false;
                }
                if (F1[k] < F2[k]) {
                    d2 = 0;
                }
                if (F1[k] > F2[k]) {
                    d3 = 0;
                }
            }

            if (class_check == true && check_equal == true) {
                classY[j] = class_count - 1;
            }
            M(i, j) = d2;
            M(j, i) = d3;
        }
    }
    if (classY[m-1] == 0){
        classY[m-1] = class_count;
    }
    List ret;
    ret["M"] = M;
    ret["ind"] = classY;
    return ret;

}



// [[Rcpp::export]]
List ecdf_list_comp_class_eps(List X, List t) {
    int m = X.size();

    NumericMatrix M(m, m);
    NumericMatrix sums_mat(m*m, 3);
    NumericMatrix vec_check(m*m, 2);

    M(m-1, m-1) = 1;
    int sd = 0;
    int count = 0;
    for (int i= 0; i < (m-1); i++) {
        M(i, i) = 1;

        NumericVector t1 = as<NumericVector>(t[i]);
        NumericVector x1 = as<NumericVector>(X[i]);

        int t_l1 = t1.size();
        for (int j = (i+1); j < m; j++) {
            NumericVector t2 = as<NumericVector>(t[j]);
            NumericVector x2 = as<NumericVector>(X[j]);

            int t_l2 = t2.size();
            int dim_thresh = t_l1 + t_l2;
            NumericVector thresholds(dim_thresh);
            for (int l = 0; l < t_l1; l++){
                thresholds[l] = t1[l];
            }
            for (int l = t_l1; l < dim_thresh; l++){
                thresholds[l] = t2[(l-t_l1)];
            }

            std::sort( thresholds.begin(), thresholds.end() );
            thresholds.erase( std::unique( thresholds.begin(), thresholds.end() ), thresholds.end() );
            int d = thresholds.size();

            NumericVector F1(d);
            NumericVector F2(d);

            int d2 = 1;
            int d3 = 1;



            for(int k = 0; k < d; k++){
                int indx1 = 0;
                int indx2 = 0;

                while (thresholds[k] >= t2[indx2] && indx2 < t_l2){
                    indx2 +=1;
                }
                if (indx2 == 0){
                    F2[k] = 0;
                } else {

                    F2[k] = x2[indx2-1];
                }

                while (thresholds[k] >= t1[indx1] && indx1 < t_l1){
                    indx1 +=1;
                }
                if (indx1 == 0){
                    F1[k] = 0;
                } else {
                    F1[k] = x1[indx1-1];
                }

                if (F1[k] < F2[k]){
                    d2 = 0;
                }
                if (F1[k] > F2[k]){
                    d3 = 0;
                }
            }

            sd = sd + d2 + d3;
            M(i, j) = d2;
            M(j, i) = d3;

            if (d2 + d3 == 0){
                double sum_s_left = 0;
                double sum_s_right = 0;
                double sum_all = 0;
                for (int k = 0; k< (d-1); k++) {
                    sum_all += std::abs(F1[k] - F2[k])*(thresholds[k+1] - thresholds[k]);
                    if (F2[k] > F1[k]) {
                        sum_s_left += (F2[k] - F1[k])*(thresholds[k+1] - thresholds[k]);
                    }
                    if (F1[k] > F2[k]) {
                        sum_s_right += (F1[k] - F2[k])*(thresholds[k+1] - thresholds[k]);
                    }
                }

                vec_check(count, 0) = i;
                vec_check(count, 1) = j;
                sums_mat(count, 0) = sum_all;

                sums_mat(count, 1) = sum_s_left;
                sums_mat(count, 2) = sum_s_right;
                count = count+1;
            }
        }
    }
    List list_d_mat;
    list_d_mat["d"] = sd;
    list_d_mat["k"] = count;
    list_d_mat["mat"] = sums_mat;
    list_d_mat["M"] = M;
    list_d_mat["indx"] = vec_check;
    return list_d_mat;
}

// [[Rcpp::export]]
NumericMatrix ecdf_comp_sd(NumericMatrix X, NumericVector t) {
    int m = X.nrow();
    int d = X.ncol();
    NumericMatrix M(m, m);

    M(m-1, m-1) = 1;
    for (int i= 0; i < (m-1); i++) {
        M(i, i) = 1;

        for (int j = (i+1); j < m; j++) {
            NumericVector F1(d);
            NumericVector F2(d);
            for(int k = 0; k < d; k++){
                F1[k] = X(i, k);
                F2[k] = X(j, k);
            }

            int test = 1;
            for (int k = 0; k< (d-1); k++) {

                if (F2[k] > F1[k]) {
                    test  = 0;
                }
            }

            if (test == 1) {
                M(i, j) = 1;
            }

            if (M(i, j) == 0) {
                int test2 = 1;
                for (int k = 0; k< (d-1); k++) {

                    if (F1[k] > F2[k]) {
                        test2 = 0;
                    }
                }
                if (test2 == 1) {
                    M(j, i) = 1;
                }
            }

        }
    }

    return M;

}

// [[Rcpp::export]]
List ecdf_comp_eps(NumericMatrix X, NumericVector t) {
    int m = X.nrow();
    int d = X.ncol();
    NumericMatrix M(m, m);
    NumericMatrix sums_mat(m*m, 3);
    NumericMatrix vec_check(m*m, 2);
    int sd = 0;
    int count = 0;
    M(m-1, m-1) = 1;
    for (int i= 0; i < (m-1); i++) {
        M(i, i) = 1;

        for (int j = (i+1); j < m; j++) {
            NumericVector F1(d);
            NumericVector F2(d);
            for(int k = 0; k < d; k++){
                F1[k] = X(i, k);
                F2[k] = X(j, k);
            }

            int test = 1;
            double sum_all = 0;
            double sum_s_left = 0;
            for (int k = 0; k< (d-1); k++) {
                sum_all += std::abs(F1[k] - F2[k])*(t[k+1] - t[k]);
                if (F2[k] > F1[k]) {
                    test  = 0;
                    sum_s_left += (F2[k] - F1[k])*(t[k+1] - t[k]);
                }
            }

            if (test == 1) {
                M(i, j) = 1;
                sd = sd + 1;
            }

            if (M(i, j) == 0) {
                int test2 = 1;
                double sum_s_right = 0;
                for (int k = 0; k< (d-1); k++) {


                    if (F1[k] > F2[k]) {
                        sum_s_right += (F1[k] - F2[k])*(t[k+1] - t[k]);
                        test2 = 0;
                    }
                }
                if (test2 == 1) {
                    M(j, i) = 1;
                    sd = sd + 1;
                } else {
                    vec_check(count, 0) = i;
                    vec_check(count, 1) = j;
                    sums_mat(count, 0) = sum_all;
                    sums_mat(count, 1) = sum_s_left;
                    sums_mat(count, 2) = sum_s_right;
                    count = count + 1;
                }
            }

        }
    }

    List list_d_mat;
    list_d_mat["d"] = sd;
    list_d_mat["k"] = count;
    list_d_mat["mat"] = sums_mat;
    list_d_mat["M"] = M;
    list_d_mat["indx"] = vec_check;
    return list_d_mat;

}


// [[Rcpp::export]]
List ecdf_comp_class_sd(NumericMatrix X, NumericVector t) {
    int m = X.nrow();
    int d = X.ncol();
    NumericMatrix M(m, m);
    NumericVector classY(m);
    int class_count = 1;

    M(m-1, m-1) = 1;
    for (int i= 0; i < (m-1); i++) {
        M(i, i) = 1;
        bool class_check = false;
        if (classY[i] == 0){
            classY[i] = class_count;
            class_count += 1;
            class_check = true;
        }
        for (int j = (i+1); j < m; j++) {
            NumericVector F1(d);
            NumericVector F2(d);
            int d2 = 1;
            int d3 = 1;
            bool check_equal = true;
            if (class_check == true){
                for(int k = 0; k < d; k++){
                    F1[k] = X(i, k);
                    F2[k] = X(j, k);
                    if (F1[k] != F2[k]){
                        check_equal = false;
                    }
                    if (F1[k] < F2[k]){
                        d2 = 0;
                    }
                    if (F1[k] > F2[k]){
                        d3 = 0;
                    }
                }
                if (check_equal == true){
                    classY[j] = class_count -1;
                }
            } else {
                for(int k = 0; k < d; k++){
                    F1[k] = X(i, k);
                    F2[k] = X(j, k);
                    if (F1[k] < F2[k]){
                        d2 = 0;
                    }
                    if (F1[k] > F2[k]){
                        d3 = 0;
                    }
                }
            }
            M(i, j) = d2;
            M(j, i) = d3;
        }
    }
    if (classY[m-1] == 0){
        classY[m-1] = class_count;
    }

    List ret;
    ret["M"] = M;
    ret["ind"] = classY;

    return ret;

}

// [[Rcpp::export]]
NumericMatrix ecdf_func(List X, List t, NumericVector thresholds) {
    int n = X.size();
    int d = thresholds.size();
    NumericMatrix F1_update(n, d);
    for (int i = 0; i < n; i++){

        NumericVector F1 = as<NumericVector>(X[i]);
        NumericVector points = as<NumericVector>(t[i]);

        for (int j = 0; j < d; j++){
            int indx = 0;
            int r = points.size();
            while(thresholds[j] >= points[indx] && indx < r){
                indx += 1;
            }
            if (indx == 0){
                F1_update(i, j) = 0;
            } else {
                F1_update(i, j) =  F1[(indx-1)];
            }

        }
    }
    return F1_update;
}

// [[Rcpp::export]]
NumericMatrix ecdf_func_mat(NumericMatrix X, NumericVector t, NumericVector thresholds) {
    int n = X.nrow();
    int d = thresholds.size();
    int r = t.size();
    NumericMatrix F1_update(n, d);
    for (int j = 0; j < d; j++){
        int indx = 0;
        while(thresholds[j] > t[indx] && indx < r){
            indx += 1;
        }
        if (indx == 0){
            for (int i = 0; i < n; i++){
                F1_update(i, j) = 0;
            }
        } else {
            for (int i = 0; i < n; i++){
                F1_update(i, j) = X(i, (indx-1));
            }
        }
    }
    return F1_update;
}

// [[Rcpp::export]]
NumericVector ecdf_comp_class_ind(NumericMatrix X) {
    int m = X.nrow();
    int d = X.ncol();
    NumericVector classY(m);
    int class_count = 1;


    for (int i= 0; i < (m-1); i++) {
        bool class_check = false;
        if (classY[i] == 0){
            classY[i] = class_count;
            class_count += 1;
            class_check = true;
        }
        for (int j = (i+1); j < m; j++) {
            NumericVector F1(d);
            NumericVector F2(d);
            bool check_equal = true;
            if (class_check == true){
                for(int k = 0; k < d; k++){
                    F1[k] = X(i, k);
                    F2[k] = X(j, k);
                    if (F1[k] != F2[k]){
                        check_equal = false;
                    }
                }
                if (check_equal == true){
                    classY[j] = class_count -1;
                }
            }
        }
    }
    if (classY[m-1] == 0){
        classY[m-1] = class_count;
    }

    return classY;

}

// [[Rcpp::export]]
NumericVector ecdf_list_comp_class_ind(List X, List t) {
    int m = X.size();
    NumericVector classY(m);
    int class_count = 1;

    for (int i= 0; i < (m-1); i++) {
        bool class_check = false;
        if (classY[i] == 0){
            classY[i] = class_count;
            class_count += 1;
            class_check = true;
        }

        NumericVector t1 = as<NumericVector>(t[i]);
        NumericVector x1 = as<NumericVector>(X[i]);
        int t_l1 = t1.size();
        for (int j = (i+1); j < m; j++) {

            NumericVector t2 = as<NumericVector>(t[j]);
            NumericVector x2 = as<NumericVector>(X[j]);

            int t_l2 = t2.size();
            int dim_thresh = t_l1 + t_l2;
            NumericVector thresholds(dim_thresh);
            for (int l = 0; l < t_l1; l++){
                thresholds[l] = t1[l];
            }
            for (int l = t_l1; l < dim_thresh; l++){
                thresholds[l] = t2[(l-t_l1)];
            }

            std::sort( thresholds.begin(), thresholds.end() );
            thresholds.erase( std::unique( thresholds.begin(), thresholds.end() ), thresholds.end() );
            int d = thresholds.size();

            NumericVector F1(d);
            NumericVector F2(d);

            bool check_equal = true;
            if (class_check == true){
                for(int k = 0; k < d; k++){
                    int indx1 = 0;
                    int indx2 = 0;
                    while (thresholds[k] >= t1[indx1] && indx1 < t_l1){
                        indx1 +=1;
                    }
                    if (indx1 == 0){
                        F1[k] = 0;
                    } else {
                        F1[k] = x1[indx1-1];
                    }
                    while (thresholds[k] >= t2[indx2] && indx2 < t_l2){
                        indx2 +=1;
                    }
                    if (indx2 == 0){
                        F2[k] = 0;
                    } else {
                        F2[k] = x2[indx2-1];
                    }
                    if (F1[k] != F2[k]){
                        check_equal = false;
                    }
                }
                if (check_equal == true){
                    classY[j] = class_count -1;
                }
            }
        }
    }
    if (classY[m-1] == 0){
        classY[m-1] = class_count;
    }
    return classY;

}

// [[Rcpp::export]]
List ecdf_comp_class_eps(NumericMatrix X, NumericVector t) {
    int m = X.nrow();
    int d = X.ncol();
    NumericMatrix M(m, m);
    NumericMatrix sums_mat(m*m, 3);
    NumericMatrix vec_check(m*m, 2);
    int sd = 0;
    M(m-1, m-1) = 1;
    int count = 0;

    for (int i= 0; i < (m-1); i++) {
        M(i, i) = 1;
        for (int j = (i+1); j < m; j++) {
            NumericVector F1(d);
            NumericVector F2(d);
            int d2 = 1;
            int d3 = 1;
            for(int k = 0; k < d; k++){
                F1[k] = X(i, k);
                F2[k] = X(j, k);
                if (F1[k] < F2[k]){
                    d2 = 0;
                }
                if (F1[k] > F2[k]){
                    d3 = 0;
                }
            }

            sd = sd + d2 + d3;
            M(i, j) = d2;
            M(j, i) = d3;

            if (d2 + d3 == 0){

                double sum_s_left = 0;
                double sum_s_right = 0;
                double sum_all = 0;
                for (int k = 0; k< (d-1); k++) {
                    sum_all += std::abs(F1[k] - F2[k])*(t[k+1] - t[k]);
                    if (F2[k] > F1[k]) {
                        sum_s_left += (F2[k] - F1[k])*(t[k+1] - t[k]);
                    }
                    if (F1[k] > F2[k]) {
                        sum_s_right += (F1[k] - F2[k])*(t[k+1] - t[k]);
                    }
                }
                vec_check(count, 0) = i;
                vec_check(count, 1) = j;
                sums_mat(count, 0) = sum_all;
                sums_mat(count, 1) = sum_s_left;
                sums_mat(count, 2) = sum_s_right;
                count = count+1;
            }
        }
    }
    List list_d_mat;
    list_d_mat["d"] = sd;
    list_d_mat["k"] = count;
    list_d_mat["mat"] = sums_mat;
    list_d_mat["M"] = M;
    list_d_mat["indx"] = vec_check;
    return list_d_mat;
}



