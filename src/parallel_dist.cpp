// colum subtraction followed by squaring elements in column
#include <Rcpp.h>
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;
using namespace Rcpp;
using namespace std;

// helper function to find difference & square it for a pair of 
// corresponding vector entries
inline double calc_dist(double v1, double v2) {
    double diff = v1 - v2;
    return diff * diff;
}

// CIDR helper function
inline double calc_cidr_dist(double v1, double v2) {
    double diff = v1 - v2;
    return diff * diff;
}

// worker function for parallel processing - euclidean
struct SquareDist: public Worker
{
    // source matrix - RMatrix provided for thread-safe parallel calculations
    RMatrix<double> input;

    // destination matrix
    RMatrix<double> output;

    // input range - the range for the columns - this worker requires
    // an iterator as input - but we want to know the column index 
    // on which we are working
    IntegerVector range;
    
    // input vector - the column with which the differece will be applied
    NumericMatrix::Column col;
    

    // initialize with source and destination
    SquareDist(NumericMatrix input, NumericMatrix output, NumericMatrix::Column col, IntegerVector range) 
        : input(input), output(output), range(range), col(col) {}
    

    // perform the "distance calculation" 
    void operator()(std::size_t begin, std::size_t end) {
        // Iterator passed in is for vector of column indexes
        int start_col_indx = (int)*(range.begin() + begin);
        // note that end will be one past the last element
        int end_col_indx = (int)*(range.begin() + end - 1);
        for (int i = start_col_indx; i <= end_col_indx; ++i) {
            RMatrix<double>::Column curr_col(input.column(i));
            // apply the calculation & save to the output matrix
            std::transform(curr_col.begin(),    // range 1 begin
                           curr_col.end(),      // range 1 end
                           col.begin(),         // range 2 begin
                           output.column(i).begin(),    // output begin
                           calc_dist);          // function to apply
        }
    
    
    }
};

// worker function for parallel processing - cidr
struct CidrDist: public Worker
{
    // source matrix - RMatrix provided for thread-safe parallel calculations
    RMatrix<double> input;

    // destination matrix
    RMatrix<double> output;

    // truth matrix
    RMatrix<int> truth;

    const double threshold;

    // input range - the range for the columns - this worker requires
    // an iterator as input - but we want to know the column index 
    // on which we are working
    IntegerVector range;
    
    // input vector - the column with which the differece will be applied
    NumericMatrix::Column col;
    
    // correpsonding column from truth matrix
    IntegerMatrix::Column truth_col;
    
    // initialize with source and destination
    CidrDist(NumericMatrix input, IntegerMatrix truth, NumericMatrix output, NumericMatrix::Column col, IntegerMatrix::Column truth_col, IntegerVector range, double threshold) 
        : input(input), truth(truth), output(output), range(range), col(col), truth_col(truth_col), threshold(threshold) {}
    
    // perform the "distance calculation" 
    void operator()(std::size_t begin, std::size_t end) {
        const double EPSILON = 0.00000001;
        // Iterator passed in is for vector of column indexes
        int start_col_indx = (int)*(range.begin() + begin);
        // note that end will be one past the last element
        int end_col_indx = (int)*(range.begin() + end - 1);
        // number of rows
        int nrows = input.nrow();
        for (int i = start_col_indx; i <= end_col_indx; ++i) {
            RMatrix<double>::Column curr_col(input.column(i));
            RMatrix<int>::Column curr_truth_col(truth.column(i));
            // apply the calculation & save to the output matrix
            for (int j=0; j<nrows; ++j) {
                double dist = 0;
                if ((!curr_truth_col[j] & !truth_col[j]) | 
                    (std::max(curr_col[j], col[j]) > threshold)) {
                    dist = curr_col[j] - col[j];
                    dist *= dist;
                }
                output(j, i) = dist;
            }
        }
    }
};

// input NumericMatrix x - the input matrix
// input int i - the R index of the "reference column"
// - foreach column in the input matrix, we subtract the 
// "reference column" & square the individual entries that result
// NOTE that R columns start at 1 while in cpp, they start from 0
// [[Rcpp::export]]
NumericMatrix parallelMatrixSquareDist(NumericMatrix x, int i) {
    // allocate the output matrix
    NumericMatrix output(x.nrow(), x.ncol());
    // create a vector from 0 to x.length()-1
    IntegerVector range = seq_len(x.cols());
    range = range - 1;
    // distance functor (pass input and output matrixes)
    SquareDist square_dist(x, output, x(_, i-1), range);
    // call parallelFor to do the work
    int items = range.length();
    parallelFor(0, items, square_dist);

    // return the output matrix
    return output;
}
//
// input NumericMatrix x - the input matrix
// input int i - the R index of the "reference column"
// - foreach column in the input matrix, we subtract the 
// "reference column" & square the individual entries that result
// NOTE that R columns start at 1 while in cpp, they start from 0
// [[Rcpp::export]]
NumericMatrix parallelMatrixCidrDist(NumericMatrix x, IntegerMatrix truth, int i, double threshold) {
    // allocate the output matrix
    NumericMatrix output(x.nrow(), x.ncol());
    // create a vector from 0 to x.length()-1
    IntegerVector range = seq_len(x.cols());
    range = range - 1;
    // distance functor (pass input and output matrixes)
    CidrDist cidr_dist(x, truth, output, x(_, i-1), truth(_, i-1), range, threshold);
    // call parallelFor to do the work
    int items = range.length();
    parallelFor(0, items, cidr_dist);

    // return the output matrix
    return output;
}
