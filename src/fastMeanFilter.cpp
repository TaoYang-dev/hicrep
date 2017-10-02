#include <Rcpp.h>
using namespace Rcpp;

// fast mean filter Rcpp by Fan Song
// Algorithm derived from Nakariyakul S. J Supercomput. 2013
//  

// [[Rcpp::export]]
NumericMatrix extendMatrix(NumericMatrix mat, int h) {
    int nrow = mat.nrow();
    NumericMatrix extendedMat(nrow + 2*h, nrow + 2*h);
    for (int i = 0; i < nrow; i++)
        for (int j = 0; j < nrow; j++)
        extendedMat(h + i, h + j) = mat(i, j);
    return extendedMat;
}

// [[Rcpp::export]]
int estimateLengthOrWidth(int x, int nrow, int h) { 
    int lw = 0;
    if (x <= h) {
        lw = x + h + 1;
    } else {
        if (x + h + 1 < nrow)
        lw = 2*h + 1;
    else
        lw = nrow + h - x;
    }
    return lw;
}

// [[Rcpp::export]]
NumericMatrix fastMeanFilter(NumericMatrix mat, int h) {
    int nrow = mat.nrow();
    NumericMatrix extendedMat = extendMatrix(mat, h);
    NumericMatrix smoothedMat(nrow, nrow);
    int nrowExt = nrow + 2*h;
    double sum = 0;
    double N = 0;
    int H = 2*h + 1;
    NumericVector colSum(nrowExt);
    
    for (int i = 0; i < nrowExt; i++)
        for (int j = 0; j < H; j++)
            colSum[i] += extendedMat(j, i);
  
    for (int i = 0; i < H; i++)
        sum += colSum[i];
  
    smoothedMat(0, 0) = sum / ((h + 1) * (h + 1));
  
    int l = 0; // length of window
    for (int col = 1; col < nrow; col++) {
        sum = sum - colSum[col-1] + colSum[col+H-1];
        l = estimateLengthOrWidth(col, nrow, h);
        N = l * (h + 1);
        smoothedMat(0, col) = sum / N;
    }
  
    int w = 0; // width of window
    for (int row = 1; row < nrow; row++) {
        sum = 0;
        for (int i = 0; i < nrowExt; i++)
            colSum[i] = colSum[i] - extendedMat(row - 1, i) + 
            extendedMat(row + H - 1, i);
    
        for (int i = 0; i < H; i++)
            sum += colSum[i];
    
        w = estimateLengthOrWidth(row, nrow, h);
        N = w * (h + 1);
        smoothedMat(row, 0) = sum / N;
    
        for (int col = 1; col < nrow; col++) {
            sum = sum - colSum[col-1] + colSum[col+H-1];
            l = estimateLengthOrWidth(col, nrow, h);
            N = l * w;
            smoothedMat(row, col) = roundf((sum / N) * 100)/100;
        }
    }
    return smoothedMat;
}