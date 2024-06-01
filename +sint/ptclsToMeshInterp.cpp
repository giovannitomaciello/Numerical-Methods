#include "mex.h"
#include <cmath>
#include <omp.h>
#include <Eigen/Dense>

using namespace Eigen;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for proper number of arguments.
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("MATLAB:rho_lr_mex:invalidNumInputs", "Five input arguments required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("MATLAB:rho_lr_mex:maxlhs", "Too many output arguments.");
    }

    // Input arguments
    double *X = mxGetPr(prhs[0]);
    double *Y = mxGetPr(prhs[1]);
    double *Z = mxGetPr(prhs[2]);
    double *q = mxGetPr(prhs[3]);
    double *ptcls_x = mxGetPr(prhs[4]);

    // Get dimensions
    mwSize numPoints = mxGetNumberOfElements(prhs[0]);
    mwSize numPtcls = mxGetN(prhs[4]);  // Assuming ptcls_x is a 3 x N matrix

    if (numPoints != mxGetNumberOfElements(prhs[1]) || numPoints != mxGetNumberOfElements(prhs[2])) {
        mexErrMsgIdAndTxt("MATLAB:rho_lr_mex:dimMismatch", "Dimensions of X, Y, and Z must match.");
    }

    // Output argument
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), mxGetN(prhs[0]), mxREAL);
    double *rho_lr = mxGetPr(plhs[0]);

    // Constants
    double rCut = 1.0; // Assuming rCut is predefined or passed as an argument
    double rCutSq = rCut * rCut;
    double H = 3.0 / (M_PI * rCutSq);

    // Initialize rho_lr to zero
    std::fill(rho_lr, rho_lr + numPoints, 0.0);

    // Create Eigen maps for efficient access
    Map<ArrayXd> X_map(X, numPoints);
    Map<ArrayXd> Y_map(Y, numPoints);
    Map<ArrayXd> Z_map(Z, numPoints);
    Map<ArrayXd> q_map(q, numPtcls);
    Map<MatrixXd> ptcls_x_map(ptcls_x, 3, numPtcls);

    // Compute rho_lr using Eigen
    #pragma omp parallel for
    for (mwSize i = 0; i < numPoints; ++i) {
        double xi = X_map(i);
        double yi = Y_map(i);
        double zi = Z_map(i);

        for (mwSize k = 0; k < numPtcls; ++k) {
            double dx = xi - ptcls_x_map(0, k);
            double dy = yi - ptcls_x_map(1, k);
            double dz = zi - ptcls_x_map(2, k);
            double rSq = dx * dx + dy * dy + dz * dz;

            if (rSq <= rCutSq) {
                double r = std::sqrt(rSq);
                #pragma omp atomic
                rho_lr[i] += q_map(k) * H * (1 - r / rCut);
            }
        }
    }
}
