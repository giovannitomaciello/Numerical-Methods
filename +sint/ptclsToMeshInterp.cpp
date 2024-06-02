#include "mex.h"
#include <cmath>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>
#include <tbb/tbb.h>

using namespace Eigen;
using namespace tbb;

struct Particle {
    double x, y, z, q;
};

// Function to find the grid index
inline int gridIndex(double coord, double minCoord, double invGridSpacing) {
    return static_cast<int>((coord - minCoord) * invGridSpacing);
}

// Custom hash function for grid indices
struct GridIndex {
    int x, y, z;

    bool operator==(const GridIndex &other) const {
        return x == other.x && y == other.y && z == other.z;
    }

    struct HashFunction {
        std::size_t operator()(const GridIndex &k) const {
            return std::hash<int>()(k.x) ^ (std::hash<int>()(k.y) << 1) ^ (std::hash<int>()(k.z) << 2);
        }
    };
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for proper number of arguments.
    if (nrhs != 6) {
        mexErrMsgIdAndTxt("MATLAB:rho_lr_mex:invalidNumInputs", "Six input arguments required.");
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
    double rCut = mxGetScalar(prhs[5]); // Extract rCut from input arguments

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
    double rCutSq = rCut * rCut;
    double H = 3.0 / (M_PI * rCutSq);

    // Initialize rho_lr to zero
    parallel_for(blocked_range<size_t>(0, numPoints), [&](const blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            rho_lr[i] = 0.0;
        }
    });

    // Create Eigen maps for efficient access
    Map<ArrayXd> X_map(X, numPoints);
    Map<ArrayXd> Y_map(Y, numPoints);
    Map<ArrayXd> Z_map(Z, numPoints);
    Map<ArrayXd> q_map(q, numPtcls);
    Map<MatrixXd> ptcls_x_map(ptcls_x, 3, numPtcls);

    // Determine the grid spacing
    double minX = X_map.minCoeff();
    double minY = Y_map.minCoeff();
    double minZ = Z_map.minCoeff();
    double maxX = X_map.maxCoeff();
    double maxY = Y_map.maxCoeff();
    double maxZ = Z_map.maxCoeff();
    double gridSpacing = rCut;
    double invGridSpacing = 1.0 / gridSpacing;

    // Create a grid to hold particles
    std::unordered_map<GridIndex, std::vector<Particle>, GridIndex::HashFunction> grid;

    for (size_t i = 0; i < numPtcls; ++i) {
        Particle p = {ptcls_x_map(0, i), ptcls_x_map(1, i), ptcls_x_map(2, i), q_map(i)};
        GridIndex index = {gridIndex(p.x, minX, invGridSpacing), gridIndex(p.y, minY, invGridSpacing), gridIndex(p.z, minZ, invGridSpacing)};
        grid[index].push_back(p);
    }

    // Precompute all the possible grid indices for neighbor searching
    std::vector<GridIndex> neighborOffsets;
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                neighborOffsets.push_back({dx, dy, dz});
            }
        }
    }


    // Compute rho_lr using the precomputed grid in parallel
    tbb::parallel_for(tbb::blocked_range<size_t>(0, numPoints), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            double xi = X_map(i);
            double yi = Y_map(i);
            double zi = Z_map(i);

            int gridX = gridIndex(xi, minX, invGridSpacing);
            int gridY = gridIndex(yi, minY, invGridSpacing);
            int gridZ = gridIndex(zi, minZ, invGridSpacing);

            double localRho = 0.0; // Local variable for accumulation to avoid atomic operation

            for (const auto &offset : neighborOffsets) {
                int nx = gridX + offset.x;
                int ny = gridY + offset.y;
                int nz = gridZ + offset.z;

                GridIndex neighborIndex = {nx, ny, nz};
                auto it = grid.find(neighborIndex);
                if (it != grid.end()) {
                    for (const auto &p : it->second) {
                        double dx = xi - p.x;
                        double dy = yi - p.y;
                        double dz = zi - p.z;
                        double rSq = dx * dx + dy * dy + dz * dz;

                        if (rSq <= rCutSq) {
                            localRho += p.q * H * (1 - rSq / rCutSq);
                        }
                    }
                }
            }

            // Accumulate the localRho to the global rho_lr array
            rho_lr[i] = localRho;
        }
    });
}
