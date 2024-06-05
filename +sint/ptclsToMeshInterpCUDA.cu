#include "mex.h"
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <gpu/mxGPUArray.h>

#define M_PI 3.14

// Define the Particle structure
struct Particle {
    double x, y, z, q;
};

// Function to find the grid index
__device__ int gridIndex(double coord, double minCoord, double invGridSpacing) {
    return static_cast<int>((coord - minCoord) * invGridSpacing);
}

// Define the GridIndex structure
struct GridIndex {
    int x, y, z;

    __host__ __device__ bool operator==(const GridIndex &other) const {
        return x == other.x && y == other.y && z == other.z;
    }

    struct HashFunction {
        __host__ __device__ std::size_t operator()(const GridIndex &k) const {
            return std::hash<int>()(k.x) ^ (std::hash<int>()(k.y) << 1) ^ (std::hash<int>()(k.z) << 2);
        }
    };
};

// CUDA kernel function to compute rho_lr
__global__ void computeRhoLrKernel(const double* X, const double* Y, const double* Z, const double* q,
                                   const Particle* particles, const GridIndex* neighborOffsets,
                                   double* rho_lr, size_t numPoints, size_t numPtcls,
                                   double minX, double minY, double minZ, double invGridSpacing,
                                   double rCutSq, double H, size_t numOffsets) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numPoints) return;

    double xi = X[i];
    double yi = Y[i];
    double zi = Z[i];

    int gridX = gridIndex(xi, minX, invGridSpacing);
    int gridY = gridIndex(yi, minY, invGridSpacing);
    int gridZ = gridIndex(zi, minZ, invGridSpacing);

    double localRho = 0.0;

    for (size_t j = 0; j < numOffsets; ++j) {
        int nx = gridX + neighborOffsets[j].x;
        int ny = gridY + neighborOffsets[j].y;
        int nz = gridZ + neighborOffsets[j].z;

        for (size_t k = 0; k < numPtcls; ++k) {
            const Particle& p = particles[k];
            int px = gridIndex(p.x, minX, invGridSpacing);
            int py = gridIndex(p.y, minY, invGridSpacing);
            int pz = gridIndex(p.z, minZ, invGridSpacing);

            if (nx == px && ny == py && nz == pz) {
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

    rho_lr[i] = localRho;
}

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
    // Variables for GPU array
    const mxGPUArray *ptcls_x_gpu = nullptr;
    const double *ptcls_x = nullptr;

    // Check if input is a GPU array
    if (mxIsGPUArray(prhs[4])) {
        ptcls_x_gpu = mxGPUCreateFromMxArray(prhs[4]);
        ptcls_x = static_cast<const double*>(mxGPUGetDataReadOnly(ptcls_x_gpu));
    } else {
        // If input is not a GPU array, handle accordingly
        ptcls_x = mxGetPr(prhs[4]);
    }
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

    // Determine the grid spacing
    double minX = X[0];
    double minY = Y[0];
    double minZ = Z[0];
    double maxX = X[0];
    double maxY = Y[0];
    double maxZ = Z[0];
    for (mwSize i = 1; i < numPoints; ++i) {
        if (X[i] < minX) minX = X[i];
        if (X[i] > maxX) maxX = X[i];
        if (Y[i] < minY) minY = Y[i];
        if (Y[i] > maxY) maxY = Y[i];
        if (Z[i] < minZ) minZ = Z[i];
        if (Z[i] > maxZ) maxZ = Z[i];
    }
    double gridSpacing = (Y[1] - Y[0]); // Assuming uniform spacing
    double invGridSpacing = 1.0 / gridSpacing;

    // Create particles vector
    thrust::host_vector<Particle> h_particles(numPtcls);
    for (size_t i = 0; i < numPtcls; ++i) {
        h_particles[i] = {ptcls_x[3*i], ptcls_x[3*i+1], ptcls_x[3*i+2], q[i]};
    }

    // Precompute all the possible grid indices for neighbor searching within cutoff radius
    std::vector<GridIndex> neighborOffsets;
    int maxOffset = std::ceil(rCut / gridSpacing);
    for (int dx = -maxOffset; dx <= maxOffset; ++dx) {
        for (int dy = -maxOffset; dy <= maxOffset; ++dy) {
            for (int dz = -maxOffset; dz <= maxOffset; ++dz) {
                double distance = std::sqrt(dx*dx + dy*dy + dz*dz) * gridSpacing;
                if (distance <= rCut) {
                    neighborOffsets.push_back({dx, dy, dz});
                }
            }
        }
    }

    // Allocate GPU memory
    double *d_X, *d_Y, *d_Z, *d_q, *d_rho_lr;
    Particle *d_particles;
    GridIndex *d_neighborOffsets;

    cudaMalloc(&d_X, numPoints * sizeof(double));
    cudaMalloc(&d_Y, numPoints * sizeof(double));
    cudaMalloc(&d_Z, numPoints * sizeof(double));
    cudaMalloc(&d_q, numPtcls * sizeof(double));
    cudaMalloc(&d_rho_lr, numPoints * sizeof(double));
    cudaMalloc(&d_particles, numPtcls * sizeof(Particle));
    cudaMalloc(&d_neighborOffsets, neighborOffsets.size() * sizeof(GridIndex));

    // Copy data to GPU
    cudaMemcpy(d_X, X, numPoints * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Y, Y, numPoints * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Z, Z, numPoints * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_q, q, numPtcls * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_neighborOffsets, neighborOffsets.data(), neighborOffsets.size() * sizeof(GridIndex), cudaMemcpyHostToDevice);

    // Copy particles data to GPU
    cudaMemcpy(d_particles, thrust::raw_pointer_cast(h_particles.data()), numPtcls * sizeof(Particle), cudaMemcpyHostToDevice);

    // Determine the number of threads and blocks
    int blockSize = 256;
    int numBlocks = (numPoints + blockSize - 1) / blockSize;

    // Launch kernel
    computeRhoLrKernel<<<numBlocks, blockSize>>>(d_X, d_Y, d_Z, d_q, d_particles, d_neighborOffsets, 
                                                 d_rho_lr, numPoints, numPtcls, minX, minY, minZ, 
                                                 invGridSpacing, rCutSq, H, neighborOffsets.size());

    // Copy the result back to the host
    cudaMemcpy(rho_lr, d_rho_lr, numPoints * sizeof(double), cudaMemcpyDeviceToHost);

    // Free GPU memory
    cudaFree(d_X);
    cudaFree(d_Y);
    cudaFree(d_Z);
    cudaFree(d_q);
    cudaFree(d_rho_lr);
    cudaFree(d_particles);
    cudaFree(d_neighborOffsets);

    // Destroy GPU array if created
    if (ptcls_x_gpu != nullptr) {
        mxGPUDestroyGPUArray(ptcls_x_gpu);
    }
}
