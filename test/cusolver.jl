
using LinearAlgebra
using CUDAKernels
using BlockPowerFlow

using CUDA.CUSPARSE
import ExaPF: LinearSolvers
import BlockPowerFlow: CUSOLVERRF

const LS = LinearSolvers

ExaPF.default_sparse_matrix(::CUDADevice) = CUSPARSE.CuSparseMatrixCSR

# Overload factorization routine to use cusolverRF
LS.exa_factorize(J::CuSparseMatrixCSR) = CUSOLVERRF.CusolverRfLU(J)
LS.exa_factorize(J::CuSparseMatrixCSC) = CUSOLVERRF.CusolverRfLU(J)

# Overload factorization for batch Hessian computation
function ExaOpt._batch_hessian_factorization(J::CuSparseMatrixCSR, nbatch)
    Jtrans = CUSPARSE.CuSparseMatrixCSC(J)
    if nbatch == 1
        lufac = CUSOLVERRF.CusolverRfLU(J)
        lufact = CUSOLVERRF.CusolverRfLU(Jtrans)
    else
        lufac = CUSOLVERRF.CusolverRfLUBatch(J, nbatch)
        lufact = CUSOLVERRF.CusolverRfLUBatch(Jtrans, nbatch)
    end
    return (lufac, lufact)
end

function ExaOpt.update_factorization!(hlag::ExaOpt.AbstractHessianStorage, J::CUSPARSE.CuSparseMatrixCSR)
    LinearAlgebra.lu!(hlag.lu, J)
    ∇gₓᵀ = CUSPARSE.CuSparseMatrixCSC(J)
    LinearAlgebra.lu!(hlag.adjlu, ∇gₓᵀ)
    return
end

