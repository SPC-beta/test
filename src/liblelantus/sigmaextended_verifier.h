#ifndef BZX_LIBLELANTUS_SIGMAEXTENDED_VERIFIER_H
#define BZX_LIBLELANTUS_SIGMAEXTENDED_VERIFIER_H

#include "lelantus_primitives.h"

namespace lelantus {

class SigmaExtendedVerifier{

public:
    SigmaExtendedVerifier(const GroupElement& g,
                      const std::vector<GroupElement>& h_gens,
                      uint64_t n, uint64_t m_);

    //gets initial double-blinded Pedersen commitments,
    //verifies proofs from single transaction, where set size and challenge are the same
    bool batchverify(const std::vector<GroupElement>& commits,
                     const Scalar& x,
                     const std::vector<Scalar>& serials,
                     const std::vector<SigmaExtendedProof>& proofs) const;
    //gets initial double-blinded Pedersen commitments
    //verifies proofs from different transactions, where set sizes and challenges are different
    bool batchverify(const std::vector<GroupElement>& commits,
                     const std::vector<Scalar>& challenges,
                     const std::vector<Scalar>& serials,
                     const std::vector<size_t>& setSizes,
                     const std::vector<SigmaExtendedProof>& proofs) const;

private:
    //auxiliary functions
    bool membership_checks(const SigmaExtendedProof& proof) const;
    bool compute_fs(
            const SigmaExtendedProof& proof,
            const Scalar& x,
            std::vector<Scalar>& f_) const;
    bool abcd_checks(
            const SigmaExtendedProof& proof,
            const Scalar& x,
            const std::vector<Scalar>& f_) const;

    void compute_fis(int j, const std::vector<Scalar>& f, std::vector<Scalar>& f_i_) const;
    void compute_fis(
            const Scalar& f_i,
            int j,
            const std::vector<Scalar>& f,
            std::vector<Scalar>::iterator& ptr,
            std::vector<Scalar>::iterator end_ptr) const;
    void compute_batch_fis(
            const Scalar& f_i,
            int j,
            const std::vector<Scalar>& f,
            const Scalar& y,
            Scalar& e,
            std::vector<Scalar>::iterator& ptr,
            std::vector<Scalar>::iterator start_ptr,
            std::vector<Scalar>::iterator end_ptr) const;

private:
    GroupElement g_;
    std::vector<GroupElement> h_;
    uint64_t n;
    uint64_t m;
};

} // namespace lelantus

#endif //BZX_LIBLELANTUS_SIGMAEXTENDED_VERIFIER_H
