#ifndef SPIN_HALF_HEISENBERG_MODEL_HEADER
#define SPIN_HALF_HEISENBERG_MODEL_HEADER

#include "../sites/spinhalf.h"
#include "../../mps/tt.h"
#include "../../mps/qtt.h"
#include "AutoMPO.h"


/*
A. Heisenberg model (with optional random fields), and MPO layout:

  H = \sum_i J/2(S+_i S-_{i+1} + S-_i S+_{i+1}) + J Sz_i Sz_{+1} + h_i S_z_i

  MPO layout of an Hamiltonian is not unique. There are many equivalent way of
  spelling out the matrices for each sites. We will discuss the one layout
  that is most appropriate for working with quantum numbers.

  e.g.
  1D: On each site, except at the boundary, the matrix looks like
    I    0    0    0    0
   hSz   I   JSz  JS-  JS+
    Sz   0    0    0    0
    S+   0    0    0    0
    S-   0    0    0    0
  On the first site (leftmost site), the matrix is just the second row of it.
  On the last site (rightmost site), the matrix is just the first column of it.

  You may try to multiply matrices for two adjacent sites together, and convince
  yourself that this in the end gives you the correct Hamiltonian.

  There are a few things sepcial about this kind of layout in general:
  (1) The (0,0) element is always Identity
  (2) The (1,1) element is always Identity
  (3) The (1,0) element always contain single body on-site operators, like Zeeman fields or Hubbard U
  (4) On site i, elements in the 1st column following (2,0) are connected to operators
  on 2nd row of site i-1.
  (5) On site i, elements in the 2nd row following (1,2) are connected to operators
  on site on 1st column of i+1.



B.  Long range interaction
  Instead of nearest neighbor interaction, assume that we have a 1D Heisenberg chain
  with interaction length up to 3 sites away, so we need NN, NNN, NNNN interactions.

             NN   NNN NNNN  NN   NNN NNNN  NN   NNN NNNN
    -----------------------------------------------------
    I    0    0    0    0    0    0    0    0    0    0
   hSz   I   JSz  JSz  JSz  JS-  JS-  JS-  JS+  JS+  JS+
    Sz        0    0    0
    0         I    0    0
    0         0    I    0
    S+                       0    0    0
    0                        I    0    0
    0                        0    I    0
    S-                                      0    0    0
    0                                       I    0    0
    0                                       0    I    0

  (Unfilled spots above are zero)
  (The 3-by-3 diagonal blocks make NNN and NNNN interactions possible.)

  Convince yourself that the added blocks passes through NNN and NNNN interactions.
  You may try to work out the MPOs for Heisenberg ladders.



C. Quantum number of MPO:
  For matrix A on site i, the flow of quantum numbers by default looks like the following.

                    top state (bra)
                         |
                         v
    vitual bond in ----> A ----> virtual bond out
                         |
                         v
                  bottom state (ket)

  We know that a vertical bond, or a physical bond, can either be spin down or up,
  so it carrys quatum number -1 or +1.
  We will try to understand the virtual bonds' quantum numbers below.

  Before we label the MPO's quantum number sector, there are a few things that can
  help our understanding.
  (1) Total Sz for this model is conserved, and its quantum number is abelian (additive).
  (2) Conservation of S^2 (if present) is not considered.
  (2) A S+ is always paired with a S-, possibly connected by some Identities, if they
  correspond to long range interaction.
  (3) Non-zero quantum numbers flow between a paired S+ and S-, but not outside of them.
  (4) There is no quantum number flow associated with I or Sz.

  Now, for the 1D nearest neighbor (NN) Heisenberg chain with magnetic fields,
  we can label the virtual bonds' quantum numbers as

      out 0    0    0   -2   +2
  in  | -----------------------
  0   |   I    0    0    0    0
  0   |  hSz   I   JSz  JS-  JS+
  0   |   Sz   0    0    0    0
  -2  |   S+   0    0    0    0
  +2  |   S-   0    0    0    0

  It is easy to understand that the 3-by-3 block (usually called "Identity block")
  at the top left corner has no in/out virtual bond quantum numbers.
  For the S+ operator in the first column, let's understand why it has in-QN=-2
  and out-QN=0.
  (1) This S+ is connected to the S- to the left of this site.
  (2) To its right, there is nothing connected to it. So for this S+, we have

                        <+| (bra)
                         | (+1)
                         v
    vitual bond in ----> S+ ----> 0
                         |
                         v (-1)
                        |-> (ket)

      In order to have the sum of all in QNs equal to the sum of all out QNs,
      vitual bond in must carry a QN=-2.

  Similarly, we can show that the S- operator in the first column has QN=2.
  You may work out the QNs for the S+/- on the 2nd row as exercise.



D. qtensor format
  This is a good place to introduce the qtensor class (quantum numbered tensor),
  since we need to store the above blocked Hamiltonian using it.

  qtensor class is based on tensor diagram like the following one

                  top bond
                      |
                      v
  vitual bond 1 ----> A ----> virtual bond 2
                      |
                      v
                bottom bond

  Using it as an example, we see that it has 4 indices (or bonds), each bond may
  carry certain numbers of quantum number, like 0,-2,+2 for the NN Heisenberg
  chain with magnetic fields; and each quantum number has a certain size,
  like 3 or 1 for the the NN Heisenberg chain with magnetic fields.

  So the basic idea of qtensor is:
  (1) A class called "qtensor_index" will spell out the index, whose attributes are:
  arrow (Inward or Outward), name, type (Link or Site), level (prime level), and QNs.
  The QNs contains a certain number of pairs of (qn, qdim), where the first element
  is the quantum number and the second element is the associated size.
  (2) Given all the "qtensor_index"s, qtensor filters out all the legal combinations
  of quantum numbers (\sum in-QNs = \sum out-QNs). Each legal block is stored as
  a separate tensor of the same rank. Each block manages its data indiidually.

  For example, this tensor below has 4 indices,
  (1) v1: Inward,  left_virtual,  Link, level 0, QNs {(0,3),  (-2,1), (+2,1)}
  (2) v2: Outward, right_virtual, Link, level 0, QNs {(0,3),  (-2,1), (+2,1)}
  (3) p1: Inward,  top_physical,  Site, level 1, QNs {(-1,1), (+1,1)}
  (3) p2: Outward, bot_physical,  Site, level 0, QNs {(-1,1), (+1,1)}

  It has 5 legal blocks:
  (1) (v1.qn=0,  v2.qn=0,  p1.qn=-1, p2.qn=-1), size = 3*3*1*1 = 9
  (2) (v1.qn=-2, v2.qn=0,  p1.qn=+1, p2.qn=-1), size = 1*3*1*1 = 3
  (3) (v1.qn=+2, v2.qn=0,  p1.qn=-1, p2.qn=+1), size = 1*3*1*1 = 3
  (4) (v1.qn=0,  v2.qn=-2, p1.qn=-1, p2.qn=+1), size = 3*1*1*1 = 3
  (5) (v1.qn=0,  v2.qn=+2, p1.qn=+1, p2.qn=-1), size = 3*1*1*1 = 3

      out 0    0    0   -2   +2
  in  | -----------------------
  0   |   I    0    0    0    0
  0   |  hSz   I   JSz  JS-  JS+
  0   |   Sz   0    0    0    0
  -2  |   S+   0    0    0    0
  +2  |   S-   0    0    0    0
  
*/

template <typename T>
class Heisenberg{
public:
  Heisenberg(spinhalf* s, double tE=0, double* dh=nullptr) {_s=s; _tE=tE; _dh=dh;}
  ~Heisenberg(){}

  void addOperators(MPO<T>& A, unsigned site, unsigned r, unsigned c, string op, double val);
  void addOperators(qMPO<T>& A, unsigned site, unsigned r, unsigned c, string op, double val, QN_t Qi, QN_t Qo);

  void buildHam(MPO<T>&  H);
  void buildHam(AutoMPO& ampo, MPO<T>&  H);
  void buildHam(qMPO<T>& H);
  void buildHam(AutoMPO& ampo, qMPO<T>& H);

private:
  double _tE;
  double* _dh;
  spinhalf* _s;
  vector< vector< vector<string> > > ops;
  vector< vector< vector<string> > > val;
};


#endif
