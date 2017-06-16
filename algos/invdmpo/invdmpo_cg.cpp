#ifndef INVDMPO_CG_SOLVER
#define INVDMPO_CG_SOLVER

// #ifdef USE_CG

#include "invdmpo_cg.h"

double dotVec(Mxd& A, Mxd& B) {
  int size = A.size();

  // double tp = 0;
  // for(int i = 0; i < size; ++i)
  //   {
  //     tp += A.data()[i] * B.data()[i];
  //   }
  // return tp;

  Eigen::Map<Eigen::VectorXd> aVec(A.data(),size);
  Eigen::Map<Eigen::VectorXd> bVec(B.data(),size);
  return aVec.dot(bVec);
}


Mxd simps_mult(MPO& H, Mxd& v, Mxd ** CRM, Mxd ** CLM, const int& site){
	int phy  = std::round(std::sqrt(H.index_size));
	int L    = H.length;
  int maxbd = *std::max_element(H.bond_dims.begin(), H.bond_dims.end());
	int r = v.rows();
	int c = v.cols()/phy;
	Mxd vec(r,c*phy);
	vec.setZero();
	////////////////////////////////////////////////
	Mxd ** TM1 = new Mxd * [phy];
	Mxd ** TM2 = new Mxd * [phy];
	for(int i = 0; i < phy; i++){
		TM1[i] = new Mxd [maxbd];
		TM2[i] = new Mxd [maxbd];
		for(int j = 0; j < maxbd; j++){
			TM1[i][j].setZero(r,c);
			TM2[i][j].setZero(r,c);
		}
	}
	////////////////////////////////////////////////
	if(site==0){
		for(int i = 0; i < phy; i++){
			for(int j = 0; j < H.bond_dims[site+1]; j++){
				TM1[i][j].noalias() = v.block(0,i*c,r,c) * CRM[site+1][j].transpose();
			}
		}
		for(int tid = 0; tid < phy; tid++){
			for(int k = 0; k < phy; k++){
				for(size_t l = 0; l < H.H[site][tid*phy+k].r.size(); l++){
					vec.block(0,tid*c,r,c) += H.H[site][phy*tid+k].v[l] * TM1[k][H.H[site][tid*phy+k].c[l]];
				}
			}
		}
	}else if(site==L-1){
		for(int i = 0; i < phy; i++){
			for(int j = 0; j < H.bond_dims[site]; j++){
				TM1[i][j].noalias() = CLM[site-1][j] * v.block(0,i*c,r,c);
			}
		}
		for(int tid = 0; tid < phy; tid++){
			for(int k = 0; k < phy; k++){
				for(size_t l = 0; l < H.H[site][phy*tid+k].r.size(); l++){
					vec.block(0,tid*c,r,c) += TM1[k][H.H[site][phy*tid+k].r[l]] * H.H[site][phy*tid+k].v[l];
				}
			}
		}
	}else{
		for(int i = 0; i < phy; i++){
			for(int j = 0; j < H.bond_dims[site+1]; j++){
				TM1[i][j].noalias() = v.block(0,i*c,r,c) * CRM[site+1][j].transpose();
			}
		}
		for(int tid = 0; tid < phy; tid++){
			for(int k = 0; k < phy; k++){
				for(size_t l = 0; l < H.H[site][tid*phy+k].r.size(); l++){
					TM2[tid][H.H[site][tid*phy+k].r[l]] += H.H[site][phy*tid+k].v[l] * TM1[k][H.H[site][tid*phy+k].c[l]];
				}
			}
		}
		for(int i = 0; i < phy; i++){
			for(int k = 0; k < H.bond_dims[site]; k++){
				vec.block(0,i*c,r,c) += CLM[site-1][k] * TM2[i][k];
			}
		}
	}
	////////////////////////////////////////////////
	for(int i = 0; i < phy; i++){
		delete [] TM1[i];
		delete [] TM2[i];
	}
	delete [] TM1;
	delete [] TM2;
	return vec;
}

// CG with Jacobi Preconditioner -- (efficient for Diagonal dominant matrix)
int CG_linearSolver(int n, MPO& H, Mxd& M, Mxd& v, Mxd& x, Mxd ** CRM, Mxd ** CLM, const int& site)
{
  int row = v.rows();
  int col = v.cols();
  Eigen::Map<Mxd> Mv(M.data(),v.rows(),v.cols());
  int max_iter = 40;
  double tol = 1e-16;
  double resid, alpha, beta, rho, rho_1, normb;
  normb = v.norm();

  Mxd r(row,col), p(row,col), z(row,col), q(row,col);
  r = v - simps_mult(H,x,CRM,CLM,site);

  if (normb == 0.0)
    normb = 1;

  if ((resid = r.norm() / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (int i = 1; i <= max_iter; i++){
    // std::cout << z.rows() << " " << z.cols() << '\n';
    // std::cout << Mv.rows() << " " << Mv.cols() << '\n';
    // std::cout << r.rows() << " " << r.cols() << '\n';
    z.array() = Mv.array() * r.array();
    // std::cout << "Ah oh" << std::endl;
    rho       = dotVec(r,z);
    if(i==1)
      p = z;
    else{
      beta = rho / rho_1;
      p    = z + beta * p;
    }

    q     = simps_mult(H,p,CRM,CLM,site);
    alpha = rho / dotVec(p,q);
    x    += alpha * p;
    r    -= alpha * q;

    if ((resid = r.norm() / normb) <= tol) {
      tol = resid;
      max_iter = i;
      v = x;
      return 0;
    }

    rho_1 = rho;
  }
  v = x;
  return 1;
}

// #endif

#endif
