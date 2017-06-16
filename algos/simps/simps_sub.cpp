#ifndef MY_SIMPS_SUB_ROUTINES
#define MY_SIMPS_SUB_ROUTINES

#ifdef USE_CG
#include "cg.h"
#endif

#include "simps_sub.h"

void inv_buildR(MPS& psi, MPO& H, Mxd ** CRM){
	int tid;
	int r1, c1;
	int phy = psi.index_size;
	int L = psi.length;
	int maxbd = *std::max_element(H.bond_dims.begin(), H.bond_dims.end());
	////////////////////////////////////////////
	Mxd ** TM1 = new Mxd * [phy];
	Mxd ** TM2 = new Mxd * [phy];
	for(int i = 0; i < phy; i++){
		TM1[i] = new Mxd [maxbd];
		TM2[i] = new Mxd [maxbd];
	}
	////////////////////////////////////////////////
	// Last site //
	// Building the FR //
	r1=psi.M[L-1][0].rows();
	c1=psi.M[L-1][0].cols();
	for(int i = 0; i < H.bond_dims[L-1]; i++){
		CRM[L-1][i].setZero(r1,r1);
	}
	for(int i = 0; i < phy; i++){
		for(int k = 0; k < phy; k++){
			for(size_t j = 0; j < H.H[L-1][phy*i+k].r.size(); j++){
				CRM[L-1][H.H[L-1][phy*i+k].r[j]].noalias() += H.H[L-1][phy*i+k].v[j]*(psi.M[L-1][i]*psi.M[L-1][k].transpose());
			}
		}
	}
	// Last-but-1 to the second site //
	for (int i = L-2; i > 0; i--) {
		// cout<<"Site "<<i<<endl;
		r1=psi.M[i][0].rows();
		c1=psi.M[i][0].cols();
		// Building the FR //
		for(int j = 0; j < H.bond_dims[i]; j++){
			CRM[i][j].setZero(r1,r1);
			for(int k = 0; k < phy; k++){
				TM1[k][j].setZero(c1,r1);
				TM2[k][j].setZero(c1,r1);
			}
		}
		for(tid = 0; tid < phy; tid++){
			for(int j = 0; j < H.bond_dims[i+1]; j++){
				TM1[tid][j].noalias() = CRM[i+1][j] * psi.M[i][tid].transpose();
			}
		}
		for(tid = 0; tid < phy; tid++){
			for(int k = 0; k < phy; k++){
				for(size_t l = 0; l < H.H[i][phy*tid+k].r.size(); l++){
					TM2[tid][H.H[i][phy*tid+k].r[l]] += H.H[i][phy*tid+k].v[l] * TM1[k][H.H[i][phy*tid+k].c[l]];
				}
			}
		}
		for(int j = 0; j < H.bond_dims[i]; j++){
			for(int k = 0; k < phy; k++){
				CRM[i][j].noalias() += psi.M[i][k] * TM2[k][j];
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
}

void buildMatrix(MPS& psi, MPO& H, Mxd ** CR, Mxd ** CL, const int& site, Mxd& A){
	int phy  = psi.index_size;
	int L    = psi.length;
	int r    = psi.bond_dims[site];
	int c    = psi.bond_dims[site+1];
	int n    = phy * r * c;
	A.setZero(n,n);
	////////////////////////////////////////////////
	// build matrix of linear equation
	if(site==0){
		for(int i = 0; i < phy; ++i){
			for(int j = 0; j < phy; ++j){
				for(size_t k = 0; k < H.H[site][i*phy+j].v.size(); ++k){
					A.block(i*c,j*c,c,c) += H.H[site][i*phy+j].v[k] * CR[site+1][H.H[site][i*phy+j].c[k]];
				}
			}
		}
	}else if(site==L-1){
		for(int i = 0; i < phy; ++i){
			for(int j = 0; j < phy; ++j){
				for(size_t k = 0; k < H.H[site][i*phy+j].v.size(); ++k){
					A.block(i*r,j*r,r,r) += H.H[site][i*phy+j].v[k] * CL[site-1][H.H[site][i*phy+j].r[k]];
				}
			}
		}
	}else{
		for(int tid = 0; tid < phy*phy; ++tid){
			int i = tid/phy;
			int j = tid%phy;
			for(int p2 = 0; p2 < c; ++p2){
				for(int p1 = 0; p1 < c; ++p1){
					for(size_t k = 0; k < H.H[site][tid].v.size(); ++k){
						A.block(i*r*c+p1*r,j*r*c+p2*r,r,r) += H.H[site][tid].v[k] * CL[site-1][H.H[site][tid].r[k]] * CR[site+1][H.H[site][tid].c[k]](p1,p2);
					}
				}
			}
		}
	}
}

void buildDiag(MPS& psi, MPO& H, Mxd ** CR, Mxd ** CL, const int& site, Mxd& A){
	int phy  = psi.index_size;
	int L    = psi.length;
	int r    = psi.bond_dims[site];
	int c    = psi.bond_dims[site+1];
	int n    = phy * r * c;
	A.setZero(n,1);
	////////////////////////////////////////////////
	// build diagonal of the matrix of linear equation
	if(site==0){
		for(int i = 0; i < phy; ++i){
			for(size_t k = 0; k < H.H[site][i*phy+i].v.size(); ++k){
				A.block(i*c,0,c,1) += H.H[site][i*phy+i].v[k] * CR[site+1][H.H[site][i*phy+i].c[k]].diagonal();
			}
		}
	}else if(site==L-1){
		for(int i = 0; i < phy; ++i){
			for(size_t k = 0; k < H.H[site][i*phy+i].v.size(); ++k){
				A.block(i*r,0,r,1) += H.H[site][i*phy+i].v[k] * CL[site-1][H.H[site][i*phy+i].r[k]].diagonal();
			}
		}
	}else{
		for(int i = 0; i < phy; ++i){
			for(int p1 = 0; p1 < c; ++p1){
				for(size_t k = 0; k < H.H[site][i*phy+i].v.size(); ++k){
					A.block(i*r*c+p1*r,0,r,1) += H.H[site][i*phy+i].v[k] * CL[site-1][H.H[site][i*phy+i].r[k]].diagonal() * CR[site+1][H.H[site][i*phy+i].c[k]](p1,p1);
				}
			}
		}
	}
}

void inv_updateSite(const char& direction, MPS& psi, MPO& H, MPO& HS, Mxd ** CR1, Mxd ** CL1, Mxd ** CR2, Mxd ** CL2, const int& site){
	int phy = psi.index_size;
	int L   = psi.length;
	int r   = psi.bond_dims[site];
	int c   = psi.bond_dims[site+1];
	int n = phy * r * c;
	////////////////////////////////////////////////
	Mxd A, B; // Linear equation problem: A*x = B*u. Close to convergence, x is close to u.
	Mxd vec, init_vec;
	////////////////////////////////////////////////
	vec.setZero(r,phy*c);
	init_vec.setZero(r,phy*c);
	////////////////////////////////////////////////
	for(int tid=0; tid<phy; tid++){
		init_vec.block(0,tid*c,r,c) = psi.M[site][tid].block(0,0,r,c);
	}
	double norm = psiHphi(psi,H,psi);
	init_vec /= norm;
	////////////////////////////////////////////////
	// buildMatrix(psi, HS, CR2, CL2, site, A);
	// buildMatrix(psi, H,  CR1, CL1, site, B);
	////////////////////////////////////////////////
	// build vector of linear equation
	if(site==0){
		for(int i = 0; i < phy; ++i){
			for(int j = 0; j < phy; ++j){
				for(size_t k = 0; k < H.H[site][i*phy+j].v.size(); ++k){
					vec.block(0,i*c,1,c) += H.H[site][i*phy+j].v[k] * (psi.M[site][j] * CR1[site+1][H.H[site][i*phy+j].c[k]].transpose());
				}
			}
		}
	}else if(site==L-1){
		for(int i = 0; i < phy; ++i){
			for(int j = 0; j < phy; ++j){
				for(size_t k = 0; k < H.H[site][i*phy+j].v.size(); ++k){
					vec.block(0,i,r,1) += H.H[site][i*phy+j].v[k] * (CL1[site-1][H.H[site][i*phy+j].r[k]] * psi.M[site][j]);
				}
			}
		}
	}else
	{
		for(int i = 0; i < phy; ++i){
			for(int j = 0; j < phy; ++j){
				for(size_t k = 0; k < H.H[site][i*phy+j].v.size(); ++k){
					vec.block(0,i*c,r,c) += H.H[site][i*phy+j].v[k] * (CL1[site-1][H.H[site][i*phy+j].r[k]] * psi.M[site][j] * CR1[site+1][H.H[site][i*phy+j].c[k]].transpose());
				}
			}
		}
	}
	////////////////////////////////////////////////
#ifdef USE_Def
	buildMatrix(psi, HS, CR2, CL2, site, A);
	def_linearSolver(n, A, vec);
#endif

#ifdef USE_sDef
	buildMatrix(psi, HS, CR2, CL2, site, A);
	sdef_linearSolver(n, A, vec);
#endif

#ifdef USE_Indef
	buildMatrix(psi, HS, CR2, CL2, site, A);
	indef_linearSolver(n, A, vec);
#endif

#ifdef USE_sIndef
	buildMatrix(psi, HS, CR2, CL2, site, A);
	sindef_linearSolver(n, A, vec);
#endif

#ifdef USE_CG
	// Jacobi preconidtioner can be applied
	// CG_linearSolver(n, HS, vec, init_vec, CR2, CL2, site);
	// buildMatrix(psi, HS, CR2, CL2, site, A);
	// Mxd dg = A.diagonal();
	Mxd dg;
	buildDiag(psi, HS, CR2, CL2, site, dg);
	for(int i = 0; i < n; ++i)
	{
		if(dg(i)!=0)
		  dg(i) = 1.0/dg(i);
		else
		  dg(i) = 1e30;
	}
	CG_linearSolver(n, HS, dg, vec, init_vec, CR2, CL2, site);
#endif
	vec.normalize();
	////////////////////////////////////////////////
	for(int tid=0; tid<phy; tid++){
		psi.M[site][tid].block(0,0,r,c) = vec.block(0,tid*c,r,c);
	}
	if(direction == 'r'){
		psi.moveRight(site);
	}else if(direction == 'l'){
		psi.moveLeft(site);
	}
}

void inv_updateLR(const char& direction, MPS& psi, MPO& H, Mxd ** CR, Mxd ** CL, const int& site){
	int r1, c1;
	int phy = psi.index_size;
	int L = psi.length;
	int maxbd = *std::max_element(H.bond_dims.begin(), H.bond_dims.end());
	////////////////////////////////////////////
	Mxd ** TM1 = new Mxd * [phy];
	Mxd ** TM2 = new Mxd * [phy];
	for(int i = 0; i < phy; i++){
		TM1[i] = new Mxd [maxbd];
		TM2[i] = new Mxd [maxbd];
	}
	////////////////////////////////////////////////
	r1=psi.M[site][0].rows();
	c1=psi.M[site][0].cols();
	if(direction == 'r'){
		if(site==0){
			for(int j = 0; j < H.bond_dims[1]; j++){
				CL[0][j].setZero(c1,c1);
			}
			for(int i = 0; i < phy; i++){
				for(int k = 0; k < phy; k++){
					for(size_t j = 0; j < H.H[site][i*phy+k].r.size(); j++){
						CL[0][H.H[site][phy*i+k].c[j]].noalias() += H.H[site][phy*i+k].v[j] * (psi.M[0][i].transpose()*psi.M[0][k]);
					}
				}
			}
		}else{
			for(int j = 0; j < H.bond_dims[site+1]; j++){
				CL[site][j].setZero(c1,c1);
				for(int k = 0; k < phy; k++){
					TM1[k][j].setZero(r1,c1);
					TM2[k][j].setZero(r1,c1);
				}
			}
			for(int tid = 0; tid < phy; tid++){
				for(int j = 0; j < H.bond_dims[site]; j++){
					TM1[tid][j].noalias() = CL[site-1][j] * psi.M[site][tid];
				}
			}
			for(int tid = 0; tid < phy; tid++){
				for(int k = 0; k < phy; k++){
					for(size_t l = 0; l < H.H[site][tid*phy+k].r.size(); l++){
						TM2[tid][H.H[site][tid*phy+k].c[l]] += TM1[k][H.H[site][tid*phy+k].r[l]] * H.H[site][phy*tid+k].v[l];
					}
				}
			}
			for(int j = 0; j < H.bond_dims[site+1]; j++){
				for(int k = 0; k < phy; k++){
					CL[site][j].noalias() += psi.M[site][k].transpose() * TM2[k][j];
				}
			}
		}
	}else if(direction == 'l'){
		if(site==L-1){
			for(int j = 0; j < H.bond_dims[L-1]; j++){
				CR[L-1][j].setZero(r1,r1);
			}
			for(int i = 0; i < phy; i++){
				for(int k = 0; k < phy; k++){
					for(size_t j = 0; j < H.H[site][i*phy+k].r.size(); j++){
						CR[L-1][H.H[site][phy*i+k].r[j]].noalias() += H.H[site][phy*i+k].v[j] * (psi.M[L-1][i]*psi.M[L-1][k].transpose());
					}
				}
			}
		}else{
			for(int j = 0; j < H.bond_dims[site]; j++){
				CR[site][j].setZero(r1,r1);
				for(int k = 0; k < phy; k++){
					TM1[k][j].setZero(c1,r1);
					TM2[k][j].setZero(c1,r1);
				}
			}
			for(int tid = 0; tid < phy; tid++){
				for(int j = 0; j < H.bond_dims[site+1]; j++){
					TM1[tid][j].noalias() = CR[site+1][j] * psi.M[site][tid].transpose();
				}
			}
			for(int tid = 0; tid < phy; tid++){
				for(int k = 0; k < phy; k++){
					for(size_t l = 0; l < H.H[site][tid*phy+k].r.size(); l++){
						TM2[tid][H.H[site][tid*phy+k].r[l]] += H.H[site][phy*tid+k].v[l] * TM1[k][H.H[site][tid*phy+k].c[l]];
					}
				}
			}
			for(int j = 0; j < H.bond_dims[site+1]; j++){
				for(int k = 0; k < phy; k++){
					CR[site][j].noalias() += psi.M[site][k] * TM2[k][j];
				}
			}
		}
	}
	////////////////////////////////////////////
	for(int i = 0; i < phy; i++){
		delete [] TM1[i];
		delete [] TM2[i];
	}
	delete [] TM1;
	delete [] TM2;
}

#endif
