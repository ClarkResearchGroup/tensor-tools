#ifndef My_DMRG_DIAG_CLASS
#define My_DMRG_DIAG_CLASS

#include "dmrg_diag.h"

void dmrg_mtv(int n, double *in, double *out, MPS& psi, MPO& H, Mxd ** CRM, Mxd ** CLM, int site){
	int phy = psi.index_size;
	int L = psi.length;
	int r = psi.M[site][0].rows();
	int c = psi.M[site][0].cols();
	int maxbd = *std::max_element(H.bond_dims.begin(), H.bond_dims.end());
	////////////////////////////////////////////////
	Mxd * inMPS  = new Mxd [phy];
	Mxd * outMPS = new Mxd [phy];
	for(int i = 0; i < phy; i++){
		inMPS[i].setZero(r,c);
		outMPS[i].setZero(r,c);
	}
	for(int i = 0; i < phy; i++){
		for(int j = 0; j < r*c; j++){
			inMPS[i](j) = in[i*r*c+j];
		}
	}
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
				TM1[i][j].noalias() = inMPS[i] * CRM[site+1][j].transpose();
			}
		}
		for(int tid = 0; tid < phy; tid++){
			for(int k = 0; k < phy; k++){
				for(size_t l = 0; l < H.H[site][tid*phy+k].r.size(); l++){
					outMPS[tid] += H.H[site][phy*tid+k].v[l] * TM1[k][H.H[site][tid*phy+k].c[l]];
				}
			}
		}
	}else if(site==L-1){
		for(int i = 0; i < phy; i++){
			for(int j = 0; j < H.bond_dims[site]; j++){
				TM1[i][j].noalias() = CLM[site-1][j] * inMPS[i];
			}
		}
		for(int tid = 0; tid < phy; tid++){
			for(int k = 0; k < phy; k++){
				for(size_t l = 0; l < H.H[site][phy*tid+k].r.size(); l++){
					outMPS[tid] += TM1[k][H.H[site][phy*tid+k].r[l]] * H.H[site][phy*tid+k].v[l];
				}
			}
		}
	}else{
		for(int i = 0; i < phy; i++){
			for(int j = 0; j < H.bond_dims[site+1]; j++){
				TM1[i][j].noalias() = inMPS[i] * CRM[site+1][j].transpose();
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
				outMPS[i].noalias() += CLM[site-1][k] * TM2[i][k];
			}
		}
	}
	////////////////////////////////////////////////
	for(int i = 0; i < phy; i++){
		for(int j = 0; j < r*c; j++){
			out[i*r*c+j] = outMPS[i](j);
		}
	}
	////////////////////////////////////////////////
	for(int i = 0; i < phy; i++){
		delete [] TM1[i];
		delete [] TM2[i];
	}
	delete [] TM1;
	delete [] TM2;
	delete [] inMPS;
	delete [] outMPS;
}

void dsaupd_LE(int n, int nev, double *Evals, double **Evecs, double *initVec, MPS& psi, MPO& H, Mxd ** CRM, Mxd ** CLM, int site){
	int ido = 0;
	char bmat[2] = "I";
	char which[3] = "SA";
	double tol = 0.0;
	// double tol = 1e-10;
	double *resid;
	resid = initVec;
	// resid = new double[n];
	int ncv = 10*nev;
	if (ncv>n) ncv = n;
	double *v;
	int ldv = n;
	v = new double[ldv*ncv];
	int *iparam;
	iparam = new int[11];
	iparam[0] = 1;
	iparam[2] = 100;
	// iparam[2] = 2*n;
	iparam[6] = 1;
	int *ipntr;
	ipntr = new int[11];
	double *workd;
	workd = new double[3*n];
	double *workl;
	workl = new double[ncv*(ncv+8)+16];
	int lworkl = ncv*(ncv+8)+16;
	int info = 1;
	int rvec = 1;  // Changed from above
	int *select;
	select = new int[ncv];
	double *d;
	d = new double[2*ncv];
	double sigma;
	int ierr;

	do {
	dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
	    &ncv, v, &ldv, iparam, ipntr, workd, workl,
	    &lworkl, &info);

	if ((ido==1)||(ido==-1)) dmrg_mtv(n, workd+ipntr[0]-1, workd+ipntr[1]-1, psi, H, CRM, CLM, site);
	} while ((ido==1)||(ido==-1));

	if (info<0) {
	   //  cout << "Error with dsaupd, info = " << info << "\n";
	//	 cout << "Check documentation in dsaupd\n\n";
	} else {
	dseupd_(&rvec, (char *) "A", select, d, v, &ldv, &sigma, bmat,
	    &n, which, &nev, &tol, resid, &ncv, v, &ldv,
	    iparam, ipntr, workd, workl, &lworkl, &ierr);

	if (ierr!=0) {
	//  cout << "Error with dseupd, info = " << ierr << "\n";
	 // cout << "Check the documentation of dseupd.\n\n";
	} else if (info==1) {
	//  cout << "Maximum number of iterations reached.\n\n";
	} else if (info==3) {
	//  cout << "No shifts could be applied during implicit\n";
	//  cout << "Arnoldi update, try increasing NCV.\n\n";
	}

	// cout<<n/10<<" "<<iparam[2]<<endl;

	int i, j;
	for (i=0; i<nev; i++) Evals[i] = d[i];
	for (i=0; i<nev; i++) for (j=0; j<n; j++) Evecs[i][j] = v[i*n+j];

	// delete resid;
	delete [] v;
	delete [] iparam;
	delete [] ipntr;
	delete [] workd;
	delete [] workl;
	delete [] select;
	delete [] d;
	}
}

void dsaupd_HE(int n, int nev, double *Evals, double **Evecs, double *initVec, MPS& psi, MPO& H, Mxd ** CRM, Mxd ** CLM, int site){
	int ido = 0;
	char bmat[2] = "I";
	char which[3] = "LA";
	double tol = 0.0;
	// double tol = 1e-10;
	double *resid;
	resid = initVec;
	// resid = new double[n];
	int ncv = 10*nev;
	if (ncv>n) ncv = n;
	double *v;
	int ldv = n;
	v = new double[ldv*ncv];
	int *iparam;
	iparam = new int[11];
	iparam[0] = 1;
	iparam[2] = 100;
	// iparam[2] = 2*n;
	iparam[6] = 1;
	int *ipntr;
	ipntr = new int[11];
	double *workd;
	workd = new double[3*n];
	double *workl;
	workl = new double[ncv*(ncv+8)+16];
	int lworkl = ncv*(ncv+8)+16;
	int info = 1;
	int rvec = 1;  // Changed from above
	int *select;
	select = new int[ncv];
	double *d;
	d = new double[2*ncv];
	double sigma;
	int ierr;

	do {
	dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
	    &ncv, v, &ldv, iparam, ipntr, workd, workl,
	    &lworkl, &info);

	if ((ido==1)||(ido==-1)) dmrg_mtv(n, workd+ipntr[0]-1, workd+ipntr[1]-1, psi, H, CRM, CLM, site);
	} while ((ido==1)||(ido==-1));

	if (info<0) {
	   //  cout << "Error with dsaupd, info = " << info << "\n";
	//	 cout << "Check documentation in dsaupd\n\n";
	} else {
	dseupd_(&rvec, (char *) "A", select, d, v, &ldv, &sigma, bmat,
	    &n, which, &nev, &tol, resid, &ncv, v, &ldv,
	    iparam, ipntr, workd, workl, &lworkl, &ierr);

	if (ierr!=0) {
	//  cout << "Error with dseupd, info = " << ierr << "\n";
	 // cout << "Check the documentation of dseupd.\n\n";
	} else if (info==1) {
	//  cout << "Maximum number of iterations reached.\n\n";
	} else if (info==3) {
	//  cout << "No shifts could be applied during implicit\n";
	//  cout << "Arnoldi update, try increasing NCV.\n\n";
	}

	// cout<<n/10<<" "<<iparam[2]<<endl;

	int i, j;
	for (i=0; i<nev; i++) Evals[i] = d[i];
	for (i=0; i<nev; i++) for (j=0; j<n; j++) Evecs[i][j] = v[i*n+j];

	// delete resid;
	delete [] v;
	delete [] iparam;
	delete [] ipntr;
	delete [] workd;
	delete [] workl;
	delete [] select;
	delete [] d;
	}
}

#endif
