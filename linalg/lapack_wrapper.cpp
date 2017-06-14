#ifndef My_LAPACK_WRAPPERS
#define My_LAPACK_WRAPPERS

#include "lapack_wrapper.h"

typedef Eigen::MatrixXd Mxd;

void diag(double* M, double* evals, int nn)
{
	int N = nn;
    char jobz = 'V';
    char uplo = 'U';
    int lwork = std::max(1,3*N-1);//std::max(1, 1+6*N+2*N*N);
    double *work = new double [lwork];
    int info;
	dsyev_(&jobz,&uplo,&N,M,&N,evals,work,&lwork,&info);
	if(info!=0)
	{
		std::cout<<"dsyev error info is "<<info<<std::endl;
	}
	delete [] work;
}

void QR(const Mxd& MT, Mxd& Q){
	///////////////////////////////////
	// INITIALIZATION
	int M = MT.rows();
	int N = MT.cols();

	Mxd tQ = MT;

	int LDA, K;
	LDA = std::max(1,M);
	K = std::min(M,N);

	double * TAU = new double [K];
	int lwork = std::max(M,N)*std::max(M,N);//The dimension of the array WORK.
	double * work = new double [lwork];
	int INFO;
	///////////////////////////////////
	dgeqrf_(&M, &N, tQ.data(), &LDA, TAU, work, &lwork, &INFO);
	if (INFO!=0) {
		std::cout<<"Illegal value at "<<INFO<<std::endl;
	}
	if (N>M) {
		N = M;
		dorgqr_(&M, &N, &K, tQ.data(), &LDA, TAU, work, &lwork, &INFO);
		if (INFO!=0) {
			std::cout<<"QR Failed!\n";
		}
		Q=tQ.block(0,0,M,N);
	} else {
		dorgqr_(&M, &N, &K, tQ.data(), &LDA, TAU, work, &lwork, &INFO);
		if (INFO!=0) {
			std::cout<<"QR Failed!\n";
		}
		Q=tQ;
	}
	///////////////////////////////////
	//Free Space
	delete [] TAU;
	delete [] work;
}

void QR(const Mxc& MT, Mxc& Q){
	///////////////////////////////////
	// INITIALIZATION
	int M = MT.rows();
	int N = MT.cols();

	Mxc tQ = MT;

	int LDA, K;
	LDA = std::max(1,M);
	K = std::min(M,N);

	std::complex<double> * TAU = new std::complex<double> [K];
	int lwork = std::max(M,N)*std::max(M,N);//The dimension of the array WORK.
	std::complex<double> * work = new std::complex<double> [lwork];
	int INFO;
	///////////////////////////////////
	zgeqrf_(&M, &N, tQ.data(), &LDA, TAU, work, &lwork, &INFO);
	if (INFO!=0) {
		std::cout<<"Illegal value at "<<INFO<<std::endl;
	}
	if (N>M) {
		N = M;
		zungqr_(&M, &N, &K, tQ.data(), &LDA, TAU, work, &lwork, &INFO);
		if (INFO!=0) {
			std::cout<<"QR Failed!\n";
		}
		Q=tQ.block(0,0,M,N);
	}else {
		zungqr_(&M, &N, &K, tQ.data(), &LDA, TAU, work, &lwork, &INFO);
		if (INFO!=0) {
			std::cout<<"QR Failed!\n";
		}
		Q=tQ;
	}
	///////////////////////////////////
	//Free Space
	delete [] TAU;
	delete [] work;
}


void SVD(const Mxd& M, int ds, double * sv, Mxd& UM, Mxd& VTM, char direction)
{
	///////////////////////////////////
	// INITIALIZATION
	char jobU; //all M columns of U are returned in array U
	char jobV; //all N rows of V^T are returned in the array VT
	if(direction=='r')
	{
		jobU = 'S'; //all M columns of U are returned in array U
		jobV = 'N'; //all N rows of V^T are returned in the array VT
	}else if(direction=='l')
	{
		jobU = 'N'; //all M columns of U are returned in array U
		jobV = 'S'; //all N rows of V^T are returned in the array VT
	}
	int rowN = M.rows();
	int colN = M.cols();
	int LDA = rowN;//The leading dimension of the array A
	int DofS = std::min(rowN, colN);
	if (ds<DofS) {
		std::cout<<"Please assign more space to the singular value container!\n";
		exit(1);
	}
	for (int i = 0; i < ds; i++) {
		sv[i]=0;
	}

	UM.setZero(rowN,DofS);
	VTM.setZero(DofS,colN);
	int LDU = rowN;//The leading dimension of the array U
	int LDV = DofS;//The leading dimension of the array VT
	int lwork = std::max(4*std::min(rowN,colN)+std::max(rowN,colN),6*std::min(rowN,colN));
	double * work = new double [lwork];
	int INFO;
	Mxd tp = M;
	///////////////////////////////////
	// CALL LAPACK
	dgesvd_(&jobU, &jobV, &rowN, &colN, tp.data(), &LDA, sv, UM.data(), &LDU, VTM.data(), &LDV, work, &lwork, &INFO);
	///////////////////////////////////
	if (INFO!=0) {
		std::cout<<"SVD Failed!\n";
	}
	///////////////////////////////////
	// No truncation
	if(direction=='r')
	{
		VTM.noalias() = UM.transpose() * M;
	}else if(direction=='l')
	{
		UM.noalias() = M * VTM.transpose();
	}
	///////////////////////////////////
	//Free Space
	delete [] work;
}

void SVD(const Mxd& M, int ds, double * sv, Mxd& UM, Mxd& VTM, int truncD, double& svd_error)
{
	///////////////////////////////////
	// INITIALIZATION
	char jobU = 'S'; //all M columns of U are returned in array U
	char jobV = 'N'; //all N rows of V^T are returned in the array VT

	int rowN = M.rows();
	int colN = M.cols();
	int LDA = rowN;//The leading dimension of the array A

	int DofS;
	if (rowN >= colN) {
		DofS = colN;
	}else {
		DofS = rowN;
	}

	if (ds<DofS) {
		std::cout<<"Please assign more space to the singular value container!\n";
		exit(1);
	}
	for (int i = 0; i < ds; i++) {
		sv[i]=0.0;
	}

	UM.setZero(rowN,DofS);
	VTM.setZero(DofS,colN);

	int LDU = rowN;//The leading dimension of the array U

	int LDV = DofS;//The leading dimension of the array VT

	int lwork = std::max(3*std::min(rowN,colN)+std::max(rowN,colN),5*std::min(rowN,colN));

	double * work = new double [lwork];

	int INFO;

	Mxd tp = M;
	Mxd temp;
	///////////////////////////////////
	// CALL LAPACK
	dgesvd_(&jobU, &jobV, &rowN, &colN, tp.data(), &LDA, sv, UM.data(), &LDU, VTM.data(), &LDV, work, &lwork, &INFO);
	///////////////////////////////////
	if (INFO!=0) {
		std::cout<<"SVD Failed!\n";
	}
	///////////////////////////////////
	// Truncate
	if(truncD>=DofS) // do nothing
	{
		VTM.noalias() = UM.transpose() * M;
	}else if(truncD>0)
	{
		svd_error = 0.0;
		for(int i = truncD; i < DofS; i++)
		{
			svd_error += sv[i]*sv[i];
		}
		temp=UM.block(0,0,M.rows(),truncD);
		UM=temp;
		VTM.noalias() = UM.transpose() * M;
	}
	///////////////////////////////////
	//Free Space
	delete [] work;
}

void SVD(const Mxd& M, std::vector<double>& sv, Mxd& UM, Mxd& VTM, double& svd_cutoff, bool dry_run)
{
  ///////////////////////////////////
  // INITIALIZATION
  char jobU = 'S'; //all M columns of U are returned in array U
  char jobV = 'N'; //all N rows of V^T are returned in the array VT

  int rowN = M.rows();
  int colN = M.cols();
  int LDA = rowN;//The leading dimension of the array A

  int DofS;
  if (rowN >= colN) {
    DofS = colN;
  }else {
    DofS = rowN;
  }

  sv.resize(DofS,0);

  UM.setZero(rowN,DofS);
  VTM.setZero(DofS,colN);

  int LDU = rowN;//The leading dimension of the array U

  int LDV = DofS;//The leading dimension of the array VT

  int lwork = std::max(3*std::min(rowN,colN)+std::max(rowN,colN),5*std::min(rowN,colN));

  double * work = new double [lwork];

  int INFO;

  Mxd tp = M;
  Mxd temp;
  ///////////////////////////////////
  // CALL LAPACK
  dgesvd_(&jobU, &jobV, &rowN, &colN, tp.data(), &LDA, sv.data(), UM.data(), &LDU, VTM.data(), &LDV, work, &lwork, &INFO);
  ///////////////////////////////////
  if (INFO!=0) {
    std::cout<<"SVD Failed!\n";
  }
  ///////////////////////////////////
  // Truncate
	double sv_all = 0.0;
	for(int i = 0; i < DofS; i++){
		sv_all += sv[i]*sv[i];
	}
  double sv_sum = 0.0;
  int new_D = 0;
  for(int i = 0; i < DofS; i++){
      sv_sum += sv[i]*sv[i]/sv_all;
      ++new_D;
      if((1.0-sv_sum)<svd_cutoff) break;
			// if(sv[i]/sv[0]<svd_tol) break;
  }
  if(dry_run){
      // std::cout<<"#(SV) would drop from "<<DofS<<" to "<<new_D<<std::endl;
      VTM.noalias() = UM.transpose() * M;
  }else{
      // std::cout<<"#(SV) droped from "<<DofS<<" to "<<new_D<<std::endl;
      temp=UM.block(0,0,M.rows(),new_D);
      UM=temp;
      VTM.noalias() = UM.transpose() * M;
			VTM /= std::sqrt(sv_sum);
  }
  ///////////////////////////////////
  //Free Space
  delete [] work;
}


void SVD(const Mxc& M, int ds, double * sv, Mxc& UM, Mxc& VTM, char direction)
{
	///////////////////////////////////
	// INITIALIZATION
	char jobU; //all M columns of U are returned in array U
	char jobV; //all N rows of V^T are returned in the array VT
	if(direction=='r')
	{
		jobU = 'S'; //all M columns of U are returned in array U
		jobV = 'N'; //all N rows of V^T are returned in the array VT
	}else if(direction=='l')
	{
		jobU = 'N'; //all M columns of U are returned in array U
		jobV = 'S'; //all N rows of V^T are returned in the array VT
	}
	int rowN = M.rows();
	int colN = M.cols();
	int LDA = rowN;//The leading dimension of the array A
	int DofS = std::min(rowN, colN);
	if (ds<DofS) {
		std::cout<<"Please assign more space to the singular value container!\n";
		exit(1);
	}
	for (int i = 0; i < ds; i++) {
		sv[i]=0;
	}

	UM.setZero(rowN,DofS);
	VTM.setZero(DofS,colN);
	int LDU = rowN;//The leading dimension of the array U
	int LDV = DofS;//The leading dimension of the array VT
	int lwork = std::max(4*std::min(rowN,colN)+std::max(rowN,colN),6*std::min(rowN,colN));
	std::complex<double> * work = new std::complex<double> [lwork];
	double * rwork = new double [5*DofS];
	int INFO;
	Mxc tp = M;
	///////////////////////////////////
	// CALL LAPACK
	zgesvd_(&jobU, &jobV, &rowN, &colN, tp.data(), &LDA, sv, UM.data(), &LDU, VTM.data(), &LDV, work, &lwork, rwork, &INFO);
	///////////////////////////////////
	if (INFO!=0) {
		std::cout<<"SVD Failed!\n";
	}
	///////////////////////////////////
	// No truncation
	if(direction=='r')
	{
		VTM.noalias() = UM.adjoint() * M;
	}else if(direction=='l')
	{
		UM.noalias() = M * VTM.adjoint();
	}
	///////////////////////////////////
	//Free Space
	delete [] work;
	delete [] rwork;
}

void SVD(const Mxc& M, int ds, double * sv, Mxc& UM, Mxc& VTM, int truncD, double& svd_error)
{
	///////////////////////////////////
	// INITIALIZATION
	char jobU = 'S'; //all M columns of U are returned in array U
	char jobV = 'N'; //all N rows of V^T are returned in the array VT

	int rowN = M.rows();
	int colN = M.cols();
	int LDA = rowN;//The leading dimension of the array A

	int DofS;
	if (rowN >= colN) {
		DofS = colN;
	}else {
		DofS = rowN;
	}

	if (ds<DofS) {
		std::cout<<"Please assign more space to the singular value container!\n";
		exit(1);
	}
	for (int i = 0; i < ds; i++) {
		sv[i]=0.0;
	}

	UM.setZero(rowN,DofS);
	VTM.setZero(DofS,colN);

	int LDU = rowN;//The leading dimension of the array U
	int LDV = DofS;//The leading dimension of the array VT
	int lwork = std::max(3*std::min(rowN,colN)+std::max(rowN,colN),5*std::min(rowN,colN));
	std::complex<double> * work = new std::complex<double> [lwork];
	double * rwork = new double [5*DofS];
	int INFO;
	Mxc tp = M;
	Mxc temp;
	///////////////////////////////////
	// CALL LAPACK
	zgesvd_(&jobU, &jobV, &rowN, &colN, tp.data(), &LDA, sv, UM.data(), &LDU, VTM.data(), &LDV, work, &lwork, rwork, &INFO);
	///////////////////////////////////
	if (INFO!=0) {
		std::cout<<"SVD Failed!\n";
	}
	///////////////////////////////////
	// Truncate
	if(truncD>=DofS) // do nothing
	{
		VTM.noalias() = UM.adjoint() * M;
	}else if(truncD>0)
	{
		svd_error = 0.0;
		for(int i = truncD; i < DofS; i++)
		{
			svd_error += sv[i]*sv[i];
		}
		temp=UM.block(0,0,M.rows(),truncD);
		UM=temp;
		VTM.noalias() = UM.adjoint() * M;
	}
	///////////////////////////////////
	//Free Space
	delete [] work;
	delete [] rwork;
}

void SVD(const Mxc& M, std::vector<double>& sv, Mxc& UM, Mxc& VTM, double& svd_cutoff, bool dry_run)
{
  ///////////////////////////////////
  // INITIALIZATION
  char jobU = 'S'; //all M columns of U are returned in array U
  char jobV = 'N'; //all N rows of V^T are returned in the array VT

  int rowN = M.rows();
  int colN = M.cols();
  int LDA = rowN;//The leading dimension of the array A

  int DofS;
  if (rowN >= colN) {
    DofS = colN;
  }else {
    DofS = rowN;
  }

  sv.resize(DofS,0);

  UM.setZero(rowN,DofS);
  VTM.setZero(DofS,colN);

  int LDU = rowN;//The leading dimension of the array U
  int LDV = DofS;//The leading dimension of the array VT
  int lwork = std::max(3*std::min(rowN,colN)+std::max(rowN,colN),5*std::min(rowN,colN));
  std::complex<double> * work = new std::complex<double> [lwork];
	double * rwork = new double [5*DofS];
  int INFO;
  Mxc tp = M;
  Mxc temp;
  ///////////////////////////////////
  // CALL LAPACK
  zgesvd_(&jobU, &jobV, &rowN, &colN, tp.data(), &LDA, sv.data(), UM.data(), &LDU, VTM.data(), &LDV, work, &lwork, rwork, &INFO);
  ///////////////////////////////////
  if (INFO!=0) {
    std::cout<<"SVD Failed!\n";
  }
  ///////////////////////////////////
  // Truncate
	double sv_all = 0.0;
	for(int i = 0; i < DofS; i++){
		sv_all += sv[i]*sv[i];
	}
  double sv_sum = 0.0;
  int new_D = 0;
  for(int i = 0; i < DofS; i++){
      sv_sum += sv[i]*sv[i]/sv_all;
      ++new_D;
      if((1.0-sv_sum)<svd_cutoff) break;
			// if(sv[i]/sv[0]<svd_tol) break;
  }
  if(dry_run){
      // std::cout<<"#(SV) would drop from "<<DofS<<" to "<<new_D<<std::endl;
      VTM.noalias() = UM.adjoint() * M;
  }else{
      // std::cout<<"#(SV) droped from "<<DofS<<" to "<<new_D<<std::endl;
      temp=UM.block(0,0,M.rows(),new_D);
      UM=temp;
      VTM.noalias() = UM.adjoint() * M;
			VTM /= std::sqrt(sv_sum);
  }
  ///////////////////////////////////
  //Free Space
  delete [] work;
	delete [] rwork;
}


void linearSolver(int ord, Mxd& A, Mxd& vec)
{
	// std::cout<<"In linear solver"<<std::endl;
	char uplo = 'U';
	int n = ord;
	int nrhs = 1;
	int * ipiv = new int [2*n](); //del
	int info;
	dsytrs_(&uplo,&n,&nrhs,A.data(),&n,ipiv,vec.data(),&n,&info);
	if(info!=0)
	{
		std::cout<<"dsytrs error! The "<<-info<<"-th argument had an illegal value!"<<std::endl;
		abort();
	}
	// for(int i = 0; i < n; ++i)
	// {
	// 	if(std::isnan(vec[i]))
	// 	{
	// 		std::cout<<"dstrs error! Invalid solution!"<<std::endl;
	// 		abort();
	// 	}
	// }
	// std::cout<<"Out linear solver"<<std::endl;
	delete [] ipiv;
}

void indef_linearSolver(int ord, Mxd& A, Mxd& vec)
{
	// std::cout<<"In linear solver"<<std::endl;
	char uplo = 'L';
	int nn = ord;
	int nrhs = 1;
	int * ipiv = new int [nn]();
	int info;
	int lwork = 16*nn;
	double* work = new double [lwork];
	dsysv_(&uplo,&nn,&nrhs,A.data(),&nn,ipiv,vec.data(),&nn,work,&lwork,&info);
	if(info!=0)
	{
		std::cout<<"dsysv error! The "<<-info<<"-th argument had an illegal value!"<<std::endl;
		abort();
	}
	// for(int i = 0; i < n; ++i)
	// {
	// 	if(std::isnan(vec[i]))
	// 	{
	// 		std::cout<<"dsysv error! Invalid solution!"<<std::endl;
	// 		abort();
	// 	}
	// }
	delete [] work;
	delete [] ipiv;
	// std::cout<<"Out linear solver"<<std::endl;
}


int def_linearSolver(int ord, Mxd& A, Mxd& vec)
{
	// std::cout<<"In linear solver"<<std::endl;
	char uplo = 'U';
	int n = ord;
	int nrhs = 1;
	int info;
	dposv_(&uplo,&n,&nrhs,A.data(),&n,vec.data(),&n,&info);
	if(info!=0)
	{
		if(info<0)
		{
			std::cout<<"dposv error! The "<<-info<<"-th argument had an illegal value!"<<std::endl;
			abort();
		}
		if(info>0)
		{
			 std::cout<<"dposv error! Matrix is exactly singular!"<<std::endl;
			 return -1;
		}
	}
	// for(int i = 0; i < n; ++i)
	// {
	// 	if(std::isnan(vec[i]))
	// 	{
	// 		std::cout<<"dstrs error! Invalid solution!"<<std::endl;
	// 		abort();
	// 	}
	// }
	// std::cout<<"Out linear solver"<<std::endl;
	return 0;
}

void geindef_linearSolver(int ord, Mxd& A, Mxd& vec)
{
	// std::cout<<"In linear solver"<<std::endl;
	int nn = ord;
	int nrhs = 1;
	int * ipiv = new int [nn]();
	int info;
	dgesv_(&nn,&nrhs,A.data(),&nn,ipiv,vec.data(),&nn,&info);
	if(info!=0)
	{
		std::cout<<"dgesv error! The "<<-info<<"-th argument had an illegal value!"<<std::endl;
		abort();
	}
	// for(int i = 0; i < n; ++i)
	// {
	// 	if(std::isnan(vec[i]))
	// 	{
	// 		std::cout<<"dsysv error! Invalid solution!"<<std::endl;
	// 		abort();
	// 	}
	// }
	delete [] ipiv;
}

int sdef_linearSolver(int ord, Mxd& A, Mxd& vec)
{
	// std::cout<<"In linear solver"<<std::endl;
	char uplo = 'U';
	int n = ord;
	int nrhs = 1;
	int info;
	// Mxd x = vec;
	double* x = new double [n];
	for(int i = 0; i < n; ++i)
	{
		x[i] = vec.data()[i];
	}
	double* work = new double [n];
	float* swork = new float [n*(n+nrhs)];
	int iter;
	dsposv_(&uplo,&n,&nrhs,A.data(),&n,x,&n,vec.data(),&n,work,swork,&iter,&info);
	if(info!=0)
	{
		if(info<0)
		{
			std::cout<<"dsposv error! The "<<-info<<"-th argument had an illegal value!"<<std::endl;
			abort();
		}
		if(info>0)
		{
			 std::cout<<"dsposv error! Matrix is exactly singular!"<<std::endl;
			 return -1;
		}
	}
	delete [] x;
	delete [] work;
	delete [] swork;
	return 0;
}

void sindef_linearSolver(int ord, Mxd& A, Mxd& vec)
{
	char fact = 'N';
	char uplo = 'U';
	int  n = ord;
	int nrhs = 1;
	double* af = new double [n*n];
	int* ipiv = new int [n];
	Mxd B = vec;
	double rcond;
	double ferr;
	double berr;
	int lwork = std::max (1,12*n);
	double* work = new double [lwork];
	int* iwork = new int [n];
	int info;

	dsysvx_(&fact,&uplo,&n,&nrhs,A.data(),&n,af,&n,ipiv,B.data(),&n,vec.data(),&n,&rcond,&ferr,&berr,work,&lwork,iwork,&info);

	if(info!=0)
	{
		if(info<0)
		{
			std::cout<<"dsysvx error! The "<<-info<<"-th argument had an illegal value!"<<std::endl;
			abort();
		}
		if(info>0&&info<n)
		{
			std::cout<<"dsysvx error! Matrix is exactly singular!"<<std::endl;
			abort();
		}
		if(info==n+1)
		{
			std::cout<<"Condition number smaller than eps!"<<std::endl;
		}
	}
	// std::cout<<"Condition Number: "<<rcond<<std::endl;
	delete [] af;
	delete [] ipiv;
	delete [] work;
	delete [] iwork;
}

void sxdef_linearSolver(int ord, Mxd& A, Mxd& vec)
{
	char fact = 'N';
	char uplo = 'U';
	int  n = ord;
	int nrhs = 1;
	double* af = new double [n*n];
	char equed;
	double* s = new double [n];
	int* ipiv = new int [n];
	Mxd B = vec;
	double rcond;
	double ferr;
	double berr;
	double* work = new double [3*n];
	int* iwork = new int [n];
	int info;

	dposvx_(&fact,&uplo,&n,&nrhs,A.data(),&n,af,&n,&equed,s,B.data(),&n,vec.data(),&n,&rcond,&ferr,&berr,work,iwork,&info);
	if(info!=0)
	{
		if(info<0)
		{
			std::cout<<"dposvx error! The "<<-info<<"-th argument had an illegal value!"<<std::endl;
			abort();
		}
		if(info>0&&info<n)
		{
			std::cout<<"dposvx error! Matrix is exactly singular!"<<std::endl;
			abort();
		}
		if(info==n+1)
		{
			std::cout<<"Condition number smaller than eps!"<<std::endl;
		}
	}
	// std::cout<<"Condition Number: "<<rcond<<std::endl;
}

#endif
