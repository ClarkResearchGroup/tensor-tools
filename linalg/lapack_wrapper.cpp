#ifndef My_LAPACK_WRAPPERS
#define My_LAPACK_WRAPPERS

#include "lapack_wrapper.h"

void MAT_VEC(int& m, int& k, int& n, double* A, double* B, double* C){
	char transa = 'N';
	char transb = 'N';
	double alpha = 1;
	double beta  = 1;
	dgemm_(&transa,&transb,&m,&n,&k,&alpha,A,&m,B,&k,&beta,C,&m);
}

void MAT_VEC(int& m, int& k, int& n, std::complex<double>* A, std::complex<double>* B, std::complex<double>* C){
	char transa = 'N';
	char transb = 'N';
	std::complex<double> alpha = 1;
	std::complex<double> beta  = 1;
	zgemm_(&transa,&transb,&m,&n,&k,&alpha,A,&m,B,&k,&beta,C,&m);
}


void DIAG(int m, double* A, double* evals){
	int N = m;
  char jobz = 'V';
  char uplo = 'U';
  int lwork = std::max(1,2*N+N*N);
  double * work = new double [lwork];
  int info;
	dsyev_(&jobz,&uplo,&N,A,&N,evals,work,&lwork,&info);
	if(info!=0){
		std::cout<<"dsyev error info is "<<info<<std::endl;
    delete [] work;
    assert(1==2);
	}
	delete [] work;
}

void DIAG(int m, std::complex<double>* A, double* evals){
	int N = m;
  char jobz = 'V';
  char uplo = 'U';
  int lwork = std::max(1,N+N*N);
  std::complex<double> * work = new std::complex<double> [lwork];
	double * rwork = new double [std::max(1,3*N-2)];
  int info;
	zheev_(&jobz,&uplo,&N,A,&N,evals,work,&lwork,rwork,&info);
	if(info!=0){
		std::cout<<"zheev error info is "<<info<<std::endl;
	}
	delete [] work;
	delete [] rwork;
}

void ORTHO(int r, int c, double* A){
	///////////////////////////////////
	// INITIALIZATION
	int M = r;
	int N = c;

	int LDA, K;
	LDA = std::max(1,M);
	K = std::min(M,N);

	double * TAU = new double [K];
	int lwork = std::max(M,N)*std::max(M,N);//The dimension of the array WORK.
	double * work = new double [lwork];
	int INFO;
	///////////////////////////////////
	dgeqrf_(&M, &N, A, &LDA, TAU, work, &lwork, &INFO);
	if (INFO!=0) {
		std::cout<<"QR illegal value at "<<INFO<<std::endl;
	}
	dorgqr_(&M, &K, &K, A, &LDA, TAU, work, &lwork, &INFO);
	if (INFO!=0) {
		std::cout<<"QR Failed!\n";
	}
	///////////////////////////////////
	//Free Space
	delete [] TAU;
	delete [] work;
}


void ORTHO(int r, int c, std::complex<double>* A){
	///////////////////////////////////
	// INITIALIZATION
	int M = r;
	int N = c;

	int LDA, K;
	LDA = std::max(1,M);
	K = std::min(M,N);

	std::complex<double> * TAU = new std::complex<double> [K];
	int lwork = std::max(M,N)*std::max(M,N);//The dimension of the array WORK.
	std::complex<double> * work = new std::complex<double> [lwork];
	int INFO;
	///////////////////////////////////
	zgeqrf_(&M, &N, A, &LDA, TAU, work, &lwork, &INFO);
	if (INFO!=0) {
		std::cout<<"QR illegal value at "<<INFO<<std::endl;
	}
	zungqr_(&M, &K, &K, A, &LDA, TAU, work, &lwork, &INFO);
	if (INFO!=0) {
		std::cout<<"QR Failed!\n";
	}
	///////////////////////////////////
	//Free Space
	delete [] TAU;
	delete [] work;
}

void QR(int r, int c, double* A, double* Q, double* R){
	///////////////////////////////////
	// INITIALIZATION
	int M = r;
	int N = c;
	double* B = new double [M*N];
	std::copy(A, A+M*N, B);

	int LDA, K;
	LDA = std::max(1,M);
	K = std::min(M,N);

	double * TAU = new double [K];
	int lwork = std::max(M,N)*std::max(M,N);//The dimension of the array WORK.
	double * work = new double [lwork];
	int INFO;
	///////////////////////////////////
	dgeqrf_(&M, &N, B, &LDA, TAU, work, &lwork, &INFO);
	if (INFO!=0) {
		std::cout<<"QR illegal value at "<<INFO<<std::endl;
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < std::min(i+1,M); j++) {
			R[i*K+j] = B[i*M+j];
		}
	}
	dorgqr_(&M, &K, &K, B, &LDA, TAU, work, &lwork, &INFO);
	if (INFO!=0) {
		std::cout<<"QR Failed!\n";
	}
	std::copy(B, B+M*K, Q);
	///////////////////////////////////
	//Free Space
	delete [] B;
	delete [] TAU;
	delete [] work;
}


void QR(int r, int c, std::complex<double>* A, std::complex<double>* Q, std::complex<double>* R){
	///////////////////////////////////
	// INITIALIZATION
	int M = r;
	int N = c;
	std::complex<double>* B = new std::complex<double> [M*N];
	std::copy(A, A+M*N, B);

	int LDA, K;
	LDA = std::max(1,M);
	K = std::min(M,N);

	std::complex<double> * TAU = new std::complex<double> [K];
	int lwork = std::max(M,N)*std::max(M,N);//The dimension of the array WORK.
	std::complex<double> * work = new std::complex<double> [lwork];
	int INFO;
	///////////////////////////////////
	zgeqrf_(&M, &N, B, &LDA, TAU, work, &lwork, &INFO);
	if (INFO!=0) {
		std::cout<<"QR illegal value at "<<INFO<<std::endl;
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < std::min(i+1,M); j++) {
			R[i*K+j] = B[i*M+j];
		}
	}
	zungqr_(&M, &K, &K, B, &LDA, TAU, work, &lwork, &INFO);
	if (INFO!=0) {
		std::cout<<"QR Failed!\n";
	}
	std::copy(B, B+M*K, Q);
	///////////////////////////////////
	//Free Space
	delete [] B;
	delete [] TAU;
	delete [] work;
}


void SVD(int r, int c, double* A, double* U, vector<double>& S, double* V)
{
	///////////////////////////////////
	// INITIALIZATION
	char jobU = 'S'; //all M columns of U are returned in array U
	char jobV = 'S'; //all N rows of V^T are returned in the array VT
	int M = r;
	int N = c;
	double* B = new double [M*N];
	std::copy(A, A+M*N, B);

	int LDA = M;//The leading dimension of the array A
	int K = std::min(M, N);
	S.resize(K);

	int LDU = M;//The leading dimension of the array U
	int LDV = K;//The leading dimension of the array VT
	int lwork = std::max(6*K+2*std::max(M,N),10*K);
	double * work = new double [lwork];
	int INFO;
	///////////////////////////////////
	// CALL LAPACK
	dgesvd_(&jobU, &jobV, &M, &N, B, &LDA, S.data(), U, &LDU, V, &LDV, work, &lwork, &INFO);
	///////////////////////////////////
	if (INFO<0) {
		std::cout<<"SVD illegal value at "<<-INFO<<"!\n";
	}else if (INFO>0){
		std::cout << "SVD did not converge!\n";
	}
	///////////////////////////////////
	//Free Space
	delete [] B;
	delete [] work;
}


void SVD(int r, int c, std::complex<double>* A, std::complex<double>* U, vector<double>& S, std::complex<double>* V)
{
	///////////////////////////////////
	// INITIALIZATION
	char jobU = 'S'; //all M columns of U are returned in array U
	char jobV = 'S'; //all N rows of V^T are returned in the array VT

	int M = r;
	int N = c;
	std::complex<double>* B = new std::complex<double> [M*N];
	std::copy(A, A+M*N, B);

	int LDA = M;//The leading dimension of the array A
	int K = std::min(M, N);
	S.resize(K);

	int LDU = M;//The leading dimension of the array U
	int LDV = K;//The leading dimension of the array VT
	int lwork = std::max(6*K+2*std::max(M,N),10*K);
	std::complex<double> * work = new std::complex<double> [lwork];
	double * rwork = new double [5*K];
	int INFO;
	///////////////////////////////////
	// CALL LAPACK
	zgesvd_(&jobU, &jobV, &M, &N, B, &LDA, S.data(), U, &LDU, V, &LDV, work, &lwork, rwork, &INFO);
	///////////////////////////////////
	if (INFO<0) {
		std::cout<<"SVD illegal value at "<<-INFO<<"!\n";
	}else if (INFO>0){
		std::cout << "SVD did not converge!\n";
	}
	///////////////////////////////////
	//Free Space
	delete [] B;
	delete [] work;
	delete [] rwork;
}

void SVD(int r, int c, double* A, double* U, vector<double>& S, double* V, char direction)
{
	///////////////////////////////////
	// INITIALIZATION
	char jobU; //all M columns of U are returned in array U
	char jobV; //all N rows of V^T are returned in the array VT
	if(direction=='L')
	{
		jobU = 'S'; //all M columns of U are returned in array U
		jobV = 'N'; //all N rows of V^T are returned in the array VT
	}else if(direction=='R')
	{
		jobU = 'N'; //all M columns of U are returned in array U
		jobV = 'S'; //all N rows of V^T are returned in the array VT
	}
	int M = r;
	int N = c;
	double* B = new double [M*N];
	std::copy(A, A+M*N, B);

	int LDA = M;//The leading dimension of the array A
	int K = std::min(M, N);
	S.resize(K);

	int LDU = M;//The leading dimension of the array U
	int LDV = K;//The leading dimension of the array VT
	int lwork = std::max(6*K+2*std::max(M,N),10*K);
	double * work = new double [lwork];
	int INFO;
	///////////////////////////////////
	// CALL LAPACK
	dgesvd_(&jobU, &jobV, &M, &N, B, &LDA, S.data(), U, &LDU, V, &LDV, work, &lwork, &INFO);
	///////////////////////////////////
	if (INFO<0) {
		std::cout<<"SVD illegal value at "<<-INFO<<"!\n";
	}else if (INFO>0){
		std::cout << "SVD did not converge!\n";
	}
	///////////////////////////////////
	//Free Space
	delete [] B;
	delete [] work;
}


void SVD(int r, int c, std::complex<double>* A, std::complex<double>* U, vector<double>& S, std::complex<double>* V, char direction)
{
	///////////////////////////////////
	// INITIALIZATION
	char jobU; //all M columns of U are returned in array U
	char jobV; //all N rows of V^T are returned in the array VT
	if(direction=='L')
	{
		jobU = 'S'; //all M columns of U are returned in array U
		jobV = 'N'; //all N rows of V^T are returned in the array VT
	}else if(direction=='R')
	{
		jobU = 'N'; //all M columns of U are returned in array U
		jobV = 'S'; //all N rows of V^T are returned in the array VT
	}
	int M = r;
	int N = c;
	std::complex<double>* B = new std::complex<double> [M*N];
	std::copy(A, A+M*N, B);

	int LDA = M;//The leading dimension of the array A
	int K = std::min(M, N);
	S.resize(K);

	int LDU = M;//The leading dimension of the array U
	int LDV = K;//The leading dimension of the array VT
	int lwork = std::max(6*K+2*std::max(M,N),10*K);
	std::complex<double> * work = new std::complex<double> [lwork];
	double * rwork = new double [5*K];
	int INFO;
	///////////////////////////////////
	// CALL LAPACK
	zgesvd_(&jobU, &jobV, &M, &N, B, &LDA, S.data(), U, &LDU, V, &LDV, work, &lwork, rwork, &INFO);
	///////////////////////////////////
	if (INFO<0) {
		std::cout<<"SVD illegal value at "<<-INFO<<"!\n";
	}else if (INFO>0){
		std::cout << "SVD did not converge!\n";
	}
	///////////////////////////////////
	//Free Space
	delete [] B;
	delete [] work;
	delete [] rwork;
}

void SVD(int r, int c, double* A, double* U, vector<double>& S, double* V, char direction, double cutoff)
{
	///////////////////////////////////
	// INITIALIZATION
	char jobU; //all M columns of U are returned in array U
	char jobV; //all N rows of V^T are returned in the array VT
	if(direction=='L')
	{
		jobU = 'S'; //all M columns of U are returned in array U
		jobV = 'N'; //all N rows of V^T are returned in the array VT
	}else if(direction=='R')
	{
		jobU = 'N'; //all M columns of U are returned in array U
		jobV = 'S'; //all N rows of V^T are returned in the array VT
	}
	int M = r;
	int N = c;
	double* B = new double [M*N];
	std::copy(A, A+M*N, B);

	int LDA = M;//The leading dimension of the array A
	int K = std::min(M, N);
	S.resize(K);

	int LDU = M;//The leading dimension of the array U
	int LDV = K;//The leading dimension of the array VT
	int lwork = std::max(6*K+2*std::max(M,N),10*K);
	double * work = new double [lwork];
	int INFO;
	///////////////////////////////////
	// CALL LAPACK
	dgesvd_(&jobU, &jobV, &M, &N, B, &LDA, S.data(), U, &LDU, V, &LDV, work, &lwork, &INFO);
	///////////////////////////////////
	if (INFO<0) {
		std::cout<<"SVD illegal value at "<<-INFO<<"!\n";
	}else if (INFO>0){
		std::cout << "SVD did not converge!\n";
	}
	///////////////////////////////////
	double norm = 0;
	for (int i = 0; i < K; i++) {
		norm += S[i]*S[i];
	}
	double running_weight = 0;
	int Kp = 0;
	for (int i = 0; i < K; i++) {
		running_weight += S[i]*S[i]/norm;
		Kp += 1;
		if(1-running_weight <= cutoff){
			break;
		}
	}
	S.resize(Kp);
	// running_weight = std::sqrt(running_weight);
	// for (int i = 0; i < Kp; i++) {
	// 	S[i] /= running_weight;
	// }
	///////////////////////////////////
	//Free Space
	delete [] B;
	delete [] work;
}


void SVD(int r, int c, std::complex<double>* A, std::complex<double>* U, vector<double>& S, std::complex<double>* V, char direction, double cutoff)
{
	///////////////////////////////////
	// INITIALIZATION
	char jobU; //all M columns of U are returned in array U
	char jobV; //all N rows of V^T are returned in the array VT
	if(direction=='L')
	{
		jobU = 'S'; //all M columns of U are returned in array U
		jobV = 'N'; //all N rows of V^T are returned in the array VT
	}else if(direction=='R')
	{
		jobU = 'N'; //all M columns of U are returned in array U
		jobV = 'S'; //all N rows of V^T are returned in the array VT
	}
	int M = r;
	int N = c;
	std::complex<double>* B = new std::complex<double> [M*N];
	std::copy(A, A+M*N, B);

	int LDA = M;//The leading dimension of the array A
	int K = std::min(M, N);
	S.resize(K);

	int LDU = M;//The leading dimension of the array U
	int LDV = K;//The leading dimension of the array VT
	int lwork = std::max(6*K+2*std::max(M,N),10*K);
	std::complex<double> * work = new std::complex<double> [lwork];
	double * rwork = new double [5*K];
	int INFO;
	///////////////////////////////////
	// CALL LAPACK
	zgesvd_(&jobU, &jobV, &M, &N, B, &LDA, S.data(), U, &LDU, V, &LDV, work, &lwork, rwork, &INFO);
	///////////////////////////////////
	if (INFO<0) {
		std::cout<<"SVD illegal value at "<<-INFO<<"!\n";
	}else if (INFO>0){
		std::cout << "SVD did not converge!\n";
	}
	///////////////////////////////////
	double norm = 0;
	for (int i = 0; i < K; i++) {
		norm += S[i]*S[i];
	}
	double running_weight = 0;
	int Kp = 0;
	for (int i = 0; i < K; i++) {
		running_weight += S[i]*S[i]/norm;
		Kp += 1;
		if(1-running_weight <= cutoff){
			break;
		}
	}
	S.resize(Kp);
	// running_weight = std::sqrt(running_weight);
	// for (int i = 0; i < Kp; i++) {
	// 	S[i] /= running_weight;
	// }
	///////////////////////////////////
	//Free Space
	delete [] B;
	delete [] work;
	delete [] rwork;
}

void SVD(int r, int c, double* A, double* U, vector<double>& S, double* V, char direction, int max_size)
{
	///////////////////////////////////
	// INITIALIZATION
	char jobU; //all M columns of U are returned in array U
	char jobV; //all N rows of V^T are returned in the array VT
	if(direction=='L')
	{
		jobU = 'S'; //all M columns of U are returned in array U
		jobV = 'N'; //all N rows of V^T are returned in the array VT
	}else if(direction=='R')
	{
		jobU = 'N'; //all M columns of U are returned in array U
		jobV = 'S'; //all N rows of V^T are returned in the array VT
	}
	int M = r;
	int N = c;
	double* B = new double [M*N];
	std::copy(A, A+M*N, B);

	int LDA = M;//The leading dimension of the array A
	int K = std::min(M, N);
	S.resize(K);

	int LDU = M;//The leading dimension of the array U
	int LDV = K;//The leading dimension of the array VT
	int lwork = std::max(6*K+2*std::max(M,N),10*K);
	double * work = new double [lwork];
	int INFO;
	///////////////////////////////////
	// CALL LAPACK
	dgesvd_(&jobU, &jobV, &M, &N, B, &LDA, S.data(), U, &LDU, V, &LDV, work, &lwork, &INFO);
	///////////////////////////////////
	if (INFO<0) {
		std::cout<<"SVD illegal value at "<<-INFO<<"!\n";
	}else if (INFO>0){
		std::cout << "SVD did not converge!\n";
	}
	if(K>max_size) S.resize(max_size);
	///////////////////////////////////
	//Free Space
	delete [] B;
	delete [] work;
}


void SVD(int r, int c, std::complex<double>* A, std::complex<double>* U, vector<double>& S, std::complex<double>* V, char direction, int max_size)
{
	///////////////////////////////////
	// INITIALIZATION
	char jobU; //all M columns of U are returned in array U
	char jobV; //all N rows of V^T are returned in the array VT
	if(direction=='L')
	{
		jobU = 'S'; //all M columns of U are returned in array U
		jobV = 'N'; //all N rows of V^T are returned in the array VT
	}else if(direction=='R')
	{
		jobU = 'N'; //all M columns of U are returned in array U
		jobV = 'S'; //all N rows of V^T are returned in the array VT
	}
	int M = r;
	int N = c;
	std::complex<double>* B = new std::complex<double> [M*N];
	std::copy(A, A+M*N, B);

	int LDA = M;//The leading dimension of the array A
	int K = std::min(M, N);
	S.resize(K);

	int LDU = M;//The leading dimension of the array U
	int LDV = K;//The leading dimension of the array VT
	int lwork = std::max(6*K+2*std::max(M,N),10*K);
	std::complex<double> * work = new std::complex<double> [lwork];
	double * rwork = new double [5*K];
	int INFO;
	///////////////////////////////////
	// CALL LAPACK
	zgesvd_(&jobU, &jobV, &M, &N, B, &LDA, S.data(), U, &LDU, V, &LDV, work, &lwork, rwork, &INFO);
	///////////////////////////////////
	if (INFO<0) {
		std::cout<<"SVD illegal value at "<<-INFO<<"!\n";
	}else if (INFO>0){
		std::cout << "SVD did not converge!\n";
	}
	if(K>max_size) S.resize(max_size);
	///////////////////////////////////
	//Free Space
	delete [] B;
	delete [] work;
	delete [] rwork;
}


#endif
