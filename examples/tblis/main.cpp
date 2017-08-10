#include <iostream>
#include <complex>
#include <Eigen/Dense>

#include "tblis.h"

int main(int argc, char const *argv[]) {
	int s1=2, s2=3, s3=4;
	Eigen::MatrixXcd mA(s1,s2);
	mA.setRandom();
	Eigen::MatrixXcd mB(s2,s3);
	mB.setRandom();
	Eigen::MatrixXcd mC(s1,s3);
	mC.setRandom();
	std::cout << "/-------------------------------------------------/" << '\n';
	std::cout << "Matrix mA:" << '\n';
	std::cout << mA << '\n' << '\n';
	std::cout << "Matrix mB:" << '\n';
	std::cout << mB << '\n' << '\n';
	std::cout << "Matrix mA*mB+mC:" << '\n';
	std::cout << mA*mB+mC << '\n' << '\n';
	std::cout << "/-------------------------------------------------/" << '\n';

	tblis::tensor_view< std::complex<double> > A({s1,s2}, mA.data());
	tblis::tensor_view< std::complex<double> > B({s2,s3}, mB.data());
	tblis::tensor_view< std::complex<double> > C({s1,s3}, mC.data());

	tblis::mult(std::complex<double>(1.0),A,"ab",B,"bc",std::complex<double>(1.0),C,"ac");

	std::cout << "tblis results:" << '\n';
	std::cout << mC << '\n' << '\n';

	std::cout << "/-------------------------------------------------/" << '\n';
  return 0;
}
