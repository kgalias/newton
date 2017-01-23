#ifndef NEWTON_H
#define NEWTON_H

#include <iostream>
#include <functional>
#include <vector>
#include <cmath>
#include <cassert>
#include <limits>
#include <stdexcept>

using namespace std;
	
template<typename scalar_t>
using matrix = vector<vector<scalar_t>>;

template<typename scalar_t> 
scalar_t max_norm(const vector<scalar_t>& x) {
	scalar_t m = 0;
	for(auto elem : x) {
		if (abs(elem) > m) m = abs(elem);
	}
	return m;
}

template<typename scalar_t>
vector<scalar_t> subtract(const vector<scalar_t>& x, const vector<scalar_t>& y) {
	assert(x.size() == y.size());
	vector<scalar_t> result;
	for(int i=0; i<x.size(); ++i) result.push_back(x[i] - y[i]);
	return result;
}

template<typename scalar_t>
vector<scalar_t> scale(const vector<scalar_t>& x, const scalar_t& a) {
	vector<scalar_t> result;
	for(auto elem : x) result.push_back(a * elem);
	return result;
}

template<typename scalar_t>
vector<scalar_t> solve_linear_system(matrix<scalar_t> A, vector<scalar_t> b) {
	unsigned int no_of_rows = A.size();
	unsigned int no_of_cols = A[0].size();

	//transform matrix into triangular form
	for(int i=0; i<no_of_cols; ++i) {
		
		//find row with maximum i_th coefficient & swap it with the first considered
		int index_max = i;
		scalar_t value_max = abs(A[i][i]);
		for(int j=i; j<no_of_rows; ++j) {
			if(abs(A[j][i]) > value_max) {
				index_max = j;
				value_max = abs(A[j][i]);
			}
		}
		scalar_t a = 10.0;
		if(value_max < a * numeric_limits<scalar_t>::epsilon()) {
			throw domain_error("Degenerate system of equations.");
		}
		swap(A[i], A[index_max]);
		swap(b[i], b[index_max]);
		
		//zero i-th coefficient in following rows
		for(int j=i+1; j<no_of_rows; ++j) {
			scalar_t coeff = A[j][i] / A[i][i];
			A[j] = subtract(A[j], scale(A[i], coeff));
			b[j] = b[j] - coeff * b[i];
		}
	}
	
	//solve for matrix in triangular form
	vector<scalar_t> result;
	result.resize(no_of_rows);
	for(int j=no_of_rows-1; j>-1; --j) {
		scalar_t partial_sum = b[j];
		for(int i=j+1; i<no_of_cols; i++) {
			partial_sum = partial_sum - (A[j][i] * result[i]);	
		}
		partial_sum = partial_sum / A[j][j];
		result[j] = partial_sum;
	}
	
	return result;
}

template<typename scalar_t>
vector<scalar_t> correction(
							const function<vector<scalar_t> (vector<scalar_t>)>& func, 
						  	const function<matrix<scalar_t> (vector<scalar_t>)>& jac, 
						  	const vector<scalar_t>& xn
						  	) {
	try {
		return solve_linear_system(jac(xn), func(xn));
	} catch(const domain_error& e) {
		throw;
	}
}

//approximate jacobian matrix
template<typename scalar_t>
matrix<scalar_t> approx_jac(
							const function<vector<scalar_t> (vector<scalar_t>)>& func, 
							const vector<scalar_t>& x,
							const scalar_t& epsilon
							) {
	scalar_t t = sqrt(epsilon);
	unsigned int n = x.size();
	matrix<scalar_t> result;
	result.resize(n);
	for(int i=0; i<n; ++i) {
		for(int j=0; j<n; ++j) {	
			vector<scalar_t> x_mod = x;
			x_mod[j] = x_mod[j] + t;
			scalar_t partial_deriv = (func(x_mod)[i] - func(x)[i]) / t;
			result[i].push_back(partial_deriv);
		}
	}
	return result;
}

//calculate zero with jacobian matrix specified
template<typename scalar_t>
vector<scalar_t> newton(
						const function<vector<scalar_t> (vector<scalar_t>)>& func, 
						const function<matrix<scalar_t> (vector<scalar_t>)>& jac, 
						const vector<scalar_t>& x0, 
						const unsigned int& max_iter_num = 100, 
						const scalar_t& epsilon = numeric_limits<scalar_t>::epsilon()
						) {
	unsigned int iter_num = 0;
	vector<scalar_t> xn = x0;
	while(true) {	
		try {
			vector<scalar_t> h = correction(func, jac, xn);
			xn = subtract(xn, h);
			++iter_num;
			if (max_norm(h) < epsilon) return xn;
			if (iter_num > max_iter_num) {
				throw domain_error("Too many iterations.");
			}  
		} catch(const domain_error& e) {
			cout << e.what() << endl;
			throw domain_error("Newton's method does not work.");
		}
	}
}

//calculate zero with no jacobian matrix specified
template<typename scalar_t>
vector<scalar_t> newton(
						const function<vector<scalar_t> (vector<scalar_t>)>& func, 
						const vector<scalar_t>& x0, 
						const unsigned int& max_iter_num = 100, 
						const scalar_t& epsilon = numeric_limits<scalar_t>::epsilon()
						) {
	function<matrix<scalar_t> (vector<scalar_t>)> approx_jac_at_point = [func, epsilon](vector<scalar_t> x) -> matrix<scalar_t> {
		return approx_jac<scalar_t>(func, x, epsilon);
	};						
	return newton(func, approx_jac_at_point, x0, max_iter_num, epsilon);
}

//approximate derivative of 1d function
template<typename scalar_t>
scalar_t approx_deriv(const function<scalar_t (scalar_t)>& func, const scalar_t& x, const scalar_t& epsilon) {
	scalar_t t = sqrt(epsilon);
	return (func(x + t) - func(x)) / t;
}

//calculate zero of 1d function with derivative specified
template<typename scalar_t>
scalar_t newton(
				const function<scalar_t (scalar_t)>& func, 
				const function<scalar_t (scalar_t)>& deriv,
				const scalar_t& x0, 
				const unsigned int& max_iter_num = 100, 
				const scalar_t& epsilon = numeric_limits<scalar_t>::epsilon()
				) {
	auto wrapped_func = [func](vector<scalar_t> x) -> vector<scalar_t> {
		vector<scalar_t> result = {func(x[0])};
		return result; 
	}; 
	auto wrapped_deriv = [deriv](vector<scalar_t> x) -> matrix<scalar_t> {
		matrix<scalar_t> result = {{}};
		result[0].push_back(deriv(x[0]));
		return result;
	};				
	vector<scalar_t> wrapped_x0 = {x0};
	return newton<scalar_t>(wrapped_func, wrapped_deriv, wrapped_x0, max_iter_num, epsilon)[0];
}

//calculate zero of 1d function with no derivative specified
template<typename scalar_t>
scalar_t newton(
				const function<scalar_t (scalar_t)>& func, 
				const scalar_t& x0, 
				const unsigned int& max_iter_num = 100, 
				const scalar_t& epsilon = numeric_limits<scalar_t>::epsilon()
				) {
	function<scalar_t (scalar_t)> approx_deriv_at_point = [func, epsilon](scalar_t x) -> scalar_t {
		return approx_deriv(func, x, epsilon);
	};						
	return newton(func, approx_deriv_at_point, x0, max_iter_num, epsilon);
}

#endif
