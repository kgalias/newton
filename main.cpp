#include "newton.h"
#include "my_double.h"
#include <iostream>
#include <iomanip>
using namespace std;

template<typename scalar_t>
struct fun {
	virtual vector<scalar_t> operator()(const vector<scalar_t>& x) = 0;
	virtual ~fun() {}
};

template<typename scalar_t>
struct jac {
	virtual matrix<scalar_t> operator()(const vector<scalar_t>& x) = 0;
	virtual ~jac() {}
};

template<typename scalar_t>
struct func4 : public fun<scalar_t>{
	vector<scalar_t> operator()(const vector<scalar_t>& x) {
		vector<scalar_t> result(2);
		scalar_t one = 1;
		scalar_t five = 5;
		result[0] = x[0] * x[0] + x[1] * x[1] - one;
		result[1] = x[0] + five * x[1];
		return result;
	}
};

template<typename scalar_t>
struct jac4 : public jac<scalar_t>{
	matrix<scalar_t> operator()(const vector<scalar_t>& x) {
		matrix<scalar_t> result(2);
		scalar_t one = 1;
		scalar_t two = 2;
		scalar_t five = 5;
		result[0].resize(2);
		result[1].resize(2);
		result[0][0] = two * x[0];
		result[0][1] = two * x[1];
		result[1][0] = one;
		result[1][1] = five;
		return result;
	}
};

template<typename scalar_t>
ostream& operator<<(ostream& s, const vector<scalar_t>& v) {
	for(auto elem : v) {
		s << elem << " ";
	}
    return s; 
}

template<typename scalar_t>
function<scalar_t(scalar_t)> make_poly(const int& n, const vector<scalar_t>& coeff) {
	auto poly = [n, coeff](scalar_t x) -> scalar_t {
		scalar_t result = coeff[0];
		for(int i=1; i<=n; ++i) {
			result = result * x + coeff[i];
		}
		return result;
	};
	return poly;
}

template<typename scalar_t>
function<scalar_t(scalar_t)> make_poly_deriv(const int& n, const vector<scalar_t>& coeff) {
	auto poly = [n, coeff](scalar_t x) -> scalar_t {
		scalar_t result = coeff[0] * n;
		for(int i=1; i<n; ++i) {
			result = result * x + coeff[i] * (n-i);
		}
		return result;
	};
	return poly;
}

template<typename scalar_t>
function<scalar_t(scalar_t)> func3(const scalar_t& a, const scalar_t& b) {
	auto func = [a, b](scalar_t x) -> scalar_t { return a * exp(x) + b * x; };
	return func;
}

template<typename scalar_t>
function<scalar_t(scalar_t)> deriv3(const scalar_t& a, const scalar_t& b) {
	auto func = [a, b](scalar_t x) -> scalar_t { return a * exp(x) + b; };
	return func;
}

template<typename scalar_t>
scalar_t newton_poly_numeric(const int& max_iter_num) {
	int n = -1;
	cout << "Input degree of polynomial: ";
	cin >> n;
	vector<scalar_t> coeff;
	coeff.resize(n+1);
	cout << "Input " << n + 1 << " coefficients of polynomial [a_n .. a_0]: ";
	for(int i=0; i<=n; ++i) cin >> coeff[i];
	scalar_t pt;
	cout << "Input starting point: ";
	cin >> pt;
	cout << endl;									
	function<scalar_t(scalar_t)> f = make_poly(n, coeff);
	if(max_iter_num == -1) return newton(f, pt);
	else return newton(f, pt, max_iter_num);								
}

template<typename scalar_t>
scalar_t newton_poly_algebraic(const int& max_iter_num) {
	int n = -1;
	cout << "Input degree of polynomial: ";
	cin >> n;
	vector<scalar_t> coeff;
	coeff.resize(n+1);
	cout << "Input " << n + 1 << " coefficients of polynomial [a_n .. a_0]: ";
	for(int i=0; i<=n; ++i) cin >> coeff[i];
	scalar_t pt;
	cout << "Input starting point: ";
	cin >> pt;
	cout << endl;									
	function<scalar_t(scalar_t)> f = make_poly(n, coeff);
	function<scalar_t(scalar_t)> df = make_poly_deriv(n, coeff);								
	if(max_iter_num == -1) return newton(f, df, pt);
	else return newton(f, df, pt, max_iter_num);
}

template<typename scalar_t>
scalar_t func3_algebraic(const int& max_iter_num) {
	scalar_t a, b;
	cout << "Input parameters a, b: ";
	cin >> a >> b;
	scalar_t pt;
	cout << "Input starting point: ";
	cin >> pt;
	cout << endl;
	function<scalar_t(scalar_t)> f = func3(a, b);
	function<scalar_t(scalar_t)> df = deriv3(a, b);
	if(max_iter_num == -1) return newton(f, df, pt);
	else return newton(f, df, pt, max_iter_num);
}

template<typename scalar_t>
scalar_t func3_numerical(const int& max_iter_num) {
	scalar_t a, b;
	cout << "Input parameters a, b: ";
	cin >> a >> b;
	scalar_t pt;
	cout << "Input starting point: ";
	cin >> pt;
	cout << endl;
	function<scalar_t(scalar_t)> f = func3(a, b);
	if(max_iter_num == -1) return newton(f, pt);
	else return newton(f, pt, max_iter_num);	
}

template<typename scalar_t>
vector<scalar_t> func4_algebraic(const int& max_iter_num) {
	func4<scalar_t> f4;
	jac4<scalar_t> df4;
	vector<scalar_t> pt;
	pt.resize(2);
	cout << "Input starting point (x0, y0): ";
	cin >> pt[0] >> pt[1];
	cout << endl;
	function<vector<scalar_t>(vector<scalar_t>)> f = f4;
	function<matrix<scalar_t>(vector<scalar_t>)> df = df4;
	if(max_iter_num == -1) return newton(f, df, pt);
	else  return newton(f, df, pt, max_iter_num);
}

template<typename scalar_t>
vector<scalar_t> func5_numerical(const int& max_iter_num) {
	scalar_t a, b;
	cout << "Input parameters a, b: ";
	cin >> a >> b;
	vector<scalar_t> pt;
	pt.resize(2);
	cout << "Input starting point (x0, y0): ";
	cin >> pt[0] >> pt[1];
	cout << endl;
	function<vector<scalar_t> (vector<scalar_t>)> f = [a, b](vector<scalar_t> x) -> vector<scalar_t> {
		vector<scalar_t> result;
		result.resize(2);
		result[0] = x[0] * x[0] * x[0] + x[0] * x[1] + a;
		result[1] = x[0] * x[0] - x[1] * x[1] + b;
		return result;
	};
	if(max_iter_num == -1) return newton(f, pt);
	else  return newton(f, pt, max_iter_num);
}

int main(int argc, char** argv) {
	
	int action1 = -1;
	int action2 = -1;
	int action3 = -1;
	int max_iter_num = -1;
	enum rep_type {dbl, flt, long_dbl, my_dbl};
	rep_type current_type = dbl;
	
	while(true) {
		cout << "1. Choice of function" << endl;
		cout << "2. Options" << endl;
		cout << "3. End" << endl;
		cout << endl;
		cin >> action1;
		cout << endl;
		switch(action1) {
			case 1:
				while(true) { 
					cout << "1. Polynomial f(x) = a_n * x^n + ... + a_0 with numeric derivative." << endl;
					cout << "2. Polynomial f(x) = a_n * x^n + ... + a_0 with algebraic derivative." << endl;
					cout << "3. Function f(x) = a * e^x + b * x with algebraic derivative." << endl;
					cout << "4. Function f(x) = a * e^x + b * x with numerical derivative." << endl;
					cout << "5. Function f(x, y) = (x^2 + y^2 - 1, x + 5 * y) with algebraic derivative." << endl;
					cout << "6. Function f(x, y) = (x^3 + x * y + a, x^2 - y^2 + b) with numeric derivative."  << endl;
					cout << "7. Back" << endl;
					cout << endl;
					cin >> action2;
					cout << endl;
					switch(action2) {
						case 1:
							try{ 
								switch(current_type) {
									case dbl:
										cout << std::setprecision(numeric_limits<double>::digits10) << "Using representation type double the zero is: " << newton_poly_numeric<double>(max_iter_num) << endl;
										break;
									case flt:
										cout << std::setprecision(numeric_limits<float>::digits10) << "Using representation type float the zero is: " << newton_poly_numeric<float>(max_iter_num) << endl;
										break;
									case long_dbl:
										cout << std::setprecision(numeric_limits<long double>::digits10) << "Using representation type long double the zero is: " << newton_poly_numeric<long double>(max_iter_num) << endl;
										break;
									case my_dbl:
										cout << std::setprecision(numeric_limits<my_double>::digits10) << "Using representation type my_double the zero is: " << newton_poly_numeric<my_double>(max_iter_num) << endl;
										break;
								}
								cout << endl;
							} catch(const domain_error& e){
								cout << e.what() << endl;
								cout << endl;
								continue;
							}
							break;
						case 2:
							try{ 
								switch(current_type) {
									case dbl:
										cout << std::setprecision(numeric_limits<double>::digits10) <<  "Using representation type double the zero is: " << newton_poly_algebraic<double>(max_iter_num) << endl;
										break;
									case flt:
										cout << std::setprecision(numeric_limits<float>::digits10) <<  "Using representation type float the zero is: " << newton_poly_algebraic<float>(max_iter_num) << endl;
										break;
									case long_dbl:
										cout << std::setprecision(numeric_limits<long double>::digits10) <<  "Using representation type long double the zero is: " << newton_poly_algebraic<long double>(max_iter_num) << endl;
										break;
									case my_dbl:
										cout << std::setprecision(numeric_limits<my_double>::digits10) <<  "Using representation type my_double the zero is: " << newton_poly_algebraic<my_double>(max_iter_num) << endl;
										break;
								}
								cout << endl;
							} catch(const domain_error& e){
								cout << e.what() << endl;
								cout << endl;
								continue;
							}
							break;
						case 3:
							try{ 
								switch(current_type) {
									case dbl:
										cout << std::setprecision(numeric_limits<double>::digits10) <<  "Using representation type double the zero is: " << func3_algebraic<double>(max_iter_num) << endl;
										break;
									case flt:
										cout << std::setprecision(numeric_limits<float>::digits10) <<  "Using representation type float the zero is: " << func3_algebraic<float>(max_iter_num) << endl;
										break;
									case long_dbl:
										cout << std::setprecision(numeric_limits<long double>::digits10) <<  "Using representation type long double the zero is: " << func3_algebraic<long double>(max_iter_num) << endl;
										break;
									case my_dbl:
										cout << std::setprecision(numeric_limits<my_double>::digits10) <<  "Using representation type my_double the zero is: " << func3_algebraic<my_double>(max_iter_num) << endl;
										break;
								}
								cout << endl;
							} catch(const domain_error& e){
								cout << e.what() << endl;
								cout << endl;
								continue;
							}
							break;
						case 4:
							try{ 
								switch(current_type) {
									case dbl:
										cout << std::setprecision(numeric_limits<double>::digits10) <<  "Using representation type double the zero is: " << func3_numerical<double>(max_iter_num) << endl;
										break;
									case flt:
										cout << std::setprecision(numeric_limits<float>::digits10) <<  "Using representation type float the zero is: " << func3_numerical<float>(max_iter_num) << endl;
										break;
									case long_dbl:
										cout << std::setprecision(numeric_limits<long double>::digits10) <<  "Using representation type long double the zero is: " << func3_numerical<long double>(max_iter_num) << endl;
										break;
									case my_dbl:
										cout << std::setprecision(numeric_limits<my_double>::digits10) <<  "Using representation type my_double the zero is: " << func3_numerical<my_double>(max_iter_num) << endl;
										break;
								}
								cout << endl;
							} catch(const domain_error& e){
								cout << e.what() << endl;
								cout << endl;
								continue;
							}
							break;
						case 5:
							try{ 
								switch(current_type) {
									case dbl:
										cout << std::setprecision(numeric_limits<double>::digits10) <<  "Using representation type double the zero is: " << func4_algebraic<double>(max_iter_num) << endl;
										break;
									case flt:
										cout << std::setprecision(numeric_limits<float>::digits10) <<  "Using representation type float the zero is: " << func4_algebraic<float>(max_iter_num) << endl;
										break;
									case long_dbl:
										cout << std::setprecision(numeric_limits<long double>::digits10) <<  "Using representation type long double the zero is: " << func4_algebraic<long double>(max_iter_num) << endl;
										break;
									case my_dbl:
										cout << std::setprecision(numeric_limits<my_double>::digits10) <<  "Using representation type my_double the zero is: " << func4_algebraic<my_double>(max_iter_num) << endl;
										break;
								}
								cout << endl;
							} catch(const domain_error& e){
								cout << e.what() << endl;
								cout << endl;
								continue;
							}
							break;
						case 6:
							try{ 
								switch(current_type) {
									case dbl:
										cout << std::setprecision(numeric_limits<double>::digits10) <<  "Using representation type double the zero is: " << func5_numerical<double>(max_iter_num) << endl;
										break;
									case flt:
										cout << std::setprecision(numeric_limits<float>::digits10) <<  "Using representation type float the zero is: " << func5_numerical<float>(max_iter_num) << endl;
										break;
									case long_dbl:
										cout << std::setprecision(numeric_limits<long double>::digits10) <<  "Using representation type long double the zero is: " << func5_numerical<long double>(max_iter_num) << endl;
										break;
									case my_dbl:
										cout << std::setprecision(numeric_limits<my_double>::digits10) <<  "Using representation type my_double the zero is: " << func5_numerical<my_double>(max_iter_num) << endl;
										break;
								}
								cout << endl;
							} catch(const domain_error& e){
								cout << e.what() << endl;
								cout << endl;
								continue;
							}
							break;	
						case 7: 
							break; 
						default:
							cout << "Improper command, try again." << endl;
							cout << endl;
							continue;
					}
					break;
				}
				break;
			case 2:
				while(true) {
					cout << "1. Change maximum number of iterations." << endl;
					cout << "2. Change representation type." << endl;
					cout << "3. Back." << endl;
					cout << endl;
					cin >> action2;
					cout << endl;
					switch(action2) {
						case 1:
							cout << "Input maximum numer of iterations: ";
							cin >> max_iter_num;
							cout << endl;
							continue;
						case 2:
							while(true) {
								cout << "1. Use double" << endl;
								cout << "2. Use float." << endl;
								cout << "3. Use long double." << endl;
								cout << "4. Use my_double." << endl;
								cin >> action3;
								cout << endl;
								switch(action3) {
									case 1:
										current_type = dbl;
										break;
									case 2:
										current_type = flt;
										break;
									case 3:
										current_type = long_dbl;
										break;
									case 4:
										current_type = my_dbl;
										break;		
									default:
										cout << "Improper command, try again." << endl;
										cout << endl;
										continue;
								}
								break;
							}
							break; 
						case 3: break; 
						default:
							cout << "Improper command, try again." << endl;
							cout << endl;
							continue; 
					}
					break;
				}
				continue;
			case 3: return 0;
			default:
				cout << "Improper command, try again." << endl;
				cout << endl;	
				continue;	
		}
	}

	///////////////////////////////////////////////////////////
	//double(*)(double)
	//vector<double> (*ptr1)(vector<double>) = &F;
	/*shared_ptr< vector<double>(*)(vector<double>) > ptr_to_func = &F;
	shared_ptr< matrix<double>*(vector<double>) > ptr_to_jac = &J;
	shared_ptr< vector<double> > ptr_to_vec = &pt;
	
	vector<double> answw = newton(ptr_to_func.get(), ptr_to_jac.get(), ptr_to_vec.get());
	cout << "newton passing smart pointers:" << endl;
	for(auto a : answw) cout << a << " ";
	cout << endl;*/
	
	return 0;
}
