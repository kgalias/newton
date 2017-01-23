#ifndef MY_DOUBLE_H
#define MY_DOUBLE_H

#include <iostream>
#include <cmath>
#include <limits>
using namespace std;

class my_double {
	private:
		double _d;
	public:
		my_double();
		my_double(const double& d);
		my_double(const my_double& dbl);
		double get_double() const;
		my_double operator+(const my_double& dbl) const;
		my_double operator-(const my_double& dbl) const;
		my_double operator*(const my_double& dbl) const;
		my_double operator/(const my_double& dbl) const;
		friend bool operator<(const my_double dbl1, const my_double dbl2);
		friend bool operator>(const my_double dbl1, const my_double dbl2);
		friend bool operator==(const my_double dbl1, const my_double dbl2);
		friend ostream& operator<<(ostream& s, const my_double& dbl);
		friend istream& operator>>(istream& s, my_double& dbl);
};

namespace std {
	template<> struct numeric_limits<my_double> {
		static constexpr int digits10 = numeric_limits<double>::digits10;
		static my_double epsilon() {
			return my_double(numeric_limits<double>::epsilon());
		}
	};
}

my_double abs(const my_double& dbl);

my_double exp(const my_double& dbl);

my_double sqrt(const my_double& dbl);

#endif
