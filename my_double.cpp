#include "my_double.h"

my_double::my_double(): _d(0.0) {}

my_double::my_double(const double& d) : _d(d) {}

my_double::my_double(const my_double& dbl) : _d(dbl.get_double()) {}

double my_double::get_double() const { return this->_d; }

my_double my_double::operator+(const my_double& dbl) const {
	return my_double(this->get_double() + dbl.get_double());
}

my_double my_double::operator-(const my_double& dbl) const {
	return my_double(this->get_double() - dbl.get_double());
}

my_double my_double::operator*(const my_double& dbl) const {
	return my_double(this->get_double() * dbl.get_double());
}

my_double my_double::operator/(const my_double& dbl) const {
	return my_double(this->get_double() / dbl.get_double());
}

bool operator<(const my_double dbl1, const my_double dbl2) {
	return dbl1.get_double() < dbl2.get_double();
}

bool operator>(const my_double dbl1, const my_double dbl2) {
	return dbl1.get_double() > dbl2.get_double();
}

bool operator==(const my_double dbl1, const my_double dbl2) {
	return dbl1.get_double() == dbl2.get_double();
}

ostream& operator<<(ostream& s, const my_double& dbl) {
	s << dbl.get_double();
    return s; 
}

istream& operator>>(istream& s, my_double& dbl) {
	s >> dbl._d;
	return s;
}

my_double abs(const my_double& dbl) {
	if(dbl < my_double(0.0)) {
		return my_double(-dbl.get_double());
	} else {
		return dbl;
	}
}

my_double exp(const my_double& dbl) {
	double d = dbl.get_double();
	double e = exp(d);
	return my_double(e);
}

my_double sqrt(const my_double& dbl) {
	double d = dbl.get_double();
	double s = sqrt(d);
	return my_double(s);
}
