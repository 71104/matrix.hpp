#include <iostream>
#include <matrix.hpp>

using namespace std;
using namespace math;

int main() {
	mat<3, 3>() * vec<3>();
	mat<5, 5>().det();
	vec<3>().dot(vec<3>());
	vec<3>().vector(vec<3>());
	vec<3>{ 1, 2, 3 };
	mat<3>() / mat<3>();
	return 0;
}
