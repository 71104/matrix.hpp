#include <iostream>
#include <matrix.hpp>

using namespace std;
using namespace math;

int main() {
	mat<3, 3>() * vec<3>();
	mat<5, 5>().det();
	return 0;
}
