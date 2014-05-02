#pragma once

#include <algorithm>
#include <cassert>


namespace math {
	template<unsigned int _n, unsigned int _m, typename _Type>
	struct Minor {};


	template<unsigned int _n, unsigned int _m, typename _Type>
	struct Determinant {};


	template<unsigned int _n, unsigned int _m = _n, typename _Type = double>
	struct mat {
		_Type m_a[_m][_n];

		mat() {}

		mat(_Type const (&a_a)[_m][_n])
			:
		m_a(a_a) {}

		mat(_Type (&&a_a)[_m][_n])
			:
		m_a(move(a_a)) {}

		struct row {
			_Type (&m_ra)[_m][_n];
			unsigned int m_i;

			inline row(_Type (&a_ra)[_m][_n], unsigned int a_i)
				:
			m_ra(a_ra),
				m_i(a_i) {}

			inline _Type &operator [] (unsigned int j) {
				assert(j < _m);
				return m_ra[j][m_i];
			}
		};

		struct const_row {
			_Type const (&m_ra)[_m][_n];
			unsigned int m_i;

			inline const_row(_Type const (&a_ra)[_m][_n], unsigned int a_i)
				:
			m_ra(a_ra),
				m_i(a_i) {}

			inline _Type const &operator [] (unsigned int j) const {
				assert(j < _m);
				return m_ra[j][m_i];
			}
		};

		inline row operator [] (unsigned int i) {
			assert(i < _n);
			return row(m_a, i);
		}

		inline const_row operator [] (unsigned int i) const {
			assert(i < _n);
			return const_row(m_a, i);
		}

		mat<_m, _n, _Type> transpose() const {
			mat<_m, _n, _Type> Result;
			for(unsigned int i = 0; i < _n; ++i) {
				for(unsigned int j = 0; j < _m; ++j) {
					Result[j][i] = m_a[j][i];
				}
			}
			return Result;
		}

		inline _Type minor(unsigned int const i, unsigned int const j) const {
			return Minor<_n, _n, _Type>::Compute(*this, i, j);
		}

		inline _Type det() const {
			return Determinant<_n, _m, _Type>::Compute(*this);
		}

		mat<_n, _m, _Type> operator + (mat<_n, _m, _Type> const &r) const {
			mat<_n, _m, _Type> Result;
			for (unsigned int i = 0; i < _n; ++i) {
				for (unsigned int j = 0; j < _m; ++j) {
					Result[i][j] += m_a[j][i] + r[i][j];
				}
			}
			return Result;
		}

		mat<_n, _m, _Type> &operator += (mat<_n, _m, _Type> const &r) {
			for (unsigned int i = 0; i < _n; ++i) {
				for (unsigned int j = 0; j < _m; ++j) {
					m_a[j][i] += r[i][j];
				}
			}
			return *this;
		}

		mat<_n, _m, _Type> operator - (mat<_n, _m, _Type> const &r) const {
			mat<_n, _m, _Type> Result;
			for (unsigned int i = 0; i < _n; ++i) {
				for (unsigned int j = 0; j < _m; ++j) {
					Result[i][j] += m_a[j][i] + r[i][j];
				}
			}
			return Result;
		}

		mat<_n, _m, _Type> &operator -= (mat<_n, _m, _Type> const &r) {
			for (unsigned int i = 0; i < _n; ++i) {
				for (unsigned int j = 0; j < _m; ++j) {
					m_a[j][i] -= r[i][j];
				}
			}
			return *this;
		}

		mat<_n, _m, _Type> operator * (_Type const &r) const {
			mat<_n, _m, _Type> Result;
			for (unsigned int i = 0; i < _n; ++i) {
				for (unsigned int j = 0; j < _m; ++j) {
					Result[i][j] = m_a[j][i] * r;
				}
			}
			return Result;
		}

		mat<_n, _m, _Type> &operator *= (_Type const &r) {
			for (unsigned int i = 0; i < _n; ++i) {
				for (unsigned int j = 0; j < _m; ++j) {
					m_a[j][i] *= r;
				}
			}
			return *this;
		}

		mat<_n, _m, _Type> operator / (_Type const &r) const {
			mat<_n, _m, _Type> Result;
			for (unsigned int i = 0; i < _n; ++i) {
				for (unsigned int j = 0; j < _m; ++j) {
					Result[i][j] = m_a[j][i] / r;
				}
			}
			return Result;
		}

		mat<_n, _m, _Type> &operator /= (_Type const &r) {
			for (unsigned int i = 0; i < _n; ++i) {
				for (unsigned int j = 0; j < _m; ++j) {
					m_a[j][i] /= r;
				}
			}
			return *this;
		}

		template<unsigned int _l>
		mat<_n, _l, _Type> operator * (mat<_m, _l, _Type> const &r) const {
			mat<_n, _l, _Type> Result;
			for (unsigned int i = 0; i < _n; ++i) {
				for (unsigned int j = 0; j < _l; ++j) {
					for (unsigned int k = 0; k < _m; ++k) {
						Result[i][j] += m_a[k][i] * r[k][j];
					}
				}
			}
			return Result;
		}
	};


	template<typename _Type>
	struct Minor<1, 1, _Type> {
		static inline _Type Compute(mat<1, 1, _Type> const &r, unsigned int const i, unsigned int const j) {
			assert((i < 1) && (j < 1));
			return 1;
		};
	};


	template<typename _Type>
	struct Minor<2, 2, _Type> {
		static inline _Type Compute(mat<2, 2, _Type> const &r, unsigned int const i, unsigned int const j) {
			assert((i < 2) && (j < 2));
			return r.m_a[1 - j][1 - i];
		};
	};


	template<unsigned int _n, typename _Type>
	struct Minor<_n, _n, _Type> {
		static _Type Compute(mat<_n, _n, _Type> const &r, unsigned int const i, unsigned int const j) {
			mat<_n - 1, _n - 1, _Type> m;
			for (unsigned int i2 = 0; i2 < i; ++i2) {
				for (unsigned int j2 = 0; j2 < j; ++j2) {
					m.m_a[j2][i2] = r.m_a[j2][i2];
				}
				for (unsigned int j2 = j + 1; j2 < _n; ++j2) {
					m.m_a[j2 - 1][i2] = r.m_a[j2][i2];
				}
			}
			for (unsigned int i2 = i + 1; i2 < _n; ++i2) {
				for (unsigned int j2 = 0; j2 < j; ++j2) {
					m.m_a[j2][i2 - 1] = r.m_a[j2][i2];
				}
				for (unsigned int j2 = j + 1; j2 < _n; ++j2) {
					m.m_a[j2 - 1][i2 - 1] = r.m_a[j2][i2];
				}
			}
			return m.det();
		}
	};


	template<typename _Type>
	struct Determinant<1, 1, _Type> {
		static inline _Type Compute(mat<1, 1, _Type> const &r) {
			return r.m_a[0][0];
		}
	};


	template<typename _Type>
	struct Determinant<2, 2, _Type> {
		static inline _Type Compute(mat<2, 2, _Type> const &r) {
			return r.m_a[0][0] * r.m_a[1][1] - r.m_a[0][1] * r.m_a[1][0];
		}
	};


	template<typename _Type>
	struct Determinant<3, 3, _Type> {
		static inline _Type Compute(mat<3, 3, _Type> const &r) {
			return
				r.m_a[0][0] * (r.m_a[1][1] * r.m_a[2][2] - r.m_a[1][2] * r.m_a[2][1]) -
				r.m_a[1][0] * (r.m_a[0][1] * r.m_a[2][2] - r.m_a[0][2] * r.m_a[2][1]) +
				r.m_a[2][0] * (r.m_a[0][1] * r.m_a[1][2] - r.m_a[0][2] * r.m_a[1][1]);
		}
	};


	template<typename _Type>
	struct Determinant<4, 4, _Type> {
		static inline _Type Compute(mat<4, 4, _Type> const &r) {
			return
				r.m_a[0][0] * (r.m_a[1][1] * (r.m_a[2][2] * r.m_a[3][3] - r.m_a[2][3] * r.m_a[3][2]) -
				               r.m_a[1][2] * (r.m_a[2][1] * r.m_a[3][3] - r.m_a[2][3] * r.m_a[3][1]) +
				               r.m_a[1][3] * (r.m_a[2][1] * r.m_a[3][2] - r.m_a[2][2] * r.m_a[3][1])) -
				r.m_a[0][1] * (r.m_a[1][0] * (r.m_a[2][2] * r.m_a[3][3] - r.m_a[2][3] * r.m_a[3][2]) -
				               r.m_a[1][2] * (r.m_a[2][0] * r.m_a[3][3] - r.m_a[2][3] * r.m_a[3][0]) +
				               r.m_a[1][3] * (r.m_a[2][0] * r.m_a[3][2] - r.m_a[2][2] * r.m_a[3][0])) +
				r.m_a[0][2] * (r.m_a[1][0] * (r.m_a[2][1] * r.m_a[3][3] - r.m_a[2][3] * r.m_a[3][1]) -
				               r.m_a[1][1] * (r.m_a[2][0] * r.m_a[3][3] - r.m_a[2][3] * r.m_a[3][0]) +
				               r.m_a[1][3] * (r.m_a[2][0] * r.m_a[3][1] - r.m_a[2][1] * r.m_a[3][0])) -
				r.m_a[0][3] * (r.m_a[1][0] * (r.m_a[2][1] * r.m_a[3][2] - r.m_a[2][2] * r.m_a[3][1]) -
				               r.m_a[1][1] * (r.m_a[2][0] * r.m_a[3][2] - r.m_a[2][2] * r.m_a[3][0]) +
				               r.m_a[1][2] * (r.m_a[2][0] * r.m_a[3][1] - r.m_a[2][1] * r.m_a[3][0]));
		}
	};


	template<unsigned int _n, typename _Type>
	struct Determinant<_n, _n, _Type> {
		static _Type Compute(mat<_n, _n, _Type> const &r) {
			_Type d(0);
			int n = 1;
			for (unsigned int j = 0; j < _n; ++j) {
				d += n * r.m_a[j][0] * r.minor(0, j);
				n = -n;
			}
			return d;
		}
	};


	typedef mat<2, 2> mat22, mat2;
	typedef mat<2, 3> mat23;
	typedef mat<2, 4> mat24;
	typedef mat<3, 2> mat32;
	typedef mat<3, 3> mat33, mat3;
	typedef mat<3, 4> mat34;
	typedef mat<4, 2> mat42;
	typedef mat<4, 3> mat43;
	typedef mat<4, 4> mat44, mat4;


	template<unsigned int _n, typename _Type = double>
	struct vec : public mat<_n, 1, _Type> {};


	typedef vec<2> vec2;
	typedef vec<3> vec3;
	typedef vec<4> vec4;
}
