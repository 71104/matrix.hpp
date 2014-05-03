#pragma once

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <cmath>


namespace math {
	using namespace std;


	template<unsigned int const _n, unsigned int const _m, typename _Type>
	struct Minor {};


	template<unsigned int const _n, unsigned int const _m, typename _Type>
	struct Determinant {};


	template<unsigned int const _n, unsigned int const _m, typename _Type>
	struct InvertedMatrix {};


	template<unsigned int const _n, unsigned int const _m = _n, typename _Type = double>
	struct mat {
		_Type m_a[_m][_n];

		mat() {}

		mat(_Type const (&a_ra)[_m][_n])
			:
		m_a(a_ra) {}

		mat(_Type (&&a_rra)[_m][_n])
			:
		m_a(move(a_rra)) {}

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

		inline mat<_m, _n, _Type> invert() const {
			return InvertedMatrix<_n, _m, _Type>::Compute(*this);
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

		template<unsigned int const _l>
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

		mat<_n, _n, _Type> operator / (mat<_m, _n, _Type> const &r) const {
			return operator * (r.invert());
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


	template<unsigned int const _n, typename _Type>
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


	template<unsigned int const _n, typename _Type>
	struct Determinant<_n, _n, _Type> {
		static _Type Compute(mat<_n, _n, _Type> const &r) {
			_Type d(0);
			int n = -1;
			for (unsigned int j = 0; j < _n; ++j) {
				d += (n = -n) * r.m_a[j][0] * r.minor(0, j);
			}
			return d;
		}
	};


	template<unsigned int const _n, typename _Type>
	struct InvertedMatrix<_n, _n, _Type> {
		static mat<_n, _n, _Type> Compute(mat<_n, _n, _Type> const &r) {
			mat<_n, _n, _Type> m;
			int n = -1;
			for (unsigned int i = 0; i < _n; ++i) {
				for (unsigned int j = 0; j < _n; ++j) {
					m.m_a[i][j] = (n = -n) * r.minor(i, j);
				}
			}
			return m / r.det();
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


	template<unsigned int const _n, typename _Type>
	struct VectorProduct {};


	template<unsigned int const _n, typename _Type = double>
	struct vec :
		public mat<_n, 1, _Type>
	{
		_Type (&m_a)[_n];

		vec()
			:
		m_a(mat<_n, 1, _Type>::m_a[0]) {}

		vec(_Type const (&a_ra)[_n])
			:
		m_a(mat<_n, 1, _Type>::m_a[0]) {
			for (unsigned int i = 0; i < _n; ++i) {
				m_a[i] = a_ra[i];
			}
		}

		vec(_Type (&&a_rra)[_n])
			:
		m_a(mat<_n, 1, _Type>::m_a[0]) {
			for (unsigned int i = 0; i < _n; ++i) {
				m_a[i] = move(a_rra[i]);
			}
		}

		vec(initializer_list<_Type> const &il)
			:
		m_a(mat<_n, 1, _Type>::m_a[0]) {
			unsigned int i = 0;
			for (auto it = il.begin(); it != il.end(); ++it) {
				m_a[i++] = *it;
			}
		}

		vec(vec<_n, _Type> const &r)
			:
		mat<_n, 1, _Type>(r),
			m_a(mat<_n, 1, _Type>::m_a[0]) {}

		vec(vec<_n, _Type> &&rr)
			:
		mat<_n, 1, _Type>(move(rr)),
			m_a(mat<_n, 1, _Type>::m_a[0]) {}

		vec(mat<_n, 1, _Type> const &r)
			:
		mat<_n, 1, _Type>(r),
			m_a(mat<_n, 1, _Type>::m_a[0]) {}

		vec(mat<_n, 1, _Type> &&rr)
			:
		mat<_n, 1, _Type>(move(rr)),
			m_a(mat<_n, 1, _Type>::m_a[0]) {}

		inline vec<_n, _Type> &operator = (vec<_n, _Type> const &r) {
			return mat<_n, 1, _Type>::operator = (r);
		}

		inline vec<_n, _Type> &operator = (vec<_n, _Type> &&rr) {
			return mat<_n, 1, _Type>::operator = (move(rr));
		}

		inline vec<_n, _Type> &operator = (mat<_n, 1, _Type> const &r) {
			return mat<_n, 1, _Type>::operator = (r);
		}

		inline vec<_n, _Type> &operator = (mat<_n, 1, _Type> &&rr) {
			return mat<_n, 1, _Type>::operator = (move(rr));
		}

		inline auto modulus() const -> decltype(sqrt(m_a[0])) {
			_Type s(0);
			for (unsigned int i = 0; i < _n; ++i) {
				s += m_a[i] * m_a[i];
			}
			return sqrt(s);
		}

		inline auto length() const -> decltype(sqrt(m_a[0])) {
			return modulus();
		}

		_Type scalar(vec<_n, _Type> const &r) const {
			_Type s(0);
			for (unsigned int i = 0; i < _n; ++i) {
				s += m_a[i] * r.m_a[i];
			}
			return s;
		}

		inline _Type inner(vec<_n, _Type> const &r) const {
			return scalar(r);
		}

		inline _Type dot(vec<_n, _Type> const &r) const {
			return scalar(r);
		}

		inline vec<_n, _Type> vector(vec<_n, _Type> const &r) const {
			return VectorProduct<_n, _Type>::Compute(*this, r);
		}

		inline vec<_n, _Type> outer(vec<_n, _Type> const &r) const {
			return vector(r);
		}

		inline vec<_n, _Type> cross(vec<_n, _Type> const &r) const {
			return vector(r);
		}
	};


	template<typename _Type>
	struct VectorProduct<3, _Type> {
		static vec<3, _Type> Compute(vec<3, _Type> const &ru, vec<3, _Type> const &rv) {
			return {
				ru.m_a[1] * rv.m_a[2] - ru.m_a[2] * rv.m_a[1],
				ru.m_a[2] * rv.m_a[0] - ru.m_a[0] * rv.m_a[2],
				ru.m_a[0] * rv.m_a[1] - ru.m_a[1] * rv.m_a[0]
			};
		}
	};


	typedef vec<2> vec2;
	typedef vec<3> vec3;
	typedef vec<4> vec4;
}
