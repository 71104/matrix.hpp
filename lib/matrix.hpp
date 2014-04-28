#pragma once

#include <algorithm>
#include <cassert>


template<unsigned int _n, unsigned int _m, typename _Type = double>
struct mat {
	_Type m_a[_m][_n];

	mat() {}

	mat(_Type (&&a_a)[_m][_n])
		:
	m_a(move(a_a)) {}

	mat(_Type const &rValue) {
		for (unsigned int i = 0; i < _n; ++i) {
			for (unsigned int j = 0; j < _m; ++j) {
				m_a[i][j] = rValue;
			}
		}
	}

	struct row {
		_Type (&m_ra)[_m][_n];
		unsigned int m_i;

		row(_Type (&a_ra)[_m][_n], unsigned int a_i)
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

		const_row(_Type const (&a_ra)[_m][_n], unsigned int a_i)
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

	mat<_n, _m, _Type> operator + (mat<_n, _m, _Type> const &r) const {
		mat<_n, _m> Result;
		for (unsigned int i = 0; i < _n; ++i) {
			for (unsigned int j = 0; j < _m; ++j) {
				Result[i][j] += m_a[j][i] + r[i][j];
			}
		}
		return Result;
	}

	mat<_n, _m, _Type> operator - (mat<_n, _m, _Type> const &r) const {
		mat<_n, _m> Result;
		for (unsigned int i = 0; i < _n; ++i) {
			for (unsigned int j = 0; j < _m; ++j) {
				Result[i][j] += m_a[j][i] + r[i][j];
			}
		}
		return Result;
	}

	template<unsigned int _l>
	mat<_n, _l> operator * (mat<_m, _l, _Type> const &r) const {
		mat<_n, _l> Result;
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


typedef mat<2, 2> mat2;
typedef mat<2, 3> mat23;
typedef mat<2, 4> mat24;
typedef mat<3, 2> mat32;
typedef mat<3, 3> mat3;
typedef mat<3, 4> mat34;
typedef mat<4, 2> mat42;
typedef mat<4, 3> mat43;
typedef mat<4, 4> mat4;


template<unsigned int _n, typename _Type = double>
struct vec : public mat<_n, 1, _Type> {};


typedef vec<2> vec2;
typedef vec<3> vec3;
typedef vec<4> vec4;
