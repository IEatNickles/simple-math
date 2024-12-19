#pragma once

#include <cstdio>
#include <ostream>
#include <string>

#include "vec2.hpp"

struct mat2 {
    constexpr mat2() : mat2(1.0f) {}
    constexpr mat2(float diag) : value {
        vec2(diag, 0.0f),
        vec2(0.0f, diag)} {}
    constexpr mat2(vec2 a, vec2 b) : value { a, b } {}

    vec2 value[2];

    constexpr vec2 const& operator[](int idx) const {
      return value[idx];
    }
    constexpr vec2& operator[](int idx) {
      return value[idx];
    }

    constexpr mat2 operator+(const mat2& other) const {
      mat2 result;
      result[0] = value[0] + other[0];
      result[1] = value[1] + other[1];
      return result;
    }
    constexpr mat2 operator-(const mat2& other) const {
      mat2 result;
      result[0] = value[0] - other[0];
      result[1] = value[1] - other[1];
      return result;
    }
    constexpr mat2 operator*(const mat2& other) const;

    constexpr mat2 operator*(float scalar) const {
      mat2 result;
      result[0] = value[0] * scalar;
      result[1] = value[1] * scalar;
      return result;
    }
    constexpr mat2 operator/(float scalar) const {
      mat2 result;
      result[0] = value[0] / scalar;
      result[1] = value[1] / scalar;
      return result;
    }

    constexpr mat2 transpose() const {
      return mat2(get_column(0), get_column(1));
    }
    constexpr mat2 inverse() const {
      return adjugate() / determinant();
    }

    constexpr float determinant() const {
      return value[0][0] * value[1][1] - value[0][1] * value[1][0];
    }
    constexpr mat2 adjugate() const {
      mat2 adj;
      for (int r = 0; r < 3; r ++) {
        for (int c = 0; c < 3; c ++) {
          float m = value[r][c];
          if (((r + c) & 1) != 0) {
            m = -m;
          }
          adj.value[r][c] = m;
        }
      }
      return adj.transpose();
    }

    constexpr vec2 get_column(int idx) const {
      return vec2(value[0][idx], value[1][idx]);
    }
};

constexpr std::ostream& operator << (std::ostream& out, const mat2& m) {
  out << m[0] << "\n";
  out << m[1];
  return out;
}
