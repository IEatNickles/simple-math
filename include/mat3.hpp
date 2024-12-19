#pragma once

#include "vec3.hpp"
#include "mat2.hpp"

struct mat3 {
    constexpr mat3() : mat3(1.0f) {}
    constexpr mat3(float diag) : value {
        vec3(diag, 0.0f, 0.0f),
        vec3(0.0f, diag, 0.0f),
        vec3(0.0f, 0.0f, diag)} {}
    constexpr mat3(vec3 a, vec3 b, vec3 c) : value { a, b, c } {}

    vec3 value[3];

    constexpr vec3 const& operator[](int idx) const {
      return value[idx];
    }
    constexpr vec3& operator[](int idx) {
      return value[idx];
    }

    constexpr mat3 operator+(const mat3& other) const{
      mat3 res;
      res[0] = value[0] + other[0];
      res[1] = value[1] + other[1];
      res[2] = value[2] + other[2];
      return res;
    }
    constexpr mat3 operator-(const mat3& other) const{
      mat3 res;
      res[0] = value[0] - other[0];
      res[1] = value[1] - other[1];
      res[2] = value[2] - other[2];
      return res;
    }
    constexpr mat3 operator*(const mat3& other) const {
      mat3 res;
      res[0][0] = vec3::dot(value[0], other.col(0));
      res[0][1] = vec3::dot(value[0], other.col(1));
      res[0][2] = vec3::dot(value[0], other.col(2));

      res[1][0] = vec3::dot(value[1], other.col(0));
      res[1][1] = vec3::dot(value[1], other.col(1));
      res[1][2] = vec3::dot(value[1], other.col(2));

      res[2][0] = vec3::dot(value[2], other.col(0));
      res[2][1] = vec3::dot(value[2], other.col(1));
      res[2][2] = vec3::dot(value[2], other.col(2));
      return res;
    }

    constexpr mat3 operator*(float scalar) const {
      mat3 res;
      res[0] = value[0] * scalar;
      res[1] = value[1] * scalar;
      res[2] = value[2] * scalar;
      return res;
    }
    constexpr mat3 operator/(float scalar) const {
      mat3 res;
      res[0] = value[0] / scalar;
      res[1] = value[1] / scalar;
      res[2] = value[2] / scalar;
      return res;
    }

    constexpr mat3 transpose() const {
      return mat3(col(0), col(1), col(2));
    }
    constexpr mat3 inverse() const {
      return adjugate() / determinant();
    }

    constexpr float determinant() const {
      float a = value[0][0];
      float b = value[0][1];
      float c = value[0][2];

      float d = value[1][0];
      float e = value[1][1];
      float f = value[1][2];

      float g = value[2][0];
      float h = value[2][1];
      float k = value[2][2];

      return a * (e*k - f*h) - b * (d*k - f*g) + c * (d*h - e*g);
    }
    constexpr mat3 adjugate() const {
      mat3 adj;
      for (int r = 0; r < 3; r ++) {
        for (int c = 0; c < 3; c ++) {
          float m = minor(r, c);
          if (((r + c) & 1) != 0) {
            m = -m;
          }
          adj.value[r][c] = m;
        }
      }
      return adj.transpose();
    }

    constexpr float minor(int i, int j) const {
      mat2 minor;
      int row = 0;
      for (int r = 0; r < 2; r++) {
        if (r == i) {
          row += 1;
        }
        int col = 0;
        for (int c = 0; c < 2; c++) {
          if (c == j) {
            col += 1;
          }
          minor.value[r][c] = value[row][col];
          col += 1;
        }
        row += 1;
      }
      return minor.determinant();
    }

    constexpr vec3 col(int idx) const {
      return vec3(value[0][idx], value[1][idx], value[2][idx]);
    }
};

constexpr std::ostream& operator << (std::ostream& out, const mat3& m) {
  out << m[0] << "\n";
  out << m[1] << "\n";
  out << m[2];
  return out;
}
