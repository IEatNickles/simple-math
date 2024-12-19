#pragma once

#include "vec4.hpp"
#include "mat3.hpp"

struct mat4 {
  constexpr mat4() : mat4(1.0f) {}
  constexpr mat4(float diag) : value {
    vec4(diag, 0.0f, 0.0f, 0.0f),
      vec4(0.0f, diag, 0.0f, 0.0f),
      vec4(0.0f, 0.0f, diag, 0.0f),
      vec4(0.0f, 0.0f, 0.0f, diag) }{ }
  constexpr mat4(vec4 a, vec4 b, vec4 c, vec4 d) : value { a, b, c, d } {}

  vec4 value[4];

  constexpr mat4 operator+(const mat4& other) const{
    mat4 res;
    res[0] = value[0] + other.value[0];
    res[0] = value[1] + other.value[1];
    res[0] = value[2] + other.value[2];
    res[0] = value[3] + other.value[3];
    return res;
  }
  constexpr mat4 operator-(const mat4& other) const{
    mat4 res;
    res[0] = value[0] - other.value[0];
    res[0] = value[1] - other.value[1];
    res[0] = value[2] - other.value[2];
    res[0] = value[3] - other.value[3];
    return res;
  }
  constexpr mat4 operator*(const mat4& other) const{
    mat4 res;
    res[0][0] = vec4::dot(value[0], other.column(0));
    res[0][1] = vec4::dot(value[0], other.column(1));
    res[0][2] = vec4::dot(value[0], other.column(2));
    res[0][3] = vec4::dot(value[0], other.column(3));

    res[1][0] = vec4::dot(value[1], other.column(0));
    res[1][1] = vec4::dot(value[1], other.column(1));
    res[1][2] = vec4::dot(value[1], other.column(2));
    res[1][3] = vec4::dot(value[1], other.column(3));

    res[2][0] = vec4::dot(value[2], other.column(0));
    res[2][1] = vec4::dot(value[2], other.column(1));
    res[2][2] = vec4::dot(value[2], other.column(2));
    res[2][3] = vec4::dot(value[2], other.column(3));

    res[3][0] = vec4::dot(value[3], other.column(0));
    res[3][1] = vec4::dot(value[3], other.column(1));
    res[3][2] = vec4::dot(value[3], other.column(2));
    res[3][3] = vec4::dot(value[3], other.column(3));
    return res;
  }

  constexpr mat4 operator*(float scalar) const {
    mat4 res;
    res[0] = value[0] * scalar;
    res[1] = value[1] * scalar;
    res[2] = value[2] * scalar;
    res[3] = value[3] * scalar;
    return res;
  }
  constexpr mat4 operator/(float scalar) const{
    mat4 res;
    res[0] = value[0] / scalar;
    res[1] = value[1] / scalar;
    res[2] = value[2] / scalar;
    res[3] = value[3] / scalar;
    return res;
  }

  constexpr vec4 const& operator[](int idx) const {
    return value[idx];
  }
  constexpr vec4& operator[](int idx) {
    return value[idx];
  }

  constexpr mat4 transpose() const {
    return mat4(column(0), column(1), column(2), column(3));
  }
  constexpr mat4 inverse() const {
    return adjugate() / determinant();
  }

  constexpr float determinant() const {
    float a = value[0][0];
    float b = value[0][1];
    float c = value[0][2];
    float d = value[0][3];

    float e = value[1][0];
    float f = value[1][1];
    float g = value[1][2];
    float h = value[1][3];

    float i = value[2][0];
    float j = value[2][1];
    float k = value[2][2];
    float l = value[2][3];

    float m = value[3][0];
    float n = value[3][1];
    float o = value[3][2];
    float p = value[3][3];

    return a * f * (k * p - o * l) - g * (j * p - n * l) + h * (j * o - n * k) -
      b * e * (k * p - o * l) - g * (i * p - m * l) + h * (i * o - m * k) +
      c * e * (j * p - n * l) - f * (i * p - m * l) + h * (i * n - m * j) -
      d * e * (j * o - n * k) - f * (i * o - m * k) + g * (i * n - m * j);
  }
  constexpr mat4 adjugate() const{
    mat4 adj;
    for (int r = 0; r < 4; r ++) {
      for (int c = 0; c < 4; c ++) {
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
    mat3 minor;
    int row = 0;
    for (int r = 0; r < 3; r++) {
      if (r == i) row += 1;
      int col = 0;
      for (int c = 0; c < 3; c++) {
        if (c == j) col += 1;
        minor.value[r][c] = value[row][col];
        col += 1;
      }
      row += 1;
    }
    return minor.determinant();
  }

  constexpr vec4 column(int idx) const {
    return vec4(value[0][idx], value[1][idx], value[2][idx], value[3][idx]);
  }

  constexpr static mat4 translation(const vec3& translation) {
    mat4 result;
    result.value[3] = vec4(translation, 1.0f);
    return result;
  }
  constexpr static mat4 rotation(const vec3& rotation){
    float sina = std::sin(rotation.x);
    float cosa = std::cos(rotation.x);
    float sinb = std::sin(rotation.y);
    float cosb = std::cos(rotation.y);
    float sinc = std::sin(rotation.z);
    float cosc = std::cos(rotation.z);

    mat4 result;
    result[0][0] =  cosc * cosb;
    result[0][1] =  cosc * sinb * sina - sinc * cosa;
    result[0][2] =  cosc * sinb * cosa + sinc * sina;

    result[1][0] =  sinc * cosb;
    result[1][1] =  sinc * sinb * sina + cosc * cosa;
    result[1][2] =  sinc * sinb * cosa - cosc * sina;

    result[2][0] = -sinb;
    result[2][1] =  cosb * sina;
    result[2][2] =  cosb * cosa;
    return result;
  }
  constexpr static mat4 scale(const vec3& scale){
    mat4 result;
    result[0][0] = scale.x;
    result[1][1] = scale.y;
    result[2][2] = scale.z;
    return result;
  }
  constexpr static mat4 shear(const vec3& shear){
    mat4 result;
    result[0][1] = shear.y;
    result[0][2] = shear.z;

    result[1][0] = shear.x;
    result[1][2] = shear.z;

    result[2][0] = shear.x;
    result[2][1] = shear.y;
    return result;
  }

  constexpr static mat4 orthographic(float right, float left, float top, float bottom, float near, float far) {
    float r = right, l = left;
    float t = top, b = bottom;
    float f = far, n = near;
    mat4 result;
    float inv_width =  1.0f / (r - l);
    float inv_height = 1.0f / (t - b);
    float inv_length = 1.0f / (f - n);

    result[0][0] = 2.0f * inv_width;
    result[1][1] = 2.0f * inv_height;
    result[2][2] = 2.0f * inv_length;

    result[0][3] = -((r + l) * inv_width);
    result[1][3] = -((t + b) * inv_height);
    result[2][3] = -((f + n) * inv_length);
    return result;
  }
  constexpr static mat4 perspective(float fov, float aspect, float near, float far) {
    float a = aspect;
    float f = 1.0f / std::tan(fov * 0.5f);
    float l = far / (far - near);

    mat4 result;
    result.value[0][0] = f * a;
    result.value[1][1] = f;
    result.value[2][2] = l;
    result.value[2][3] = 1.0f;
    result.value[3][2] = -l * near;
    return result;
  }
};

constexpr std::ostream& operator << (std::ostream& out, const mat4& m) {
  out << m[0] << "\n";
  out << m[1] << "\n";
  out << m[2] << "\n";
  out << m[3];
  return out;
}
