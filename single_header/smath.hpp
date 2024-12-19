#ifndef SMATH_HPP
#define SMATH_HPP

#include <ostream>
#include <cassert>
#include <cmath>

struct vec2 {
  float x, y;

  constexpr vec2() : x(0.0f), y(0.0f) {}
  constexpr vec2(float x, float y) : x(x), y(y) {}
  constexpr vec2(float xy) : x(xy), y(xy) {}

  constexpr float const& operator[](int idx) const {
    assert(idx < 2);
    switch (idx) {
      default:
      case 0: return x;
      case 1: return y;
    }
  }
  constexpr float& operator[](int idx) {
    assert(idx < 2);
    switch (idx) {
      default:
      case 0: return x;
      case 1: return y;
    }
  }

  constexpr vec2 operator+(const vec2& b) const {
    return vec2(x + b.x, y + b.y);
  }
  constexpr vec2 operator-(const vec2& b) const {
    return vec2(x - b.x, y - b.y);
  }
  constexpr vec2 operator*(const vec2& b) const {
    return vec2(x * b.x, y * b.y);
  }

  constexpr vec2 operator*(float scalar) const {
    return vec2(x * scalar, y * scalar);
  }
  constexpr vec2 operator/(float scalar) const {
    return vec2(x / scalar, y / scalar);
  }

  constexpr static float dot(vec2 a, vec2 b) {
    return a.x * b.x + a.y * b.y;
  }

};

constexpr static vec2 operator*(float scalar, vec2 vec) {
  return vec2(scalar * vec.x, scalar * vec.y);
}
static vec2 operator/(float scalar, vec2 vec) {
  return vec2(scalar / vec.x, scalar / vec.y);
}

constexpr std::ostream& operator << (std::ostream& out, const vec2& v) {
  out << "[ " << v.x << ", " << v.y << " ]";
  return out;
}

struct vec3 {
  float x, y, z;

  constexpr vec3() : x(0.0f), y(0.0f), z(0.0f) {}
  constexpr vec3(float x, float y, float z) : x(x), y(y), z(z) {}
  constexpr vec3(float x, float y) : x(x), y(y), z(0.0f) {}
  constexpr vec3(float xyz) : x(xyz), y(xyz), z(xyz) {}
  constexpr vec3(vec2 xy, float z) : x(xy.x), y(xy.y), z(z) {}

  constexpr float const& operator[](int idx) const {
    assert(idx < 3);
    switch (idx) {
      default:
      case 0: return x;
      case 1: return y;
      case 2: return z;
    }
  }
  constexpr float& operator[](int idx) {
    assert(idx < 3);
    switch (idx) {
      default:
      case 0: return x;
      case 1: return y;
      case 2: return z;
    }
  }

  constexpr vec3 operator+(const vec3& b) const {
    return vec3(x + b.x, y + b.y, z + b.z);
  }
  constexpr vec3 operator-(const vec3& b) const {
    return vec3(x - b.x, y - b.y, z - b.z);
  }
  constexpr vec3 operator*(const vec3& b) const {
    return vec3(x * b.x, y * b.y, z * b.z);
  }

  constexpr vec3 operator*(float scalar) const {
    return vec3(x * scalar, y * scalar, z * scalar);
  }
  constexpr vec3 operator/(float scalar) const {
    return vec3(x / scalar, y / scalar, z / scalar);
  }

  constexpr static float dot(const vec3& a, const vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }
  constexpr static vec3 cross(const vec3& a, const vec3& b) {
    return vec3(a.y*b.z - a.z*b.y, -(a.x*b.z - a.z*b.x), a.x*b.y - a.y*b.x);
  }
};

constexpr static vec3 operator*(float scalar, vec3 vec) {
  return vec3(scalar * vec.x, scalar * vec.y, scalar * vec.z);
}
constexpr static vec3 operator/(float scalar, vec3 vec) {
  return vec3(scalar / vec.x, scalar / vec.y, scalar / vec.z);
}

constexpr std::ostream& operator << (std::ostream& out, const vec3& v) {
  out << "[ " << v.x << ", " << v.y << v.z << " ]";
  return out;
}

struct vec4 {
  float x, y, z, w;

  constexpr vec4() : x(0.0f), y(0.0f), z(0.0f), w(0.0f) {}
  constexpr vec4(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}
  constexpr vec4(float xyzw) : x(xyzw), y(xyzw), z(xyzw), w(xyzw) {}
  constexpr vec4(const vec3& xyz, float w) : x(xyz.x), y(xyz.y), z(xyz.z), w(w) {}
  constexpr vec4(const vec2& xy, const vec2& zw) : x(xy.x), y(xy.y), z(zw.x), w(zw.y) {}
  constexpr vec4(const vec2& xy, float z, float w) : x(xy.x), y(xy.y), z(z), w(w) {}

  constexpr float const& operator[](int idx) const {
    assert(idx < 4);
    switch (idx) {
      default:
      case 0: return x;
      case 1: return y;
      case 2: return z;
      case 3: return w;
    }
  }
  constexpr float& operator[](int idx) {
    assert(idx < 4);
    switch (idx) {
      default:
      case 0: return x;
      case 1: return y;
      case 2: return z;
      case 3: return w;
    }
  }

  constexpr vec4 operator+(const vec4& other) const {
    return vec4(x + other.x, y + other.y, z + other.z, w + other.w);
  }
  constexpr vec4 operator-(const vec4& other) const {
    return vec4(x - other.x, y - other.y, z - other.z, w - other.w);
  }
  constexpr vec4 operator*(const vec4& other) const {
    return vec4(x * other.x, y * other.y, z * other.z, w * other.w);
  }

  constexpr vec4 operator*(float scalar) const {
    return vec4(x * scalar, y * scalar, z * scalar, w * scalar);
  }
  constexpr vec4 operator/(float scalar) const {
    return vec4(x / scalar, y / scalar, z / scalar, w / scalar);
  }

  constexpr vec4 normalized() const {
    return *this / length();
  }

  constexpr float length() const {
    return std::sqrt(x * x + y * y + z * z + w * w);
  }

  constexpr float length_sq() const {
    return x * x + y * y + z * z + w * w;
  }

  constexpr static float dot(const vec4& a, const vec4& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
  }
};

constexpr vec4 operator*(float scalar, vec4 vec) {
  return vec4(scalar * vec.x, scalar * vec.y, scalar * vec.z, scalar * vec.w);
}

constexpr vec4 operator/(float scalar, vec4 vec) {
  return vec4(scalar / vec.x, scalar / vec.y, scalar / vec.z, scalar / vec.w);
}

constexpr std::ostream& operator << (std::ostream& out, const vec4& v) {
  out << "[ " << v.x << ", " << v.y << v.z << v.w << " ]";
  return out;
}
struct mat2 {
  vec2 value[2];

  constexpr mat2() : mat2(1.0f) {}
  constexpr mat2(float diag) : value {
    vec2(diag, 0.0f),
      vec2(0.0f, diag)} {}
  constexpr mat2(vec2 a, vec2 b) : value { a, b } {}

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

struct mat3 {
  vec3 value[3];

  constexpr mat3() : mat3(1.0f) {}
  constexpr mat3(float diag) : value {
    vec3(diag, 0.0f, 0.0f),
      vec3(0.0f, diag, 0.0f),
      vec3(0.0f, 0.0f, diag)} {}
  constexpr mat3(vec3 a, vec3 b, vec3 c) : value { a, b, c } {}

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

struct mat4 {
  vec4 value[4];

  constexpr mat4() : mat4(1.0f) {}
  constexpr mat4(float diag) : value {
    vec4(diag, 0.0f, 0.0f, 0.0f),
      vec4(0.0f, diag, 0.0f, 0.0f),
      vec4(0.0f, 0.0f, diag, 0.0f),
      vec4(0.0f, 0.0f, 0.0f, diag) }{ }
  constexpr mat4(vec4 a, vec4 b, vec4 c, vec4 d) : value { a, b, c, d } {}

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

#endif
