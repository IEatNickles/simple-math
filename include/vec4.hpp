#pragma once

#include "vec2.hpp"
#include "vec3.hpp"

#include <cassert>
#include <cmath>

struct vec4 {
  constexpr vec4() : x(0.0f), y(0.0f), z(0.0f), w(0.0f) {}
  constexpr vec4(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}
  constexpr vec4(float xyzw) : x(xyzw), y(xyzw), z(xyzw), w(xyzw) {}
  constexpr vec4(const vec3& xyz, float w) : x(xyz.x), y(xyz.y), z(xyz.z), w(w) {}
  constexpr vec4(const vec2& xy, const vec2& zw) : x(xy.x), y(xy.y), z(zw.x), w(zw.y) {}
  constexpr vec4(const vec2& xy, float z, float w) : x(xy.x), y(xy.y), z(z), w(w) {}

  float x, y, z, w;

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
