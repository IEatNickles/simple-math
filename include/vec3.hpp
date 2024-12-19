#pragma once

#include "vec2.hpp"

#include <cassert>
#include <cmath>

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
