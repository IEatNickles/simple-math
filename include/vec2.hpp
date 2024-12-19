#pragma once

#include <cassert>
#include <cmath>
#include <ostream>

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
