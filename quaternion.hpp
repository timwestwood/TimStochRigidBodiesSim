// quaternion.hpp

// =============================================================================
// Include guard
#ifndef MY_QUATERNION_HEADER_INCLUDED
#define MY_QUATERNION_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies

// =============================================================================
// Included dependencies
#include <armadillo>
using namespace arma;

class quaternion {

public:

  double scalar_part;
  vec3 vector_part;

  ~quaternion();
  quaternion();
  quaternion(const double q0, const double q1, const double q2, const double q3);
  quaternion(const double q0, const vec& q);
  quaternion(const quaternion& quaternion_to_copy);

  double norm() const;
  void normalise_in_place();
  quaternion conj() const;
  mat33 matrix() const;
  vec rotate(const vec& v) const;
  void sqrt_in_place();

  // Operator overloading

  // q*s
  quaternion operator *(const double& s) const {

    quaternion out(s*scalar_part, s*vector_part);

    return out;

  };

  // s*q
  friend quaternion operator *(const double& s, const quaternion& q){

    return q*s;

  };

  // q*p
  quaternion operator *(const quaternion& p) const {

    double q0 = scalar_part*p.scalar_part - dot(vector_part, p.vector_part);
    vec q = scalar_part*p.vector_part + p.scalar_part*vector_part + cross(vector_part, p.vector_part);

    return quaternion(q0, q);

  };

  // q/s
  quaternion operator /(const double& s) const {

    return quaternion(scalar_part/s, vector_part/s);

  };

  // q+p
  quaternion operator +(const quaternion& p) const {

    return quaternion(scalar_part+p.scalar_part, vector_part+p.vector_part);

  };

  quaternion& operator +=(const quaternion& p) {

    scalar_part += p.scalar_part;
    vector_part += p.vector_part;

    return *this;

  };

  //q-p
  quaternion operator -(const quaternion& p) const {

    return quaternion(scalar_part-p.scalar_part, vector_part-p.vector_part);

  };

  quaternion& operator -=(const quaternion& p) {

    scalar_part -= p.scalar_part;
    vector_part -= p.vector_part;

    return *this;

  };

}; // End of quaternion class

quaternion lie_exp(const vec& u);

#endif // MY_QUATERNION_HEADER_INCLUDED
