
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#ifndef PBRT_CORE_TRANSFORM_H
#define PBRT_CORE_TRANSFORM_H

#include "../core/geometry.h"
#include "quaternion.h"

// Matrix4x4 Declarations
struct Matrix4x4 {
    // Matrix4x4 Public Methods
    Matrix4x4() {
        m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.f;
        m[0][1] = m[0][2] = m[0][3] = m[1][0] =
             m[1][2] = m[1][3] = m[2][0] = m[2][1] = m[2][3] =
             m[3][0] = m[3][1] = m[3][2] = 0.0;
    }
    Matrix4x4(double mat[4][4]);
    Matrix4x4(double t00, double t01, double t02, double t03,
              double t10, double t11, double t12, double t13,
              double t20, double t21, double t22, double t23,
              double t30, double t31, double t32, double t33);
    bool operator==(const Matrix4x4 &m2) const {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (m[i][j] != m2.m[i][j]) return false;
        return true;
    }
    bool operator!=(const Matrix4x4 &m2) const {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (m[i][j] != m2.m[i][j]) return true;
        return false;
    }
    friend Matrix4x4 Transpose(const Matrix4x4 &);
    void Print(FILE *f) const {
        fprintf(f, "[ ");
        for (int i = 0; i < 4; ++i) {
            fprintf(f, "  [ ");
            for (int j = 0; j < 4; ++j)  {
                fprintf(f, "%f", m[i][j]);
                if (j != 3) fprintf(f, ", ");
            }
            fprintf(f, " ]\n");
        }
        fprintf(f, " ] ");
    }
    static Matrix4x4 Mul(const Matrix4x4 &m1, const Matrix4x4 &m2) {
        Matrix4x4 r;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                r.m[i][j] = m1.m[i][0] * m2.m[0][j] +
                            m1.m[i][1] * m2.m[1][j] +
                            m1.m[i][2] * m2.m[2][j] +
                            m1.m[i][3] * m2.m[3][j];
        return r;
    }
    friend Matrix4x4 Inverse(const Matrix4x4 &);
    double m[4][4];
};



// Transform Declarations
class Transform {
public:
    // Transform Public Methods
    Transform() { }
    Transform(const double mat[4][4]) {
        m = Matrix4x4(mat[0][0], mat[0][1], mat[0][2], mat[0][3],
                      mat[1][0], mat[1][1], mat[1][2], mat[1][3],
                      mat[2][0], mat[2][1], mat[2][2], mat[2][3],
                      mat[3][0], mat[3][1], mat[3][2], mat[3][3]);
        mInv = Inverse(m);
    }
    Transform(const Matrix4x4 &mat)
        : m(mat), mInv(Inverse(mat)) {
    }
    Transform(const Matrix4x4 &mat, const Matrix4x4 &minv)
       : m(mat), mInv(minv) {
    }
    void Print(FILE *f) const;
    friend Transform Inverse(const Transform &t) {
        return Transform(t.mInv, t.m);
    }
    friend Transform Transpose(const Transform &t) {
        return Transform(Transpose(t.m), Transpose(t.mInv));
    }
    bool operator==(const Transform &t) const {
        return t.m == m && t.mInv == mInv;
    }
    bool operator!=(const Transform &t) const {
        return t.m != m || t.mInv != mInv;
    }
    bool operator<(const Transform &t2) const {
        for (uint32_t i = 0; i < 4; ++i)
            for (uint32_t j = 0; j < 4; ++j) {
                if (m.m[i][j] < t2.m.m[i][j]) return true;
                if (m.m[i][j] > t2.m.m[i][j]) return false;
            }
        return false;
    }
    bool IsIdentity() const {
        return (m.m[0][0] == 1.f && m.m[0][1] == 0.0 &&
                m.m[0][2] == 0.0 && m.m[0][3] == 0.0 &&
                m.m[1][0] == 0.0 && m.m[1][1] == 1.f &&
                m.m[1][2] == 0.0 && m.m[1][3] == 0.0 &&
                m.m[2][0] == 0.0 && m.m[2][1] == 0.0 &&
                m.m[2][2] == 1.f && m.m[2][3] == 0.0 &&
                m.m[3][0] == 0.0 && m.m[3][1] == 0.0 &&
                m.m[3][2] == 0.0 && m.m[3][3] == 1.f);
    }
    const Matrix4x4 &GetMatrix() const { return m; }
    const Matrix4x4 &GetInverseMatrix() const { return mInv; }
    bool HasScale() const {
        double la2 = (*this)(Vector3D(1,0,0)).LengthSquared();
        double lb2 = (*this)(Vector3D(0,1,0)).LengthSquared();
        double lc2 = (*this)(Vector3D(0,0,1)).LengthSquared();
#define NOT_ONE(x) ((x) < .999f || (x) > 1.001f)
        return (NOT_ONE(la2) || NOT_ONE(lb2) || NOT_ONE(lc2));
#undef NOT_ONE
    }
    inline Point3D operator()(const Point3D &pt) const;
    inline void operator()(const Point3D &pt, Point3D *ptrans) const;
    inline Vector3D operator()(const Vector3D &v) const;
    inline void operator()(const Vector3D &v, Vector3D *vt) const;
    inline Normal operator()(const Normal &) const;
    inline void operator()(const Normal &, Normal *nt) const;
    inline Ray operator()(const Ray &r) const;
    inline void operator()(const Ray &r, Ray *rt) const;
    inline RayDifferential operator()(const RayDifferential &r) const;
    inline void operator()(const RayDifferential &r, RayDifferential *rt) const;
    BBox operator()(const BBox &b) const;
    Transform operator*(const Transform &t2) const;
    bool SwapsHandedness() const;
private:
    // Transform Private Data
    Matrix4x4 m, mInv;
    friend class AnimatedTransform;
    friend struct Quaternion;
};


Transform Translate(const Vector3D &delta);
Transform Scale(double x, double y, double z);
Transform RotateX(double angle);
Transform RotateY(double angle);
Transform RotateZ(double angle);
Transform Rotate(double angle, const Vector3D &axis);
Transform Rotatec(double cosAngle, double sinAngle, const Vector3D &axis);
Transform LookAt(const Point3D &pos, const Point3D &look, const Vector3D &up);
bool SolveLinearSystem2x2(const double A[2][2], const double B[2],
    double *x0, double *x1);
Transform Orthographic(double znear, double zfar);
Transform Perspective(double fov, double znear, double zfar);

// Transform Inline Functions
inline Point3D Transform::operator()(const Point3D &pt) const {
    double x = pt.x, y = pt.y, z = pt.z;
    double xp = m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z + m.m[0][3];
    double yp = m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z + m.m[1][3];
    double zp = m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z + m.m[2][3];
    double wp = m.m[3][0]*x + m.m[3][1]*y + m.m[3][2]*z + m.m[3][3];
    Assert(wp != 0);
    if (wp == 1.) return Point3D(xp, yp, zp);
    else          return Point3D(xp, yp, zp)/wp;
}


inline void Transform::operator()(const Point3D &pt,
                                  Point3D *ptrans) const {
    double x = pt.x, y = pt.y, z = pt.z;
    ptrans->x = m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z + m.m[0][3];
    ptrans->y = m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z + m.m[1][3];
    ptrans->z = m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z + m.m[2][3];
    double w   = m.m[3][0]*x + m.m[3][1]*y + m.m[3][2]*z + m.m[3][3];
    if (w != 1.) *ptrans /= w;
}


inline Vector3D Transform::operator()(const Vector3D &v) const {
  double x = v.x, y = v.y, z = v.z;
  return Vector3D(m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z,
                m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z,
                m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z);
}


inline void Transform::operator()(const Vector3D &v,
        Vector3D *vt) const {
  double x = v.x, y = v.y, z = v.z;
  vt->x = m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z;
  vt->y = m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z;
  vt->z = m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z;
}


inline Normal Transform::operator()(const Normal &n) const {
    double x = n.x, y = n.y, z = n.z;
    return Normal(mInv.m[0][0]*x + mInv.m[1][0]*y + mInv.m[2][0]*z,
                  mInv.m[0][1]*x + mInv.m[1][1]*y + mInv.m[2][1]*z,
                  mInv.m[0][2]*x + mInv.m[1][2]*y + mInv.m[2][2]*z);
}


inline void Transform::operator()(const Normal &n,
        Normal *nt) const {
    double x = n.x, y = n.y, z = n.z;
    nt->x = mInv.m[0][0] * x + mInv.m[1][0] * y +
        mInv.m[2][0] * z;
    nt->y = mInv.m[0][1] * x + mInv.m[1][1] * y +
        mInv.m[2][1] * z;
    nt->z = mInv.m[0][2] * x + mInv.m[1][2] * y +
        mInv.m[2][2] * z;
}


inline Ray Transform::operator()(const Ray &r) const {
    Ray ret = r;
    (*this)(ret.m_origin, &ret.m_origin);
    (*this)(ret.m_direction, &ret.m_direction);
    return ret;
}


inline void Transform::operator()(const Ray &r, Ray *rt) const {
    (*this)(r.m_origin, &rt->m_origin);
    (*this)(r.m_direction, &rt->m_direction);
    if (rt != &r) {
        rt->m_mint = r.m_mint;
        rt->m_maxt = r.m_maxt;
        rt->m_time = r.m_time;
        rt->m_depth = r.m_depth;
    }
}



inline void Transform::operator()(const RayDifferential &r, RayDifferential *rt) const {
    (*this)(Ray(r), rt);
    rt->hasDifferentials = r.hasDifferentials;
    (*this)(r.rxOrigin, &rt->rxOrigin);
    (*this)(r.ryOrigin, &rt->ryOrigin);
    (*this)(r.rxDirection, &rt->rxDirection);
    (*this)(r.ryDirection, &rt->ryDirection);
}



inline RayDifferential Transform::operator()(const RayDifferential &r) const {
    RayDifferential ret;
    (*this)(Ray(r), &ret);
    ret.hasDifferentials = r.hasDifferentials;
    (*this)(r.rxOrigin, &ret.rxOrigin);
    (*this)(r.ryOrigin, &ret.ryOrigin);
    (*this)(r.rxDirection, &ret.rxDirection);
    (*this)(r.ryDirection, &ret.ryDirection);
    return ret;
}




// AnimatedTransform Declarations
class AnimatedTransform {
public:
    // AnimatedTransform Public Methods
    AnimatedTransform(const Transform *transform1, double time1,
                      const Transform *transform2, double time2)
        : startTime(time1), endTime(time2),
          startTransform(transform1), endTransform(transform2),
          actuallyAnimated(*startTransform != *endTransform) {
        Decompose(startTransform->m, &T[0], &R[0], &S[0]);
        Decompose(endTransform->m, &T[1], &R[1], &S[1]);
    }
    static void Decompose(const Matrix4x4 &m, Vector3D *T, Quaternion *R, Matrix4x4 *S);
    void Interpolate(double time, Transform *t) const;
    void operator()(const Ray &r, Ray *tr) const;
    void operator()(const RayDifferential &r, RayDifferential *tr) const;
    Point3D operator()(double time, const Point3D &p) const;
    Vector3D operator()(double time, const Vector3D &v) const;
    Ray operator()(const Ray &r) const;
    BBox MotionBounds(const BBox &b, bool useInverse) const;
    bool HasScale() const { return startTransform->HasScale() || endTransform->HasScale(); }
private:
    // AnimatedTransform Private Data
    const double startTime, endTime;
    const Transform *startTransform, *endTransform;
    const bool actuallyAnimated;
    Vector3D T[2];
    Quaternion R[2];
    Matrix4x4 S[2];
};



#endif // PBRT_CORE_TRANSFORM_H
