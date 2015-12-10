
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


#include "transform.h"
#include "../core/common.h"
#include <memory.h>

// Matrix4x4 Method Definitions
bool solveLinearSystem2x2(const double A[2][2],
        const double B[2], double *x0, double *x1) {
    double det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    if (fabsf(det) < 1e-10f)
        return false;
    *x0 = (A[1][1]*B[0] - A[0][1]*B[1]) / det;
    *x1 = (A[0][0]*B[1] - A[1][0]*B[0]) / det;
    if (isnan(*x0) || isnan(*x1))
        return false;
    return true;
}


Matrix4x4::Matrix4x4(double mat[4][4]) {
    memcpy(m, mat, 16*sizeof(double));
}


Matrix4x4::Matrix4x4(double t00, double t01, double t02, double t03,
                     double t10, double t11, double t12, double t13,
                     double t20, double t21, double t22, double t23,
                     double t30, double t31, double t32, double t33) {
    m[0][0] = t00; m[0][1] = t01; m[0][2] = t02; m[0][3] = t03;
    m[1][0] = t10; m[1][1] = t11; m[1][2] = t12; m[1][3] = t13;
    m[2][0] = t20; m[2][1] = t21; m[2][2] = t22; m[2][3] = t23;
    m[3][0] = t30; m[3][1] = t31; m[3][2] = t32; m[3][3] = t33;
}


Matrix4x4 Transpose(const Matrix4x4 &m) {
   return Matrix4x4(m.m[0][0], m.m[1][0], m.m[2][0], m.m[3][0],
                    m.m[0][1], m.m[1][1], m.m[2][1], m.m[3][1],
                    m.m[0][2], m.m[1][2], m.m[2][2], m.m[3][2],
                    m.m[0][3], m.m[1][3], m.m[2][3], m.m[3][3]);
}


Matrix4x4 Inverse(const Matrix4x4 &m) {
    int indxc[4], indxr[4];
    int ipiv[4] = { 0, 0, 0, 0 };
    double minv[4][4];
    memcpy(minv, m.m, 4*4*sizeof(double));
    for (int i = 0; i < 4; i++) {
        int irow = -1, icol = -1;
        double big = 0.;
        // Choose pivot
        for (int j = 0; j < 4; j++) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < 4; k++) {
                    if (ipiv[k] == 0) {
                        if (fabsf(minv[j][k]) >= big) {
                            big = double(fabsf(minv[j][k]));
                            irow = j;
                            icol = k;
                        }
                    }
                    else if (ipiv[k] > 1)
                        Error("Singular matrix in MatrixInvert");
                }
            }
        }
        ++ipiv[icol];
        // Swap rows _irow_ and _icol_ for pivot
        if (irow != icol) {
            for (int k = 0; k < 4; ++k)
                swap(minv[irow][k], minv[icol][k]);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (minv[icol][icol] == 0.)
            Error("Singular matrix in MatrixInvert");

        // Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
        double pivinv = 1.f / minv[icol][icol];
        minv[icol][icol] = 1.f;
        for (int j = 0; j < 4; j++)
            minv[icol][j] *= pivinv;

        // Subtract this row from others to zero out their columns
        for (int j = 0; j < 4; j++) {
            if (j != icol) {
                double save = minv[j][icol];
                minv[j][icol] = 0;
                for (int k = 0; k < 4; k++)
                    minv[j][k] -= minv[icol][k]*save;
            }
        }
    }
    // Swap columns to reflect permutation
    for (int j = 3; j >= 0; j--) {
        if (indxr[j] != indxc[j]) {
            for (int k = 0; k < 4; k++)
                swap(minv[k][indxr[j]], minv[k][indxc[j]]);
        }
    }
    return Matrix4x4(minv);
}



// Transform Method Definitions
//void Transform::print(FILE *f) const {
//    m.Print(f);
//}


Transform translate(const Vector3D &delta) {
    Matrix4x4 m(1, 0, 0, delta.x.value(),
                0, 1, 0, delta.y.value(),
                0, 0, 1, delta.z.value(),
                0, 0, 0,       1);
    Matrix4x4 minv(1, 0, 0, -delta.x.value(),
                   0, 1, 0, -delta.y.value(),
                   0, 0, 1, -delta.z.value(),
                   0, 0, 0,        1);
    return Transform(m, minv);
}


Transform scale(double x, double y, double z) {
    Matrix4x4 m(x, 0, 0, 0,
                0, y, 0, 0,
                0, 0, z, 0,
                0, 0, 0, 1);
    Matrix4x4 minv(1.f/x,     0,     0,     0,
                   0,     1.f/y,     0,     0,
                   0,         0,     1.f/z, 0,
                   0,         0,     0,     1);
    return Transform(m, minv);
}


Transform rotateX(double angle) {
    double sin_t = sinf((angle));
    double cos_t = cosf((angle));
    Matrix4x4 m(1,     0,      0, 0,
                0, cos_t, -sin_t, 0,
                0, sin_t,  cos_t, 0,
                0,     0,      0, 1);
    return Transform(m, Transpose(m));
}


Transform rotateY(double angle) {
    double sin_t = sinf((angle));
    double cos_t = cosf((angle));
    Matrix4x4 m( cos_t,   0,  sin_t, 0,
                 0,   1,      0, 0,
                -sin_t,   0,  cos_t, 0,
                 0,   0,   0,    1);
    return Transform(m, Transpose(m));
}



Transform rotateZ(double angle) {
    double sin_t = sinf((angle));
    double cos_t = cosf((angle));
    Matrix4x4 m(cos_t, -sin_t, 0, 0,
                sin_t,  cos_t, 0, 0,
                0,      0, 1, 0,
                0,      0, 0, 1);
    return Transform(m, Transpose(m));
}


Transform rotate(double angle, const Vector3D &axis) {
    Vector3D a = normalize(axis);
    double s = sinf((angle));
    double c = cosf((angle));
    double m[4][4];

    // TODO: Verify use of units here

    m[0][0] = a.x.value() * a.x.value() + (1.f - a.x.value() * a.x.value()) * c;
    m[0][1] = a.x.value() * a.y.value() * (1.f - c) - a.z.value() * s;
    m[0][2] = a.x.value() * a.z.value() * (1.f - c) + a.y.value() * s;
    m[0][3] = 0;

    m[1][0] = a.x.value() * a.y.value() * (1.f - c) + a.z.value() * s;
    m[1][1] = a.y.value() * a.y.value() + (1.f - a.y.value() * a.y.value()) * c;
    m[1][2] = a.y.value() * a.z.value() * (1.f - c) - a.x.value() * s;
    m[1][3] = 0;

    m[2][0] = a.x.value() * a.z.value() * (1.f - c) - a.y.value() * s;
    m[2][1] = a.y.value() * a.z.value() * (1.f - c) + a.x.value() * s;
    m[2][2] = a.z.value() * a.z.value() + (1.f - a.z.value() * a.z.value()) * c;
    m[2][3] = 0;

    m[3][0] = 0;
    m[3][1] = 0;
    m[3][2] = 0;
    m[3][3] = 1;

    Matrix4x4 mat(m);
    return Transform(mat, Transpose(mat));
}

Transform rotatec(double cosAngle, double sinAngle, const Vector3D &axis) {
    Vector3D a = normalize(axis);
    double s = sinAngle;
    double c = cosAngle;
    double m[4][4];

    m[0][0] = a.x.value() * a.x.value() + (1.f - a.x.value() * a.x.value()) * c;
    m[0][1] = a.x.value() * a.y.value() * (1.f - c) - a.z.value() * s;
    m[0][2] = a.x.value() * a.z.value() * (1.f - c) + a.y.value() * s;
    m[0][3] = 0;

    m[1][0] = a.x.value() * a.y.value() * (1.f - c) + a.z.value() * s;
    m[1][1] = a.y.value() * a.y.value() + (1.f - a.y.value() * a.y.value()) * c;
    m[1][2] = a.y.value() * a.z.value() * (1.f - c) - a.x.value() * s;
    m[1][3] = 0;

    m[2][0] = a.x.value() * a.z.value() * (1.f - c) - a.y.value() * s;
    m[2][1] = a.y.value() * a.z.value() * (1.f - c) + a.x.value() * s;
    m[2][2] = a.z.value() * a.z.value() + (1.f - a.z.value() * a.z.value()) * c;
    m[2][3] = 0;

    m[3][0] = 0;
    m[3][1] = 0;
    m[3][2] = 0;
    m[3][3] = 1;

    Matrix4x4 mat(m);
    return Transform(mat, Transpose(mat));
}


//Transform lookAt(const Point3D &pos, const Point3D &look, const Vector3D &up) {
//    double m[4][4];
//    // Initialize fourth column of viewing matrix
//    m[0][3] = pos.x;
//    m[1][3] = pos.y;
//    m[2][3] = pos.z;
//    m[3][3] = 1;

//    // Initialize first three columns of viewing matrix
//    Vector3D dir = normalize(look - pos);
//    if (cross(normalize(up), dir).length() == 0) {
//        Error("\"up\" vector (%f, %f, %f) and viewing direction (%f, %f, %f) "
//              "passed to LookAt are pointing in the same direction.  Using "
//              "the identity transformation.", up.x, up.y, up.z, dir.x, dir.y,
//              dir.z);
//        return Transform();
//    }
//    Vector3D left = normalize(cross(normalize(up), dir));
//    Vector3D newUp = cross(dir, left);
//    m[0][0] = left.x;
//    m[1][0] = left.y;
//    m[2][0] = left.z;
//    m[3][0] = 0.;
//    m[0][1] = newUp.x;
//    m[1][1] = newUp.y;
//    m[2][1] = newUp.z;
//    m[3][1] = 0.;
//    m[0][2] = dir.x;
//    m[1][2] = dir.y;
//    m[2][2] = dir.z;
//    m[3][2] = 0.;
//    Matrix4x4 camToWorld(m);
//    return Transform(Inverse(camToWorld), camToWorld);
//}


BBox Transform::operator()(const BBox &b) const {
    const Transform &M = *this;
    BBox ret(        M(Point3D(b.pMin.x, b.pMin.y, b.pMin.z)));
    ret = makeUnion(ret, M(Point3D(b.pMax.x, b.pMin.y, b.pMin.z)));
    ret = makeUnion(ret, M(Point3D(b.pMin.x, b.pMax.y, b.pMin.z)));
    ret = makeUnion(ret, M(Point3D(b.pMin.x, b.pMin.y, b.pMax.z)));
    ret = makeUnion(ret, M(Point3D(b.pMin.x, b.pMax.y, b.pMax.z)));
    ret = makeUnion(ret, M(Point3D(b.pMax.x, b.pMax.y, b.pMin.z)));
    ret = makeUnion(ret, M(Point3D(b.pMax.x, b.pMin.y, b.pMax.z)));
    ret = makeUnion(ret, M(Point3D(b.pMax.x, b.pMax.y, b.pMax.z)));
    return ret;
}


Transform Transform::operator*(const Transform &t2) const {
    Matrix4x4 m1 = Matrix4x4::Mul(m, t2.m);
    Matrix4x4 m2 = Matrix4x4::Mul(t2.mInv, mInv);
    return Transform(m1, m2);
}


bool Transform::swapsHandedness() const {
    double det = ((m.m[0][0] *
                  (m.m[1][1] * m.m[2][2] -
                   m.m[1][2] * m.m[2][1])) -
                 (m.m[0][1] *
                  (m.m[1][0] * m.m[2][2] -
                   m.m[1][2] * m.m[2][0])) +
                 (m.m[0][2] *
                  (m.m[1][0] * m.m[2][1] -
                   m.m[1][1] * m.m[2][0])));
    return det < 0.0;
}


Transform orthographic(double znear, double zfar) {
    return scale(1.f, 1.f, 1.f / (zfar-znear)) *
           translate(Vector3D(0.0_um, 0.0_um, -znear * 1.0_um));
}


Transform perspective(double fov, double n, double f) {
    // Perform projective divide
    Matrix4x4 persp = Matrix4x4(1, 0,           0,              0,
                                0, 1,           0,              0,
                                0, 0, f / (f - n), -f*n / (f - n),
                                0, 0,           1,              0);

    // Scale to canonical viewing volume
    double invTanAng = 1.f / tanf(Radians(fov) / 2.f);
    return scale(invTanAng, invTanAng, 1) * Transform(persp);
}



// AnimatedTransform Method Definitions
//void AnimatedTransform::decompose(const Matrix4x4 &m, Vector3D *T,
//                                  Quaternion *Rquat, Matrix4x4 *S) {
//    // Extract translation _T_ from transformation matrix
//    T->x = m.m[0][3];
//    T->y = m.m[1][3];
//    T->z = m.m[2][3];

//    // Compute new transformation matrix _M_ without translation
//    Matrix4x4 M = m;
//    for (int i = 0; i < 3; ++i)
//        M.m[i][3] = M.m[3][i] = 0.0;
//    M.m[3][3] = 1.f;

//    // Extract rotation _R_ from transformation matrix
//    double norm;
//    int count = 0;
//    Matrix4x4 R = M;
//    do {
//        // Compute next matrix _Rnext_ in series
//        Matrix4x4 Rnext;
//        Matrix4x4 Rit = Inverse(Transpose(R));
//        for (int i = 0; i < 4; ++i)
//            for (int j = 0; j < 4; ++j)
//                Rnext.m[i][j] = 0.5f * (R.m[i][j] + Rit.m[i][j]);

//        // Compute norm of difference between _R_ and _Rnext_
//        norm = 0.0;
//        for (int i = 0; i < 3; ++i) {
//            double n = fabsf(R.m[i][0] - Rnext.m[i][0]) +
//                      fabsf(R.m[i][1] - Rnext.m[i][1]) +
//                      fabsf(R.m[i][2] - Rnext.m[i][2]);
//            norm = max(norm, n);
//        }
//        R = Rnext;
//    } while (++count < 100 && norm > .0001f);
//    // XXX TODO FIXME deal with flip...
//    *Rquat = Quaternion(R);

//    // Compute scale _S_ using rotation and original matrix
//    *S = Matrix4x4::Mul(Inverse(R), M);
//}


//void AnimatedTransform::interpolate(double time, Transform *t) const {
//    // Handle boundary conditions for matrix interpolation
//    if (!actuallyAnimated || time <= startTime) {
//        *t = *startTransform;
//        return;
//    }
//    if (time >= endTime) {
//        *t = *endTransform;
//        return;
//    }
//    double dt = (time - startTime) / (endTime - startTime);
//    // Interpolate translation at _dt_
//    Vector3D trans = (1.f - dt) * T[0] + dt * T[1];

//    // Interpolate rotation at _dt_
////    Quaternion rotate = slerp(dt, R[0], R[1]);

//    // Interpolate scale at _dt_
//    Matrix4x4 scale;
//    for (int i = 0; i < 3; ++i)
//        for (int j = 0; j < 3; ++j)
//            scale.m[i][j] = lerp(dt, S[0].m[i][j], S[1].m[i][j]);

//    // Compute interpolated matrix as product of interpolated components
//    *t = translate(trans) * rotate.toTransform() * Transform(scale);
//}


//BBox AnimatedTransform::motionBounds(const BBox &b,
//                                     bool useInverse) const {
//    if (!actuallyAnimated)
//      return useInverse ? Inverse(*startTransform)(b) : (*startTransform)(b);
//    BBox ret;
//    const int nSteps = 128;
//    for (int i = 0; i < nSteps; ++i) {
//        Transform t;
//        double time = lerp(double(i)/double(nSteps-1), startTime, endTime);
//        interpolate(time, &t);
//        if (useInverse) t = Inverse(t);
//        ret = makeUnion(ret, t(b));
//    }
//    return ret;
//}


//void AnimatedTransform::operator()(const Ray &r, Ray *tr) const {
//    if (!actuallyAnimated || r.m_time <= startTime)
//        (*startTransform)(r, tr);
//    else if (r.m_time >= endTime)
//        (*endTransform)(r, tr);
//    else {
//        Transform t;
//        interpolate(r.m_time, &t);
//        t(r, tr);
//    }
//    tr->m_time = r.m_time;
//}


//void AnimatedTransform::operator()(const RayDifferential &r,
//    RayDifferential *tr) const {
//    if (!actuallyAnimated || r.m_time <= startTime)
//        (*startTransform)(r, tr);
//    else if (r.m_time >= endTime)
//        (*endTransform)(r, tr);
//    else {
//        Transform t;
//        interpolate(r.m_time, &t);
//        t(r, tr);
//    }
//    tr->m_time = r.m_time;
//}


//Point3D AnimatedTransform::operator()(double time, const Point3D &p) const {
//    if (!actuallyAnimated || time <= startTime)
//        return (*startTransform)(p);
//    else if (time >= endTime)
//        return (*endTransform)(p);
//    Transform t;
//    interpolate(time, &t);
//    return t(p);
//}


//Vector3D AnimatedTransform::operator()(double time, const Vector3D &v) const {
//    if (!actuallyAnimated || time <= startTime)
//        return (*startTransform)(v);
//    else if (time >= endTime)
//        return (*endTransform)(v);
//    Transform t;
//    interpolate(time, &t);
//    return t(v);
//}


//Ray AnimatedTransform::operator()(const Ray &r) const {
//    Ray ret;
//    (*this)(r, &ret);
//    return ret;
//}


