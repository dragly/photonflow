#ifndef COMMON
#define COMMON

#include "../core/error.h"
#include "../core/units.h"

#include <assert.h>
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>

using std::min;
using std::max;
using std::swap;
using std::sort;
using std::vector;
using std::string;

class LightSample;
class VisibilityTester;
class Scene;
class Renderer;
class Sample;
class RNG;
class LightSampleOffsets;
template <typename T, int logBlockSize = 2> class BlockedArray;

#define UNUSED(x) (void)(x)

#define ALLOCA(TYPE, COUNT) (TYPE *)alloca((COUNT) * sizeof(TYPE))

#ifdef M_PI
#undef M_PI
#endif
#define M_PI       3.14159265358979323846f
#define INV_PI     0.31830988618379067154f
#define INV_TWOPI  0.15915494309189533577f
#define INV_FOURPI 0.07957747154594766788f

#ifndef PBRT_L1_CACHE_LINE_SIZE
#define PBRT_L1_CACHE_LINE_SIZE 64
#endif

// Global Inline Functions
inline auto lerp(double t, auto v1, auto v2) {
    return (1.f - t) * v1 + t * v2;
}


inline auto clamp(auto val, auto low, auto high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}


inline int clamp(int val, int low, int high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}


inline int Mod(int a, int b) {
    int n = int(a/b);
    a -= n*b;
    if (a < 0) a += b;
    return a;
}


inline double Radians(double deg) {
    return ((double)M_PI/180.0) * deg;
}


inline double Degrees(double rad) {
    return (180.0/(double)M_PI) * rad;
}


inline double Log2(double x) {
    static double invLog2 = 1.f / logf(2.f);
    return logf(x) * invLog2;
}


inline int Floor2Int(double val);
inline int Log2Int(double v) {
    return Floor2Int(Log2(v));
}


inline bool IsPowerOf2(int v) {
    return v && !(v & (v - 1));
}


inline uint32_t RoundUpPow2(uint32_t v) {
    v--;
    v |= v >> 1;    v |= v >> 2;
    v |= v >> 4;    v |= v >> 8;
    v |= v >> 16;
    return v+1;
}


inline int Floor2Int(double val) {
    return (int)floorf(val);
}


inline int Round2Int(double val) {
    return Floor2Int(val + 0.5f);
}


inline int Float2Int(double val) {
    return (int)val;
}


inline int Ceil2Int(double val) {
    return (int)ceilf(val);
}

//template<typename U>
//bool isNaN(U val) {
//    return isnan(val.value());
//}

template<typename T>
bool isnan(boost::units::quantity<T> val) {
    return isnan(val.value());
}

//template<typename T>
//auto sqrt(boost::units::quantity<T> val) {
//    return sqrt(val.value());
//}

template<typename T>
auto min(boost::units::quantity<T> val) {
    return min(val.value());
}

template<typename T>
auto max(boost::units::quantity<T> val) {
    return max(val.value());
}

#ifdef NDEBUG
#define photonFlowAssert(expr) ((void)0)
#else
#define photonFlowAssert(expr) \
    ((expr) ? (void)0 : \
        Severe("Assertion \"%s\" failed in %s, line %d", \
               #expr, __FILE__, __LINE__))
#endif // NDEBUG
inline bool Quadratic(double A, double B, double C, double *t0, double *t1) {
    // Find quadratic discriminant
    double discrim = B * B - 4.f * A * C;
    if (discrim < 0.) return false;
    double rootDiscrim = sqrtf(discrim);

    // Compute quadratic _t_ values
    double q;
    if (B < 0) q = -.5f * (B - rootDiscrim);
    else       q = -.5f * (B + rootDiscrim);
    *t0 = q / A;
    *t1 = C / q;
    if (*t0 > *t1) swap(*t0, *t1);
    return true;
}


#endif // COMMON

