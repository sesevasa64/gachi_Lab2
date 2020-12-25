from math import *
from linalg import *


class Mat4d:
    @staticmethod
    def ones() -> Matrix:
        return Matrix.ones(4)

    @staticmethod
    def saf(p: Vector):
        a, _e, i = p[0], p[1], p[2]
        return Mat([[a,  0, 0, 0],
                    [0, _e, 0, 0],
                    [0,  0, i, 0],
                    [0,  0, 0, 1]])

    @staticmethod
    def distortion(b, c, d, f, g, h):
        return Mat([[1, b, c, 0],
                    [d, 1, f, 0],
                    [g, h, 1, 0],
                    [0, 0, 0, 1]])

    @staticmethod
    def parallel(p: Vector) -> Matrix:
        k, _l, m = p[0], p[1], p[2]
        return Mat([[1,  0, 0, 0],
                    [0,  1, 0, 0],
                    [0,  0, 1, 0],
                    [k, _l, m, 1]])

    @staticmethod
    def scale(s) -> Matrix:
        res = Mat4d.ones()
        res[3][3] = s
        return res

    @staticmethod
    def rotate_0x(a, vec: Vector = Vec([0, 0, 0])):
        y, z = vec[1], vec[2]
        _l = y * (1-cos(a)) + z * sin(a)
        m = z * (1-cos(a)) - y * sin(a)
        return Mat([[1,       0,      0, 0],
                    [0,  cos(a), sin(a), 0],
                    [0, -sin(a), cos(a), 0],
                    [0,      _l,      m, 1]])

    @staticmethod
    def rotate_0y(b, vec: Vector = Vec([0, 0, 0])):
        x, z = vec[0], vec[2]
        k = x * (1-cos(b)) - z * sin(b)
        m = z * (1-cos(b)) + x * sin(b)
        return Mat([[cos(b), 0, -sin(b), 0],
                    [     0, 1,       0, 0],
                    [sin(b), 0,  cos(b), 0],
                    [     k, 0,       m, 1]])

    @staticmethod
    def rotate_0z(g, vec: Vector = Vec([0, 0, 0])):
        x, y = vec[0], vec[1]
        k = x * (1-cos(g)) + y * sin(g)
        _l = y * (1-cos(g)) - x * sin(g)
        return Mat([[ cos(g),  sin(g), 0, 0],
                    [-sin(g),  cos(g), 0, 0],
                    [      0,       0, 1, 0],
                    [      k,      _l, 0, 1]])

    @staticmethod
    def reflection(x=1, y=1, z=1, k=0, _l=0, m=0):
        return Mat([[x,      0,   0, 0],
                    [0,      y,   0, 0],
                    [0,      0,   z, 0],
                    [2*k, 2*_l, 2*m, 1]])

    @staticmethod
    def r_0x():
        return Mat4d.reflection(x=-1)

    @staticmethod
    def r_0y():
        return Mat4d.reflection(y=-1)

    @staticmethod
    def r_0z():
        return Mat4d.reflection(z=-1)

    @staticmethod
    def r_kx(k):
        return Mat4d.reflection(x=-1, k=k)

    @staticmethod
    def r_ky(_l):
        return Mat4d.reflection(y=-1, _l=_l)

    @staticmethod
    def r_kz(m):
        return Mat4d.reflection(z=-1, m=m)

    @staticmethod
    def scale_from_vector(p: Vector, vec: Vector):
        a, _e, i = p[0], p[1], p[2]
        x, y, z = vec[0], vec[1], vec[2]
        return Mat([[a,               0,       0, 0],
                    [0,              _e,       0, 0],
                    [0,               0,       i, 0],
                    [x*(1-a),  y*(1-_e), z*(1-i), 1]])

    @staticmethod
    def p_0x():
        return Mat4d.reflection(x=0)

    @staticmethod
    def p_0y():
        return Mat4d.reflection(y=0)

    @staticmethod
    def p_0z():
        return Mat4d.reflection(z=0)

    @staticmethod
    def p_kx(k):
        return Mat4d.reflection(x=0, k=k)

    @staticmethod
    def p_ky(_l):
        return Mat4d.reflection(y=0, _l=_l)

    @staticmethod
    def p_kz(m):
        return Mat4d.reflection(z=0, m=m)

    @staticmethod
    def ref_plane_2_points(vecA: Vector, vecB : Vector):
        vx = vecB[0] - vecA[0]
        vy = vecB[1] - vecA[1]
        vz = vecB[2] - vecA[2]
        lx = vx/sqrt(vx**2+vy**2+vz**2)
        ly = vy/sqrt(vx**2+vy**2+vz**2)
        lz = vz/sqrt(vx**2+vy**2+vz**2)
        h = lx*vecA[0]+ly*vecA[1]+lz*vecA[2]
        return Mat([[1 - 2*lx**2,               -2*lx*ly,       -2*lx*lz, 0],
                    [-2*lx*ly,                  1-2*ly**2,      -2*ly*lz, 0],
                    [-2*lx*lz,                  -2*ly*lz,       1-2*lz**2, 0],
                    [2*lx*h,                    2*ly*h,         2*lz*h, 1]])

    @staticmethod
    def ref_plane_3_points(pA: Vector, pB : Vector, pC: Vector):
        vx = (pB[1] - pA[1])*(pC[2]-pA[2])-(pB[2]-pA[2])*(pC[1]-pA[1])
        vy = (pB[2] - pA[2])*(pC[0]-pA[0])-(pB[0]-pA[0])*(pC[2]-pA[2])
        vz = (pB[0] - pA[0])*(pC[1]-pA[1])-(pB[1]-pA[1])*(pC[0]-pA[0])
        lx = vx/sqrt(vx**2+vy**2+vz**2)
        ly = vy/sqrt(vx**2+vy**2+vz**2)
        lz = vz/sqrt(vx**2+vy**2+vz**2)
        h = lx*pA[0]+ly*pA[1]+lz*pA[2]
        return Mat([[1 - 2*lx**2,               -2*lx*ly,       -2*lx*lz, 0],
                    [-2*lx*ly,                  1-2*ly**2,      -2*ly*lz, 0],
                    [-2*lx*lz,                  -2*ly*lz,       1-2*lz**2, 0],
                    [2*lx*h,                    2*ly*h,         2*lz*h, 1]])
                    
    @staticmethod
    def rotate_straight_line(pA: Vector, pB : Vector, angle):
        ABx = pB[0] - pA[0]
        ABy = pB[1] - pA[1]
        ABz = pB[2] - pA[2]
        lx = ABx/sqrt(ABx**2+ABy**2+ABz**2)
        ly = ABy/sqrt(ABx**2+ABy**2+ABz**2)
        lz = ABz/sqrt(ABx**2+ABy**2+ABz**2)
        a = lx**2+cos(angle)*(1-lx**2)
        b = lx*ly*(1-cos(angle))+lz*sin(angle)
        c = lx*lz*(1-cos(angle))-ly*sin(angle)
        d = lx*ly*(1-cos(angle))-lz*sin(angle)
        e = ly**2+cos(angle)*(1-ly**2)
        f = lz*ly*(1-cos(angle))+lx*sin(angle)
        g = lz*lx*(1-cos(angle))+ly*sin(angle)
        h = lz*ly*(1-cos(angle))-lx*sin(angle)
        i = lz**2+cos(angle)*(1-lz**2)
        k = pA[0]*(1-a)-pA[1]*d-pA[2]*g
        l = -pA[0]*b+pA[1]*(1-e)-pA[2]*h
        m = -pA[0]*c-pA[1]*f+pA[2]*(1-i)
        return Mat([[a,b,c,0],
                    [d,e,f,0],
                    [g,h,i,0],
                    [k,l,m,1]])

    @staticmethod
    def ortho(a: Vector, b: Vector, c: Vector):
        ABx = b[0] - a[0]
        ABy = b[1] - a[1]
        ABz = b[2] - a[2]
        ACx = c[0] - a[0]
        ACy = c[1] - a[1]
        ACz = c[2] - a[2]
        vx = ABy * ACz - ABz * ACy
        vy = ABz * ACx - ABx * ACz
        vz = ABx * ACy - ABy * ACx
        lx = vx / sqrt(vx ** 2 + vy ** 2 + vz ** 2)
        ly = vy / sqrt(vx ** 2 + vy ** 2 + vz ** 2)
        lz = vz / sqrt(vx ** 2 + vy ** 2 + vz ** 2)
        lam = sqrt(ly ** 2 + lz ** 2)
        k = (-lam ** 2 * a[0] + lx * (ly * a[1] + lz * a[2])) / lam
        l = (ly * a[2] - lz * a[1]) / lam
        return Mat([[lam, 0, 0, 0],
                    [-lx * ly / lam, lz / lam, 0, 0],
                    [-lx * lz / lam, -ly / lam, 0, 0],
                    [k, l, 0, 1]])

    @staticmethod
    def kosoug(a: Vector, b: Vector, c: Vector, fz = 1/2, angle = pi/6):
        ABx = b[0] - a[0]
        ABy = b[1] - a[1]
        ABz = b[2] - a[2]
        ACx = c[0] - a[0]
        ACy = c[1] - a[1]
        ACz = c[2] - a[2]
        vx = ABy * ACz - ABz * ACy
        vy = ABz * ACx - ABx * ACz
        vz = ABx * ACy - ABy * ACx
        lx = vx / sqrt(vx ** 2 + vy ** 2 + vz ** 2)
        ly = vy / sqrt(vx ** 2 + vy ** 2 + vz ** 2)
        lz = vz / sqrt(vx ** 2 + vy ** 2 + vz ** 2)
        lam = sqrt(ly ** 2 + lz ** 2)
        a1 = lam - lx * fz * cos(angle)
        b1 = -lx * fz * sin(angle)
        d1 = -ly * (lx / lam + fz * cos(angle))
        e1 = lz / lam - ly * fz * sin(angle)
        g1 = -lz * (lx / lam + fz * cos(angle))
        h1 = -ly / lam - lz * fz * sin(angle)
        l1 = -a1 * a[0] - d1 * a[1] - g1 * a[2]
        m1 = -b1 * a[0] - e1 * a[1] - h1 * a[2]
        return Mat([[a1, b1, 0, 0],
                    [d1, e1, 0, 0],
                    [g1, h1, 0, 0],
                    [l1, m1, 0, 1]])

    @staticmethod
    def kosougOXY(fz, angle):
        return Mat([[1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [-fz * cos(angle), -fz * sin(angle), 0, 0],
                    [0, 0, 0, 1]])

    @staticmethod
    def kabine():
        return Mat4d.kosougOXY(0.5, pi/6)
    
    @staticmethod
    def kavalie():
        return Mat4d.kosougOXY(1, pi/3)

    @staticmethod
    def aksoug(fx, fy, fz):
        a, b = None, None
        if fx == fy and fy != fz:
            b = asin(fz / sqrt(2 - fz ** 2))
            a = asin(fz / sqrt(2))
        if fx == fz and fz != fy:
            b = pi / 4
            a = 228
        if fy == fz and fz != fx:
            b = asin(sqrt((2 * fx ** 2 - 2)/(fx ** 2 - 2)))
            a = asin(fx / sqrt(2))
        if fx == fy == fz:
            b = pi / 4
            a = acos(sqrt(2 / 3))
        return Mat([[cos(b), sin(b) * sin(a), 0, 0],
                    [0, cos(a), 0, 0],
                    [sin(b), -cos(b) * sin(a), 0, 0],
                    [0, 0, 0, 1]])
