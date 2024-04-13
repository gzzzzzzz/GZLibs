#pragma once
#include <string>
#include <Runtime\Math\Vector3.h>
#include <Runtime\Math\Matrix4x4.h>
#include <Runtime\Math\Vector4.h>
#include <Runtime\Math\Quaternion.h>
#include <cassert>
#include <sstream>
#include <iostream>
using namespace std;

class COMP
{
public:
    double rmz;
    double imz;
    string ToString(string fmt)
    {
        return "[" + to_string(rmz) + "," + to_string(imz) + "i]";
    }
    friend COMP operator + (double y, COMP x)
    {
        COMP ret;
        ret.rmz = x.rmz + y;
        ret.imz = x.imz;
        return ret;
    }
    friend COMP operator + (COMP x, double y)
    {
        COMP ret;
        ret.rmz = x.rmz + y;
        ret.imz = x.imz;
        return ret;
    }
    friend COMP operator + (COMP x, COMP y)
    {
        COMP ret;
        ret.rmz = x.rmz + y.rmz;
        ret.imz = x.imz + y.imz;
        return ret;
    }
    friend COMP operator * (COMP x, COMP y)
    {
        COMP ret;
        ret.rmz = x.rmz * y.rmz - x.imz * y.imz;
        ret.imz = x.rmz * y.imz + x.imz * y.rmz;
        return ret;
    }
    friend COMP operator * (double b, COMP x)
    {
        COMP ret;
        ret.rmz = x.rmz * b;
        ret.imz = x.imz * b;
        return ret;
    }
    friend COMP operator * (COMP x, double b)
    {
        COMP ret;
        ret.rmz = x.rmz * b;
        ret.imz = x.imz * b;
        return ret;
    }
};

class Vec
{
public:
    double* m = nullptr;
    int nRow = 0;

    ~Vec()
    {
        delete m;
    }
    Vec()
    {
    }
    Vec(int _size)
    {
        Resize(_size);
    }
    Vec(int _size, double val)
    {
        Resize(_size);
        for (int i = 0; i < _size; i++)
            m[i] = val;
    }
    Vec(const Vec& b)
    {
        Resize(b.Size());
        for (int i = 0; i < nRow; i++)
            m[i] = b.m[i];
    }
    Vec(Vec&& b) noexcept
    {
        m = b.m;
        b.m = nullptr;
        nRow = b.nRow;
    }
    Vec& operator = (const Vec& b)
    {
        Resize(b.Size());
        for (int i = 0; i < nRow; i++)
            m[i] = b.m[i];
        return *this;
    }
    Vec& operator = (Vec&& b) noexcept
    {
        m = b.m;
        b.m = nullptr;
        nRow = b.nRow;
        return *this;
    }
    Vec(const Vector3f& b)
    {
        Resize(3);
        m[0] = (double)b.x;
        m[1] = (double)b.y;
        m[2] = (double)b.z;
    }
    int Size() const { return nRow; }
    void Resize(int size)
    {
        nRow = size;
        if (m != nullptr)
            delete m;
        m = new double[size];
        memset(m, 0, sizeof(double) * size);
    }
    friend Vec operator +(const Vec& A, const Vec& B)
    {
        assert(B.Size() == A.Size());
        Vec ret(A.Size());
        for (int i = 0; i < A.Size(); i++)
            ret.m[i] = A.m[i] + B.m[i];
        return ret;
    }
    friend Vec operator -(const Vec& A, const Vec& B)
    {
        assert(B.Size() == A.Size());
        Vec ret(A.Size());
        for (int i = 0; i < A.Size(); i++)
            ret.m[i] = A.m[i] - B.m[i];
        return ret;
    }
    friend Vec operator * (const Vec& A, double f)
    {
        Vec ret(A);
        for (int i = 0; i < ret.nRow; i++)
            ret.m[i] = ret.m[i] * f;
        return ret;
    }
    friend Vec operator / (const Vec& A, double f)
    {
        return A * (1.0 / f);
    }
    double Dot(Vec A)
    {
        assert(Size() == A.Size());
        double sum = 0;
        for (int i = 0; i < nRow; i++)
            sum += m[i] * A.m[i];
        return sum;
    }
    string ToString(string fmt = "f3")
    {
        ostringstream str;
        str << "[";
        for (int i = 0; i < (nRow - 1); i++)
            str << to_string(m[i]) << ",";
        str << to_string(m[nRow - 1]) << "]";
        return str.str();
    }
    Vector3f ToUnityVec3()
    {
        Vector3f uvec3;
        uvec3.x = (float)m[0];
        uvec3.y = (float)m[1];
        uvec3.z = (float)m[2];
        return uvec3;
    }
    void SetByTCrossR(double tx, double ty, double tz, Vec R)
    {
        assert(R.Size() == 9);
        Resize(9);
        m[0] = -tz * R.m[3] + ty * R.m[6];
        m[1] = -tz * R.m[4] + ty * R.m[7];
        m[2] = -tz * R.m[5] + ty * R.m[8];

        m[3] = tz * R.m[0] - tx * R.m[6];
        m[4] = tz * R.m[1] - tx * R.m[7];
        m[5] = tz * R.m[2] - tx * R.m[8];

        m[6] = -ty * R.m[0] + tx * R.m[3];
        m[7] = -ty * R.m[1] + tx * R.m[4];
        m[8] = -ty * R.m[2] + tx * R.m[5];
    }
    void SetByMatrix44(const Matrix4x4f& mat)
    {
        Resize(9);
        Vector4f v;
        v = mat.GetRow(0); m[0] = v.x; m[1] = v.y; m[2] = v.z;//GetRow, GetColumn
        v = mat.GetRow(1); m[3] = v.x; m[4] = v.y; m[5] = v.z;
        v = mat.GetRow(2); m[6] = v.x; m[7] = v.y; m[8] = v.z;
    }
    double Length()
    {
        double ret = 0;
        for (int i = 0; i < nRow; i++)
            ret += m[i] * m[i];
        return sqrt(ret);
    }
    
};
class MyMat
{
    //=======================
    using DP = double*;
public:
    double** m = nullptr;
private:
    int nRow = 0;
    int nColumn = 0;
    long nTotal = 0;
public:
    MyMat()
    {
    }
    ~MyMat()
    {
        if (m != nullptr)
        {
            for (int i = 0; i < nRow; i++)
            {
                delete m[i];
            }
            delete m;
        }
    }
    void Init(int rows, int columns)
    {
        if (m != nullptr)
        {
            for (int i = 0; i < nRow; i++)
            {
                delete m[i];
            }
            delete m;
        }
        nRow = rows;
        nColumn = columns;
        nTotal = nRow * nColumn;

        m = new double* [nRow];
        for (int i = 0; i < nRow; i++)
        {
            m[i] = new double[nColumn];
        }
    }
    MyMat(int rows, int columns, double x = 0)
    {
        Init(rows, columns);
        SetByScalar(x);
    }
    MyMat(const MyMat& A)
    {
        Init(A.nRow, A.nColumn);
        CopyFromArray(A.m);
    }
    MyMat(MyMat&& A) noexcept
    {
        m = A.m;
        A.m = nullptr;
        nRow = A.nRow;
        nColumn = A.nColumn;
        nTotal = A.nTotal;
    }
    MyMat& operator = (const MyMat& A)
    {
        Init(A.nRow, A.nColumn);
        CopyFromArray(A.m);

        return *this;
    }
    MyMat& operator = (MyMat&& A) noexcept
    {
        m = A.m;
        A.m = nullptr;
        nRow = A.nRow;
        nColumn = A.nColumn;
        nTotal = A.nTotal;
        return *this;
    }
    void CopyFromArray(const vector<double>& arr)
    {
        assert(arr.size() == nTotal);
        int id = 0;
        for (int i = 0; i < nRow; ++i)
            for (int j = 0; j < nColumn; ++j)
                m[i][j] = arr[id++];
    }
    void CopyFromArray(const double* arr, int count)
    {
        assert(count == nTotal);
        int id = 0;
        for (int i = 0; i < nRow; ++i)
            for (int j = 0; j < nColumn; ++j)
                m[i][j] = arr[id++];
    }
    
    void CopyFromArray(const DP* arr)
    {
        for (int i = 0; i < nRow; ++i)
            for (int j = 0; j < nColumn; ++j)
                m[i][j] = arr[i][j];
    }
    void SetByScalar(double x)
    {
        for (int i = 0; i < nRow; ++i)
            for (int j = 0; j < nColumn; ++j)
                m[i][j] = x;
    }
    static MyMat SetByTCross(double tx, double ty, double tz)
    {
        MyMat A(3, 3);
        A.m[0][0] = 0; A.m[0][1] = -tz; A.m[0][2] = ty;
        A.m[1][0] = tz; A.m[1][1] = 0; A.m[1][2] = -tx;
        A.m[2][0] = -ty; A.m[2][1] = tx; A.m[2][2] = 0;
        return A;
    }
    

    long Size() const { return nTotal; }
    int RowSize() const { return nRow; }
    int ColSize() const { return nColumn; }
    MyMat Resize(int rows, int columns)
    {
        if (rows == nRow && columns == nColumn)
            return *this;

        Init(rows, columns);
        return *this;
    }
    Vec GetRow(int row)
    {
#if BOUNDS_CHECK
        assert(row >= 0);
        assert(row < nRow);
#endif
        Vec tmp(nColumn);
        for (int j = 0; j < nColumn; ++j)
            tmp.m[j] = m[row][j];
        return tmp;
    }
    Vec GetColumn(int column)
    {
#if BOUNDS_CHECK
        assert(column >= 0);
        assert(column < nColumn);
#endif
        Vec tmp(nRow);
        for (int i = 0; i < nRow; ++i)
            tmp.m[i] = m[i][column];
        return tmp;
    }
    void SetRow(const Vec& v, int row)
    {
#if BOUNDS_CHECK
        assert(row >= 0);
        assert(row < nRow);
        assert(v.dim() == nColumn);
#endif

        for (int j = 0; j < nColumn; ++j)
            m[row][j] = v.m[j];
    }
    void SetColumn(const Vec& v, int column)
    {
#if  BOUNDS_CHECK
        assert(column >= 0);
        assert(column < nColumn);
        assert(v.dim() == nRow);
#endif
        assert(true);
        for (int i = 0; i < nRow; ++i)
            m[i][column] = v.m[i];
    }
    friend MyMat operator *(const MyMat& A, double x)
    {
        MyMat ret(A);
        for (int i = 0; i < A.nRow; ++i)
        {
            for (int j = 0; j < A.nColumn; ++j)
                ret.m[i][j] *= x;
        }
        return ret;
    }
    friend MyMat operator /(const MyMat& A, double x)
    {
        MyMat ret(A);
        x = 1.0 / x;
        for (int i = 0; i < A.nRow; ++i)
        {
            for (int j = 0; j < A.nColumn; ++j)
                ret.m[i][j] *= x;
        }
        return ret;
    }
    friend MyMat operator +(const MyMat& A, const MyMat& B)
    {
        assert(A.RowSize() == B.RowSize());
        assert(A.ColSize() == B.ColSize());
        MyMat C(A.RowSize(), A.ColSize());
        for (int i = 0; i < A.RowSize(); ++i)
        {
            for (int j = 0; j < A.ColSize(); ++j)
                C.m[i][j] = A.m[i][j] + B.m[i][j];
        }
        return C;
    }
    friend MyMat operator -(const MyMat& A, const MyMat& B)
    {
        assert(A.RowSize() == B.RowSize());
        assert(A.ColSize() == B.ColSize());
        MyMat C(A.RowSize(), A.ColSize());
        for (int i = 0; i < A.RowSize(); ++i)
        {
            for (int j = 0; j < A.ColSize(); ++j)
                C.m[i][j] = A.m[i][j] - B.m[i][j];
        }
        return C;
    }
    friend MyMat operator *(const MyMat& A, const MyMat& B)
    {
        int M = A.RowSize();
        int N = B.ColSize();
        int K = A.ColSize();
        assert(B.RowSize() == K);
        MyMat C(M, N);
        double sum;
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; ++j)
            {
                sum = 0;
                for (int k = 0; k < K; ++k)
                    sum += A.m[i][k] * B.m[k][j];
                C.m[i][j] = sum;
            }
        return C;
    }
    friend Vec operator *(const MyMat& A, const Vec& b)
    {
        int M = A.RowSize();
        int N = A.ColSize();
        assert(b.Size() == N);
        Vec c(M);
        double sum;
        for (int i = 0; i < M; i++)
        {
            sum = 0;
            for (int j = 0; j < N; ++j)
                sum += A.m[i][j] * b.m[j];
            c.m[i] = sum;
        }
        return c;
    }
    MyMat MulEbyE(const MyMat& A)
    {
        assert(A.RowSize() == nRow);
        assert(A.ColSize() == nColumn);
        for (int i = 0; i < nRow; ++i)
        {
            for (int j = 0; j < nColumn; ++j)
                m[i][j] *= A.m[i][j];
        }
        return *this;
    }
    MyMat T() const
    {
        int rows = nColumn;
        int clumns = nRow;
        MyMat tmp(rows, clumns);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < clumns; ++j)
                tmp.m[i][j] = m[j][i];
        return tmp;
    }
    double Det3()
    {
        assert(nRow == 3 && nColumn == 3);
        double ret = m[0][0] * m[1][1] * m[2][2] + m[0][1] * m[1][2] * m[2][0] + m[0][2] * m[1][0] * m[2][1];
        ret -= m[0][0] * m[1][2] * m[2][1] + m[0][1] * m[1][0] * m[2][2] + m[0][2] * m[1][1] * m[2][0];
        return ret;
    }
    static double DetIJ(const MyMat& mat44, int i, int j)
    {
        assert(mat44.nRow == 4 && mat44.nColumn == 4);
        int x, y, ii, jj;
        double ret = 0;
        MyMat mat(3, 3, 0);
        x = 0;

        for (ii = 0; ii < 4; ii++)
        {
            if (ii == i) continue;
            y = 0;
            for (jj = 0; jj < 4; jj++)
            {
                if (jj == j) continue;
                mat.m[x][y] = mat44.m[ii][jj];
                y++;
            }
            x++;
        }
        ret += mat.m[0][0] * (mat.m[1][1] * mat.m[2][2] - mat.m[2][1] * mat.m[1][2]);
        ret -= mat.m[0][1] * (mat.m[1][0] * mat.m[2][2] - mat.m[2][0] * mat.m[1][2]);
        ret += mat.m[0][2] * (mat.m[1][0] * mat.m[2][1] - mat.m[2][0] * mat.m[1][1]);

        return ret;
    }
    static MyMat Inv44(const MyMat& mat44)
    {
        assert(mat44.nRow == 4 && mat44.nColumn == 4);
        MyMat mInverse(4, 4, 0);
        int i, j;
        double det, detij;

        // calculate 4x4 determinant
        det = 0.0;
        for (i = 0; i < 4; i++)
        {
            //(i & 0x1)与操作代替-1……（i+j）次幂
            det += (i & 0x1) != 0 ? (-mat44.m[0][i] * DetIJ(mat44, 0, i)) : (mat44.m[0][i] * DetIJ(mat44, 0, i));
        }
        if (abs(det) < 1e-10)
        {
            cout << ("function Inv44(m) : det(m) == 0, return Mat44(4,4,0)\n");
            return mInverse;
        }
        det = 1.0 / det;

        // calculate inverse
        for (i = 0; i < 4; i++)
        {
            for (j = 0; j < 4; j++)
            {
                //取伴随矩阵的值
                detij = DetIJ(mat44, j, i);
                mInverse.m[i][j] = ((i + j) & 0x1) != 0 ? (-detij * det) : (detij * det);
            }
        }
        return mInverse;
    }
    MyMat Inv33()
    {
        assert(nRow == 3 && nColumn == 3);
        double _d = Det3();
        MyMat ret(3, 3, 0);
        if (_d == 0.0)
            return ret;
        _d = 1.0 / _d;
        double a = m[0][0];
        double b = m[0][1];
        double c = m[0][2];
        double d = m[1][0];
        double e = m[1][1];
        double f = m[1][2];
        double g = m[2][0];
        double h = m[2][1];
        double i = m[2][2];
        ret.m[0][0] = _d * (e * i - h * f);
        ret.m[0][1] = -_d * (b * i - h * c);
        ret.m[0][2] = _d * (b * f - c * e);

        ret.m[1][0] = _d * (f * g - i * d);
        ret.m[1][1] = -_d * (c * g - i * a);
        ret.m[1][2] = _d * (c * d - a * f);

        ret.m[2][0] = _d * (d * h - g * e);
        ret.m[2][1] = -_d * (a * h - g * b);
        ret.m[2][2] = _d * (a * e - b * d);
        return ret;
    }
    Vec ToQuaternionVec()
    {
        assert(nRow >= 3 && nColumn >= 3);
        Vec q(4);
        double tr = m[0][0] + m[1][1] + m[2][2];
        double qw, qx, qy, qz;
        if (tr > 0)
        {
            double S = sqrt(tr + 1.0) * 2; // S=4*qw 
            qw = 0.25 * S;
            qx = (m[2][1] - m[1][2]) / S;
            qy = (m[0][2] - m[2][0]) / S;
            qz = (m[1][0] - m[0][1]) / S;
        }
        else if ((m[0][0] > m[1][1]) & (m[0][0] > m[2][2]))
        {
            double S = sqrt(1.0 + m[0][0] - m[1][1] - m[2][2]) * 2; // S=4*qx 
            qw = (m[2][1] - m[1][2]) / S;
            qx = 0.25 * S;
            qy = (m[0][1] + m[1][0]) / S;
            qz = (m[0][2] + m[2][0]) / S;
        }
        else if (m[1][1] > m[2][2])
        {
            double S = sqrt(1.0 + m[1][1] - m[0][0] - m[2][2]) * 2; // S=4*qy
            qw = (m[0][2] - m[2][0]) / S;
            qx = (m[0][1] + m[1][0]) / S;
            qy = 0.25 * S;
            qz = (m[1][2] + m[2][1]) / S;
        }
        else
        {
            double S = sqrt(1.0 + m[2][2] - m[0][0] - m[1][1]) * 2; // S=4*qz
            qw = (m[1][0] - m[0][1]) / S;
            qx = (m[0][2] + m[2][0]) / S;
            qy = (m[1][2] + m[2][1]) / S;
            qz = 0.25 * S;
        }
        q.m[0] = qw; q.m[1] = qx; q.m[2] = qy; q.m[3] = qz;
        return q;
    }
    Quaternionf ToQuaternion()
    {
        Vec vq = ToQuaternionVec();
        return Quaternionf((float)vq.m[1], (float)vq.m[2], (float)vq.m[3], (float)vq.m[0]);
    }
    string ToString(const string& fmt = "f4")
    {
        string str = "size: " + to_string(nRow) + " by " + to_string(nColumn) + "\n";
        for (int i = 0; i < nRow; ++i)
        {
            for (int j = 0; j < nColumn; ++j)
                str += (to_string(m[i][j]) + "\t");
            str += "\n";
        }
        return str;
    }
    
};

class SVD
{
public:
    SVD()
    {

    }
    void dec(const MyMat& A)
    {
        int m = A.RowSize(),
            n = A.ColSize(),
            p = min(m, n);

        U = MyMat(m, p);
        V = MyMat(n, p);
        S = Vec(p);
        if (m >= n)
        {
            MyMat B(A);
            decomposition(B, U, S, V);
        }
        else
        {
            MyMat B = A.T();
            decomposition(B, V, S, U);
        }
    }
    static double hypot(double x, double y)
    {
        return sqrt(x * x + y * y);
    }
    static void swap(double& x, double& y)
    {
        double temp = x;
        x = y;
        y = temp;
    }
    static void decomposition(const MyMat& B, MyMat& U, Vec& S, MyMat& V)
    {
        int m = B.RowSize(),
            n = B.ColSize();

        Vec e(n);
        Vec work(m);

        // boolean
        int wantu = 1;
        int wantv = 1;

        // Reduce A to bidiagonal form, storing the diagonal elements
        // in s and the super-diagonal elements in e.
        int nct = min(m - 1, n);
        int nrt = max(0, n - 2);
        int i = 0,
            j = 0,
            k = 0;

        for (k = 0; k < max(nct, nrt); ++k)
        {
            if (k < nct)
            {
                // Compute the transformation for the k-th column and
                // place the k-th diagonal in s[k].
                // Compute 2-norm of k-th column without under/overflow.
                S.m[k] = 0;
                for (i = k; i < m; ++i)
                    S.m[k] = hypot(S.m[k], B.m[i][k]);

                if (S.m[k] != 0)
                {
                    if (B.m[k][k] < 0)
                        S.m[k] = -S.m[k];

                    for (i = k; i < m; ++i)
                        B.m[i][k] /= S.m[k];
                    B.m[k][k] += 1;
                }
                S.m[k] = -S.m[k];
            }

            for (j = k + 1; j < n; ++j)
            {
                if ((k < nct) && (S.m[k] != 0))
                {
                    // apply the transformation
                    double t = 0;
                    for (i = k; i < m; ++i)
                        t += B.m[i][k] * B.m[i][j];

                    t = -t / B.m[k][k];
                    for (i = k; i < m; ++i)
                        B.m[i][j] += t * B.m[i][k];
                }
                e.m[j] = B.m[k][j];
            }

            // Place the transformation in U for subsequent back
            // multiplication.
            if ((wantu & ((k < nct) ? 1 : 0)) != 0)
                for (i = k; i < m; ++i)
                    U.m[i][k] = B.m[i][k];

            if (k < nrt)
            {
                // Compute the k-th row transformation and place the
                // k-th super-diagonal in e[k].
                // Compute 2-norm without under/overflow.
                e.m[k] = 0;
                for (i = k + 1; i < n; ++i)
                    e.m[k] = hypot(e.m[k], e.m[i]);

                if (e.m[k] != 0)
                {
                    if (e.m[k + 1] < 0)
                        e.m[k] = -e.m[k];

                    for (i = k + 1; i < n; ++i)
                        e.m[i] /= e.m[k];
                    e.m[k + 1] += 1;
                }
                e.m[k] = -e.m[k];

                if ((k + 1 < m) && (e.m[k] != 0))
                {
                    // apply the transformation
                    for (i = k + 1; i < m; ++i)
                        work.m[i] = 0;

                    for (j = k + 1; j < n; ++j)
                        for (i = k + 1; i < m; ++i)
                            work.m[i] += e.m[j] * B.m[i][j];

                    for (j = k + 1; j < n; ++j)
                    {
                        double t = -e.m[j] / e.m[k + 1];
                        for (i = k + 1; i < m; ++i)
                            B.m[i][j] += t * work.m[i];
                    }
                }

                // Place the transformation in V for subsequent
                // back multiplication.
                if (wantv != 0)
                    for (i = k + 1; i < n; ++i)
                        V.m[i][k] = e.m[i];
            }
        }

        // Set up the final bidiagonal matrix or order p.
        int p = n;

        if (nct < n)
            S.m[nct] = B.m[nct][nct];
        if (m < p)
            S.m[p - 1] = 0;

        if (nrt + 1 < p)
            e.m[nrt] = B.m[nrt][p - 1];
        e.m[p - 1] = 0;

        // if required, generate U
        if (wantu != 0)
        {
            for (j = nct; j < n; ++j)
            {
                for (i = 0; i < m; ++i)
                    U.m[i][j] = 0;
                U.m[j][j] = 1;
            }

            for (k = nct - 1; k >= 0; --k)
                if (S.m[k] != 0)
                {
                    for (j = k + 1; j < n; ++j)
                    {
                        double t = 0;
                        for (i = k; i < m; ++i)
                            t += U.m[i][k] * U.m[i][j];
                        t = -t / U.m[k][k];

                        for (i = k; i < m; ++i)
                            U.m[i][j] += t * U.m[i][k];
                    }

                    for (i = k; i < m; ++i)
                        U.m[i][k] = -U.m[i][k];
                    U.m[k][k] = 1 + U.m[k][k];

                    for (i = 0; i < k - 1; ++i)
                        U.m[i][k] = 0;
                }
                else
                {
                    for (i = 0; i < m; ++i)
                        U.m[i][k] = 0;
                    U.m[k][k] = 1;
                }
        }

        // if required, generate V
        if (wantv != 0)
            for (k = n - 1; k >= 0; --k)
            {
                if ((k < nrt) && (e.m[k] != 0))
                    for (j = k + 1; j < n; ++j)
                    {
                        double t = 0;
                        for (i = k + 1; i < n; ++i)
                            t += V.m[i][k] * V.m[i][j];
                        t = -t / V.m[k + 1][k];

                        for (i = k + 1; i < n; ++i)
                            V.m[i][j] += t * V.m[i][k];
                    }

                for (i = 0; i < n; ++i)
                    V.m[i][k] = 0;
                V.m[k][k] = 1;
            }

        // main iteration loop for the singular values
        int pp = p - 1;
        int iter = 0;
        double eps = pow(2.0, -52.0);

        while (p > 0)
        {
            int kk = 0;
            int kase = 0;

            // Here is where a test for too many iterations would go.
            // This section of the program inspects for negligible
            // elements in the s and e arrays. On completion the
            // variables kase and kk are set as follows.
            // kase = 1     if s(p) and e[kk-1] are negligible and kk<p
            // kase = 2     if s(kk) is negligible and kk<p
            // kase = 3     if e[kk-1] is negligible, kk<p, and
            //				s(kk), ..., s(p) are not negligible
            // kase = 4     if e(p-1) is negligible (convergence).
            for (kk = p - 2; kk >= -1; --kk)
            {
                if (kk == -1)
                    break;

                if (abs(e.m[kk]) <= eps * (abs(S.m[kk]) + abs(S.m[kk + 1])))
                {
                    e.m[kk] = 0;
                    break;
                }
            }

            if (kk == p - 2)
                kase = 4;
            else
            {
                int ks;
                for (ks = p - 1; ks >= kk; --ks)
                {
                    if (ks == kk)
                        break;

                    double t = ((ks != p) ? abs(e.m[ks]) : 0) +
                        ((ks != kk + 1) ? abs(e.m[ks - 1]) : 0);

                    if (abs(S.m[ks]) <= eps * t)
                    {
                        S.m[ks] = 0;
                        break;
                    }
                }

                if (ks == kk)
                    kase = 3;
                else if (ks == p - 1)
                    kase = 1;
                else
                {
                    kase = 2;
                    kk = ks;
                }
            }
            kk++;

            // Perform the task indicated by kase.
            switch (kase)
            {
                // deflate negligible s(p)
            case 1:
            {
                double f = e.m[p - 2];
                e.m[p - 2] = 0;

                for (j = p - 2; j >= kk; --j)
                {
                    double t = hypot(S.m[j], f);
                    double cs = S.m[j] / t;
                    double sn = f / t;
                    S.m[j] = t;

                    if (j != kk)
                    {
                        f = -sn * e.m[j - 1];
                        e.m[j - 1] = cs * e.m[j - 1];
                    }

                    if (wantv != 0)
                        for (i = 0; i < n; ++i)
                        {
                            t = cs * V.m[i][j] + sn * V.m[i][p - 1];
                            V.m[i][p - 1] = -sn * V.m[i][j] + cs * V.m[i][p - 1];
                            V.m[i][j] = t;
                        }
                }
            }
            break;

            // split at negligible s(kk)
            case 2:
            {
                double f = e.m[kk - 1];
                e.m[kk - 1] = 0;

                for (j = kk; j < p; ++j)
                {
                    double t = hypot(S.m[j], f);
                    double cs = S.m[j] / t;
                    double sn = f / t;
                    S.m[j] = t;
                    f = -sn * e.m[j];
                    e.m[j] = cs * e.m[j];

                    if (wantu != 0)
                        for (i = 0; i < m; ++i)
                        {
                            t = cs * U.m[i][j] + sn * U.m[i][kk - 1];
                            U.m[i][kk - 1] = -sn * U.m[i][j] + cs * U.m[i][kk - 1];
                            U.m[i][j] = t;
                        }
                }
            }
            break;

            // perform one qr step
            case 3:
            {
                // calculate the shift
                double scale = max(max(max(max(
                    abs(S.m[p - 1]), abs(S.m[p - 2])), abs(e.m[p - 2])),
                    abs(S.m[kk])), abs(e.m[kk]));
                double sp = S.m[p - 1] / scale;
                double spm1 = S.m[p - 2] / scale;
                double epm1 = e.m[p - 2] / scale;
                double sk = S.m[kk] / scale;
                double ek = e.m[kk] / scale;
                double b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
                double c = (sp * epm1) * (sp * epm1);
                double shift = 0;

                if ((b != 0) || (c != 0))
                {
                    shift = sqrt(b * b + c);
                    if (b < 0)
                        shift = -shift;
                    shift = c / (b + shift);
                }
                double f = (sk + sp) * (sk - sp) + shift;
                double g = sk * ek;

                // chase zeros
                for (j = kk; j < p - 1; ++j)
                {
                    double t = hypot(f, g);
                    double cs = f / t;
                    double sn = g / t;
                    if (j != kk)
                        e.m[j - 1] = t;

                    f = cs * S.m[j] + sn * e.m[j];
                    e.m[j] = cs * e.m[j] - sn * S.m[j];
                    g = sn * S.m[j + 1];
                    S.m[j + 1] = cs * S.m[j + 1];

                    if (wantv != 0)
                        for (i = 0; i < n; ++i)
                        {
                            t = cs * V.m[i][j] + sn * V.m[i][j + 1];
                            V.m[i][j + 1] = -sn * V.m[i][j] + cs * V.m[i][j + 1];
                            V.m[i][j] = t;
                        }

                    t = hypot(f, g);
                    cs = f / t;
                    sn = g / t;
                    S.m[j] = t;
                    f = cs * e.m[j] + sn * S.m[j + 1];
                    S.m[j + 1] = -sn * e.m[j] + cs * S.m[j + 1];
                    g = sn * e.m[j + 1];
                    e.m[j + 1] = cs * e.m[j + 1];

                    if (wantu != 0 && (j < m - 1))
                        for (i = 0; i < m; ++i)
                        {
                            t = cs * U.m[i][j] + sn * U.m[i][j + 1];
                            U.m[i][j + 1] = -sn * U.m[i][j] + cs * U.m[i][j + 1];
                            U.m[i][j] = t;
                        }
                }
                e.m[p - 2] = f;
                iter = iter + 1;
            }
            break;

            // convergence
            case 4:
            {
                // Make the singular values positive.
                if (S.m[kk] <= 0)
                {
                    S.m[kk] = (S.m[kk] < 0) ? -S.m[kk] : 0;
                    if (wantv != 0)
                        for (i = 0; i <= pp; ++i)
                            V.m[i][kk] = -V.m[i][kk];
                }

                // Order the singular values.
                while (kk < pp)
                {
                    if (S.m[kk] >= S.m[kk + 1])
                        break;

                    double t = S.m[kk];
                    S.m[kk] = S.m[kk + 1];
                    S.m[kk + 1] = t;

                    if (wantv != 0 && (kk < n - 1))
                        for (i = 0; i < n; ++i)
                            swap(V.m[i][kk], V.m[i][kk + 1]);

                    if (wantu != 0 && (kk < m - 1))
                        for (i = 0; i < m; ++i)
                            swap(U.m[i][kk], U.m[i][kk + 1]);
                    kk++;
                }
                iter = 0;
                p--;
            }
            break;
            }
        }
    }

    MyMat U;
    MyMat V;
    Vec S;
    MyMat getSM()
    {
        int N = S.Size();
        MyMat tmp(N, N);
        for (int i = 0; i < N; ++i)
            tmp.m[i][i] = S.m[i];
        return tmp;
    }
    /*
        public Matrix getU();
        public Matrix getV();
        public Matrix getSM();
        public Vector getSV();
        public double norm2();
        public double cond();
        public int  rank();
        // the orthogonal matrix and singular value vector


        // docomposition for matrix with rows >= columns
        private void decomposition(ref Matrix B, ref Matrix U, ref Vector S, ref Matrix V);
    //*/
};

class MathUtil
{
    struct KeyValuePair
    {
        double Key;
        int Value;
    };
    //http://blog.csdn.net/zhouxuguang236/article/details/40212143
public:
    static bool JacbiCor(double* pMatrix, int nDim, double* pdblVects, double* pdbEigenValues, double dbEps, int nJt)
    {
        for (int i = 0; i < nDim; i++)
        {
            pdblVects[i * nDim + i] = 1.0f;
            for (int j = 0; j < nDim; j++)
            {
                if (i != j)
                    pdblVects[i * nDim + j] = 0.0f;
            }
        }

        int nCount = 0;     //迭代次数  
        while (true)
        {
            //在pMatrix的非对角线上找到最大元素  
            double dbMax = pMatrix[1];
            int nRow = 0;
            int nCol = 1;
            for (int i = 0; i < nDim; i++)          //行  
            {
                for (int j = 0; j < nDim; j++)      //列  
                {
                    double d = abs(pMatrix[i * nDim + j]);

                    if ((i != j) && (d > dbMax))
                    {
                        dbMax = d;
                        nRow = i;
                        nCol = j;
                    }
                }
            }

            if (dbMax < dbEps)     //精度符合要求   
                break;

            if (nCount > nJt)       //迭代次数超过限制  
                break;

            nCount++;

            double dbApp = pMatrix[nRow * nDim + nRow];
            double dbApq = pMatrix[nRow * nDim + nCol];
            double dbAqq = pMatrix[nCol * nDim + nCol];

            //计算旋转角度  
            double dbAngle = 0.5 * atan2(-2 * dbApq, dbAqq - dbApp);
            double dbSinTheta = sin(dbAngle);
            double dbCosTheta = cos(dbAngle);
            double dbSin2Theta = sin(2 * dbAngle);
            double dbCos2Theta = cos(2 * dbAngle);

            pMatrix[nRow * nDim + nRow] = dbApp * dbCosTheta * dbCosTheta + dbAqq * dbSinTheta * dbSinTheta + 2 * dbApq * dbCosTheta * dbSinTheta;
            pMatrix[nCol * nDim + nCol] = dbApp * dbSinTheta * dbSinTheta + dbAqq * dbCosTheta * dbCosTheta - 2 * dbApq * dbCosTheta * dbSinTheta;
            pMatrix[nRow * nDim + nCol] = 0.5 * (dbAqq - dbApp) * dbSin2Theta + dbApq * dbCos2Theta;
            pMatrix[nCol * nDim + nRow] = pMatrix[nRow * nDim + nCol];

            for (int i = 0; i < nDim; i++)
            {
                if ((i != nCol) && (i != nRow))
                {
                    int u = i * nDim + nRow;  //p    
                    int w = i * nDim + nCol;  //q  
                    dbMax = pMatrix[u];
                    pMatrix[u] = pMatrix[w] * dbSinTheta + dbMax * dbCosTheta;
                    pMatrix[w] = pMatrix[w] * dbCosTheta - dbMax * dbSinTheta;
                }
            }

            for (int j = 0; j < nDim; j++)
            {
                if ((j != nCol) && (j != nRow))
                {
                    int u = nRow * nDim + j;  //p  
                    int w = nCol * nDim + j;  //q  
                    dbMax = pMatrix[u];
                    pMatrix[u] = pMatrix[w] * dbSinTheta + dbMax * dbCosTheta;
                    pMatrix[w] = pMatrix[w] * dbCosTheta - dbMax * dbSinTheta;
                }
            }

            //计算特征向量  
            for (int i = 0; i < nDim; i++)
            {
                int u = i * nDim + nRow;      //p     
                int w = i * nDim + nCol;      //q  
                dbMax = pdblVects[u];
                pdblVects[u] = pdblVects[w] * dbSinTheta + dbMax * dbCosTheta;
                pdblVects[w] = pdblVects[w] * dbCosTheta - dbMax * dbSinTheta;
            }

        }

        //对特征值进行排序以及重新排列特征向量,特征值即pMatrix主对角线上的元素  
        vector<KeyValuePair> eigenList;
        for (int i = 0; i < nDim; i++)
        {
            pdbEigenValues[i] = pMatrix[i * nDim + i];
            eigenList.push_back(KeyValuePair{ pdbEigenValues[i], i });
        }
        //eigenList.Sort((p1, p2) = > p1.Key.CompareTo(p2.Key));
        sort(eigenList.begin(), eigenList.end(), [](const KeyValuePair& p1, const KeyValuePair& p2)
        {
            return p1.Key < p2.Key;
        });

        shared_ptr<double[]> pdbTmpVec(new double[nDim * nDim]);
        for (int i = 0; i < nDim; i++)
        {
            for (int j = 0; j < nDim; j++)
                pdbTmpVec[j * nDim + i] = pdblVects[j * nDim + eigenList[i].Value];

            //特征值重新排列  
            pdbEigenValues[i] = eigenList[i].Key;
        }

        //设定正负号  
        for (int i = 0; i < nDim; i++)
        {
            double dSumVec = 0;
            for (int j = 0; j < nDim; j++)
                dSumVec += pdbTmpVec[j * nDim + i];
            if (dSumVec < 0)
            {
                for (int j = 0; j < nDim; j++)
                    pdbTmpVec[j * nDim + i] *= -1;
            }
        }

        for (int i = nDim * nDim - 1; i >= 0; i--)
            pdblVects[i] = pdbTmpVec[i];

        return true;
    }

    static Vector3f MinDistPos2NRay(const Vector3f* _pos, const Vector3f* _d_input, float* _w, int size)
    {
        if (size < 2)
            return Vector3f::zero;

        float* self_w = nullptr;
        if (_w == nullptr)
        {
            self_w = new float[size];
            _w = self_w;
            for (int i = 0; i < size; i++)
                _w[i] = (1.0f / size);
        }
        shared_ptr<Vector3f[]> _d(new Vector3f[size]);
        for (int i = 0; i < size; i++)
            _d[i] = Normalize(_d_input[i]);

        Vector3f ret = Vector3f::zero;
        MyMat A(3, 3, 0);
        Vec B(3, 0);
        for (int i = 0; i < size; i++)
        {
            Vector3f p = _pos[i];
            Vector3f d = _d[i];
            double wi = _w[i];
            double d01 = (double)(d.x * d.y);
            double d02 = (double)(d.x * d.z);
            double d12 = (double)(d.y * d.z);
            A.m[0][0] += wi * (double)(d.y * d.y + d.z * d.z);
            A.m[0][1] += wi * -d01;
            A.m[0][2] += wi * -d02;
            A.m[1][0] += wi * -d01;
            A.m[1][1] += wi * (double)(d.x * d.x + d.z * d.z);
            A.m[1][2] += wi * -d12;
            A.m[2][0] += wi * -d02;
            A.m[2][1] += wi * -d12;
            A.m[2][2] += wi * (double)(d.x * d.x + d.y * d.y);

            B.m[0] += wi * (d.z * (d.z * p.x - d.x * p.z) - d.y * (d.x * p.y - d.y * p.x));
            B.m[1] += wi * (d.x * (d.x * p.y - d.y * p.x) - d.z * (d.y * p.z - d.z * p.y));
            B.m[2] += wi * (d.y * (d.y * p.z - d.z * p.y) - d.x * (d.z * p.x - d.x * p.z));
        }
        //Debug.Log("A:(det="+A.Det3()+")" + A.ToString("F3"));
        //Debug.Log("B:" + B.ToString("F3"));
        MyMat Ainv = A.Inv33();
        Vec ans = Ainv * B;
        //Debug.Log("A*Ainv:" + (A*Ainv).ToString("F2"));
        ret.x = (float)ans.m[0];
        ret.y = (float)ans.m[1];
        ret.z = (float)ans.m[2];

        if (self_w != nullptr)
            delete self_w;
        return ret;
    }
};