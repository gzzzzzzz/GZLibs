#pragma once
#include "MathClass.h"
#include <vector>

using namespace std;
class RTAns
{
public:
    RTAns(MyMat _r, Vec _t, double _resi);
    MyMat E;
    MyMat R;
    Vec T;
    double Resi = 1e10;
};

typedef shared_ptr<RTAns> RTAnsPtr;

class EPSolver
{
    //Epipolar Constraint
public: 
    EPSolver()
        : E(9), EM(3, 3), R0(3, 3), R1(3, 3), T0(3), T1(3)
    {
    }

    vector<Vec> AiArr;
    double ATA[81];
    Vec E;
    MyMat EM;
    MyMat R0;
    MyMat R1;
    Vec T0;
    Vec T1;
    vector<RTAnsPtr> candicateRTarr;

    void AddAi(float hori0, float vert0, float hori1, float vert1)
    {
        double u0 = tan(hori0);
        double v0 = tan(vert0);
        double u1 = tan(hori1);
        double v1 = tan(vert1);

        Vec v9(9);
        v9.m[0] = u1 * u0;
        v9.m[1] = u1 * v0;
        v9.m[2] = u1 * 1.0;

        v9.m[3] = v1 * u0;
        v9.m[4] = v1 * v0;
        v9.m[5] = v1 * 1.0;

        v9.m[6] = u0;
        v9.m[7] = v0;
        v9.m[8] = 1.0;
        AiArr.push_back(v9);
    }
    void AddAi(Vector3f m, Vector3f n)
    {
        Vec v9(9);
        double u = (double)m.x;
        double v = (double)m.y;
        v9.m[0] = n.x * u;
        v9.m[1] = n.x * v;
        v9.m[2] = n.x * 1.0;

        v9.m[3] = n.y * u;
        v9.m[4] = n.y * v;
        v9.m[5] = n.y * 1.0;

        v9.m[6] = u;
        v9.m[7] = v;
        v9.m[8] = 1.0;
        AiArr.push_back(v9);
    }
    void Do()
    {
        candicateRTarr.clear();
        if (AiArr.size() < 9)
        {
            //Debug.LogWarning("AiArr.Count < 9! return!");
            //return;
        }
        //CheckE(owner.realE, "real E");
        CaclATA();
        for(int i = 0; i < 9; i++)
        //for (int i = 0; i < 3; i++)
        {
            CaclE(i);
            //CheckE(E, "E");
            CaclRT(E);
        }
        //candicateRTarr.Sort((p1, p2) = > p1.Resi.CompareTo(p2.Resi));
        sort(candicateRTarr.begin(), candicateRTarr.end(), [](const RTAnsPtr& p1, const RTAnsPtr& p2)
        {
            if (p1->Resi > p2->Resi)
                return true;
            else
                return false;
        });
    }
    void CaclATA()
    {
        int n = static_cast<int>(AiArr.size());
        for (int i = 0; i < 9; i++)
            for (int j = 0; j < 9; j++)
            {
                double sum = 0.0;
                for (int k = 0; k < n; k++)
                    sum += AiArr[k].m[i] * AiArr[k].m[j];
                ATA[i * 9 + j] = sum;
            }
    }
    void CaclE(int eid = 0)
    {
        double eigenVec9[81] = { 0 };
        double eigenVal9[9] = { 0 };
        double copyATA[81] = { 0 };
        for (int i = 0; i < 81; i++)
            copyATA[i] = ATA[i];
        MathUtil::JacbiCor(copyATA, 9, eigenVec9, eigenVal9, 1e-50, 500);
        //CheckEigen(ATA, eigenVec9, eigenVal9);
        for (int i = 0; i < 9; i++)
            E.m[i] = eigenVec9[i * 9 + eid];
        double lenE = E.Length();
        E = E * (sqrt(2.0) / lenE);
    }
    bool CheckRT(const MyMat& r, Vec t)
    {
        if (AiArr.size() == 0) return false;
        float invalid = 0;
        float valid = 0;
        float reciCnt = 1.0f / AiArr.size();
        for (int i = 0; i < AiArr.size(); i++)
        {
            vector<Vec> MM = CaclPoint(r, t, i);

            if (MM[0].m[2] > 0 && MM[1].m[2] > 0)
                valid++;
            else
                invalid++;

        }
        if (valid * reciCnt >= 0.8f)
            return true;
        if (invalid * reciCnt > 0.2f)
            return false;
        return true;
    }
    void CaclRT(const Vec& vec)
    {
        auto dm = vec.m;
        SVD svd;
        EM.CopyFromArray(dm, vec.nRow);
        svd.dec(EM);
        MyMat W(3, 3);
        W.m[0][1] = -1; W.m[1][0] = 1; W.m[2][2] = 1;
        R0 = svd.U * W * svd.V.T();
        R1 = svd.U * W.T() * svd.V.T();

        if (R0.Det3() < 0) R0 = R0 * -1;
        if (R1.Det3() < 0) R1 = R1 * -1;
        T0 = svd.U.GetColumn(2);
        T1 = T0 * -1;

        if (CheckRT(R0, T0)) candicateRTarr.push_back(make_shared<RTAns>(R0, T0, GetResidual(R0, T0)));
        if (CheckRT(R0, T1)) candicateRTarr.push_back(make_shared<RTAns>(R0, T1, GetResidual(R0, T1)));
        if (CheckRT(R1, T0)) candicateRTarr.push_back(make_shared<RTAns>(R1, T0, GetResidual(R1, T0)));
        if (CheckRT(R1, T1)) candicateRTarr.push_back(make_shared<RTAns>(R1, T1, GetResidual(R1, T1)));
    }
    Vec CaclPostionByAngle(MyMat r, Vec t, float hori0, float vert0, float hori1, float vert1)
    {
        double u0 = tan(hori0);
        double v0 = tan(vert0);
        double u1 = tan(hori1);
        double v1 = tan(vert1);
        //Debug.Log(t.ToString() + ",uv0=[" + u0 + "," + v0 + ", uv1=[" + u1 + "," + v1);
        return CaclPostionByUV(r, t, u0, v0, u1, v1)[0];
    }
    vector<Vec> CaclPostionByUV(const MyMat& r, const Vec& t, double u0, double v0, double u1, double v1)
    {
        Vec m0(3);
        Vec m1(3);
        Vec M0(3);//M0 = lamda0 * m0
        Vec M1(3);//M1 = lamda1 * m1
        m0.m[0] = u0; m0.m[1] = v0; m0.m[2] = 1.0;
        m1.m[0] = u1; m1.m[1] = v1; m1.m[2] = 1.0;

        //M1 = R*M0 + T
        Vec Rm0 = r * m0;
        //(lamda*Rm0[0]+vT.x)/(lamda*Rm0[2]+vT.z) = u1/1.0
        double lamda0 = (t.m[0] - t.m[2] * u1) / (Rm0.m[2] * u1 - Rm0.m[0]);
        M0 = m0 * lamda0;
        Vec RM0 = Rm0 * lamda0;
        M1 = RM0 + t;

        vector<Vec> MM(2);
        MM[0] = M0;
        MM[1] = M1;
        return MM;
    }
    vector<Vec> CaclPoint(const MyMat& r, const Vec& t, int AiId)
    {
        Vec ai = AiArr[AiId];
        double u0 = ai.m[6];
        double v0 = ai.m[7];
        double u1 = ai.m[2];
        double v1 = ai.m[5];
        return CaclPostionByUV(r, t, u0, v0, u1, v1);
    }
    /*
    void CheckEigen(double[] ATA, double[] eigenVec9, double[] eigenVal9)
    {
        Vec eVec = new Vec(9);
        Vec f = new Vec(9);
        for (int k = 0; k < 9; k++)
        {
            double eVal = eigenVal9[k];
            for (int j = 0; j < 9; j++)
            {
                eVec.m[j] = eigenVec9[j * 9 + k];
                f.m[j] = 0.0;
            }
            //f = ATA*eVec
            for (int i = 0; i < 9; i++)
            {
                for (int j = 0; j < 9; j++)
                    f.m[i] += ATA[i * 9 + j] * eVec.m[j];
            }

            double diff = 0;
            for (int i = 0; i < 9; i++)
                diff += abs(eVal * eVec.m[i] - f.m[i]);
            Debug.Log(k + ", eigen val = " + eVal + ", |Ax-lamdax| = " + diff);
        }
    }
    */
    void CheckE(const Vec& e, string eStr = "E")
    {
        int errCnt = 0;
        double Ex[3];
        for (int i = 0; i < AiArr.size(); i++)
        {
            double* m = AiArr[i].m;
            double u0 = m[6];
            double v0 = m[7];
            double u1 = m[2];
            double v1 = m[5];
            Ex[0] = e.m[0] * u0 + e.m[1] * v0 + e.m[2];
            Ex[1] = e.m[3] * u0 + e.m[4] * v0 + e.m[5];
            Ex[2] = e.m[6] * u0 + e.m[7] * v0 + e.m[8];
            double xEx = Ex[0] * u1 + Ex[1] * v1 + Ex[2];

            if (abs(xEx) > 1e8)
            {
                cout << "Ai[" << i << "]: xEx=" << xEx << ",        [" << to_string(u0) << "," << to_string(v0)
                    << "]->[" << to_string(u1) << "," << to_string(v1) << "]"
                << endl;
                errCnt++;
            }
        }
        cout << "Check " << eStr + ", xp*E*x < eps, error cnt = " << errCnt + "/" << AiArr.size();
    }

    //e = R^-1 * t
    //e' = t
    //li = ([t]xR)^-1 * mi'
    //li'= ([t]xR)    * mi
    double GetResidual(MyMat R, Vec T)
    {
        assert(R.RowSize() == 3 && R.ColSize() == 3 && T.Size() == 3);
        double reciTlen = 1.0 / T.Length();
        MyMat tx = MyMat::SetByTCross(T.m[0] * reciTlen, T.m[1] * reciTlen, T.m[2] * reciTlen);
        MyMat txR = tx * R;
        //txR = txR.T();
        double resi = 0;
        Vec mi(3);
        Vec mj(3);
        for (int i = 0; i < AiArr.size(); i++)
        {
            Vec ai = AiArr[i];
            mi.m[0] = ai.m[6];
            mi.m[1] = ai.m[7];
            mi.m[2] = 1.0;
            mj.m[0] = ai.m[2];
            mj.m[1] = ai.m[5];
            mj.m[2] = 1.0;
            double diff = mj.Dot(txR * mi);
            diff *= diff;
            Vec lj = txR * mi;
            Vec li = txR.T() * mj;
            lj = lj * (1.0 / lj.m[2]);
            li = li * (1.0 / li.m[2]);
            double s0 = 1.0 / (li.m[0] * li.m[0] + li.m[1] * li.m[1]);
            double s1 = 1.0 / (lj.m[0] * lj.m[0] + lj.m[1] * lj.m[1]);
            resi += diff * (s0 + s1);
        }

        return resi;
    }

};

