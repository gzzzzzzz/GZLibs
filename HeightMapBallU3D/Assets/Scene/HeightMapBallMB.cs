//#define def_hm2
using System.Collections;
using System.Collections.Generic;
using System.Security.Cryptography;
using System.Text;
using UnityEngine;
public class MathUtils
{
    public static double Det2(double[,] m)
    {
        double ret = m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0];
        return ret;
    }
    public static bool Inv22(double[,] m)
    {
        Debug.Assert(m.Rank == 2 && m.GetLength(0) == 2 && m.GetLength(1) == 2);
        double _d = Det2(m);
        if (_d == 0.0)
            return false;
        _d = 1.0 / _d;
        double a = m[0, 0];
        double b = m[0, 1];
        double c = m[1, 0];
        double d = m[1, 1];
        m[0, 0] = _d * d;
        m[0, 1] = -_d * b;
        m[1, 0] = -_d * c;
        m[1, 1] = _d * a;
        return true;
    }
    public static double Det3(double[,] m)
    {
        double ret = m[0, 0] * m[1, 1] * m[2, 2] + m[0, 1] * m[1, 2] * m[2, 0] + m[0, 2] * m[1, 0] * m[2, 1];
        ret -= m[0, 0] * m[1, 2] * m[2, 1] + m[0, 1] * m[1, 0] * m[2, 2] + m[0, 2] * m[1, 1] * m[2, 0];
        return ret;
    }
    public static bool Inv33(double[,] m)
    {
        Debug.Assert(m.Rank == 2 && m.GetLength(0) == 3 && m.GetLength(1) == 3);
        double _d = Det3(m);
        if (_d == 0.0)
            return false;
        _d = 1.0 / _d;
        double a = m[0, 0];
        double b = m[0, 1];
        double c = m[0, 2];
        double d = m[1, 0];
        double e = m[1, 1];
        double f = m[1, 2];
        double g = m[2, 0];
        double h = m[2, 1];
        double i = m[2, 2];
        m[0, 0] = _d * (e * i - h * f);
        m[0, 1] = -_d * (b * i - h * c);
        m[0, 2] = _d * (b * f - c * e);

        m[1, 0] = _d * (f * g - i * d);
        m[1, 1] = -_d * (c * g - i * a);
        m[1, 2] = _d * (c * d - a * f);

        m[2, 0] = _d * (d * h - g * e);
        m[2, 1] = -_d * (a * h - g * b);
        m[2, 2] = _d * (a * e - b * d);
        return true;
    }
    public static double DetIJ(double[,] mat44, int i, int j)
    {
        Debug.Assert(mat44.GetLength(0) == 4 && mat44.GetLength(1) == 4);
        int x, y, ii, jj;
        double ret = 0;
        double[,] mat = new double[3, 3];
        x = 0;

        for (ii = 0; ii < 4; ii++)
        {
            if (ii == i) continue;
            y = 0;
            for (jj = 0; jj < 4; jj++)
            {
                if (jj == j) continue;
                mat[x, y] = mat44[ii, jj];
                y++;
            }
            x++;
        }
        ret += mat[0, 0] * (mat[1, 1] * mat[2, 2] - mat[2, 1] * mat[1, 2]);
        ret -= mat[0, 1] * (mat[1, 0] * mat[2, 2] - mat[2, 0] * mat[1, 2]);
        ret += mat[0, 2] * (mat[1, 0] * mat[2, 1] - mat[2, 0] * mat[1, 1]);

        return ret;
    }
    public static double[,] Inv44(double[,] mat44)
    {
        Debug.Assert(mat44.GetLength(0) == 4 && mat44.GetLength(1) == 4);
        double[,] mInverse = new double[4, 4];
        int i, j;
        double det, detij;

        // calculate 4x4 determinant
        det = 0.0;
        for (i = 0; i < 4; i++)
        {
            //(i & 0x1)与操作代替-1……（i+j）次幂
            det += (i & 0x1) != 0 ? (-mat44[0, i] * DetIJ(mat44, 0, i)) : (mat44[0, i] * DetIJ(mat44, 0, i));
        }
        if (Mathf.Abs((float)det) < 1e-10)
        {
            Debug.LogWarning("function Inv44(m) : det(m) == 0, return Mat44(4,4,0)");
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
                mInverse[i, j] = ((i + j) & 0x1) != 0 ? (-detij * det) : (detij * det);
            }
        }
        return mInverse;
    }
    public static double[,] MatMul(double[,] a, double[,] b)
    {
        Debug.Assert(a.GetLength(1) == b.GetLength(0));
        int row = a.GetLength(0);
        int col = b.GetLength(1);
        int m = a.GetLength(1);
        double[,] ret = new double[row, col];
        for (int i = 0; i < row; i++)
            for (int j = 0; j < col; j++)
            {
                ret[i, j] = 0;
                for (int k = 0; k < m; k++)
                    ret[i, j] += a[i, k] * b[k, j];
            }
        return ret;
    }
    public static List<Vector2> Calc_CM(List<Vector2> xy_list, ref Vector2 xy_cm)
    {
        List<Vector2> ret = new List<Vector2>();
        xy_cm = Vector2.zero;
        for (int i = 0; i < xy_list.Count; i++)
            xy_cm += xy_list[i];
        float reci_cnt = 1.0f / xy_list.Count;
        xy_cm *= reci_cnt;
        for (int i = 0; i < xy_list.Count; i++)
            ret.Add(xy_list[i] - xy_cm);
        return ret;
    }
    public static bool Calc_abcd_xy_E_uv(List<Vector2> xy_list, List<Vector2> uv_list,
        ref Vector4 abcd, ref Vector2 xy_cm, ref Vector2 uv_cm)
    {
        if (xy_list.Count != uv_list.Count || xy_list.Count < 4)
            return false;
        int row = xy_list.Count;
        xy_cm = Vector2.zero;
        uv_cm = Vector2.zero;
        List<Vector2> norm_xy_list = Calc_CM(xy_list, ref xy_cm);
        List<Vector2> norm_uv_list = Calc_CM(uv_list, ref uv_cm);

        int n = row * 2;
        double[,] A = new double[n, 4];
        double[,] At = new double[4, n];
        double[,] y = new double[n, 1];
        for (int i = 0; i < row; i++)
        {
            double cx = (double)norm_xy_list[i].x;
            double cy = (double)norm_xy_list[i].y;
            double cu = (double)norm_uv_list[i].x;
            double cv = (double)norm_uv_list[i].y;
            int j = i * 2;
            A[j, 0] = At[0, j] = cx;
            A[j, 1] = At[1, j] = cy;
            A[j, 2] = At[2, j] = 0;
            A[j, 3] = At[3, j] = 0;
            y[j, 0] = cu;
            j++;
            A[j, 0] = At[0, j] = 0;
            A[j, 1] = At[1, j] = 0;
            A[j, 2] = At[2, j] = cx;
            A[j, 3] = At[3, j] = cy;
            y[j, 0] = cv;
        }
        double[,] AtA = MatMul(At, A);//4x4
        double[,] inv_AtA = Inv44(AtA);
        double[,] AtA_inv_At = MatMul(inv_AtA, At);//4xN
        double[,] ans = MatMul(AtA_inv_At, y);
        abcd.x = (float)ans[0, 0];
        abcd.y = (float)ans[1, 0];
        abcd.z = (float)ans[2, 0];
        abcd.w = (float)ans[3, 0];
        return true;
    }
    public static float Calc_abcd_xy_E_uv2(List<Vector2> xy_list, List<Vector2> uv_list,
        ref Vector4 abcd, float speed = 0.01f, int max_itr = 30)
    {
        Vector2 xy_cm = Vector2.zero;
        Vector2 uv_cm = Vector2.zero;
        List<Vector2> norm_xy_list = Calc_CM(xy_list, ref xy_cm);
        List<Vector2> norm_uv_list = Calc_CM(uv_list, ref uv_cm);

        abcd = new Vector4(1, -1, 1, 1);
        float err = 0;
        for (int itr = 0; itr < max_itr; itr++)
        {
            err = 0;
            for (int i = 0; i < xy_list.Count; i++)
            {
                Vector2 xy = norm_xy_list[i];
                Vector2 uv = norm_uv_list[i];
                float cu = abcd.x * xy.x + abcd.y * xy.y;
                float cv = abcd.z * xy.x + abcd.w * xy.y;
                float du = uv.x - cu;
                float dv = uv.y - cv;
                err += du * du + dv * dv;
                if (Mathf.Abs(xy.x) > 1e-1f)
                {
                    abcd.x += speed * (du / xy.x);
                    abcd.z += speed * (dv / xy.x);
                }
                if (Mathf.Abs(xy.y) > 1e-1f)
                {
                    abcd.y += speed * (du / xy.y);
                    abcd.w += speed * (dv / xy.y);
                }
            }
        }
        return err;
    }
    public static Vector2 ApplyLT(Vector4 abcd, Vector3 pos, Vector2 cm)
    {
        Vector2 dir = new Vector2(pos.x - cm.x, pos.y - cm.y);
        float map_col = abcd.x * dir.x + abcd.y * dir.y;
        float map_row = abcd.z * dir.x + abcd.w * dir.y;
        return new Vector2(map_col, map_row);
    }
    public static Vector4 InvABCD(Vector4 abcd)
    {
        float a = abcd.x;
        float b = abcd.y;
        float c = abcd.z;
        float d = abcd.w;
        float det = a * d - b * c;
        if (Mathf.Abs(det) < 1e-5f)
            return Vector4.zero;
        det = 1.0f / det;
        return new Vector4(d, -b, -c, a) * det;
    }
    public static Vector4 MatLT22(Vector4 v0, Vector4 v1)
    {
        return new Vector4(
            v0.x * v1.x + v0.y * v1.z,
            v0.x * v1.y + v0.y * v1.w,
            v0.z * v1.x + v0.w * v1.z,
            v0.z * v1.y + v0.w * v1.w
            );
    }
    public static void CheckCalcABCD()
    {
        //float T = (Random.value * 2 - 1) * Mathf.PI / 4.0f;
        float T = (Random.value * 90 - 45);
        T = 45;
        T *= Mathf.Deg2Rad;
        float sa = Random.value + 0.5f;
        float sb = Random.value + 0.5f;
        float sc = Random.value + 0.5f;
        float sd = Random.value + 0.5f;
        sa = 1.29f;sb = 1.23f;sc = 1.22f;sd = 1.04f;
        Debug.Log("T=" + (T * Mathf.Rad2Deg)
            + ", [a,b]:" + sa.ToString("F2") + "/" + sb.ToString("F2")
            + ", [c,d]:" + sc.ToString("F2") + "/" + sd.ToString("F2"));

        float cosT = Mathf.Cos(T);
        float sinT = Mathf.Sin(T);
        float m00 = sa * cosT;
        float m01 = sb * -sinT;
        float m10 = sa * sinT;
        float m11 = sb * cosT;

        float n00 = sc * m00;
        float n01 = sc * m01;
        float n10 = sd * m10;
        float n11 = sd * m11;
        Vector4 gt_ans = new Vector4(n00, n01, n10, n11);
        Vector4 gt_inv = InvABCD(gt_ans);
        Vector4 check_gt0 = MatLT22(gt_ans, gt_inv);
        Vector4 check_gt1 = MatLT22(gt_inv, gt_ans);
        Debug.Log("gt_ans=" + gt_ans.ToString("F4") + ", gt_inv=" + gt_inv.ToString("F4")
            + "\n\t check0=" + check_gt0.ToString("F2") + ", check1=" + check_gt1.ToString("F2"));

        List<Vector2> xy_list = new List<Vector2>();
        List<Vector2> uv_list = new List<Vector2>();

        for (int i = 0; i < 20; i++)
        {
            Vector2 xy = new Vector2(Random.value, Random.value) * 2 - Vector2.one;
            xy *= 15;
            float rx = n00 * xy.x + n01 * xy.y;
            float ry = n10 * xy.x + n11 * xy.y;
            Vector2 uv = Vector2.one*5 + new Vector2(rx, ry);
            uv += (new Vector2(Random.value, Random.value) * 2 - Vector2.one) * 0.5f;
            xy_list.Add(xy);
            uv_list.Add(uv);
        }
        Vector4 ans = Vector4.zero;
        Vector4 ans_inv = Vector4.zero;
        Vector2 cm0 = Vector2.zero;
        Vector2 cm1 = Vector2.zero;
        MathUtils.Calc_abcd_xy_E_uv(xy_list, uv_list, ref ans, ref cm0, ref cm1);
        MathUtils.Calc_abcd_xy_E_uv(uv_list, xy_list, ref ans_inv, ref cm0, ref cm1);
        check_gt0 = MatLT22(ans, ans_inv);
        check_gt1 = MatLT22(ans_inv, ans);
        Debug.Log("calc_ans=" + ans.ToString("F4") + ", calc_ans_inv=" + ans_inv.ToString("F4")
            + "\n\t check0=" + check_gt0.ToString("F2") + ", check1=" + check_gt1.ToString("F2"));

         //ans = Vector4.zero;
        //float final_err = MathUtils.Calc_abcd_xy_E_uv2(xy_list, uv_list, ref ans, 0.005f, 20);
        //Debug.Log("calc_ans2=" + ans.ToString("F4") + ",\n final err="+final_err);
    }
}
public class MatchLoss
{
#if def_hm2
        public float loss_match_back = 0;
#endif
    public float loss_fg_no_match = 0, loss_fg2bg2fg_back = 0, loss_fg_break_sca = 0;
    public float loss_bg_no_match = 0, loss_bg2fg2bg_back = 0, loss_bg_break_sca = 0;
    public float loss_fg_keep_shape = 0, loss_sub_fg_keep_shape = 0;

    public float w_cover_match = 1;
    public float w_match_back = 1;
    public float w_break_sca = 1;
    public float w_keep_shape = 1;

    public int fg_no_match_val_cnt = 0, fg_can_not_match_back_cnt = 0, fg_break_sca_cnt = 0;
    public int bg_no_match_val_cnt = 0, bg_can_not_match_back_cnt = 0, bg_break_sca_cnt = 0;
    public int fg_break_shape_cnt = 0;
    public MatchLoss()
    {
        Clear();
        SetW(1, 1, 1);
    }
    public void DivFgBgCnt(int fg_cnt, int bg_cnt)
    {
        float reci_cnt = 1.0f / bg_cnt;
        loss_bg_no_match *= (bg_no_match_val_cnt * reci_cnt);
        loss_bg_break_sca *= (bg_break_sca_cnt * reci_cnt);
        loss_bg2fg2bg_back *= (bg_can_not_match_back_cnt * reci_cnt);

        reci_cnt = 1.0f / fg_cnt;
        loss_fg_no_match *= (fg_no_match_val_cnt * reci_cnt);
        loss_fg_break_sca *= (fg_break_sca_cnt * reci_cnt);
        loss_fg2bg2fg_back *= (fg_can_not_match_back_cnt * reci_cnt);
        loss_fg_keep_shape *= (fg_break_shape_cnt * reci_cnt);
    }
    public void SetW(float _match_back_w, float _no_match_w, float _break_sca_w)
    {
        //w_match_back = _match_back_w;
        //w_fg_no_match = _no_match_w;
        //w_bg_no_match = _no_match_w;
        //w_break_sca = _break_sca_w;
    }
    public void Clear()
    {
        loss_fg_no_match = loss_fg2bg2fg_back = loss_fg_break_sca = 0;
        loss_bg_no_match = loss_bg2fg2bg_back = loss_bg_break_sca = 0;
        loss_fg_keep_shape = 0;
        loss_sub_fg_keep_shape = 1e10f;

        fg_no_match_val_cnt = fg_can_not_match_back_cnt = fg_break_sca_cnt = 0;
        bg_no_match_val_cnt = bg_can_not_match_back_cnt = bg_break_sca_cnt = 0;
        fg_break_shape_cnt = 0;
    }
    public float GetMatchVal()
    {
        float keep_shape = Mathf.Min(loss_sub_fg_keep_shape, loss_fg_keep_shape);
        return loss_fg_no_match * w_cover_match
                + loss_bg_no_match * w_cover_match
                + loss_fg2bg2fg_back * w_match_back
                + loss_bg2fg2bg_back * w_match_back
                + loss_fg_break_sca * w_break_sca
                + loss_bg_break_sca * w_break_sca
            + keep_shape * w_keep_shape;
    }
    public string GetMatchStr(string fmt = "F2", bool split2 = false)
    {
        return "keep(" + (loss_fg_keep_shape * w_keep_shape).ToString(fmt)
            + "/" + (loss_sub_fg_keep_shape > 1e9f ? "null": (loss_sub_fg_keep_shape*w_keep_shape).ToString(fmt))
            + (split2 ? ")\n" : "), ")
            + "cov(" + (loss_fg_no_match * w_cover_match).ToString(fmt) + "/" + (loss_bg_no_match * w_cover_match).ToString(fmt)
            + "), mbk(" + (loss_fg2bg2fg_back * w_match_back).ToString(fmt) + "/" + (loss_bg2fg2bg_back * w_match_back).ToString(fmt)
            + "), brk(" + (loss_fg_break_sca * w_break_sca).ToString(fmt) + "/" + (loss_bg_break_sca * w_break_sca).ToString(fmt)
            + ")";
    }
}
public class HBNeibor
{
    //public int cnt = 0;
    HBall self_hb;
    public List<HBall> neibors = new List<HBall>();
    public List<float> neibors_w;
    //public Vector2 neibors_cw;
    //public float neibors_sum_w = 0;
    public List<Vector3> init_cw2b_dir = new List<Vector3>();
    public List<Vector3> rt_cw2b_dir = new List<Vector3>();
    public Vector3 init_cw, rt_cw;
    public float init_sz, rt_sz, lcl_cw_sz_sca = 1, lcl_cw_theta = 0, rt_cw_cnt = 1;
    public Vector3 GetTarGridPos(HBall[,] rel_hb = null)
    {
        UpdateRTCW(rel_hb);
        lcl_cw_sz_sca = CalcCWSca(rt_cw2b_dir, rt_cw, init_cw2b_dir, rel_hb);
        lcl_cw_theta = CalcCWRot();
        Quaternion rot_z = Quaternion.AngleAxis(lcl_cw_theta, Vector3.forward);
        return rt_cw + rot_z * (lcl_cw_sz_sca * init_cw2b_dir[0]);
    }
    public Vector3 UpdateRTCW(HBall[,] rel_hb = null)
    {
        rt_cw = self_hb.SamplePos(rel_hb);// tfm.localPosition;
        rt_cw_cnt = 1;
        //foreach (var n in neibors)
        for (int i = 0; i < neibors.Count; i++)
        {
            if (!neibors[i].IsValidVal())
                continue;
            rt_cw += neibors[i].SamplePos(rel_hb);// tfm.localPosition;
            rt_cw_cnt++;
        }
        rt_cw /= rt_cw_cnt;
        return rt_cw;
    }
    
    public float CalcCWSca(List<Vector3> _cw_dir, Vector3 _cw, List<Vector3> _ini_dir, HBall[,] rel_hb = null)
    {
        if (_cw_dir.Count == 0) return 1;
        Vector3 cw2b = self_hb.SamplePos(rel_hb) - _cw;// tfm.localPosition - _cw;
        float _sz_rt = cw2b.magnitude + 1e-5f;
        _cw_dir[0] = cw2b;
        float _sz_ini = _ini_dir[0].magnitude + 1e-5f;

        for (int i = 0; i < neibors.Count; i++)
        {
            if (!neibors[i].IsValidVal()) continue;
            cw2b = neibors[i].SamplePos(rel_hb) - _cw;// tfm.localPosition - _cw;
            cw2b.z = 0;
            _cw_dir[i + 1] = cw2b;
            _sz_rt += cw2b.magnitude;
            _sz_ini += _ini_dir[i + 1].magnitude;
        }
        rt_sz = _sz_rt;
        float ret = Mathf.Clamp(rt_sz / _sz_ini, 0.1f, 10f);
        return ret;
    }
    public float CalcCWRot()
    {
        float ret = 0;
        if (rt_cw2b_dir.Count != init_cw2b_dir.Count || rt_cw2b_dir.Count == 0)
        {
            Debug.Log("calc cw rot err, rt_cw2b_dir cnt=" + rt_cw2b_dir.Count);
            return 0;
        }
        ret = HBall.GetAngleBetween(init_cw2b_dir[0], rt_cw2b_dir[0]);
        float _cnt = 1;
        for (int i = 0; i < neibors.Count; i++)
        {
            if (!neibors[i].IsValidVal())
                continue;
            ret += HBall.GetAngleBetween(init_cw2b_dir[i + 1], rt_cw2b_dir[i + 1]);
            _cnt++;
        }
        return ret / _cnt;
    }
    //public float UpdateCW2BDir(List<Vector3> _cw_dir, Vector3 _cw)
    //{
    //    Vector3 cw2b = self_hb.tfm_pos - _cw;// tfm.localPosition - _cw;
    //    float _sz = cw2b.magnitude + 1e-5f;
    //    _cw_dir[0] = cw2b;

    //    for (int i = 0; i < neibors.Count; i++)
    //    {
    //        cw2b = neibors[i].tfm_pos - _cw;// tfm.localPosition - _cw;
    //        cw2b.z = 0;
    //        _cw_dir[i + 1] = cw2b;
    //        _sz += cw2b.magnitude;
    //    }
    //    return _sz;
    //}
    public float UpdateIniCW2BDir(ref List<Vector3> _cw_dir, Vector3 _ini_cw)
    {
        Vector3 cw2b = self_hb.GetInitPos() - _ini_cw;// tfm.localPosition - _cw;
        cw2b.z = 0;
        float _sz = cw2b.magnitude + 1e-5f;
        _cw_dir[0] = cw2b;

        for (int i = 0; i < neibors.Count; i++)
        {
            if (!neibors[i].IsValidVal())
                continue;
            cw2b = neibors[i].GetInitPos() - _ini_cw;// tfm.localPosition - _cw;
            cw2b.z = 0;
            _cw_dir[i + 1] = cw2b;
            _sz += cw2b.magnitude;
        }
        return _sz;
    }
    void BuildCWByNeibors()
    {
        init_cw = self_hb.GetInitPos();
        float cw_cnt = 1;
        init_cw2b_dir = new List<Vector3>();
        rt_cw2b_dir = new List<Vector3>();
        init_cw2b_dir.Add(Vector3.zero);
        rt_cw2b_dir.Add(Vector3.zero);
        for (int i = 0; i < neibors.Count; i++)
        {
            cw_cnt++;
            init_cw += neibors[i].GetInitPos();
            init_cw2b_dir.Add(Vector3.zero);
            rt_cw2b_dir.Add(Vector3.zero);
        }
        init_cw /= cw_cnt;

        init_sz = UpdateIniCW2BDir(ref init_cw2b_dir, init_cw);
        rt_sz = init_sz;
        lcl_cw_sz_sca = 1;
    }
    //public void BFSBuildNeibor(HBall hb, float max_dist = 6, int max_cnt = 40)
    public void BFSBuildNeibor(HBall hb, float max_dist = 5, int max_cnt = 30)
    {
        self_hb = hb;

        neibors = new List<HBall>();
        neibors_w = new List<float>();
        //neibors_sum_w = 0;
        //neibors_cw = Vector2.zero;
        int head = -1;
        HBall chb = hb;
        HashSet<int> hs = new HashSet<int>();
        hs.Add(chb.UID());
        while(neibors.Count < max_cnt && chb != null)
        {
            foreach(var nb8 in chb.neibor8)
                if (nb8 != null && nb8.IsValidVal())
                {
                    Vector2 dist_rc = new Vector2(nb8.col - self_hb.col, nb8.row - self_hb.row);
                    int hs_id = nb8.UID();
                    float _wi = dist_rc.magnitude;
                    if (_wi <= max_dist && neibors.Count < max_cnt && !hs.Contains(hs_id))
                    {
                        neibors.Add(nb8);
                        _wi = 1.0f / _wi;
                        neibors_w.Add(_wi);
                        //neibors_sum_w += _wi;
                        //neibors_cw += new Vector2(nb8.col, nb8.row) * _wi;
                        hs.Add(hs_id);
                    }
                }
            head++;
            if (neibors.Count <= head)
                break;
            chb = neibors[head];
        }

        BuildCWByNeibors();
    }
    public void BuildNeibor(HBall hb, List<HBall> candidates, float max_dist = 3, float min_dist = 0)
    {
        self_hb = hb;
        //if (!hb.IsValidVal()) return;

        neibors = new List<HBall>();
        neibors_w = new List<float>();
        //float max_dist2 = max_dist * max_dist;
        foreach (var c in candidates)
        {
            if (!c.IsValidVal())
                continue;
            Vector3 dir = new Vector3(c.col - self_hb.col, c.row - self_hb.row, 0);
            float dist = dir.magnitude;// sqrMagnitude;//magnitude;
            if (dist > min_dist && dist <= max_dist)//don't include self
            {
                neibors.Add(c);
                neibors_w.Add(1.0f / dist);
            }
        }
        BuildCWByNeibors();
    }
}
public class LineOrder
{

    public int Count = 0;
    //public List<int> id_list = new List<int>();
    //public List<int> pid_list = new List<int>();
    public List<HBall> id_list = new List<HBall>();
    public List<HBall> pid_list = new List<HBall>();
    public List<int> depth_list = new List<int>();
    public HashSet<int> hs = new HashSet<int>();
    public LineOrder()
    {
        Clear();
    }
    public void Clear()
    {
        Count = 0;
        hs.Clear();
        pid_list.Clear();
        depth_list.Clear();
        id_list.Clear();
    }
    public void Add(HBall new_hb, HBall par_hb, int dep, int reverses = 0)
    //public void Add(int new_hb, int par_hb, int dep)
    {
        if (hs.Contains(new_hb.UID()))
            return;
        id_list.Add(new_hb);
        pid_list.Add(par_hb);
        depth_list.Add(dep);
        hs.Add(new_hb.UID());
        Count++;
        new_hb.dep_id = dep;
        if (reverses == 0)
            new_hb.copy_dep_id = dep;
        else
            new_hb.rvs_dep_id = dep;
    }
    public int FirstId()
    {
        if (id_list.Count == 0)
            return -1;
        return id_list[0].UID();
    }
    public int LastId()
    {
        if (id_list.Count == 0)
            return -1;
        return id_list[id_list.Count - 1].UID();
    }
    public int LastDepth()
    {
        if (depth_list.Count == 0)
            return -1;
        return depth_list[depth_list.Count - 1];
    }
}
public class HBall
{
    public bool fg_or_bg = true;
    public int sub_id = 0,expand_sid=0, dep_id = 0, copy_dep_id=0, rvs_dep_id=0, temp_id = 0;
    private int owner_sub_id = 0;
    public Vector3 tfm_pos, init_pos, thin_pos, prev_tfm_pos, prev_move_dir, temp_pos;
    public float val;
    public bool is_active = true;
    public readonly int row, col;
    public int bfs_root_row, bfs_root_col;
    public bool b_no_match, b_fg_break_sca, b_bg_break_sca, b_bg_covered, b_can_not_match_back, b_fg_keep_shape;
    public float dist_fg2bg, dist_bg2fg, dist_fg2bg2fg, dist_bg2fg2bg, fg_break_loss, bg_break_loss, dist_fg_keep_shape;
    public static float ms_break_sca_thresh = 0.8f;
    public static float ms_match_err_thresh = 1.01f;
    public Vector3 ini_dir2global_cw, ini_dir2sub_active_cw;
    float max_sca0 = 0;
    float max_sca1 = 0;
    public Vector2 dir_n = Vector2.up, dir_t = Vector2.right, dir_s = Vector2.right, dir_t_mean = Vector2.right;
    public float curvature = 0;
    public Vector3 GetPos(int pt = 0)
    {
        if (pt == 0) return init_pos;
        else if (pt == 1) return tfm_pos;
        else return thin_pos;
    }
    public Vector2 PosV2(bool init_or_rt)
    {
        return init_or_rt ? new Vector2(col, row) : new Vector2(tfm_pos.x, tfm_pos.y);
    }
    public Vector3 PosV3(bool init_or_rt)
    {
        return init_or_rt ? new Vector3(col, row, 0) : new Vector3(tfm_pos.x, tfm_pos.y, 0);
    }
    public void ClearLineInfo()
    {
        dir_n = Vector2.right;
        dir_t = Vector2.up;
        dir_t_mean = Vector2.up;
        dir_s = Vector2.zero;
        curvature = 0;
    }

    public float[] b_match_bg_dist4 = new float[4];//min b_match_bg_dist
    public HBall[] match_self_bg4 = new HBall[4];
    public HBall match_self_bg = null;
    //public HBall[] match_self_fg4 = new HBall[4];
    public HBall match_self_fg = null;
    //------- neibor grid ------
    public HBNeibor neibor_lod0 = new HBNeibor();
    public HBNeibor neibor_lod1 = new HBNeibor();
    public HBNeibor neibor_lod2 = new HBNeibor();
#if def_hm2
    HBall[] match_b4 = new HBall[4];
#endif
    public HBall[] neibor4 = new HBall[4];
    public HBall[] neibor8 = new HBall[8];
    public HBall[] neibor4_ignore = new HBall[4];
    public HBall[] neibor8_ignore = new HBall[8];
    public int neibor8_null_cnt = 0;
    public HBall(int _row, int _col, float _val, int _owner_hm_sub_id = 0)
    {
        row = _row;
        col = _col;
        init_pos = new Vector3(col, row, 0);
        owner_sub_id = _owner_hm_sub_id;
        sub_id = 0;
        Reset(_val);
    }
    public void Reset(float _val)
    {
        val = _val;
        bfs_root_row = row;
        bfs_root_col = col;
        tfm_pos = new Vector3(col, row, 0);
        prev_tfm_pos = tfm_pos;
        prev_move_dir = Vector3.zero;
        match_self_fg = null;
        match_self_bg = null;

        dir_n = dir_t = dir_s = dir_t_mean = Vector2.zero;
    }
    public static HBall GetHBSameRC(HBall[,] hb_arr, HBall tar_p)
    {
        return hb_arr[tar_p.row, tar_p.col];
    }
    public Vector3 SamplePos(HBall[,] hb_arr)
    {
        if (hb_arr == null) return tfm_pos;
        return HBall.GetHBSameRC(hb_arr, this).tfm_pos;
    }
    public static int RC2I(int row, int col)
    {
        return col | (row << 16);
    }
    public static void I2RC(int i, ref int row, ref int col)
    {
        row = i >> 16;
        col = i & 0xffff;
    }
    public int UID()
    {
        return RC2I(row, col);
    }
    public static Vector2 Rot90b(Vector2 a)
    {
        return new Vector2(-a.y, a.x);
    }
    public static Vector2 Rot90(Vector2 a)
    {
        return new Vector2(a.y, -a.x);
    }
    public static Vector2 Rot45(Vector2 a)
    {
        return new Vector2(0.7071f * (a.x + a.y), 0.7071f * (a.y - a.x));
    }
    public static Vector2 RotN(Vector2 a, float angle)
    {
        angle *= Mathf.Deg2Rad;
        float cosn = Mathf.Cos(angle);
        float sinn = Mathf.Sin(angle);
        return new Vector2(cosn * a.x + sinn * a.y, -sinn * a.x + cosn * a.y);
    }
    public static Vector2 RotNb(Vector2 a, float angle)
    {
        angle *= -Mathf.Deg2Rad;
        float cosn = Mathf.Cos(angle);
        float sinn = Mathf.Sin(angle);
        return new Vector2(cosn * a.x + sinn * a.y, -sinn * a.x + cosn * a.y);
    }
    public static Vector2 Rot45b(Vector2 a)
    {
        return new Vector2(0.7071f * (a.x - a.y), 0.7071f * (a.y + a.x));
    }
    public static Vector2 Rot25(Vector2 a)
    {
        return new Vector2(0.9063f * a.x + 0.4226f * a.y, -0.4226f * a.x + 0.9063f * a.y);
    }
    public static Vector2 Rot25b(Vector2 a)
    {
        return new Vector2(0.9063f * a.x - 0.4226f * a.y, 0.4226f * a.x + 0.9063f * a.y);
    }
    public static Vector2 Rot15(Vector2 a)
    {
        return new Vector2(0.9659f * a.x + 0.2588f * a.y, -0.2588f * a.x + 0.9659f * a.y);
    }
    public static Vector2 Rot15b(Vector2 a)
    {
        return new Vector2(0.9659f * a.x - 0.2588f * a.y, 0.2588f * a.x + 0.9659f * a.y);
    }
    public static Vector2 Rot5(Vector2 a)
    {
        return new Vector2(0.9962f * a.x + 0.0872f * a.y, -0.0872f * a.x + 0.9962f * a.y);
    }
    public static Vector2 Rot5b(Vector2 a)
    {
        return new Vector2(0.9962f * a.x - 0.0872f * a.y, 0.0872f * a.x + 0.9962f * a.y);
    }
    public const float ms_active_thresh = 0.1f;
    public bool IsValGreater(float _thresh = ms_active_thresh)
    {
        return (val > _thresh);
    }
    public bool IsValidVal()
    {
        is_active = IsActive();
        return is_active;
    }
    private bool IsActive()
    {
        if (val < ms_active_thresh) return false;
        if (owner_sub_id == 0) return true;
        return (owner_sub_id & sub_id) != 0;
    }
    public void UpdateNeibor4(HBall[,] _hb_grid)
    {
        int ni4 = 0;
        int ni8 = 0;
        neibor8_null_cnt = 0;
        for (int oi = -1; oi < 2; oi++)
            for (int oj = -1; oj < 2; oj++)
            {
                int oij = oi * oi + oj * oj;
                if (oij == 0) continue;
                int nr = row + oi;
                int nc = col + oj;
                HBall nrc_hb_ignore_sub_id = null;
                HBall nrc_hb = null;
                if (HeightMap.ValidRC(nr, nc))
                {
                    nrc_hb = _hb_grid[nr, nc];
                    nrc_hb_ignore_sub_id = nrc_hb;
                    if (nrc_hb_ignore_sub_id.val < HBall.ms_active_thresh)
                        nrc_hb_ignore_sub_id = null;
                    if (!nrc_hb.IsValidVal())
                        nrc_hb = null;
                }
                if (oij == 1)
                {
                    neibor4[ni4] = nrc_hb;
                    neibor4_ignore[ni4] = nrc_hb_ignore_sub_id;
                    ni4++;
                }
                neibor8[ni8] = nrc_hb;
                neibor8_ignore[ni8] = nrc_hb_ignore_sub_id;
                ni8++;
                if (nrc_hb == null)
                    neibor8_null_cnt++;
            }
    }
    public float NeiborBreakScaRate()
    {
        HBall[] _neibors = neibor8;
        float break_cnt = 0;
        float all_cnt = 0;
        foreach (var _n in _neibors)
        {
            if (_n == null) continue;
            all_cnt++;
            if (_n.b_fg_break_sca)
                break_cnt++;
        }
        if (all_cnt == 0)
            return 0;
        return break_cnt;// / all_cnt;
    }
    
    //public float NeiborCoveredRate()
    //{
    //    HBall[] _neibors = neibor8;
    //    float covered_cnt = 0;
    //    float all_cnt = 0;
    //    foreach(var _n in _neibors)
    //    {
    //        if (_n == null) continue;
    //        all_cnt++;
    //        if (_n.b_bg_covered)
    //            covered_cnt++;
    //    }
    //    if (all_cnt == 0)
    //        return 0;
    //    return covered_cnt / all_cnt;
    //}
    public Vector3 GetInitPos()
    {
        return new Vector3(col, row, 0);
    }
    public Vector2 GetXY2()
    {
        return new Vector2(col, row);
    }
    public string RCStr()
    {
        return (new Vector2(row, col)).ToString("F0");
    }
    public string GetDebugStr(int fg_bg = 0)
    {
        string fmt = "F2";
        StringBuilder sbd = new StringBuilder();
        sbd.Append("RC[" + row + ", " + col + "], tfm_pos:" + tfm_pos.ToString(fmt));
        if (fg_bg == 0)
        {
            sbd.Append("\n\t match_fg:" + (match_self_fg == null ? "NULL" : ("[" + match_self_fg.row + "," + match_self_fg.col + "]")));
            sbd.Append("\n\t match_bg:" + (match_self_bg == null ? "NULL" : ("[" + match_self_bg.row + "," + match_self_bg.col + "]")));
            sbd.Append("\n\t cover:" + dist_fg2bg.ToString(fmt) + "/" + dist_bg2fg.ToString(fmt));
            sbd.Append("\n\t m-back:" + dist_fg2bg2fg.ToString(fmt) + "/" + dist_bg2fg2bg.ToString(fmt));
            sbd.Append("\n\t break:" + fg_break_loss.ToString(fmt) + "/" + bg_break_loss.ToString(fmt));
        }
        else if (fg_bg < 0)//bg ball
        {
            sbd.Append("\n\t match_fg:" + (match_self_fg == null ? "NULL" : ("[" + match_self_fg.row + "," + match_self_fg.col + "]")));
            sbd.Append("\n\t bg-covered:" + dist_bg2fg.ToString(fmt));
            sbd.Append("\n\t bg-m-back:" + dist_bg2fg2bg.ToString(fmt));
            sbd.Append("\n\t bg-break:" + bg_break_loss.ToString(fmt));
            sbd.Append(", bg-mean-sca(" + bg_mean_sca.ToString(fmt)
                + "/" + bg_mean_sca2.ToString(fmt) + ", " + bg_sca_cnt
                + "), dx_sca=" + (bg_dx_sca_sqrt).ToString(fmt));
        }
        else//fg ball
        {
            sbd.Append("\n\t match_bg:" + (match_self_bg == null ? "NULL" : ("[" + match_self_bg.row + "," + match_self_bg.col + "]")));
            sbd.Append("\n\t fg-cover:" + dist_fg2bg.ToString(fmt));
            sbd.Append("\n\t fg-m-back:" + dist_fg2bg2fg.ToString(fmt));
            sbd.Append("\n\t fg-break:" + fg_break_loss.ToString(fmt));
            sbd.Append(", fg-mean-sca(" + fg_mean_sca.ToString(fmt)
                + "/" + fg_mean_sca2.ToString(fmt) + ", " + fg_sca_cnt
                + "), dx_sca=" + (fg_dx_sca_sqrt).ToString(fmt));
        }
        sbd.Append("\n\t max_sca(" + max_sca0.ToString(fmt) + "/" + max_sca1.ToString(fmt)+")");
        sbd.Append("\n\t sub-id:" + sub_id + "/" + owner_sub_id);
        return sbd.ToString();
    }
    
    
    public void UnmatchFGBDragByCW(float drag_unmatch_lbd = 0.3f)
    {
        if (is_active && b_no_match && match_self_bg != null)
        {
            tfm_pos = Vector3.Lerp(tfm_pos, match_self_bg.GetInitPos(), drag_unmatch_lbd);
            /*
            for (int i = 0; i < neibor_lod0.neibors.Count; i++)
            {
                var nbb = neibor_lod0.neibors[i];
                {
                    if (!nbb.b_no_match && nbb.is_active)
                    {
                        tfm_pos = Vector3.Lerp(tfm_pos, nbb.tfm_pos, drag_unmatch_lbd);
                    }
                }
            }//*/
        }
    }
    public void SortGrid(float _grid_lbd = 0.5f)
    {
        Vector3 tar_pos_lod0 = neibor_lod0.GetTarGridPos();
        //Vector3 tar_pos_lod1 = neibor_lod1.GetTarGridPos();
        //Vector3 tar_pos_lod2 = neibor_lod2.GetTarGridPos();
        //Vector3 tar_pos = 0.5f * tar_pos_lod0 + 0.3f * tar_pos_lod1 + 0.2f * tar_pos_lod2;
        Vector3 tar_pos = tar_pos_lod0;
        tfm_pos = Vector3.Lerp(tfm_pos, tar_pos, _grid_lbd);
        //if (row == 11 && col == 10)
        //{
        //    Debug.Log("sz init=" + init_sz + ", rt sz =" + rt_sz + ", sz_sca =" + sz_sca);
        //    foreach (var dd in rt_cw2b_dir)
        //        Debug.DrawLine(dd + rt_cw, rt_cw, Color.red, 2);
        //}
    }

    public void AccMoveByMoveDir()
    {
        Vector3 cur_move_dir = tfm_pos - prev_tfm_pos;
        if (Vector3.Dot(cur_move_dir, prev_move_dir) > 0)
            tfm_pos = prev_tfm_pos + (cur_move_dir + 0.5f * prev_move_dir);
        prev_tfm_pos = tfm_pos;
        prev_move_dir = tfm_pos - prev_tfm_pos;
    }
    public static float X01Y01(float val, float x0, float x1, float y0, float y1)
    {
        if (x0 == x1)
            x0 = x1 - 1e-6f;
        val = Mathf.Clamp(val, Mathf.Min(x0, x1), Mathf.Max(x0, x1));
        return ((y0 * x1 - y1 * x0) + val * (y1 - y0)) / (x1 - x0);
    }
    public static float GetAngleBetween(Vector3 a, Vector3 b)
    {
        if (a.magnitude < 1e-5f || b.magnitude < 1e-5f)
            return 0;
        Vector3 c = Vector3.Cross(a.normalized, b.normalized);
        float neg = c.z < 0 ? -1 : 1;
        float ret = Mathf.Asin(Mathf.Clamp(neg * c.z, -1, 1)) * Mathf.Rad2Deg;
        if (Vector3.Dot(a, b) < 0)
            ret = 180 - ret;
        return ret * neg;
    }
    public static void TestGetAB()
    {
        float a0 = HBall.GetAngleBetween(new Vector3(1, 0, 0), new Vector3(1, 1, 0));
        float a1 = HBall.GetAngleBetween(new Vector3(1, 0, 0), new Vector3(0, 1, 0));
        float a2 = HBall.GetAngleBetween(new Vector3(1, 0, 0), new Vector3(-1, 1, 0));
        float b0 = HBall.GetAngleBetween(new Vector3(1, 0, 0), new Vector3(1, -1, 0));
        float b1 = HBall.GetAngleBetween(new Vector3(1, 0, 0), new Vector3(0, -1, 0));
        float b2 = HBall.GetAngleBetween(new Vector3(1, 0, 0), new Vector3(-1, -1, 0));
        Debug.Log("a:" + a0 + ", " + a1 + ", " + a2);
        Debug.Log("b:" + b0 + ", " + b1 + ", " + b2);
    }
    static void Max2(float val, ref float max0, ref float max1)
    {
        if (val > max0)
        {
            max1 = max0;
            max0 = val;
        }
        else if (val > max1)
        {
            max1 = val;
        }
    }
    float bg_mean_sca = 0;
    float bg_mean_sca2 = 0;
    float bg_dx_sca_sqrt = 0;
    int bg_sca_cnt = 0;
    public float BGMaxBreakSca()
    {
        bg_break_loss = 0;
        b_bg_break_sca = false;
        if (match_self_fg == null)
        {
            //Debug.Log("bg break sca, match self.fg == null");
            return bg_break_loss;
        }

        bg_mean_sca = 0;
        bg_mean_sca2 = 0;
        bg_sca_cnt = 0;

        Vector3 self_init_pos = GetInitPos();
        max_sca0 = max_sca1 = 0;
        //for (int i = 0; i < neibor8.Length; i++)
        for (int i = 0; i < neibor_lod0.neibors.Count; i++)
        {
            //HBall _nhb = neibor8[i];
            HBall _nhb = neibor_lod0.neibors[i];
            if (_nhb == null || !_nhb.IsValidVal() || _nhb.match_self_fg == null) continue;

            Vector3 dir0 = _nhb.GetInitPos() - self_init_pos;
            //Vector3 dir1 = _nhb.match_self_fg.tfm_pos - match_self_fg.tfm_pos;
            Vector3 dir1 = _nhb.match_self_fg.GetInitPos() - match_self_fg.GetInitPos();
            float _sca = AddStatSampleBreakSca(dir0, dir1, ref bg_sca_cnt,
                ref max_sca0, ref max_sca1,
                ref bg_mean_sca, ref bg_mean_sca2);
        }
        bg_break_loss = CaclBreakSca(bg_sca_cnt, ref max_sca0, ref max_sca1, 
            ref bg_mean_sca, ref bg_mean_sca2, ref bg_dx_sca_sqrt);
        b_bg_break_sca = bg_break_loss > ms_break_sca_thresh;
        return bg_break_loss;
    }
    float fg_mean_sca = 0;
    float fg_mean_sca2 = 0;
    float fg_dx_sca_sqrt = 0;
    int fg_sca_cnt = 0;
    public float FGMaxBreakSca()
    {
        fg_break_loss = 0;
        b_fg_break_sca = false;
        fg_mean_sca = 0;
        fg_mean_sca2 = 0;
        fg_sca_cnt = 0;

        Vector3 self_init_pos = GetInitPos();
        max_sca0 = max_sca1 = 0;
        for (int i = 0; i < neibor_lod0.neibors.Count; i++)
        //for (int i = 0; i < neibor8.Length; i++)
        {
            //HBall _nhb = neibor8[i];
            HBall _nhb = neibor_lod0.neibors[i];
            if (_nhb == null || !_nhb.IsValidVal()) continue;
            Vector3 dir0 = _nhb.GetInitPos() - self_init_pos;
            Vector3 dir1 = _nhb.tfm_pos - tfm_pos;
            AddStatSampleBreakSca(dir0, dir1, ref fg_sca_cnt,
                ref max_sca0, ref max_sca1,
                ref fg_mean_sca, ref fg_mean_sca2);

        }
        fg_break_loss = CaclBreakSca(fg_sca_cnt, ref max_sca0, ref max_sca1, 
            ref fg_mean_sca, ref fg_mean_sca2, ref fg_dx_sca_sqrt);
        b_fg_break_sca = fg_break_loss > ms_break_sca_thresh;

        return fg_break_loss;
    }
    float CaclBreakSca(int sca_cnt,  ref float max_sca0,  ref float max_sca1, 
        ref float mean_sca, ref float mean_sca2, ref float dx_sca_sqrt)
    {
        float ret = 0;
        if (max_sca1 == 0) max_sca1 = max_sca0;
        if (sca_cnt > 0)
        {
            mean_sca /= sca_cnt;
            mean_sca2 /= sca_cnt;
            float dx_sca = mean_sca2 - mean_sca * mean_sca;
            dx_sca_sqrt = Mathf.Sqrt(dx_sca);
            ret = dx_sca_sqrt / (1e-5f + mean_sca);
            ret *= max_sca0;
        }
        return ret;
    }
    float AddStatSampleBreakSca(Vector3 dir0, Vector3 dir1, ref int sca_cnt
        , ref float max_sca0, ref float max_sca1
        , ref float mean_sca, ref float mean_sca2)
    {
        dir0.z = dir1.z = 0;
        float _sca = 1;
        if (dir0.magnitude > 0)
        {
            sca_cnt++;
            _sca = Mathf.Clamp(dir1.magnitude / dir0.magnitude, 0.25f, 4);
            Max2(_sca, ref max_sca0, ref max_sca1);
            mean_sca += _sca;
            mean_sca2 += _sca * _sca;
        }
        return _sca;
    }
#region sample-unsample-write-read
    public static float SampleI(float[,] _h, int ix, int iy)
    {
        if (ix < 0) ix = 0;
        if (ix >= HeightMap.iw1) ix = HeightMap.iw1;
        if (iy < 0) iy = 0;
        if (iy >= HeightMap.ih1) iy = HeightMap.ih1;
        return _h[iy, ix];
    }
    public static float SampleF(float[,] _h, Vector3 _pxy)
    {
        return SampleF(_h, _pxy.x, _pxy.y);
    }

    public static void Fxy24Ixy(float fx, float fy,
        out int x0, out int x1, out int y0, out int y1,
        out float lbd_x, out float lbd_y)
    {
        float _px = Mathf.Clamp(fx, 0, HeightMap.iw1 - 0.0001f);
        float _py = Mathf.Clamp(fy, 0, HeightMap.ih1 - 0.0001f);
        x0 = (int)(_px);
        y0 = (int)(_py);
        x1 = x0 + 1;
        y1 = y0 + 1;
        lbd_x = _px - x0;
        lbd_y = _py - y0;
    }
    public static float SampleF(float[,] _h, float fx, float fy)
    {
        int x0, x1, y0, y1;
        float lbd_x, lbd_y;
        Fxy24Ixy(fx, fy, out x0, out x1, out y0, out y1, out lbd_x, out lbd_y);

        float h00 = _h[y0, x0];
        float h01 = _h[y1, x0];
        float h0 = Mathf.Lerp(h00, h01, lbd_y);
        float h10 = _h[y0, x1];
        float h11 = _h[y1, x1];
        float h1 = Mathf.Lerp(h10, h11, lbd_y);
        float hm = Mathf.Lerp(h0, h1, lbd_x);
        return hm;
    }
    public static void UnSampleF(float[,] _h, float fx, float fy, float _val)
    {
        int x0, x1, y0, y1;
        float lbd_x, lbd_y;
        Fxy24Ixy(fx, fy, out x0, out x1, out y0, out y1, out lbd_x, out lbd_y);

        float hy0 = 0, hy1 = 0;//y0+y1=_val
        UnLerp(_val, lbd_y, ref hy0, ref hy1);
        float hx0y0 = 0, hx1y0 = 0;
        UnLerp(hy0, lbd_x, ref hx0y0, ref hx1y0);
        float hx0y1 = 0, hx1y1 = 0;
        UnLerp(hy1, lbd_x, ref hx0y1, ref hx1y1);

        _h[y0, x0] += hx0y0;
        _h[y0, x1] += hx1y0;
        _h[y1, x0] += hx0y1;
        _h[y1, x1] += hx1y1;
    }
    static void UnLerp(float _val, float _lbd, ref float _v0, ref float _v1)
    {
        _v0 = _val * (1 - _lbd);
        _v1 = _val * _lbd;
    }
    public void WriteVal2HByPos(float[,] _h)
    {
        //if (tfm == null) return;
        

        Vector3 lpos = tfm_pos;// tfm.localPosition;
        UnSampleF(_h, tfm_pos.x, tfm_pos.y, val);
    }
    public void TfmSampleH(float[,] _h, float zbias = 1)
    {
        //if (tfm == null) return;
        //if (!IsValidVal()) return;

        Vector3 _pos = tfm_pos;// tfm.localPosition;
        float hm = SampleF(_h, _pos) + zbias;
        tfm_pos = new Vector3(_pos.x, _pos.y, -hm);
        //tfm.localPosition = new Vector3(_pos.x, _pos.y, -hm);
    }
#endregion


    public void FindFgMatchBg4(HBall[,] self_bg_ball, float sp_err = -1)
    {
        if (sp_err < 0)
            sp_err = ms_match_err_thresh;
        int x0, x1, y0, y1;
        float lbd_x, lbd_y;
        //Vector2 self_pos = new Vector2(col, row);
        Fxy24Ixy(tfm_pos.x, tfm_pos.y, out x0, out x1, out y0, out y1, out lbd_x, out lbd_y);

        int row_bfs = 0, col_bfs = 0;
        ToBFSRoot(y0, x0, self_bg_ball, ref row_bfs, ref col_bfs);
        b_match_bg_dist4[0]= TfmPosToBFSRootDist(tfm_pos, row_bfs, col_bfs);
        match_self_bg4[0] = self_bg_ball[row_bfs, col_bfs];


        ToBFSRoot(y0, x1, self_bg_ball, ref row_bfs, ref col_bfs);
        b_match_bg_dist4[1] = TfmPosToBFSRootDist(tfm_pos, row_bfs, col_bfs);
        match_self_bg4[1] = self_bg_ball[row_bfs, col_bfs];


        ToBFSRoot(y1, x0, self_bg_ball, ref row_bfs, ref col_bfs);
        b_match_bg_dist4[2] = TfmPosToBFSRootDist(tfm_pos, row_bfs, col_bfs);
        match_self_bg4[2] = self_bg_ball[row_bfs, col_bfs];


        ToBFSRoot(y1, x1, self_bg_ball, ref row_bfs, ref col_bfs);
        b_match_bg_dist4[3] = TfmPosToBFSRootDist(tfm_pos, row_bfs, col_bfs);
        match_self_bg4[3] = self_bg_ball[row_bfs, col_bfs];


        dist_fg2bg = FgbAddMatchBgb();
        dist_fg2bg = Mathf.Max(0, dist_fg2bg - sp_err);
        b_no_match = dist_fg2bg > 0;

    }
    public float FgbAddMatchBgb()
    {
        float ret = 200;
        match_self_bg = null;
        for(int i =0; i < b_match_bg_dist4.Length; i++)
        {
            if(ret > b_match_bg_dist4[i])
            {
                ret = b_match_bg_dist4[i];
                match_self_bg = match_self_bg4[i];
            }
            match_self_bg4[i].BgbAddMatchFgb(this);
        }
        dist_fg2bg = ret;
        return ret;
    }
    public void BgbAddMatchFgb(HBall fgb)
    {
        if (match_self_fg == null)
        {
            match_self_fg = fgb;
            return;
        }
        float dy0 = (match_self_fg.tfm_pos.y - row);
        float dx0 = (match_self_fg.tfm_pos.x - col);
        float dy1 = (fgb.tfm_pos.y - row);
        float dx1 = (fgb.tfm_pos.x - col);
        if ((dy0 * dy0 + dx0 * dx0) > (dx1 * dx1 + dy1 * dy1))
            match_self_fg = fgb;
    }
    public bool KeepShapeInRange(float sp_err)
    {
        var p1 = neibor_lod1.GetTarGridPos();
        var p2 = neibor_lod2.GetTarGridPos();
        var pm = Vector3.Lerp(p1, p2, 0.5f);
        pm = pm - tfm_pos;
        pm.z = 0;
        dist_fg_keep_shape = pm.sqrMagnitude;
        dist_fg_keep_shape = Mathf.Max(0, dist_fg_keep_shape - sp_err);
        b_fg_keep_shape = (dist_fg_keep_shape == 0);
        return b_fg_keep_shape;
    }
    public bool CanFg2Bg2FgbMatchback(float sp_err = -1)
    {
        if (sp_err < 0)
            sp_err = ms_match_err_thresh;
        dist_fg2bg2fg = 50;
        if (match_self_bg != null && match_self_bg.match_self_fg != null)
        {
            var _fr = GetInitPos();
            var _to = match_self_bg.match_self_fg.GetInitPos();
            dist_fg2bg2fg = (_fr - _to).sqrMagnitude;
            dist_fg2bg2fg = Mathf.Max(0, dist_fg2bg2fg - sp_err);
        }
        b_can_not_match_back = (dist_fg2bg2fg > 0);
        return !b_can_not_match_back;
    }
    public bool CanBg2Fg2BgMatchback(float sp_err = -1)
    {
        if (sp_err < 0)
            sp_err = ms_match_err_thresh;
        dist_bg2fg2bg = 50;
        if (match_self_fg != null && match_self_fg.match_self_bg != null)
        {
            var _fr = GetInitPos();
            var _to = match_self_fg.match_self_bg.GetInitPos();
            dist_bg2fg2bg = (_fr - _to).sqrMagnitude;
            dist_bg2fg2bg = Mathf.Max(0, dist_bg2fg2bg - sp_err);
        }
        b_can_not_match_back = (dist_bg2fg2bg > 0);
        return !b_can_not_match_back;
    }
    public HBall FindBgBNearestFbgInNeibors(float sp_err = -1)
    {
        dist_bg2fg = 50;
        if (sp_err < 0)
            sp_err = ms_match_err_thresh * 2;
        var self_bg_rc = GetInitPos();
        if (match_self_fg != null)
        {
            Vector3 dir_xy = self_bg_rc - match_self_fg.tfm_pos;
            dir_xy.z = 0;
            dist_bg2fg = (dir_xy).sqrMagnitude;
            if (dist_bg2fg < sp_err)
            {
                b_bg_covered = true;
                return match_self_fg;
            }
        }
        //if (match_self_fg == null)
        foreach (var nbb in neibor_lod0.neibors)
            if (nbb.match_self_fg != null)
            {
                Vector3 dir_xy = self_bg_rc - nbb.match_self_fg.tfm_pos;
                dir_xy.z = 0;
                float curr_fg_bg_dist = (dir_xy).sqrMagnitude;
                if (dist_bg2fg > curr_fg_bg_dist)
                {
                    match_self_fg = nbb.match_self_fg;
                    dist_bg2fg = curr_fg_bg_dist;
                }
            }

        dist_bg2fg = Mathf.Max(0, dist_bg2fg - sp_err);

        b_bg_covered = (dist_bg2fg == 0);
        return match_self_fg;
    }
    
    static float TfmPosToBFSRootDist(Vector3 _pos, int row_root, int col_root)
    {
        float dr = row_root - _pos.y;
        float dc = col_root - _pos.x;
        return dr * dr + dc * dc;
    }
    static float ToBFSRoot(int row, int col, HBall[,] self_bg_ball, ref int row_root, ref int col_root)
    {
        row_root = self_bg_ball[row, col].bfs_root_row;
        col_root = self_bg_ball[row, col].bfs_root_col;
        int dr = row - row_root;
        int dl = col - col_root;//or self_bg_ball[row, col].isvalid()
        return (float)(dr * dr + dl * dl);
    }
    static float CmpMatchBackLoss(float ret, HBall b00, int self_x, int self_y)
    {
        if (b00.IsValidVal())
        {
            float dx = b00.tfm_pos.x - self_x;
            float dy = b00.tfm_pos.y - self_y;
            float loss = dx * dx + dy * dy;
            ret = Mathf.Min(ret, loss);
        }
        return ret;
    }
}
public class HeightMap
{
    public static readonly int sz = 28;
    public static readonly int iw = sz;
    public static readonly int ih = sz;
    public static readonly int iw1 = iw - 1;
    public static readonly int ih1 = ih - 1;
    const float reci_255 = 1.0f / 255.0f;

    public HeightMap(int _sub_hm_id = 0)
    {
        sub_hm_id = _sub_hm_id;
        CreateHB();
    }

    private bool b_init_hb = false;
    public int sub_hm_id = 0;
    public float[,] bgH = new float[ih, iw];
    public float[,] bgH2 = new float[ih, iw];
    public float[,] ini_bgH = new float[ih, iw];
    public bool[,] bgH_checked = new bool[ih, iw];
    public float[,] fgH = new float[ih, iw];
    public float[,] ini_fgH = new float[ih, iw];
    byte[,] ini_bytes_fgH = new byte[ih, iw];
    public float[,] H = new float[ih, iw];
    float[,] dxH = new float[ih, iw];
    float[,] dyH = new float[ih, iw];
    Vector3[,] dxyH = new Vector3[ih, iw];

    float g_rot_ang = 0;
    //Vector3 sub_hm_line_top, sub_hm_line_bot;
    public List<HeightMap> sub_hm_list = new List<HeightMap>();
    public MatchLoss hm_match_loss = new MatchLoss();
    
    public HBall[,] bg_ball = new HBall[ih, iw];
    public HBall[,] fg_ball = new HBall[ih, iw];
    public List<HBall> bg_active_ball = new List<HBall>();
    //List<HBall> bg_active_ball_sub0 = new List<HBall>();
    public List<HBall> fg_active_ball = new List<HBall>();
    public List<HBall> fg_active_ball_lod1 = new List<HBall>();
    public List<HBall> fg_active_ball_lod2 = new List<HBall>();
    public float bg_sca = -3;
    public float fg_sca = 3;
    Vector3 init_bg_active_cw = new Vector3(iw >> 1, ih >> 1, 0);
    Vector3 init_fg_active_cw = new Vector3(iw >> 1, ih >> 1, 0);
    Vector3 thin_fg_active_cw = new Vector3(iw >> 1, ih >> 1, 0);
    float bg_size_mean_dist2bgcw = 0, sub_fg_size_mean_dist2bgcw = 0;
    //Vector3 sub_fg_rt_cw = new Vector3(iw >> 1, ih >> 1, 0);
    Vector3 sub_fg_ini_cw = new Vector3(iw >> 1, ih >> 1, 0);
    public int sub_fg_active_cnt = 0, prev_sub_fg_active_cnt = 0;//exclude fg.sub_id!=owner_sub_id
    public float GetMatchVal()
    {
        return hm_match_loss.GetMatchVal();
    }
    public static bool ValidRC(int _row, int _col)
    {
        return _row >= 0 && _row < ih && _col >= 0 && _col < iw;
    }
    void CreateHB()
    {
        if (b_init_hb)
            return;
        for (int row = 0; row < ih; row++)
            for (int col = 0; col < iw; col++)
            {
                fg_ball[row, col] = new HBall(row, col, 0, sub_hm_id);
                bg_ball[row, col] = new HBall(row, col, 0, sub_hm_id);
            }
        b_init_hb = true;
    }
    //===================
    public class LineHM
    {
        public int top_id, bot_id;
        public LineOrder lod0 = new LineOrder();
        public LineOrder lod1 = new LineOrder();
        //public LineOrder lod_new = new LineOrder();

        public void SwapBgLineOrder()
        {
            int tmp = top_id;
            top_id = bot_id;
            bot_id = tmp;
            LineOrder lo = lod0;
            lod0 = lod1;
            lod1 = lo;
        }
    }
    public LineHM fg_lhm = new LineHM();
    public LineHM bg_lhm = new LineHM();
    static float CalcCurvature(ref Vector2 dt, Vector2 dt_s)
    {
        float dot_v = Vector2.Dot(dt_s, dt);
        if (dot_v < 0)
        {
            dt = -dt;
            dot_v = -dot_v;
        }
        float dt_len2 = dt.sqrMagnitude;
        if (dt_len2 > 0)
            return dot_v / dt_len2;
        return 0;
    }
    public static Vector2 ReverseTanKeepNorm(Vector2 base_dir, Vector2 new_dir)
    {
        Vector2 dir = base_dir.normalized;
        float dot_v = Vector2.Dot(dir, new_dir);
        Vector2 dtan = dir * dot_v;
        Vector2 dnor = new_dir - dtan;
        return dnor - dtan;
    }
    static void CalcTan(HBall hb, bool init_or_rt)
    {
        hb.dir_t = AdjTan(hb, Vector2.up, init_or_rt, false);
        hb.dir_t = AdjTan(hb, hb.dir_t, init_or_rt, true);
        //if (hb.fg_or_bg && hb.row == 21 && hb.col == 11)
        //{
        //    Debug.Log(hb.RCStr() + ", dot = " + dot_t + "/" + dot_t2 + "\n\t tan="+ tan.ToString("F2")+"/"+tan2.ToString("F2"));
        //}
        
        hb.dir_n = HBall.Rot90(hb.dir_t);
        if (Vector2.Dot(hb.dir_n, hb.dir_s) < 0)
            hb.dir_n = -hb.dir_n;
    }
    static Vector2 AdjTan(HBall hb, Vector2 prev_tan, bool init_or_rt, bool refine)
    {
        Vector2 tan = prev_tan.normalized;
        List<Vector2> candidates_tan = new List<Vector2>();
        candidates_tan.Add(tan);
        int N = 1;
        if(refine)
        {
            candidates_tan.Add(HBall.RotN(tan, 4));
            candidates_tan.Add(HBall.RotN(tan, 8));
            candidates_tan.Add(HBall.RotN(tan, 12));
            candidates_tan.Add(HBall.RotNb(tan, 4));
            candidates_tan.Add(HBall.RotNb(tan, 8));
            candidates_tan.Add(HBall.RotNb(tan, 12));
        }
        else//coarse
        {
            candidates_tan.Add(HBall.RotN(tan, 18));
            candidates_tan.Add(HBall.RotN(tan, 36));
            candidates_tan.Add(HBall.RotNb(tan, 18));
            candidates_tan.Add(HBall.RotNb(tan, 36));
            N = candidates_tan.Count;
            for (int i = 0; i < N; i++)
                candidates_tan.Add(HBall.Rot90(candidates_tan[i]));
        }
        List<float> sum_t = new List<float>();
        List<float> neg_dot = new List<float>();
        List<float> pos_dot = new List<float>();
        N = candidates_tan.Count;
        for (int i = 0; i < N; i++)
        {
            sum_t.Add(0);
            neg_dot.Add(0);
            pos_dot.Add(0);
        }
        
        Vector2 cpos = hb.PosV2(init_or_rt);
        hb.dir_s = Vector2.zero;
        foreach (var nhb in hb.neibor_lod0.neibors)
        {
            if (!nhb.IsValidVal())
                continue;
            Vector2 npos = nhb.PosV2(init_or_rt);
            Vector2 dir = npos - cpos;
            hb.dir_s += dir;
            for (int i = 0; i < N; i++)
            {
                float dv = Vector2.Dot(dir, candidates_tan[i]);
                if (dv > 0)
                {
                    pos_dot[i] += dv;
                    sum_t[i] += dv;
                }
                else
                {
                    neg_dot[i] -= dv;
                    sum_t[i] -= dv;
                }
            }
        }
        float max_d = 0;
        int max_i = 0;
        for (int i = 0; i < N; i++)
        {
            float cmp_v = (neg_dot[i] + 1e-5f) / (pos_dot[i] + 1e-5f);
            if (cmp_v > 1) cmp_v = 1.0f / cmp_v;
            float sca_t = HBall.X01Y01(cmp_v, 0.3f, 1.0f, 1, 3);
            float proj_v = sca_t * sum_t[i];
            if(max_d < proj_v)//sum_t[i])
            {
                max_d = proj_v;// sum_t[i];
                max_i = i;
            }
        }
        return candidates_tan[max_i].normalized;
    }
    public int CalcAllPointTanNorm(List<HBall> _active_ball, bool init_or_rt = true)
    {
        int max_id = -1;
        float max_curve = 0;
        foreach (var hb in _active_ball)
        {
            if (!hb.IsValidVal())
                continue;
            CalcTan(hb, init_or_rt);
            hb.curvature = Mathf.Abs(Vector2.Dot(hb.dir_t, hb.dir_s));
            if (hb.curvature > max_curve)
            {
                max_curve = hb.curvature;
                max_id = hb.UID();
            }
        }
        return max_id;
    }
    void SmoothAllPointTanNorm(List<HBall> _active_ball)
    {
        int self_w = 2;
        foreach (var hb in _active_ball)
        {
            if (!hb.IsValidVal())
                //if (!hb.IsValGreater())
                continue;
            Vector2 sum_dt = hb.dir_t * self_w;
            float dt_cnt = self_w;
            for(int i = 0; i < hb.neibor_lod0.neibors.Count; i++)
            {
                var nhb = hb.neibor_lod0.neibors[i];
                float _wi = hb.neibor_lod0.neibors_w[i];
                if (!nhb.IsValidVal())
                    continue;
                if (Vector2.Dot(sum_dt, nhb.dir_t) < 0)
                //if (Vector2.Dot(hb.dir_t, nhb.dir_t) < 0)
                    sum_dt -= nhb.dir_t * _wi;
                else
                    sum_dt += nhb.dir_t * _wi;
                dt_cnt +=_wi;
            }
            //if (dt_cnt > 0)
            hb.dir_t_mean = (sum_dt / dt_cnt).normalized;

        }
        foreach (var hb in _active_ball)
        {
            if (!hb.IsValidVal())
                //if (!hb.IsValGreater())
                continue;
            //bool flag = sub_hm_id == 0 && hb.fg_or_bg && hb.row == 21 && hb.col == 11;
            //if(flag) Debug.Log("smooth "+ hb.RCStr() +", dir_t:" + hb.dir_t_mean.ToString("F2") + "/" + hb.dir_t.ToString("F2"));
            hb.dir_t = hb.dir_t_mean;
            hb.dir_n = HBall.Rot90(hb.dir_t);
            if (Vector2.Dot(hb.dir_n, hb.dir_s) < 0)
                hb.dir_n = -hb.dir_n;
        }
    }
    public void BuildLineOrder(HBall[,] _balls, List<HBall> _active_ball, ref LineHM lhm)
    {
        lhm.bot_id = -1;
        lhm.top_id = CalcAllPointTanNorm(_active_ball);
        if (lhm.top_id < 0)
            return;
        SmoothAllPointTanNorm(_active_ball);

        BFSSortLineOrder(_balls, ref lhm.lod0, lhm.top_id, 0);
        //TODO check line_cnt > 0.5*total_cnt
        lhm.bot_id = lhm.lod0.LastId();
        BFSSortLineOrder(_balls, ref lhm.lod1, lhm.bot_id, 1);
        lhm.top_id = lhm.lod1.LastId();
        BFSSortLineOrder(_balls, ref lhm.lod0, lhm.top_id, 0);
    }
    void BFSSortLineOrder(HBall[,] _balls, ref LineOrder lod, int start_id, int reverses = 0)
    {
        lod.Clear();
        int row = 0, col = 0;
        HBall.I2RC(start_id, ref row, ref col);
        lod.Add(_balls[row, col], null, 0, reverses);
        BFSSortLineOrderContinue(ref lod, reverses);
    }
    static Vector2 HBNeiborCW(HBall chb, int max_dist = 3, int max_cnt = 12)
    {
        Vector2 ret = Vector2.zero;
        HashSet<int> hs = new HashSet<int>();
        Vector2 self_xy2 = chb.GetXY2();
        ret += self_xy2;
        int add_cnt = 1;
        hs.Add(chb.UID());
        List<HBall> hb_list = new List<HBall>();
        hb_list.Add(chb);
        int head = 0;
        while (head < hb_list.Count)
        {
            foreach (var nhb in hb_list[head].neibor4_ignore)
            {
                if (nhb != null && !hs.Contains(nhb.UID()))
                {
                    ret += nhb.GetXY2();
                    add_cnt++;
                    hs.Add(nhb.UID());
                    hb_list.Add(nhb);
                    if (hb_list.Count >= max_cnt)
                        break;
                }
            }
            head++;
        }
        return ret / add_cnt;
    }
    void BFSSortLineOrderLimitByBGCurve(ref LineOrder lo, int pre_cnt = 10, float dist_limit = 3, float angle_limit = 20)
    {
        int ini_cnt = lo.Count;
        lo.id_list[lo.Count-1].temp_id = 1001;
        lo.id_list[0].temp_id = 1000;
        int head = 0;
        while (head < lo.Count)
        {
            HBall chb = lo.id_list[head];
            int cur_dep = lo.depth_list[head];
            foreach (var nhb in chb.neibor4_ignore)
            {
                if (nhb == null || nhb.val < HBall.ms_active_thresh || nhb.match_self_bg == null)
                    continue;
                int nxt_id = nhb.UID();// RC2I(nxt_r, nxt_c);
                if (lo.hs.Contains(nxt_id))
                    continue;
                Vector2 cw_tan = nhb.dir_t;
                int fr = head - ini_cnt;
                if(fr > 0)
                {
                    Vector2 nhb_cw = HBNeiborCW(nhb);
                    Vector2 fr_cw = HBNeiborCW(lo.id_list[fr]);
                    cw_tan = (fr_cw - nhb_cw).normalized;
                }
                //int to = lo.Count;
                //int fr = to - ini_cnt;
                //Vector2 cw_prev10 = Vector2.zero;
                //for (int k = fr; k < to; k++)
                //    cw_prev10 += lo.id_list[k].GetXY2();
                //cw_prev10 /= (to - fr);
                //Vector2 cw_tan = (nhb.GetXY2() - cw_prev10).normalized;
                
                Vector2 bg_tan = nhb.match_self_bg.dir_t;

                float dot_tan = Mathf.Abs(Vector2.Dot(cw_tan, bg_tan));
                float ang_btw = Mathf.Acos(dot_tan) * Mathf.Rad2Deg;
                if (ang_btw > 30)
                {
                    //if (sub_hm_id == 1)
                    //    Debug.Log("push [" + nhb.row + ", " + nhb.col 
                    //        + "], angle-btw:" + ang_btw.ToString("F1")
                    //        + ", cw_tan:" + cw_tan.ToString("F3")
                    //        + ", bg_tan:" + bg_tan.ToString("F3")
                    //                                    );
                    continue;
                }
                //if (sub_hm_id == 2 && nhb.row == 12 && nhb.col == 16)
                //{
                //    Debug.Log("push [" + nhb.row + ", " + nhb.col 
                //        + "], angle=" + ang_btw.ToString("F1")
                //        + ", cw_tan:"+cw_tan.ToString("F3")
                //        + ", bg_tan:"+bg_tan.ToString("F3")
                //        +", bg.xy2="+ nhb.match_self_bg.RCStr()
                //        );
                //    for (int k = fr; k < to; k++)
                //        lo.id_list[k].temp_id = 22;
                //}
                lo.Add(nhb, chb, cur_dep + 1);
            }


            head++;
        }
    }
    void BFSSortLineOrderContinue(ref LineOrder lo, int reverses )
    {
        //*
        int head = lo.id_list.Count - 1;
        while (head < lo.id_list.Count)
        {
            //if (limit_cnt > 0 && limit_cnt <= lo.Count)
            //    break;
            HBall chb = lo.id_list[head];
            int cur_dep = lo.depth_list[head];
            Vector2 ctan = chb.dir_t;//[row, col];

            foreach (var nhb in chb.neibor4_ignore)
            {
                if (nhb == null || !nhb.IsValGreater())
                    continue;
                if (!nhb.IsValidVal())
                //if (!ignore_sub_id && !nhb.IsValidVal())
                    continue;
                int nxt_id = nhb.UID();// RC2I(nxt_r, nxt_c);
                if (!lo.hs.Contains(nxt_id))
                {
                    float tan_dot = Vector2.Dot(ctan, nhb.dir_t);//[nxt_r, nxt_c]);
                    if (tan_dot < 0)
                        nhb.dir_t = -nhb.dir_t;
                    //nhb.dir_t = Vector2.Lerp(nhb.dir_t, ctan, 0.3f).normalized;
                    nhb.dir_n = HBall.Rot90(nhb.dir_t);
                    if (Vector2.Dot(nhb.dir_n, nhb.dir_s) < 0)//make norm toward line-mid
                        nhb.dir_n = -nhb.dir_n;
                    lo.Add(nhb, chb, cur_dep + 1, reverses);
                }
            }
            head++;
        }//*/
    }

    //===================
    public void UpdateFgSubId(int _row, int _col, int _sub_id, Vector3 _tfm_pos)
    {
        HBall hb = fg_ball[_row, _col];
        hb.sub_id = _sub_id;
        hb.tfm_pos = _tfm_pos;
        if (hb.IsValidVal())
        {
            sub_fg_active_cnt++;
        }
    }
    void FgActiveLod0GenLod12()
    {
        fg_active_ball_lod1.Clear();
        fg_active_ball_lod2.Clear();
        foreach (var fgb in fg_active_ball)
        {
            if (fgb.IsValidVal())
            {
                if ((1 == (fgb.row & 0x1)) && (1 == (fgb.col & 0x1)))
                    fg_active_ball_lod1.Add(fgb);
                if (((fgb.row % 3) == 2) && (2 == (fgb.col % 3)))
                    fg_active_ball_lod2.Add(fgb);
            }
        }
    }
    void BuildFgbNeibor(HBall fgb)
    {
        fgb.neibor_lod0.BFSBuildNeibor(fgb);
        //fgb.neibor_lod0.BuildNeibor(fgb, fg_active_ball, 3);
        fgb.neibor_lod1.BuildNeibor(fgb, fg_active_ball_lod1, 9, 1);
        fgb.neibor_lod2.BuildNeibor(fgb, fg_active_ball_lod2, 16, 2);
    }
    public void SubHMClearActiveFgb()
    {
        sub_fg_active_cnt = 0;
        foreach (var fgb in fg_active_ball)
        {
            fgb.sub_id = 0;
            fgb.is_active = false;
        }
    }
    //*
    public void ExpandSubHMLineHead()
    {
        foreach (var sub_hm in sub_hm_list)
        {
            sub_hm.SubHMReBuildFg();
            sub_hm.ExpandFgLine();
            sub_hm.SubHMReBuildFg();
        }
    }
    public void ExpandFgLine(int expand_cnt = 2)
    {
        LineOrder lo_new = new LineOrder();
        if (fg_lhm.lod0.Count > 0)
        {
            lo_new.Add(fg_lhm.lod0.id_list[0], null, 0);
            lo_new.Add(fg_lhm.lod0.id_list[fg_lhm.lod0.Count-1], null, 0);
        }
        if (fg_lhm.lod1.Count > 0)
        {
            lo_new.Add(fg_lhm.lod1.id_list[0], null, 0);
            lo_new.Add(fg_lhm.lod1.id_list[fg_lhm.lod1.Count - 1], null, 0);
        }
        int head = 0;
        while(lo_new.Count > head)
        {
            HBall chb = lo_new.id_list[head];
            int cur_dep = lo_new.depth_list[head];
            if (cur_dep >= expand_cnt)
                break;
            foreach(var nhb in chb.neibor8_ignore)
            {
                if (nhb == null || lo_new.hs.Contains(nhb.UID()) || nhb.sub_id != 0)
                    continue;
                if (CheckInLineRange(nhb, chb, 3, -0.5f, 2))
                {
                    lo_new.Add(nhb, chb, cur_dep + 1);
                    nhb.dir_t = chb.dir_t;
                    nhb.dir_n = chb.dir_n;
                    nhb.sub_id = sub_hm_id;
                }
            }
            head++;
            
        }
        //foreach (var fgb in fg_active_ball)
        //{
        //    if (fgb.IsValGreater() && !fgb.IsValidVal())
        //}
    }//*/
    public void SubHMReBuildFg()
    {
        FgActiveLod0GenLod12();
        foreach (var fgb in fg_active_ball)
        {
            fgb.UpdateNeibor4(fg_ball);
            fgb.dep_id = fgb.copy_dep_id = fgb.rvs_dep_id = -1;
        }
        foreach (var fgb in fg_active_ball)
            if (fgb.IsValidVal())
                BuildFgbNeibor(fgb);
        BuildLineOrder(fg_ball, fg_active_ball, ref fg_lhm);
        BuildFGMapBG();
    }
    
    public void SetFgH(byte[,] _b)
    {
        fg_active_ball.Clear();

        init_fg_active_cw = Vector3.zero;
        sub_fg_active_cnt = 0;
        float sum_val = 0;
        fg_sca = 1;
        for (int row = 0; row < ih; row++)
            for (int col = 0; col < iw; col++)
            {
                ini_bytes_fgH[row, col] = _b[row, col];
                float _val = _b[row, col] * reci_255;
                ini_fgH[row, col] = fgH[row, col] = _val;
                HBall rchb = fg_ball[row, col];
                rchb.Reset(_val);
                rchb.fg_or_bg = true;
                rchb.sub_id = 0;
                if (_val > HBall.ms_active_thresh)
                {
                    sum_val += _val;
                    fg_active_ball.Add(rchb);
                    init_fg_active_cw += new Vector3(col, row, 0);
                    sub_fg_active_cnt++;
                }
            }
        if (fg_active_ball.Count > 0)
        {
            //fg_sca = 100.0f / sum_val;
            init_fg_active_cw /= fg_active_ball.Count;
            FgActiveLod0GenLod12();
        }
        //if (sub_hm_id > 0)
        //    Debug.Log("sub id=" + sub_hm_id.ToString()
        //        + ", fg ball cnt :" + fg_active_ball.Count + ",(" + fg_active_ball_lod1.Count + ", " + fg_active_ball_lod2.Count + ")"
        //        + ", sca = " + fg_sca + ", sum_val =" + sum_val.ToString("F1")
        //        + ", init fg active cw="+ init_fg_active_cw.ToString("F1")
        //        );
        foreach (var fgb in fg_active_ball)
            fgb.UpdateNeibor4(fg_ball);
        foreach (var fgb in fg_active_ball)
        {
            fgb.ini_dir2global_cw = fgb.GetInitPos() - init_fg_active_cw;
            BuildFgbNeibor(fgb);
        }
        //BuildLineOrder(fg_ball, fg_active_ball, ref fg_lhm);
        ThinHBList(fg_active_ball, true);
        //Debug.Log("set fgh, thin");
        //AlignInitCW();
            BuildFGMapBG();
        GetH();
        for (int i = 0; i < sub_hm_list.Count; i++)
        {
            sub_hm_list[i].SetFgH(ini_bytes_fgH);
        }
    }
    
    public void SetBgH(byte[,] _b, List<byte[,]> sub_bytes = null)
    {
        CreateSubBg(sub_bytes);
        bg_active_ball.Clear();
        init_bg_active_cw = Vector3.zero;

        float sum_val = 0;
        bg_sca = -1;
        for (int row = 0; row < ih; row++)
            for (int col = 0; col < iw; col++)
            {
                float _val = _b[row, col] * reci_255;
                ini_bgH[row, col] = bgH2[row, col] = bgH[row, col] = _val;
                bgH_checked[row, col] = false;
                HBall rchb = bg_ball[row, col];
                rchb.fg_or_bg = false;
                rchb.Reset(_val);
                if (sub_bytes != null)
                {
                    for (int k = 0; k < sub_bytes.Count; k++)
                    {
                        if (sub_bytes[k][row, col] > 25)
                            rchb.sub_id |= (1 << k);
                    }
                }
                else
                    rchb.sub_id = sub_hm_id;
                if (rchb.IsValidVal())
                {
                    bgH_checked[row, col] = true;
                    sum_val += _val;
                    bg_active_ball.Add(rchb);
                    init_bg_active_cw += new Vector3(col, row, 0);
                }
            }
        BlurModifyBGH(sum_val);
        BFSModifyBGH();
        if (bg_active_ball.Count > 0)
        {
            //bg_sca = -100.0f / sum_val;// / (bg_active_ball.Count);
            init_bg_active_cw /= bg_active_ball.Count;

            bg_size_mean_dist2bgcw = 0;
            foreach (var bgb in bg_active_ball)
            {
                bgb.ini_dir2global_cw = bgb.GetInitPos() - init_bg_active_cw;
                bg_size_mean_dist2bgcw += bgb.ini_dir2global_cw.magnitude;
                bgb.UpdateNeibor4(bg_ball);
            }
            foreach (var bgb in bg_active_ball)
                bgb.neibor_lod0.BFSBuildNeibor(bgb);
                //bgb.neibor_lod0.BuildNeibor(bgb, bg_active_ball);
            BuildLineOrder(bg_ball, bg_active_ball, ref bg_lhm);
            
            bg_size_mean_dist2bgcw *= (1.0f / bg_active_ball.Count);
        }
        //Debug.Log("bg ball cnt :" + bg_active_ball.Count 
        //    + ", sca = " + bg_sca 
        //    + ", sum_val =" + sum_val.ToString("F1")
        //    +"\n\t [12,13].val="+bg_ball[12,13].val.ToString("F2")+", bfs:["+bg_ball[12,13].bfs_root_row+","+bg_ball[12,13].bfs_root_col+"]"
        //    +"\n\t [13,13].val="+bg_ball[13,13].val.ToString("F2")+", bfs:["+bg_ball[13,13].bfs_root_row+","+bg_ball[13,13].bfs_root_col+"]"
        //    +"\n\t [12,14].val="+bg_ball[12,14].val.ToString("F2")+", bfs:["+bg_ball[12,14].bfs_root_row+","+bg_ball[12,14].bfs_root_col+"]"
        //    +"\n\t [13,14].val="+bg_ball[13,14].val.ToString("F2")+", bfs:["+bg_ball[13,14].bfs_root_row+","+bg_ball[13,14].bfs_root_col+"]"
        //    );
        AlignInitCW();
        GetH();
    }
    public void SubHMFindLine()
    {
        if (bg_active_ball.Count == 0)
            return;
        foreach(var bgb in bg_active_ball)
        {

        }
    }
    void CreateSubBg(List<byte[,]> sub_bytes)
    {
        sub_hm_list.Clear();
        if (sub_bytes == null || sub_bytes.Count == 0)
            return;
        for (int i = 0; i < sub_bytes.Count; i++)
        {
            HeightMap sub_hm = new HeightMap(1 << i);
            sub_hm.SetBgH(sub_bytes[i], null);
            sub_hm.SubHMFindLine();
            //sub_hm.SetFgH(ini_bytes_fgH);
            sub_hm_list.Add(sub_hm);
        }
    }
    void BlurModifyBGH(float sum_h)
    {
        //bgH-->bgH2
        //for (int row = 1; row < ih1; row++)
        //for (int col = 1; col < iw1; col++)
        float new_sum = 0;
        foreach (var bgb in bg_active_ball)
        {
            //if (!bg_ball[row, col].IsValidVal())
            //continue;
            float sval = 0;// bgH[row, col];
            for (int oi = -1; oi < 2; oi++)
                for (int oj = -1; oj < 2; oj++)
                {
                    if (!ValidRC(bgb.row + oi, bgb.col + oj))
                        continue;
                    HBall b0 = bg_ball[bgb.row + oi, bgb.col + oj];
                    if (b0.IsValidVal())
                    {
                        float dist = oi * oi + oj * oj;
                        float wdist = 1.0f / (1 + dist * 2);//1,1/3,1/5
                        sval += bgH[b0.row, b0.col] * wdist;
                    }
                }

            bgH2[bgb.row, bgb.col] = sval;
            new_sum += sval;
        }
        if (new_sum == 0)
            return;
        float vsca = 1.1f * sum_h / new_sum;
        foreach (var bgb in bg_active_ball)
            bgH[bgb.row, bgb.col] = vsca * bgH2[bgb.row, bgb.col];
        
    }
    void BFSModifyBGH(float decayH = 0.3f)
    {
        List<HBall> temp_hb = new List<HBall>();
        foreach (var b in bg_active_ball)
            temp_hb.Add(b);
        int ihead = 0;
        while (ihead < temp_hb.Count)
        {
            HBall hb = temp_hb[ihead];
            int cx = hb.col;
            int cy = hb.row;
            float hv = bgH[cy, cx];
            float aug_dep = ihead < bg_active_ball.Count ? 8 * decayH : decayH;
            if (cy > 0)
            {
                int nx = cx;
                int ny = cy - 1;
                BFSAddNewHB(temp_hb, hb, nx, ny, hv, aug_dep);
            }
            if (cx > 0)
            {
                int nx = cx - 1;
                int ny = cy;
                BFSAddNewHB(temp_hb, hb, nx, ny, hv, aug_dep);
            }
            if (cy < ih1)
            {
                int nx = cx;
                int ny = cy + 1;
                BFSAddNewHB(temp_hb, hb, nx, ny, hv, aug_dep);
            }
            if (cx < iw1)
            {
                int nx = cx + 1;
                int ny = cy;
                BFSAddNewHB(temp_hb, hb, nx, ny, hv, aug_dep);
            }
            ihead++;
        }
        //Debug.Log("bfs cnt=" + temp_hb.Count + ", ihead=" + ihead);
    }
    void BFSAddNewHB(List<HBall> temp_hb, HBall hb, int nx, int ny, float hv, float decayH = 0.3f)
    {
        if (!bgH_checked[ny, nx])
        {
            bgH_checked[ny, nx] = true;
            bgH[ny, nx] = hv - decayH;
            HBall nb = bg_ball[ny, nx];
            nb.bfs_root_col = hb.bfs_root_col;
            nb.bfs_root_row = hb.bfs_root_row;
            temp_hb.Add(nb);
        }
    }
    public bool GuessGlobalRotAngle()
    {
        if (bg_active_ball.Count == 0 || fg_active_ball.Count == 0)
            return false;
        float mean_ang = 0;
        float rot_cnt = 0;
        float cx = init_fg_active_cw.x;
        float cy = init_fg_active_cw.y;
        for (int i = 0; i < fg_active_ball.Count; i++)
        {
            HBall fgb = fg_active_ball[i];
            if (fgb.b_no_match || fgb.b_fg_break_sca) continue;
            Vector3 dir0 = new Vector3(fgb.col - cx, fgb.row - cy, 0);
            Vector3 dir1 = new Vector3(fgb.tfm_pos.x - cx, fgb.tfm_pos.y - cy, 0);
            mean_ang += HBall.GetAngleBetween(dir0, dir1);
            rot_cnt++;
        }
        if (rot_cnt == 0)
            return false;

        mean_ang /= rot_cnt;
        float delta_ang = mean_ang - g_rot_ang;
        if (Mathf.Abs(delta_ang) < 5)
            return false;

        Debug.Log("global rot angle:" + g_rot_ang.ToString("F1") + "---->" + mean_ang.ToString("F1"));
        g_rot_ang = mean_ang;
        Quaternion dq = Quaternion.AngleAxis(mean_ang, Vector3.forward);
        for (int i = 0; i < fg_active_ball.Count; i++)
        {
            var b = fg_active_ball[i];
            b.b_no_match = false;
            b.tfm_pos = init_bg_active_cw + dq * b.ini_dir2global_cw;
        }
        return true;
    }
    public void AlignInitCW()
    {
        g_rot_ang = 0;
        for (int row = 0; row < ih; row++)
            for (int col = 0; col < iw; col++)
            {
                fgH[row, col] = ini_bgH[row, col];

            }
        foreach (var bgb in bg_active_ball)
            bgb.match_self_fg = null;
        if (bg_active_ball.Count > 0 && fg_active_ball.Count > 0)
        {
            for (int i = 0; i < fg_active_ball.Count; i++)
            {
                var b = fg_active_ball[i];
                b.b_no_match = false;
                b.tfm_pos = init_bg_active_cw + b.ini_dir2global_cw;
            }
            BuildFGMapBG();
        }
    }
    public float FindBestBiasFGBByInitCW(List<HBall> _fg_list, int search_rng = 2)
    {
        float max_cover_bg_val = 0;
        int max_oi = 0, max_oj = 0;
        
        for (int oi = -search_rng; oi <= search_rng; oi+=2)
            for (int oj = -search_rng; oj <= search_rng; oj+=2)
            {
                float cover_bg_val = ApplyFGBByBias(_fg_list,  oi, oj);
                if (cover_bg_val > max_cover_bg_val)
                {
                    max_cover_bg_val = cover_bg_val;
                    max_oi = oi;
                    max_oj = oj;
                }
            }
        //if (sub_hm_id == 2)
        //Debug.Log("sub["+sub_hm_id +"], max-cover-cnt=" + max_cover_bg_val + "/"+fg_active_ball.Count);
        float ret = ApplyFGBByBias(_fg_list, max_oi, max_oj);
        return ret;
    }
    public float ApplyFGBByBias(List<HBall> _tar_fg_list,int oi, int oj, float sca = 1.0f)
    {
        float cover_bg_val = 0;
        Vector3 _bias = new Vector3(oi, oj, 0);
        for (int i = 0; i < _tar_fg_list.Count; i++)
        {
            var tar_b = _tar_fg_list[i];
            var b = fg_active_ball[i];
            //if (b.val < HBall.ms_active_thresh)
            //    continue;
            b.tfm_pos = _bias + tar_b.tfm_pos;
            b.FindFgMatchBg4(bg_ball);
            if (!b.b_no_match)
                cover_bg_val++;
        }
        return cover_bg_val;
    }
    public float[,] GetH()
    {
        for (int row = 0; row < ih; row++)
            for (int col = 0; col < iw; col++)
            {
                float _hval = fg_sca* fgH[row, col] + bg_sca * bgH[row, col];
                H[row, col] = _hval;
            }
        return H;
    }
    void UpdateDxy(float[,] harr)
    {
        for (int row = 0; row < ih; row++)
            for (int col = 0; col < iw; col++)
            {
                float dx = HBall.SampleI(harr, col + 1, row) - HBall.SampleI(harr, col - 1, row);
                float dy = HBall.SampleI(harr, col, row + 1) - HBall.SampleI(harr, col, row - 1);
                dxH[row, col] = dx;
                dyH[row, col] = dy;
                dxyH[row, col] = new Vector3(dx, dy, 0);
            }
    }
    public void ApplyBallByH(List<HBall> blist, float[,] _h)
    {
        UpdateHbyBall();                //Ball mass/pos drive--->fgH
        foreach (var b in blist)
            //if(b.IsValidVal())
                b.TfmSampleH(_h);
    }
    public void BallFallDown(float fall_speed = 0.5f)
    {
        UpdateDxy(H);
        float speed_sca = -fall_speed;
        foreach (var b in fg_active_ball)
        {
            if (!b.IsValidVal()) continue;
            Vector3 _lclp = b.tfm_pos;// tfm.localPosition;
            float dx = HBall.SampleF(dxH, _lclp);
            float dy = HBall.SampleF(dyH, _lclp);

            b.tfm_pos += new Vector3(dx, dy, 0) * speed_sca;
            //b.tfm.localPosition += new Vector3(dx, dy, 0) * speed_sca;
        }
    }
    public void UncoveredBGDragFGB(float itr_lbd)
    {
        foreach (var bgb in bg_active_ball)
        {
            if (bgb.b_bg_covered)
                continue;
            
            var tar_fgb = bgb.match_self_fg;
            if (tar_fgb != null)
            {
                tar_fgb.tfm_pos = Vector3.Lerp(tar_fgb.tfm_pos, bgb.GetInitPos(), 0.25f* itr_lbd);
                foreach(var tar_fgb_nb in tar_fgb.neibor4)
                    if(tar_fgb_nb != null)
                        tar_fgb_nb.tfm_pos = Vector3.Lerp(tar_fgb_nb.tfm_pos, bgb.GetInitPos(), 0.15f* itr_lbd);
            }
        }
    }
    public void UnmatchFGBDragByCW(float drag_unmatch_lbd = 0.3f)
    {
        foreach (var b in fg_active_ball)
            b.UnmatchFGBDragByCW(drag_unmatch_lbd);
    }
    public void BallSortGrid(float sort_grid_speed = 0.5f)
    {
        foreach (var b in fg_active_ball)
            if(b.IsValidVal())
                b.SortGrid(sort_grid_speed);
    }
    public void UpdateHbyBall()
    {
        for (int row = 0; row < ih; row++)
            for (int col = 0; col < iw; col++)
                fgH[row, col] = 0;//clear all 2 zero
        foreach (var b in fg_active_ball)
            if(b.IsValidVal())
                b.WriteVal2HByPos(fgH);
        GetH();
        //foreach (var sum_hm in sub_hm_list)
            //sum_hm.UpdateHbyBall(); ;
    }
    public void DoFindBgBNearestFbgInNeibors()
    {
        for (int k = 0; k < 2; k++)
            foreach (var bgb in bg_active_ball)
                bgb.FindBgBNearestFbgInNeibors();
    }
    public void BuildFGMapBG(bool skip_check_valid = false)
    {
        foreach (var bgb in bg_active_ball)
        {
            bgb.match_self_fg = null;
            bgb.b_bg_covered = false;
        }

        for (int i = 0; i < fg_active_ball.Count; i++)
        {
            HBall cfgb = fg_active_ball[i];
            if (skip_check_valid || cfgb.is_active)
                cfgb.FindFgMatchBg4(bg_ball);

        }
        DoFindBgBNearestFbgInNeibors();
    }
    public void DetSubIdByTfmPos(float dist_err = -1)
    {
        foreach (var fgb in fg_active_ball)
        {
            fgb.FindFgMatchBg4(bg_ball, dist_err);
            fgb.sub_id = (fgb.b_no_match) ? 0 : sub_hm_id;
        }
    }
    void SubBg2Fb(Vector2 fg_cm, Vector2 bg_cm, Vector4 abcd)
    {
        foreach (var bgb in bg_active_ball)
        {
            Vector3 bgb_pos = new Vector3(bgb.col, bgb.row, 0);
            float dcol = (bgb.col - bg_cm.x);
            float drow = (bgb.row - bg_cm.y);
            float map_col = abcd.x * dcol + abcd.y * drow + fg_cm.x;
            float map_row = abcd.z * dcol + abcd.w * drow + fg_cm.y;
            int r0 = Mathf.RoundToInt(map_row);
            int c0 = Mathf.RoundToInt(map_col);
            if (r0 < 0) r0 = 0;
            if (c0 < 0) c0 = 0;
            int r1 = r0 + 1;
            int c1 = c0 + 1;
            if (r0 > ih1) r0 = ih1;
            if (c0 > iw1) c0 = iw1;
            if (r1 > ih1) r1 = ih1;
            if (c1 > iw1) c1 = iw1;
        }
    }
    bool CheckInLineRange(HBall fgb, HBall bgb, float l_tan = 4, float l_nor_min = -0.6f, float l_nom_max=2)
    {
        Vector2 dir = new Vector2(fgb.tfm_pos.x - bgb.tfm_pos.x, fgb.tfm_pos.y - bgb.tfm_pos.y);
        //if (Mathf.Abs(bgb.dir_t.sqrMagnitude - 1) > 0.01f)
        //    Debug.Log(bgb.RCStr() + ", dir_t:" + bgb.dir_t.ToString("F3"));
        //if (Mathf.Abs(bgb.dir_n.sqrMagnitude - 1) > 0.01f)
        //    Debug.Log(bgb.RCStr() + ", dir_n:" + bgb.dir_n.ToString("F3"));
        float dist = dir.magnitude;
        float dot_tan = Vector2.Dot(dir, bgb.dir_t);
        float dot_nor = Vector2.Dot(dir, bgb.dir_n);
        //bool ret = Mathf.Abs(dot_tan) < l_tan && dot_nor > l_nor_min && dot_nor < l_nom_max;
        //if ((sub_hm_id == 4 && fgb.row == 9 && fgb.col == 9)//[8,8]->bgb[5,6]
        // //|| (sub_hm_id == 2 && fgb.row == 6 && fgb.col == 5)
        // )
        //    Debug.Log(fgb.RCStr()
        //        + " --> " + bgb.RCStr() + (ret ? " InRange" : " OutRange") + ", dist=" + dist
        //        + ",\n\t dot-tan:" + dot_tan.ToString("F2") + ", bgb.dir_t:" + bgb.dir_t.ToString("F3")
        //        + "\n\t dot-nor:" + dot_nor.ToString("F2") + ", bgb.dir_n:" + bgb.dir_n.ToString("F3")
        //        //+ "\n\t bg." + bgb.neibor4[2].RCStr() + ".dir_t=" + bgb.neibor4[2].dir_t.ToString("F3")
        //        );
        //Vector2 dir_line_tan = bgb.dir_t * dot_tan;
        //Vector2 dir_line_nor = dir - dir_line_tan;
        return Mathf.Abs(dot_tan) < l_tan && dot_nor > l_nor_min && dot_nor < l_nom_max;
    }
    bool SubCalcFg2BgLT(ref Vector4 fg2bg_trfm, ref Vector2 fg_cm, ref Vector2 bg_cm, float err_thresh=4.0f)
    {
        List<Vector2> bg_list = new List<Vector2>();
        List<Vector2> fg_list = new List<Vector2>();
        foreach (var fgb in fg_active_ball)
            if (fgb.match_self_bg != null && fgb.dist_fg2bg < err_thresh)
            {
                fg_list.Add(new Vector2(fgb.col, fgb.row));
                bg_list.Add(new Vector2(fgb.match_self_bg.col, fgb.match_self_bg.row));
            }

        MathUtils.Calc_abcd_xy_E_uv(fg_list, bg_list, ref fg2bg_trfm, ref fg_cm, ref bg_cm);
        return true;
    }
    bool SubCalcBg2FgLT(ref Vector4 bg2fg_trfm, ref Vector2 fg_cm, ref Vector2 bg_cm)
    {
        List<Vector2> bg_list = new List<Vector2>();
        List<Vector2> fg_list = new List<Vector2>();
        foreach (var bgb in bg_active_ball)
            if (bgb.match_self_fg != null)
            {
                bg_list.Add(new Vector2(bgb.col, bgb.row));
                fg_list.Add(new Vector2(bgb.match_self_fg.col, bgb.match_self_fg.row));
            }

        MathUtils.Calc_abcd_xy_E_uv(bg_list, fg_list, ref bg2fg_trfm, ref bg_cm, ref fg_cm);
        return true;
    }
    public void SubFitFg2BgLT(float dist_err = 4.01f)//int itr, float fall_speed, float sort_grid_speed, float drag_unmatch_lbd)
    {
        //bool debug_log = false;
        //if (sub_hm_id == 2)
            //debug_log = true;
        BuildFGMapBG(true);
        Vector4 bg2fg_trfm = Vector4.zero;
        Vector4 fg2bg_trfm = Vector4.zero;
        Vector2 bg_cm = Vector2.zero;
        Vector2 fg_cm = Vector2.zero;
        Vector2 bg_cm2 = Vector2.zero;
        Vector2 fg_cm2 = Vector2.zero;
        SubCalcBg2FgLT(ref bg2fg_trfm, ref fg_cm2, ref bg_cm2);
        SubCalcFg2BgLT(ref fg2bg_trfm, ref fg_cm, ref bg_cm, dist_err);
        Vector4 inv_bg2fg = MathUtils.InvABCD(bg2fg_trfm);
        Vector4 mean_fg2bg = 0.5f * (inv_bg2fg + fg2bg_trfm);
        //Vector4 inv_fg2bg = MathUtils.InvABCD(fg2bg_trfm);
        //Vector4 mean_bg2fg = 0.5f * (inv_fg2bg + bg2fg_trfm);

        foreach (var fgb in fg_active_ball)
        {
            Vector3 trfm_pos = bg_cm + MathUtils.ApplyLT(mean_fg2bg, fgb.GetInitPos(), fg_cm);
            fgb.tfm_pos = trfm_pos;// Vector3.Lerp(fgb.tfm_pos, trfm_pos, 0.9f);//calc fg.tfm_pos by Linear-Transform
            //fgb.sub_id = 0;
            fgb.FindFgMatchBg4(bg_ball);
            fgb.sub_id = (fgb.dist_fg2bg < dist_err) ? sub_hm_id : 0;
        }
    }
    void CopySubId2Root()
    {
        for (int i = 0; i < fg_active_ball.Count; i++)
        {
            HBall fgb = fg_active_ball[i];
            fgb.sub_id = 0;
            foreach (var sub_hm in sub_hm_list)
                fgb.sub_id |= sub_hm.fg_active_ball[i].sub_id;
        }
    }
    void ExpandZeroSubId()
    {
        for (int i = 0; i < fg_active_ball.Count; i++)
        {
            HBall fgb = fg_active_ball[i];
            fgb.expand_sid = fgb.sub_id;
            if (fgb.sub_id == 0)
            {
                fgb.expand_sid = 0;
                foreach (var nhb in fgb.neibor8)
                    if (nhb != null)
                        fgb.expand_sid |= nhb.sub_id;
            }
        }
        for (int i = 0; i < fg_active_ball.Count; i++)
        {
            HBall fgb = fg_active_ball[i];
            fgb.sub_id = fgb.expand_sid;
        }
    }
    public void ExpandSubIdIfZero(int expand_cnt = 2)
    {
        for (int kk = 0; kk < expand_cnt; kk++)
        {
            //copy sub.sub_id to self.sub_id
            CopySubId2Root();
            //if (sub_id==0)expand neibor.sub_id 
            ExpandZeroSubId();
            for (int i = 0; i < fg_active_ball.Count; i++)
            {
                HBall fgb = fg_active_ball[i];
                fgb.sub_id = fgb.expand_sid;
                foreach (var sub_hm in sub_hm_list)
                    sub_hm.fg_active_ball[i].sub_id = (fgb.sub_id & sub_hm.sub_hm_id);
            }
        }
        foreach (var sub_hm in sub_hm_list)
            sub_hm.SubHMReBuildFg();
    }
    public void SortInactiveGrid(HBall[,] hb_root, float sort_grid_speed)
    {
        foreach (var hb in fg_active_ball)
            if (!hb.IsValidVal())
            {
                //hb.SortGrid(sort_grid_speed);
                var rel_hb = HBall.GetHBSameRC(hb_root, hb);
                Vector3 tar_pos_lod0 = rel_hb.neibor_lod0.GetTarGridPos(fg_ball);
                hb.tfm_pos = Vector3.Lerp(hb.tfm_pos, tar_pos_lod0, sort_grid_speed);
                //hb.FindFgMatchBg4(bg_ball, 0.5f);
                //if (!hb.b_no_match && hb.sub_id != sub_hm_id)
                //    hb.sub_id = sub_hm_id;
            }
    }
    
    public void NextStep(int itr, float fall_speed, float sort_grid_speed, float drag_unmatch_lbd)
    {
        BuildFGMapBG();
        BallFallDown(fall_speed);       //H gradient drive--->Ball
        UnmatchFGBDragByCW(drag_unmatch_lbd);
        UncoveredBGDragFGB(HBall.X01Y01(itr, 0, 10, 0, 1));
        float grid_r = 1;// HBall.X01Y01(itr, 0, 20, 0, 1);
        BallSortGrid(sort_grid_speed * grid_r);  //Ball grid-net drive--->Ball
        
        AccMoveByMoveDir();
        ApplyBallByH(fg_active_ball, H);//H height ----->Ball
    }
   
    public float EstimateAllLoss()
    {
        foreach (var sub_hm in sub_hm_list)
            sub_hm.EstimateAllLoss();
        MatchLoss mloss = hm_match_loss;

        mloss.Clear();

        int fg_cnt = fg_active_ball.Count;
        int bg_cnt = bg_active_ball.Count;
        if (fg_cnt == 0 || bg_cnt == 0)
            return 0;

        BuildFGMapBG();
        
        for (int i = 0; i < bg_cnt; i++)
        {
            HBall bgb = bg_active_ball[i];
            bgb.BGMaxBreakSca();
            if (!bgb.CanBg2Fg2BgMatchback(2))
            {
                mloss.bg_can_not_match_back_cnt++;
                mloss.loss_bg2fg2bg_back += bgb.dist_bg2fg2bg;
            }
            if (bgb.b_bg_break_sca)
            {
                mloss.bg_break_sca_cnt++;
                mloss.loss_bg_break_sca += bgb.bg_break_loss;
            }
            if (!bgb.b_bg_covered)
            {
                mloss.bg_no_match_val_cnt++;
                mloss.loss_bg_no_match += bgb.dist_bg2fg;
            }
        }
        //bool temp_debug = false;
        //if (temp_debug) Debug.Log("========================== sub-hm-id=" + sub_hm_id);
        for (int i = 0; i < fg_cnt; i++)
        {
            HBall cfgb = fg_active_ball[i];
            if (!cfgb.is_active)
                continue;
            cfgb.FGMaxBreakSca();
            if (!cfgb.KeepShapeInRange(2))
            {
                mloss.fg_break_shape_cnt++;
                mloss.loss_fg_keep_shape += cfgb.dist_fg_keep_shape;
            }
            if (!cfgb.CanFg2Bg2FgbMatchback(2))
            {
                mloss.fg_can_not_match_back_cnt++;
                mloss.loss_fg2bg2fg_back += cfgb.dist_fg2bg2fg;
            }
            if (cfgb.b_fg_break_sca)
            {
                //float nbr = HBall.X01Y01(cfgb.NeiborBreakScaRate(), 2, 8, 1, 5);
                //if (cfgb.b_no_match)
                //    nbr *= 2;
                mloss.fg_break_sca_cnt++;
                mloss.loss_fg_break_sca += cfgb.fg_break_loss;// * nbr;
            }
            if (cfgb.b_no_match)
            {
                mloss.fg_no_match_val_cnt++;
                mloss.loss_fg_no_match += cfgb.dist_fg2bg;
            }
        }
        if (sub_hm_list.Count > 0)
        {
            mloss.loss_sub_fg_keep_shape = 0;
            foreach (var sub_hm in sub_hm_list)
                mloss.loss_sub_fg_keep_shape += sub_hm.hm_match_loss.loss_fg_keep_shape;
        }

        mloss.DivFgBgCnt(fg_cnt, bg_cnt);
        
        return mloss.GetMatchVal();
    }
    public void RunThin(bool init_or_rt)
    {

        ThinHBList(fg_active_ball, init_or_rt);
        ThinHBList(bg_active_ball, init_or_rt);
    }
    public void ThinHBList(List<HBall> _hb_list, bool init_or_rt)
    {
        //return;
        //if(!init_or_rt)
            CalcAllPointTanNorm(_hb_list, init_or_rt);
        SmoothAllPointTanNorm(_hb_list);

        for (int i = 0; i < _hb_list.Count; i++)
        {
            var _hb = _hb_list[i];
        //bool flag = sub_hm_id == 0 && _hb.fg_or_bg && _hb.row == 21 && _hb.col == 11;
            _hb.temp_id = -1;
            Vector2 neibor_cw = init_or_rt ? _hb.neibor_lod0.init_cw : _hb.neibor_lod0.UpdateRTCW();
            //Vector3 fgb2cw = _hb.neibor_lod0.init_cw - _hb.GetInitPos();
            Vector2 fgb2cw = neibor_cw - _hb.PosV2(init_or_rt);
            Vector2 line_tan = _hb.dir_t * Vector2.Dot(fgb2cw, _hb.dir_t);
            Vector2 line_nor = fgb2cw - line_tan;
            _hb.temp_pos = _hb.PosV2(init_or_rt) + new Vector2(line_nor.x, line_nor.y);
            //_hb.tfm_pos = _hb.GetInitPos() + line_nor;
            //if (flag)
            //    Debug.Log(_hb.RCStr() 
            //        + ".line_nor=" + line_nor.ToString("F1")
            //        + ", line_tan=" + line_tan.ToString("F1")
            //        + ",\n\t fgb2cw=" + fgb2cw.ToString("F1")
            //        + ", neibor_cw=" + neibor_cw.ToString("F1")
            //        + ", tfm_pos=" + _hb.tfm_pos.ToString("F1")
            //        + ",\n\t dir_t="+ _hb.dir_t.ToString("F3")
            //        + ",.dir_n="+ _hb.dir_n.ToString("F3")
            //        );

        }
        for (int i = 0; i < _hb_list.Count; i++)
            _hb_list[i].tfm_pos = _hb_list[i].temp_pos;
        for (int i = 0; i < _hb_list.Count; i++)
        {
            var _hb = _hb_list[i];
            float max_wid = 0;
            foreach(var nhb in _hb.neibor8)
                if(nhb!=null && nhb.IsValidVal())
                {
                    Vector2 d0 = nhb.tfm_pos - _hb.tfm_pos;
                    float _wid = Vector2.Dot(d0, _hb.dir_n);
                    float _abs = Mathf.Abs(_wid);
                    if (max_wid < _abs)
                        max_wid = _abs;
                }
            if (max_wid < 0.2f)
                _hb.temp_id = 4;
            else if (max_wid < 0.4f)
                _hb.temp_id = 3;
            else if (max_wid < 0.8f)
                _hb.temp_id = 2;
            else if (max_wid < 1.2f)
                _hb.temp_id = 1;
            else if (max_wid < 1.6f)
                _hb.temp_id = 0;
            else
                _hb.temp_id = 6;

        }

    }
    public void RunItr0123()
    {
        RunItr0();      //to bg-bfs-root, sort-grid X 3, search sub-hm in range-4, expand-zero/head
        RunItr1();      //sub-hm NextStep-X3(include inactive hb), determinate sub-id by pos
        RunItr2();      //exclude sub-id by depth - reverse-depth
        RunItr3(3); //NextStep-X3, merge sub-hm pos to root
    }
    //to bg-bfs-root, sort-grid X 3, search sub-hm in range-4
    public void RunItr0(int sort_grid_times = 3)
    {
        for (int i = 0; i < fg_active_ball.Count; i++)
        {
            HBall chb = fg_active_ball[i];
            //HBall bgb = GetBGRCRoot(chb.row, chb.col);
            //chb.tfm_pos = Vector3.Lerp(chb.GetRC(), bgb.GetRC(), 0.9f);
            int nr = (int)(Mathf.Clamp(0.5f + chb.tfm_pos.y, 0, HeightMap.ih1));
            int nc = (int)(Mathf.Clamp(0.5f + chb.tfm_pos.x, 0, HeightMap.iw1));
            HBall bgb = GetBGRCRoot(nr, nc);
            chb.tfm_pos = Vector3.Lerp(chb.tfm_pos, bgb.GetInitPos(), 0.8f);
        }
        for (int i = 0; i < sort_grid_times; i++)
            BallSortGrid(1);
        SearchBestInitMatchSubBg();
    }//itr0
    public void RunItr2(int dep_err_thresh = 2)//exclude by dep-reverse err, TODO, find all hb with same dep_id
    {
        //var self_fg_arr = fg_active_ball;
        foreach (var sub_hm in sub_hm_list)
        {
            var sub_fgb_arr = sub_hm.fg_active_ball;
            int l1cnt = sub_hm.fg_lhm.lod1.Count;
            if (l1cnt < 1)
                break;
            int l1_max_dep = sub_hm.fg_lhm.lod1.depth_list[l1cnt - 1];
            for (int i = 0; i < sub_fgb_arr.Count; i++)
            {
                //var self_fg = self_fg_arr[i];
                var sub_fg = sub_fgb_arr[i];
                sub_fg.sub_id = 0;
                sub_fg.temp_id = -3;

                if (sub_fg.copy_dep_id >= 0 && sub_fg.rvs_dep_id >= 0)
                {
                    //int dep_err = sub_fg.rvs_dep_id - (l0_max_dep - sub_fg.copy_dep_id);
                    int dep_err = sub_fg.copy_dep_id - (l1_max_dep - sub_fg.rvs_dep_id);
                    if (dep_err < 0) dep_err = -dep_err;
                    if (dep_err < dep_err_thresh)
                        sub_fg.sub_id = sub_hm.sub_hm_id;

                    dep_err -= 1;
                    if (dep_err < 0) dep_err = 0;
                    if (dep_err > 6) dep_err = 6;
                    sub_fg.temp_id = dep_err;// sub_fg.rvs_dep_id;// (dep_err>>2);
                }
            }
            //sub_hm.ExpandFgLine();
        }
        ExpandSubIdIfZero();
    }//RunItr2
    //sub-hm NextStep-X3(include inactive hb), determinate sub-id by pos
    public void RunItr1(int run_times = 3)
    {
        foreach (var sub_hm in sub_hm_list)
        {
            //sub_hm.SubFitFg2BgLT(3.1f);
            sub_hm.SubFitFg2BgLT(1.5f);
        }
        for (int k = 0; k < run_times; k++)
            foreach (var sub_hm in sub_hm_list)
            {
                sub_hm.NextStep(100, 0.5f, 0.5f, 0.5f);
                sub_hm.SortInactiveGrid(fg_ball, 0.99f);
            }
        foreach (var sub_hm in sub_hm_list)
        {
            sub_hm.DetSubIdByTfmPos();
        }
        ExpandSubIdIfZero();

    }
    public void RunItr1_LT()//
    {
        //return;
        foreach (var sub_hm in sub_hm_list)
        {
            //if(update_tfm_pos)
            //    for (int i = 0; i < hm1.fg_active_ball.Count; i++)
            //        sub_hm.fg_active_ball[i].tfm_pos = hm1.fg_active_ball[i].tfm_pos;
            sub_hm.SubFitFg2BgLT(3.1f);
            sub_hm.SubFitFg2BgLT(1.5f);
            //sub_hm.SubFitFg2BgLT(1.1f);
        }
        ExpandSubIdIfZero();
    }
    public void RunItr3(int run_times, float fall_speed = 0.5f, float sort_grid_speed = 0.5f, float drag_unmatch_lbd = 0.5f)
    {
        for (int ii = 0; ii < run_times; ii++)
        {
            foreach (var _sub in sub_hm_list)
                _sub.NextStep(100, fall_speed, sort_grid_speed, drag_unmatch_lbd);
        }
        foreach (var fgb in fg_active_ball)
            MergeSubFgbPos(fgb);
        for (int ii = 0; ii < run_times; ii++)
            //NextStep(100, fall_speed, sort_grid_speed, drag_unmatch_lbd);
            BallSortGrid(1);  //Ball grid-net drive--->Ball
    }
    void MergeSubFgbPos(HBall fgb)
    {
        Vector3 mean_tfm_pos = fgb.tfm_pos;
        float merge_cnt = 1;
        float sub_w = 3;
        int _r = fgb.row;
        int _c = fgb.col;
        foreach (var _sub in sub_hm_list)
        {
            HBall _sub_fg = _sub.fg_ball[_r, _c];
            if (_sub_fg.is_active && !_sub_fg.b_no_match)
            {
                mean_tfm_pos += _sub_fg.tfm_pos * sub_w;
                merge_cnt += sub_w;
            }
        }
        mean_tfm_pos /= merge_cnt;
        fgb.tfm_pos = mean_tfm_pos;
        foreach (var _sub in sub_hm_list)
        {
            HBall _sub_fg = _sub.fg_ball[_r, _c];
            if (_sub_fg.is_active)// && !_sub_fg.b_no_match)
                _sub_fg.tfm_pos = mean_tfm_pos;
        }
    }
    //*
    public void SearchBestInitMatchSubBg(int search_rng = 2)
    {
        foreach (var sub_hm in sub_hm_list)
        {
            sub_hm.SubHMClearActiveFgb();
            sub_hm.FindBestBiasFGBByInitCW(fg_active_ball, search_rng);
            for (int i = 0; i < fg_active_ball.Count; i++)//record tfm_pos
            {
                var sub_fgb = sub_hm.fg_active_ball[i];
                sub_fgb.sub_id = sub_fgb.b_no_match ? 0 : sub_hm.sub_hm_id;
            }

            //sub_hm.SubHMBuildActiveFgb(false);// true);
            //sub_hm.BuildLineOrder(sub_hm.fg_ball, sub_hm.fg_active_ball, ref sub_hm.fg_lhm);
            //Debug.Log("sub_hm["+ sub_hm.sub_hm_id+"].fg_active_ball.cnt=" + sub_hm.fg_active_ball.Count);
        }
        ExpandSubIdIfZero();
        //sub_hm.ExpandFgLine(2);
    }//*/
    void AccMoveByMoveDir()
    {
        foreach (var fgb in fg_active_ball)
            fgb.AccMoveByMoveDir();
    }
    //public void SetHBSubId(int row, int col, int sid, bool active, Vector3 _tar_pos)
    //{
    //    HBall fgb = fg_ball[row, col];
    //    if (fgb.val < HBall.ms_active_thresh) return;
    //    fgb.sub_id = sid;
    //    fgb.is_active = active;
    //    fgb.tfm_pos.x = Mathf.Lerp(fgb.tfm_pos.x, _tar_pos.x, 0.8f);
    //    fgb.tfm_pos.y = Mathf.Lerp(fgb.tfm_pos.y, _tar_pos.y, 0.8f);
    //}
    //public void CalcSubHMLinearTransBg()
    //{
    //    List<Vector2> xy_list = new List<Vector2>();
    //    List<Vector2> uv_list = new List<Vector2>();
    //    foreach (var bgb in bg_active_ball)
    //        if (bgb.match_self_fg != null)
    //        {
    //            xy_list.Add(new Vector2(bgb.col, bgb.row));
    //            uv_list.Add(new Vector2(bgb.match_self_fg.col, bgb.match_self_fg.row));
    //        }
    //    Vector4 abcd = Vector4.zero;
    //    Vector2 xy_cm = Vector2.zero;
    //    Vector2 uv_cm = Vector2.zero;
    //    if (!MathUtils.Calc_abcd_xy_E_uv(xy_list, uv_list, ref abcd, ref xy_cm, ref uv_cm))
    //        return;

    //    //clear sub-id
    //    foreach (var fgb in fg_active_ball)
    //    {
    //        fgb.sub_id = 0;
    //        fgb.is_active = false;
    //        fgb.tfm_pos = fgb.GetInitPos();
    //    }
    //    foreach (var bgb in bg_active_ball)
    //    {
    //        Vector3 bgb_pos = new Vector3(bgb.col, bgb.row, 0);
    //        float dcol = (bgb.col - xy_cm.x);
    //        float drow = (bgb.row - xy_cm.y);
    //        float map_col = abcd.x * dcol + abcd.y * drow + uv_cm.x;
    //        float map_row = abcd.z * dcol + abcd.w * drow + uv_cm.y;
    //        int r0 = Mathf.RoundToInt(map_row);
    //        int c0 = Mathf.RoundToInt(map_col);
    //        if (r0 < 0) r0 = 0;
    //        if (c0 < 0) c0 = 0;
    //        int r1 = r0 + 1;
    //        int c1 = c0 + 1;
    //        if (r0 > ih1) r0 = ih1;
    //        if (c0 > iw1) c0 = iw1;
    //        if (r1 > ih1) r1 = ih1;
    //        if (c1 > iw1) c1 = iw1;
    //        SetHBSubId(r0, c0, sub_hm_id, true, bgb_pos);
    //        SetHBSubId(r1, c0, sub_hm_id, true, bgb_pos);
    //        SetHBSubId(r0, c1, sub_hm_id, true, bgb_pos);
    //        SetHBSubId(r1, c1, sub_hm_id, true, bgb_pos);
    //    }
    //}
    //public void SubHMUpdateHMActiveFg()
    //{
    //    foreach (var sub_hm in sub_hm_list)
    //    {
    //        sub_hm.CalcSubHMLinearTransBg();
    //    }
    //    foreach (var fgb in fg_active_ball)
    //    {
    //        int _r = fgb.row;
    //        int _c = fgb.col;
    //        fgb.sub_id = 0;
    //        foreach (var sub_hm in sub_hm_list)
    //            fgb.sub_id |= sub_hm.fg_ball[_r, _c].sub_id;
    //        if(false)// (fgb.sub_id == 0)
    //        {
    //            foreach (var nb8 in fgb.neibor8)
    //                if (nb8 != null)
    //                {
    //                    int _r8 = nb8.row;
    //                    int _c8 = nb8.col;
    //                    foreach (var sub_hm in sub_hm_list)
    //                    {
    //                        int csid = sub_hm.fg_ball[_r8, _c8].sub_id;
    //                        fgb.sub_id |= csid;
    //                        if (csid == sub_hm.sub_hm_id)
    //                        {
    //                            sub_hm.fg_ball[_r, _c].sub_id = csid;
    //                            sub_hm.fg_ball[_r, _c].is_active = true;
    //                            sub_hm.fg_ball[_r, _c].tfm_pos = Vector3.Lerp(sub_hm.fg_ball[_r, _c].tfm_pos, sub_hm.fg_ball[_r8, _c8].tfm_pos, 0.8f);
    //                        }
    //                    }
    //                }
    //        }
    //    }
    //    foreach (var sub_hm in sub_hm_list)
    //    {
    //        sub_hm.SubHMBuildActiveFgb(false);
    //        sub_hm.BuildFGMapBG();
    //    }
    //}
    //public void HMUpdateSubHMActiveFg()
    //{
    //    //update sub_hm.fg

    //    foreach (var fgb in fg_active_ball)
    //    {
    //        if (!fgb.IsValidVal())
    //            continue;
    //        fgb.FindFgMatchBg4(bg_ball);
    //        int fgb_sub_id = 0;
    //        for (int k = 0; k < 4; k++)
    //        {
    //            float _dist = fgb.b_match_bg_dist4[k];
    //            HBall _mbg4 = fgb.match_self_bg4[k];
    //            //if (fgb.row == 15 && fgb.col == 12)
    //            //Debug.Log("fgb[15,12] _dist=" + _dist.ToString("F1"));
    //            if (_dist < 2)
    //                fgb_sub_id |= (_mbg4.sub_id);
    //        }
    //        fgb.sub_id = fgb_sub_id;

    //    }
    //    //expand
    //    //for(int e = 0; e < 0; e++)
    //    //foreach (var fgb in fg_active_ball)
    //    //if (fgb.is_active && fgb.sub_id == 0)
    //    //{
    //    //    foreach (var nb in fgb.neibor4)
    //    //        if (nb != null)
    //    //            fgb.sub_id |= nb.sub_id;
    //    //}
    //    //apply fg 2 sub-hm
    //    foreach (var _sub in sub_hm_list)
    //    {
    //        _sub.SubHMClearActiveFgb();
    //        foreach (var fgb in fg_active_ball)
    //            if (fgb.is_active)
    //            {
    //                _sub.UpdateFgSubId(fgb.row, fgb.col, fgb.sub_id, fgb.tfm_pos);
    //            }
    //        _sub.SubHMBuildActiveFgb(false);
    //        //_sub.NextStep(0, 0, 0, 0);
    //        _sub.BuildLineOrder(_sub.fg_ball, _sub.fg_active_ball, ref _sub.fg_lhm);
    //        _sub.EstimateAllLoss();
    //    }
    //}
    public void Refresh()
    {
        //Debug.Log("refresh");
        sub_fg_active_cnt = prev_sub_fg_active_cnt = 0;
        AlignInitCW();
        ApplyBallByH(fg_active_ball, H);
        //foreach(var fgb in fg_active_ball)
            //fgb.ClearLineInfo();
        foreach (var bgb in bg_active_ball)
        {
            //bgb.ClearLineInfo();
            bgb.tfm_pos = bgb.GetInitPos();
        }
        foreach(var fgb in fg_active_ball)
        {

        }
        foreach (var _sub in sub_hm_list)
        {
            foreach(var fgb in _sub.fg_active_ball)
            {
                fgb.sub_id = 0;
                fgb.is_active = false;
                fgb.tfm_pos = fgb.GetInitPos();

            }
            _sub.Refresh();
        }
    }
    public HBall GetBGRCRoot(int row, int col)
    {
        int row_root = bg_ball[row, col].bfs_root_row;
        int col_root = bg_ball[row, col].bfs_root_col;
        return bg_ball[row_root, col_root];
    }
}
public class HeightMapPair
{
    public bool failed_to_best = false;
    public int i_run_itr = 0;
    public HeightMap hm1;
#if def_hm2
    public HeightMap hm2;
#endif
    //public float match_val = 0, fg_match = 0, bg_match = 0, no_match_val = 0, break_sca_val = 0;

    public MatchLoss[] fg_match_loss_arr = new MatchLoss[5 * 5 * 5];
#if def_hm2
    public MatchLoss bg_match_loss = new MatchLoss();
    public MatchLoss[] bg_match_loss_arr = new MatchLoss[5*5*5];
#endif
    static float[] w5 = new float[5] { 0.1f, 0.3f, 1, 3f, 10f };
    public HeightMapPair()
    {
        hm1 = new HeightMap();
#if def_hm2
        hm2 = new HeightMap();
#endif

        //fg_match_loss.SetW(1, 0.1f, 0.1f);
        //bg_match_loss.SetW(1, 0.1f, 0.1f);
        int mi = 0;
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
                for (int k = 0; k < 5; k++)
                {
                    float wi = w5[i];
                    float wj = w5[j];
                    float wk = w5[k];
                    fg_match_loss_arr[mi] = new MatchLoss();
                    fg_match_loss_arr[mi].SetW(wi, wj, wk);
#if def_hm2
                    bg_match_loss_arr[mi] = new MatchLoss();
                    bg_match_loss_arr[mi].SetW(wi, wj, wk);
#endif
                    mi++;
                }
    }
    public float GetMatchVal()
    {
#if def_hm2
        return hm1.GetMatchVal() + hm2.GetMatchVal();
#else
        return hm1.GetMatchVal();
#endif
    }
    public float GetMatchValAt(int i)
    {
#if def_hm2
        return fg_match_loss_arr[i].GetMatchVal() + bg_match_loss_arr[i].GetMatchVal();
#else
        return fg_match_loss_arr[i].GetMatchVal();
#endif
    }
    public void SetDefaultByte()
    {
        byte[,] bg_byt = GetZerosBytes();
        byte[,] fg_byt = GetZerosBytes();
        FillBytesLines(bg_byt, 5, 6, 10, 25);
        FillBytesLines(bg_byt, 5, 16, 10, 11);
        FillBytesLines(bg_byt, 15, 16, 10, 20);
        FillBytesLines(bg_byt, 5, 13, 24, 25);
        FillBytesLines(bg_byt, 10, 10, 6, 9);

        FillBytesLines(fg_byt, 10, 11, 7, 15);
        FillBytesLines(fg_byt, 10, 17, 7, 8);
        FillBytesLines(fg_byt, 17, 17, 7, 15);
        FillBytesLines(fg_byt, 10, 17, 15, 15);
        SetRawByte(fg_byt, bg_byt, null);
        //SetRawByte(bg_byt, fg_byt);
    }
    public void SetRawByte(byte[,] bgb, byte[,] fgb, List<byte[,]> sub_bgb_list = null)
    {
        hm1.SetBgH(bgb, sub_bgb_list);
        hm1.SetFgH(fgb);

#if def_hm2
        hm2.SetFgH(bgb, sub_bgb_list);
        hm2.SetBgH(fgb);
#endif
        Refresh();
    }
    public void Refresh()
    {
        failed_to_best = false;
        i_run_itr = 0;
        hm1.Refresh();
#if def_hm2
        hm2.Refresh();
#endif
    }
    public void NextStep(float fall_speed, float sort_grid_speed, float drag_unmatch_lbd)
    {
        if (failed_to_best)//skip next-step because match too bad
            return;
        if (i_run_itr == 0)
        {
            //hm1.RunItr0123();
        }
        else
        {
            hm1.NextStep(i_run_itr, fall_speed, sort_grid_speed, drag_unmatch_lbd);
        }
        EstimateMatchValue();
        i_run_itr++;
    }
    
    void SplitSubHMCrossPart()
    {
        //find mix-sub-id HB, set max-sub-id
        /*
        for (int i = 0; i < hm1.fg_active_ball.Count; i++)
        {
            HBall fgb = hm1.fg_active_ball[i];
            int fgb_sid = fgb.sub_id;
            fgb.expand_sid = fgb_sid;
            //if (fgb_sid == 0 || fgb_sid == 1 || fgb_sid == 2 || fgb_sid == 4 || fgb_sid == 8)
            //    continue;
            int sub_cnt = hm1.sub_hm_list.Count;
            int sub_valid_cnt = 0;
            int[] nei_sid_cnt = new int[sub_cnt];
            for (int k = 0; k < sub_cnt; k++)
            {
                if ((fgb_sid | (1 << k)) != 0)
                    nei_sid_cnt[k] = 1;
                else
                    nei_sid_cnt[k] = 0;
            }
            foreach (var nhb in fgb.neibor8)
                if (nhb != null)
                {
                    int nhb_sid = nhb.sub_id;
                    for (int k = 0; k < sub_cnt; k++)
                    {
                        if (((1 << k) & nhb_sid) != 0)
                            nei_sid_cnt[k]++;
                    }
                    sub_valid_cnt++;
                }
            int max_k = 0;
            int max_c = nei_sid_cnt[0];
            for (int k = 1; k < sub_cnt; k++)
                if (max_c < nei_sid_cnt[k])
                {
                    max_c = nei_sid_cnt[k];
                    max_k = k;
                }
            fgb.expand_sid = (1 << max_k);
            for (int k = 0; k < sub_cnt; k++)
                if (k != max_k && nei_sid_cnt[k] == max_c)
                    fgb.expand_sid |= (1 << k);
            //if (fgb.row == 12 && fgb.col == 15)
            //    Debug.Log("fgb[].sub-id=" + fgb.sub_id + ", exp-id=" + fgb.expand_sid
            //        + ", max k:" + max_k + ", c=" + max_c
            //        + ", nsid[" + nei_sid_cnt[0] + ", " + nei_sid_cnt[1] + ", " + nei_sid_cnt[2]);
        }
        //*/
    }

    float EstimateMatchValue()
    {
#if def_hm2
        EstimateMatchValue(ref hm1.fg_match_loss, hm1, hm2);
        EstimateMatchValue(ref hm2.fg_match_loss, hm2, hm1);
#else
        hm1.EstimateAllLoss();

#endif

        return GetMatchVal();
    }

    public void BallSortGrid(float sort_grid_speed = 0.5f)
    {
        hm1.BallSortGrid(sort_grid_speed);
#if def_hm2
        hm2.BallSortGrid(sort_grid_speed);
#endif
    }
    public void BallFallDown(float fall_speed = 0.5f)
    {
        hm1.BallFallDown(fall_speed);
#if def_hm2
        hm2.BallFallDown(fall_speed);
#endif
    }


    public static byte[,] GetZerosBytes()
    {
        byte[,] ret = new byte[HeightMap.ih, HeightMap.iw];
        for (int i = 0; i < HeightMap.ih; i++)
            for (int j = 0; j < HeightMap.iw; j++)
                ret[i, j] = 0;
        return ret;
    }
    public void FillBytesLines(byte[,] _byt, int r0, int r1, int c0, int c1)
    {
        for (int i = r0; i <= r1; i++)
            for (int j = c0; j <= c1; j++)
                _byt[i, j] = byte.MaxValue;
    }
}
public class ViewHeightMap
{
    public HeightMap hm = null;
    Mesh bg_mesh, fg_mesh, h_mesh;
    public Transform quat_fg, quat_bg, quat_h;
    public GameObject view_root = null;
    GameObject fg_ball_root = null;
    GameObject bg_ball_root = null;
    List<Transform> show_fg_active_ball = new List<Transform>();
    List<Transform> show_bg_active_ball = new List<Transform>();
    static Material _black_mat, _white_mat, _red_mat, _green_mat, _blue_mat, _gray_mat;
    static Material _yellow_mat, _cyan_mat, _magenta_mat, _orange_mat, _purple_mat, _black_gray_mat;
    static Material[] mats7 = new Material[7];

    public List<ViewHeightMap> sub_hmv_list = new List<ViewHeightMap>();
    public void CreateHeightMapTfm(HeightMap _hm, string desc)
    {
        hm = _hm;
        string _go_n = "view-Hmap-" + desc;
        //Debug.Log("create :" + desc);
        if (view_root == null)
        {
            if (GameObject.Find(_go_n) != null)
                GameObject.DestroyImmediate(GameObject.Find(_go_n));
            view_root = new GameObject(_go_n);
            if (fg_ball_root == null)
            {
                fg_ball_root = new GameObject("fg_active_ball");
                fg_ball_root.transform.parent = view_root.transform;
            }
            if (bg_ball_root == null)
            {
                bg_ball_root = new GameObject("bg_active_ball");
                bg_ball_root.transform.parent = view_root.transform;
            }
            bg_mesh = CreateMesh27();
            fg_mesh = CreateMesh27();
            h_mesh = CreateMesh27();
            quat_bg = CreateQuatBG(bg_mesh, "quat_bg");
            quat_fg = CreateQuatBG(fg_mesh, "quat_fg");
            quat_h = CreateQuatBG(h_mesh, "quat_h");
            quat_bg.parent = view_root.transform;
            quat_fg.parent = view_root.transform;
            quat_h.parent = view_root.transform;
            quat_bg.localPosition += new Vector3(29, 0, 0);
            quat_fg.localPosition -= new Vector3(29, 0, 0);
            bg_ball_root.transform.localPosition = quat_bg.localPosition;
        }
        view_root.name = _go_n;
        if (_black_mat == null)
        {
            _black_mat = U3DUtils.CreateMaterial("");
            _white_mat = U3DUtils.CreateMaterial("");
            _gray_mat = U3DUtils.CreateMaterial("");
            _black_gray_mat = U3DUtils.CreateMaterial("");

            _red_mat = U3DUtils.CreateMaterial("");
            _green_mat = U3DUtils.CreateMaterial("");
            _blue_mat = U3DUtils.CreateMaterial("");

            _yellow_mat = U3DUtils.CreateMaterial("");
            _cyan_mat = U3DUtils.CreateMaterial("");
            _magenta_mat = U3DUtils.CreateMaterial("");
            _orange_mat = U3DUtils.CreateMaterial("");
            _purple_mat = U3DUtils.CreateMaterial("");
            //, , , , ;
            _black_mat.color = Color.black;//fg_cube_color
            _white_mat.color = Color.white;
            _gray_mat.color = Color.gray;
            _black_gray_mat.color = Color.Lerp(Color.black, Color.gray, 0.5f);

            _red_mat.color = Color.red;
            _green_mat.color = Color.green;
            _blue_mat.color = Color.blue;

            _yellow_mat.color = Color.yellow;
            _cyan_mat.color = Color.cyan;
            _magenta_mat.color = Color.magenta;

            _orange_mat.color = new Color(1, 0.5f, 0);
            _purple_mat.color = new Color(0.5f, 0.1f, 0.8f);
            mats7 = new Material[7] { _red_mat, _orange_mat, _yellow_mat, _green_mat,
                _cyan_mat, _blue_mat, _purple_mat };
        }
        if (hm == null) return;

        UpdateShowTfm(ref show_fg_active_ball, hm.fg_active_ball, fg_ball_root);
        UpdateShowTfm(ref show_bg_active_ball, hm.bg_active_ball, bg_ball_root);
        //add sub hmv================
        if (hm.sub_hm_list.Count > 0)
        {
            for (int k = sub_hmv_list.Count; k < hm.sub_hm_list.Count; k++)
            {
                ViewHeightMap hmv = new ViewHeightMap();
                sub_hmv_list.Add(hmv);
            }
            for (int k = hm.sub_hm_list.Count; k < sub_hmv_list.Count; k++)
            {
                if (sub_hmv_list[k].view_root != null)
                    sub_hmv_list[k].view_root.SetActive(false);
            }
            for (int k = 0; k < hm.sub_hm_list.Count && k < sub_hmv_list.Count; k++)
            {
                ViewHeightMap chmv = sub_hmv_list[k];
                chmv.CreateHeightMapTfm(hm.sub_hm_list[k], desc + "_sub_" + k);
                chmv.view_root.transform.parent = view_root.transform;
                chmv.view_root.transform.localPosition = new Vector3(0, (k + 1) * 20, (k + 1) * 20);
            }
        }
    }
    void UpdateShowTfm(ref List<Transform> _show_hb_tfm, List<HBall> _active_hb_list, GameObject _root)
    {
        //if (hm.sub_hm_id > 0)
        //Debug.Log("sub hm[" + hm.sub_hm_id + "], fg active cnt=" + _active_hb_list.Count + ", show ball cnt:" + show_fg_active_ball.Count);
        //_show_hb_tfm.Clear();
        if (_show_hb_tfm.Count < _active_hb_list.Count)
        {
            int dcnt = _active_hb_list.Count - _show_hb_tfm.Count;
            for (int kk = 0; kk < dcnt; kk++)
            {
                var tfm = GameObject.CreatePrimitive(PrimitiveType.Cube).transform;
                tfm.gameObject.GetComponent<MeshRenderer>().material = _black_mat;
                tfm.parent = _root.transform;
                tfm.localScale = Vector3.one * 0.76f;
                _show_hb_tfm.Add(tfm);
            }
        }

        for (int i = 0; i < _active_hb_list.Count; i++)
        {
            var fg_ball = _active_hb_list[i];
            //var tfm = GameObject.CreatePrimitive(PrimitiveType.Cube).transform;
            //tfm.gameObject.GetComponent<MeshRenderer>().material = _black_mat;
            //tfm.parent = fg_ball_root.transform;
            var tfm = _show_hb_tfm[i];
            tfm.localPosition = fg_ball.tfm_pos;// new Vector3(col, row, 0);
            tfm.gameObject.name = "r[" + fg_ball.row + "]c[" + fg_ball.col + "]";
            tfm.gameObject.SetActive(true);
        }
        for (int i = _active_hb_list.Count; i < _show_hb_tfm.Count; i++)
            _show_hb_tfm[i].gameObject.SetActive(false);
    }
    Material FGGetMatByFlag(HBall chb)
    {
        Material _m = _black_mat;
        if (chb.b_no_match)
            _m = _red_mat;
        if (chb.b_fg_break_sca)
            _m = _blue_mat;
        if (chb.b_no_match && chb.b_fg_break_sca)
            _m = _yellow_mat;
        if (!chb.IsValidVal())
            _m = _gray_mat;
        return _m;
    }
    public void ShowLineInfo(float z=-3)
    {
        ShowLineInfo(hm.fg_active_ball, show_fg_active_ball, z);
        ShowLineInfo(hm.bg_active_ball, show_bg_active_ball, z);
    }
    public void ShowLineInfo(List<HBall> _hb_list, List<Transform> _tfm_list, float z= -3)
    {
        for (int i = 0; i < _hb_list.Count && i < _tfm_list.Count; i++)
        {
            var tfm = _tfm_list[i];
            var hb = _hb_list[i];
            tfm.localPosition = new Vector3(tfm.localPosition.x, tfm.localPosition.y, z);
            tfm.rotation = Quaternion.LookRotation(hb.dir_n, hb.dir_t);
            tfm.GetComponent<MeshRenderer>().material = mats7[(hb.temp_id + mats7.Length) % mats7.Length];
        }
    }
    public void FlatTfmZ(float z=-3)
    {
        if (sub_hmv_list.Count == 0)
        {
            for (int i = 0; i < show_bg_active_ball.Count && i < hm.bg_active_ball.Count; i++)
            {
                var tfm = show_bg_active_ball[i];
                var bgb = hm.bg_active_ball[i];
                tfm.localPosition = new Vector3(tfm.localPosition.x, tfm.localPosition.y, z);
                //tfm.rotation = Quaternion.FromToRotation(Vector3.right, bgb.dir_t);
                tfm.rotation = Quaternion.LookRotation(bgb.dir_n, bgb.dir_t);
            }
            foreach (var tfm in show_fg_active_ball)
                tfm.localPosition = new Vector3(tfm.localPosition.x, tfm.localPosition.y, z);
        }
        foreach (var sub_hmv in sub_hmv_list)
            sub_hmv.FlatTfmZ(z);
    }
    public void ShowTempID()
    {
        for (int i = 0; i < hm.fg_active_ball.Count && i < show_fg_active_ball.Count; i++)
        {
            HBall chb = hm.fg_active_ball[i];
            int tid = chb.temp_id;
            var drd = show_fg_active_ball[i].GetComponent<MeshRenderer>();
            if (tid == -3)
                drd.material = _black_mat;
            else if (tid == -2)
                drd.material = _black_gray_mat;
            else if (tid == -1)
                drd.material = _gray_mat;
            else
            {
                if (tid >= 0)
                    drd.material = mats7[tid % mats7.Length];
                else
                    drd.material = _white_mat;
            }
        }
    }
    public void ShowHeightMap()
    {
        if (hm == null) return;
        for (int i = 0; i < hm.fg_active_ball.Count && i < show_fg_active_ball.Count; i++)
        {
            HBall chb = hm.fg_active_ball[i];
            show_fg_active_ball[i].GetComponent<MeshRenderer>().material = FGGetMatByFlag(chb);
        }
        SetActiveBallByTfmPos(hm.fg_active_ball, show_fg_active_ball);
        SetActiveBallByTfmPos(hm.bg_active_ball, show_bg_active_ball);
        ApplyMeshByH(fg_mesh, hm.fgH, hm.fg_ball, hm.fg_sca, 1);
        ApplyMeshByH(bg_mesh, hm.bgH, hm.bg_ball, hm.bg_sca, 2);
        ApplyMeshByH(h_mesh, hm.H, hm.fg_ball, 1, 3);
        foreach (var sub_hmv in sub_hmv_list)
            sub_hmv.ShowHeightMap();
        //if (hm.sub_hm_id==0)
        //{
        //    int rr = 11;
        //    int cc = 24;
        //    var fgb = hm.fg_ball[rr, cc];
        //    int ii = -1;
        //    for (int i = 0; i < hm.fg_active_ball.Count; i++)
        //        if (hm.fg_active_ball[i].UID() == fgb.UID())
        //            ii = i;
        //    if (ii > 0)
        //    {
        //        var shb = show_fg_active_ball[ii];
        //        Debug.Log("main hm.fg: " + fgb.GetDebugStr() + ", \n\t hm.H=" + hm.H[rr, cc] + ", show_tfm=" + shb.localPosition.ToString("F2"));
        //    }
        //}
        FlatTfmZ();
    }

    void SetActiveBallByTfmPos(List<HBall> _active_ball, List<Transform> _tfm_ball)
    {
        //if (_active_ball.Count != _tfm_ball.Count)
        //Debug.Log("apply ball by H, count :" + _active_ball.Count + "/" + _tfm_ball.Count);
        for (int i = 0; i < _active_ball.Count && i < _tfm_ball.Count; i++)
        {
            int _r = _active_ball[i].row;
            int _c = _active_ball[i].col;
            float z = -hm.H[_r, _c];
            _tfm_ball[i].localPosition = _active_ball[i].tfm_pos;// new Vector3(_active_ball[i].tfm_pos.x, _active_ball[i].tfm_pos.y, z);
            _tfm_ball[i].gameObject.SetActive(_active_ball[i].val > HBall.ms_active_thresh);
            //_tfm_ball[i].gameObject.SetActive(_active_ball[i].IsValidVal());
        }
    }
    Transform CreateQuatBG(Mesh mesh, string _name)
    {
        Transform quat = new GameObject(_name).transform;
        MeshFilter mf = quat.gameObject.AddComponent<MeshFilter>();
        mf.mesh = mesh;

        MeshRenderer mr = quat.gameObject.AddComponent<MeshRenderer>();
        Material mat = U3DUtils.CreateMaterial("Mats/mat_vertlit");
        mr.material = mat;
        return quat;
    }
    static Color SubId2Color(int sid)
    {
        if (sid == 0) return Color.black;
        else if (sid == 1) return Color.red;
        else if (sid == 2) return Color.green;
        else if (sid == 4) return Color.blue;
        else if (sid == 3) return Color.yellow;
        else if (sid == 5) return Color.magenta;
        else if (sid == 6) return Color.Lerp(Color.cyan, Color.blue, 0.3f);
        return Color.Lerp(Color.gray, Color.black, 0.5f);

    }
    void ApplyMeshByH(Mesh _m, float[,] _h = null, HBall[,] _hb_arr = null, float _z_fwd = 1.0f, int color_p = -1)
    {
        //if (color_p > 99)            Debug.Log("apply mesh , color_pol=" + color_p);
        int rows = HeightMap.ih;
        int columns = HeightMap.iw;

        // 计算顶点数和三角形数
        //int vertexCount = rows * columns;
        Vector3[] vertices = new Vector3[vertexCount * 2];
        Color[] _colors = _m.colors;
        for (int row = 0; row < rows; row++)
        {
            for (int column = 0; column < columns; column++)
            {
                // 计算顶点的位置
                int x = column;
                int y = row;
                float z = _h == null ? 0 : -_z_fwd * _h[row, column];//z.forward(far), z.backward(camera)

                int vid = row * columns + column;
                vertices[vid] = vertices[vid + vertexCount] = new Vector3(x, y, z);
                _colors[vid] = _colors[vid + vertexCount] = GetColorByPolicy(_hb_arr, color_p, x, y, z);

            }
        }
        _m.vertices = vertices;
        _m.colors = _colors;
        _m.RecalculateNormals();
    }
    Color GetColorByPolicy(HBall[,] _hb_arr, int _color_p, int col, int row, float z)
    {
        if (_color_p == 1) return FGGetMatByFlag(_hb_arr[row, col]).color;
        else if (_color_p == 2) return GetColorByBgB(_hb_arr, row, col);
        else if (_color_p == 3) return Z2Color(z);
        else if (_color_p == 100) return FgKeepShape2Color(_hb_arr, row, col);
        else if (_color_p == 101) return FgCover2Color(_hb_arr, row, col);
        else if (_color_p == 104) return GetColorByFgSubId(_hb_arr, row, col);
        else if (_color_p == 103) return FgBreakSca2Color(_hb_arr, row, col);
        else if (_color_p == 201) return BgCover2Color(_hb_arr, row, col);
        else if (_color_p == 203) return BgBreakSca2Color(_hb_arr, row, col);
        else if (_color_p == 205) return BgLineInfo2Color(_hb_arr, row, col);
        return Color.gray;
    }
    static Color FgKeepShape2Color(HBall[,] _balls, int row, int column)
    {
        if (_balls != null && _balls[row, column] != null && _balls[row, column].IsValidVal())
        {
            HBall fgb = _balls[row, column];
            float dist = fgb.dist_fg_keep_shape;
            if (!fgb.b_fg_keep_shape)
            {
                if (dist < 1) return Color.green;
                if (dist < 3) return Color.cyan;
                if (dist < 5) return Color.blue;
                if (dist < 7) return Color.yellow;
                if (dist < 9) return Color.magenta;
                return Color.red;
            }
            return Color.black;
        }
        return Color.gray;
    }
    static Color FgBreakSca2Color(HBall[,] _balls, int row, int column)
    {
        if (_balls != null && _balls[row, column] != null && _balls[row, column].IsValidVal())
        {
            HBall fgb = _balls[row, column];
            if (fgb.b_fg_break_sca)
            {
                return Color.red;
            }
            return Color.black;
        }
        return Color.gray;
    }
    static Color BgBreakSca2Color(HBall[,] _balls, int row, int column)
    {
        if (_balls != null && _balls[row, column] != null && _balls[row, column].IsValidVal())
        {
            HBall bgb = _balls[row, column];
            if (bgb.b_bg_break_sca)
                return Color.yellow;
            if (bgb.match_self_fg == null)
                return Color.red;
            return Color.black;
        }
        return Color.gray;
    }
    static Color BgLineInfo2Color(HBall[,] _balls, int row, int column)
    {
        if (_balls != null && _balls[row, column] != null && _balls[row, column].IsValidVal())
        {
            HBall bgb = _balls[row, column];
            if (bgb.dep_id <= 2)
                return Color.black;
            else if (bgb.dep_id > 0)
            {
                float dep_f = bgb.dep_id * 0.5f;
                int lbdi = Mathf.RoundToInt(dep_f);
                int lbdi1 = lbdi + 1;
                float lbd = dep_f - lbdi;
                return Color.Lerp(mats7[lbdi % mats7.Length].color, mats7[lbdi1 % mats7.Length].color, lbd);
            }

            return Color.white;
        }
        return Color.gray;
    }
    static Color FgCover2Color(HBall[,] _balls, int row, int column)
    {
        if (_balls != null && _balls[row, column] != null && _balls[row, column].IsValidVal())
        {
            HBall fgb = _balls[row, column];
            if (fgb.b_no_match)
                return Color.red;
            return Color.black;
        }
        return Color.gray;
    }
    static Color BgCover2Color(HBall[,] _balls, int row, int column)
    {
        if (_balls != null && _balls[row, column] != null && _balls[row, column].IsValidVal())
        {
            HBall bgb = _balls[row, column];
            if (!bgb.b_bg_covered)
                return Color.red;
            return Color.black;
        }
        return Color.gray;
    }
    static Color GetColorByFgSubId(HBall[,] _balls, int row, int column)
    {
        if (_balls != null && _balls[row, column] != null)
        {
            if (_balls[row, column].IsValidVal())
                return SubId2Color(_balls[row, column].sub_id);
            else if (_balls[row, column].IsValGreater())
                return Color.white;
        }
        return Color.gray;
    }
    static Color GetColorByBgB(HBall[,] _balls, int row, int column)
    {
        Color ret = Color.gray;
        if (_balls != null && _balls[row, column] != null && _balls[row, column].IsValidVal())
        {
            var cbgb = _balls[row, column];
            ret = Color.black;
            if (!_balls[row, column].b_bg_covered)
            {
                ret = Color.red;
            }
            else if (cbgb.match_self_fg == null)
            {
                ret = Color.magenta;
            }
            else
            {
                float dist_mfg = (cbgb.match_self_fg.tfm_pos - cbgb.GetInitPos()).sqrMagnitude;
                if (dist_mfg > 3.0f)
                    ret = Color.yellow;
                else if (dist_mfg > 1.5f)
                    ret = Color.green;
            }
        }
        return ret;
    }
    int vertexCount = 0;// rows * columns;
    Mesh CreateMesh27()
    {
        int rows = HeightMap.ih;
        int columns = HeightMap.iw;

        // 计算顶点数和三角形数
        vertexCount = rows * columns;  // 每个方块需要4个顶点
        Vector3[] vertices = new Vector3[vertexCount * 2];
        Color[] colors = new Color[vertexCount * 2];
        for (int row = 0; row < rows; row++)
        {
            for (int column = 0; column < columns; column++)
            {
                // 计算顶点的位置
                float x = column;
                float y = row;
                float z = 0;
                int vid = row * columns + column;
                vertices[vid] = new Vector3(x, y, z);
                colors[vid] = XY2Color(column, row);
                vertices[vid + vertexCount] = vertices[vid];
                colors[vid + vertexCount] = colors[vid];
            }
        }

        int triangleCount = (rows - 1) * (columns - 1) * 2;  // 每个方块需要2个三角形
        int[] triangles = new int[triangleCount * 3 * 2];
        //int[] triangles2 = new int[triangleCount * 3];
        for (int row = 0; row < rows - 1; row++)
        {
            for (int column = 0; column < columns - 1; column++)
            {
                int tid = (row * (columns - 1) + column) * 6;
                int v0 = row * columns + column;
                int v1 = v0 + 1;
                int v2 = v0 + rows;
                int v3 = v2 + 1;
                triangles[tid + 0] = v0;
                triangles[tid + 1] = v2;
                triangles[tid + 2] = v3;
                triangles[tid + 3] = v3;
                triangles[tid + 4] = v1;
                triangles[tid + 5] = v0;

                tid += triangleCount * 3;
                triangles[tid + 0] = v0 + vertexCount;
                triangles[tid + 1] = v3 + vertexCount;
                triangles[tid + 2] = v2 + vertexCount;
                triangles[tid + 3] = v1 + vertexCount;
                triangles[tid + 4] = v3 + vertexCount;
                triangles[tid + 5] = v0 + vertexCount;
            }
        }

        Mesh mesh = new Mesh();
        //mesh.subMeshCount = 2;
        mesh.vertices = vertices;
        mesh.colors = colors;
        mesh.triangles = triangles;
        //mesh.normals = normals;
        //mesh.SetTriangles(triangles, 0);
        //mesh.SetTriangles(triangles2, 1);
        mesh.RecalculateNormals();
        return mesh;
    }
    static Color XY2Color(int x, int y)
    {
        if (x == 0)
        {
            if (y == 5)
                return Color.red;
            else if (y == 10)
                return Color.green;
            else if (y == 15)
                return Color.blue;
            else if (y == 20)
                return Color.yellow;
            else if (y == 25)
                return Color.cyan;
        }
        else if (y == 0)
        {
            if (x == 5)
                return Color.red;
            else if (x == 10)
                return Color.green;
            else if (x == 15)
                return Color.blue;
            else if (x == 20)
                return Color.yellow;
            else if (x == 25)
                return Color.cyan;
        }
        return Color.gray;
    }
    static Color Z2Color(float z)
    {
        float rng = 1;
        float zz = Mathf.Clamp(z, -5, 3) / rng;
        Color ret = Color.gray;
        return ret;
        if (zz < -4)
            ret = Color.Lerp(Color.gray, Color.cyan, 5 + zz);
        else if (zz < -3)
            ret = Color.Lerp(Color.cyan, Color.yellow + Color.gray, 4 + zz);
        else if (zz < -2)
            ret = Color.Lerp(Color.yellow + Color.gray, Color.red + Color.gray, 3 + zz);
        else if (zz < -1)
            ret = Color.Lerp(Color.red + Color.gray, Color.magenta + Color.gray, 2 + zz);
        else if (zz < 0)
            ret = Color.Lerp(Color.magenta + Color.gray, Color.gray, 1 + zz);
        else if (zz < 1)
            ret = Color.Lerp(Color.gray, Color.green + Color.gray, zz);
        else if (zz < 2)
            ret = Color.Lerp(Color.green + Color.gray, Color.blue, zz - 1);
        else// if (zz <= 3)
            ret = Color.Lerp(Color.blue, Color.blue + Color.gray, zz - 2);
        return ret;
    }
    public void DebugNeiborRC(int _r, int _c)
    {
        DebugSpHMNeiborRC(hm, show_fg_active_ball, _r, _c);
        for (int i = 0; i < hm.sub_hm_list.Count && i < sub_hmv_list.Count; i++)
            DebugSpHMNeiborRC(hm.sub_hm_list[i], sub_hmv_list[i].show_fg_active_ball, _r, _c);
    }
    static void DebugSpHMNeiborRC(HeightMap _hm, List<Transform> _show_fg_tfms, int _r, int _c)
    {
        if (_hm == null) return;
        HBall hb = _hm.fg_ball[_r, _c];
        if (hb == null || !hb.is_active) return;
        Debug.Log("hm.sub_id=" + _hm.sub_hm_id + ", hb[" + _r + ", " + _c + "].neibor cnt:"
            + hb.neibor_lod0.neibors.Count + ", " + hb.neibor_lod1.neibors.Count + ", " + hb.neibor_lod2.neibors.Count
            + "\n\t neibor_lod0: init_cw=" + hb.neibor_lod0.init_cw.ToString("F1") + ", dir2cw:"
            );
        //Transform 
        for (int i = 0; i < _hm.fg_active_ball.Count && i < _show_fg_tfms.Count; i++)
        {
            HBall chb = _hm.fg_active_ball[i];
            _show_fg_tfms[i].GetComponent<MeshRenderer>().material = _black_mat;
            if (chb == hb)
                _show_fg_tfms[i].GetComponent<MeshRenderer>().material = _red_mat;
        }
        for (int ni = 0; ni < hb.neibor_lod0.neibors.Count; ni++)
        {
            HBall nhb = hb.neibor_lod0.neibors[ni];
            for (int i = 0; i < _hm.fg_active_ball.Count && i < _show_fg_tfms.Count; i++)
            {
                HBall chb = _hm.fg_active_ball[i];
                if (chb == nhb)
                {
                    _show_fg_tfms[i].GetComponent<MeshRenderer>().material = _green_mat;
                    break;
                }
            }
        }
        for (int ni = 0; ni < hb.neibor_lod1.neibors.Count; ni++)
        {
            HBall nhb = hb.neibor_lod1.neibors[ni];
            for (int i = 0; i < _hm.fg_active_ball.Count && i < _show_fg_tfms.Count; i++)
            {
                HBall chb = _hm.fg_active_ball[i];
                if (chb == nhb)
                {
                    _show_fg_tfms[i].GetComponent<MeshRenderer>().material = _blue_mat;
                    break;
                }
            }
        }
        for (int ni = 0; ni < hb.neibor_lod2.neibors.Count; ni++)
        {
            HBall nhb = hb.neibor_lod2.neibors[ni];
            for (int i = 0; i < _hm.fg_active_ball.Count && i < _show_fg_tfms.Count; i++)
            {
                HBall chb = _hm.fg_active_ball[i];
                if (chb == nhb)
                {
                    _show_fg_tfms[i].GetComponent<MeshRenderer>().material = _yellow_mat;
                    break;
                }
            }
        }
    }

    //"0-keep-shape", "1-covered", "2-match-back", "3-break-sca", "4-sub-hm-id", "5-line-info"
    public void DebugBgColor(int color_p, int drow, int dcol)
    {
        if (hm == null) return;
        foreach (var sub_hmv in sub_hmv_list)
            sub_hmv.DebugBgColor(color_p, drow, dcol);
        ApplyMeshByH(bg_mesh, hm.bgH, hm.bg_ball, hm.bg_sca, color_p);
        if (color_p == 205)
        {

        }
        else if (color_p == 203)
        {
            HBall bgb = hm.bg_ball[drow, dcol];
            if (bgb == null || !bgb.is_active)
                return;
            StringBuilder sbd = new StringBuilder();
            sbd.Append(bgb.GetDebugStr(-1));

            Vector3 self_init_pos = bgb.GetInitPos();
            int nhb_cnt = bgb.neibor8.Length;
            for (int i = 0; i < nhb_cnt; i++)
            {
                HBall _nhb = bgb.neibor8[i];
                //HBall _nhb = bgb.neibor_lod0.neibors[i];
                if (_nhb == null || !_nhb.IsValidVal() || _nhb.match_self_fg == null) continue;

                Vector3 dir0 = _nhb.GetInitPos() - self_init_pos;
                Vector3 dir1 = _nhb.match_self_fg.tfm_pos - bgb.match_self_fg.tfm_pos;
                dir0.z = dir1.z = 0;
                if (dir0.magnitude > 0)
                {
                    float sca = Mathf.Clamp(dir1.magnitude / dir0.magnitude, 0.25f, 4);
                    sbd.Append("\n\t sca[" + i + "/" + nhb_cnt + "] :" + sca.ToString("F2"));
                    sbd.Append(",\t dir0/1:" + dir0.ToString("F2") + ", " + dir1.ToString("F2"));
                    sbd.Append(",\t nhb.m-fg:[" + _nhb.match_self_fg.row + ", " + _nhb.match_self_fg.col + "].tfm_pos=" + _nhb.match_self_fg.tfm_pos.ToString("F2"));
                }
            }
            Debug.Log(sbd.ToString());
        }
    }
    public void DebugFgColor(int color_p, int drow, int dcol)
    {
        if (hm == null) return;
        foreach (var sub_hmv in sub_hmv_list)
            sub_hmv.DebugFgColor(color_p, drow, dcol);

        if (color_p == 105)//fg line-info
        {
            for (int i = 0; i < hm.fg_active_ball.Count; i++)
            {
                int cuid = hm.fg_active_ball[i].UID();
                MeshRenderer mrdr = show_fg_active_ball[i].GetComponent<MeshRenderer>();
                if (hm.fg_active_ball[i].IsValidVal())
                {
                    int lbdi = hm.fg_active_ball[i].dep_id;
                    mrdr.material = _white_mat;// mats7[lbdi % mats7.Length];
                }
                else
                {
                    mrdr.material = _black_gray_mat;
                }

                switch (hm.fg_active_ball[i].temp_id)
                {
                    case 1: mrdr.material = _black_mat; break;
                    case 2: mrdr.material = _cyan_mat; break;
                    case 3: mrdr.material = _magenta_mat; break;
                    //case 1000: mrdr.material = _red_mat; break;
                    //case 1001: mrdr.material = _green_mat; break;
                    default: mrdr.material = _white_mat; break;
                }

                //if (hm.fg_active_ball[i].dir_t.magnitude < 0.1f)
                //    mrdr.material = _orange_mat;
            }
            SetActiveBallByTfmPos(hm.fg_active_ball, show_fg_active_ball);
        }
        else if (color_p == 104)//sub-hm-id
        {
            //hm.HMUpdateSubHMActiveFg();
            //hm.NextStep(0, 0, 0, 0);
            hm.EstimateAllLoss();

            //ShowHeightMap();
        }
        else if (color_p == 100)//0-keep-shape
        {
            for (int i = 0; i < hm.fg_active_ball.Count; i++)
            {
                var fgb = hm.fg_active_ball[i];
                var p1 = fgb.neibor_lod1.GetTarGridPos();
                var p2 = fgb.neibor_lod2.GetTarGridPos();
                var pm = Vector3.Lerp(p1, p2, 0.5f);
                show_fg_active_ball[i].localPosition = pm;
            }
            float sum_sub_keep_shape_loss = 0;
            foreach (var sub_hm in hm.sub_hm_list)
            {
                sum_sub_keep_shape_loss += sub_hm.hm_match_loss.loss_fg_keep_shape;
                Debug.Log("sub-hm[" + sub_hm.sub_hm_id + "]"
                    + sub_hm.hm_match_loss.GetMatchStr());
            }
            Debug.Log("self keep-shape=" + hm.hm_match_loss.loss_fg_keep_shape.ToString("F3")
                + ", check-sum-sub-hm-keep-shape = " + sum_sub_keep_shape_loss.ToString("F3")
                + ", self.sum-sub-hm-keep-shape = " + hm.hm_match_loss.loss_sub_fg_keep_shape.ToString("F3")
                );

        }

        ApplyMeshByH(fg_mesh, hm.fgH, hm.fg_ball, hm.fg_sca, color_p);
    }
}
public class ViewHeightMapPair
{
    HeightMapPair hmp;
    public ViewHeightMap hmv1;
#if def_hm2
    public ViewHeightMap hmv2;
#endif
    public void CreateHeightMapTfm(HeightMapPair _hmp, string desc = "Default", float ox = 0, float oy = 0)
    {
        hmp = _hmp;
        if (hmv1 == null)
            hmv1 = new ViewHeightMap();
#if def_hm2
        if (hmv2 == null)
            hmv2 = new ViewHeightMap();
#endif
        hmv1.CreateHeightMapTfm(_hmp.hm1, desc + "_fwd");
#if def_hm2
        hmv2.CreateHeightMapTfm(_hmp.hm2, desc + "_bwd");
        hmv2.view_root.transform.parent = hmv1.view_root.transform;
#endif

        root_pos = new Vector3(ox, oy, 0) * 30;
        hmv1.view_root.transform.position = root_pos;
#if def_hm2
        hmv2.view_root.transform.position = root_pos - new Vector3(0, HeightMap.ih1, 0);
#endif
    }
    public Vector3 root_pos = Vector3.zero;
    public bool skip_show_hm = false;
    public void ShowHeightMap()
    {
        if (hmp == null
#if def_hm2
            || hmv2 == null
#endif
            || hmv1 == null || skip_show_hm) return;
        if (hmp.failed_to_best)
        {
            hmv1.view_root.SetActive(false);
#if def_hm2
            hmv2.view_root.SetActive(false);
#endif
            return;
        }
        else
        {
            hmv1.view_root.SetActive(true);
            hmv1.ShowHeightMap();
#if def_hm2
            hmv2.view_root.SetActive(true);
            hmv2.ShowHeightMap();
#endif
        }
    }
}
public class MNIST_HMP
{
    public int gt_lbl, font_id;
    public HeightMapPair hmp = null;
    public ViewHeightMapPair hmpv = null;
    string desc = "font_hm";
    public MNIST_HMP(int _gt_lbl, int _font_id)
    {
        gt_lbl = _gt_lbl;
        font_id = _font_id;
        desc = "font_hm_" + gt_lbl.ToString() + "_" + font_id.ToString();
    }

    public void Create()
    {
        if (hmp == null)
            hmp = new HeightMapPair();
        if (hmpv == null)
            hmpv = new ViewHeightMapPair();
#if def_hm2
        hmpv.CreateHeightMapTfm(hmp, desc, font_id * 3.1f, gt_lbl * 2.1f);
#else
        hmpv.CreateHeightMapTfm(hmp, desc, font_id * 3.1f, gt_lbl * 1.05f);
#endif
    }
    public void SetBGFont(Pix28x28 _bg_font)
    {
        //bg_font = _bg_font;
        if (hmp == null)
            hmp = new HeightMapPair();
        hmp.hm1.SetBgH(_bg_font.raw_byte, _bg_font.sub_raw_byte);
#if def_hm2
        hmp.hm2.SetFgH(_bg_font.raw_byte, _bg_font.sub_raw_byte);
        hmpv.CreateHeightMapTfm(hmp, desc, font_id * 3.1f, gt_lbl * 2.1f);
#else
        hmpv.CreateHeightMapTfm(hmp, desc, font_id * 3.1f, gt_lbl * 1.05f);
#endif
        hmpv.ShowHeightMap();
    }
    public void SetFGHW(Pix28x28 _fg_hw)
    {
        if (hmp == null)
            hmp = new HeightMapPair();
        hmp.hm1.SetFgH(_fg_hw.raw_byte);
#if def_hm2
        hmp.hm2.SetBgH(_fg_hw.raw_byte, null);
        hmpv.CreateHeightMapTfm(hmp, desc, font_id * 3.1f, gt_lbl * 2.1f);
#else
        hmpv.CreateHeightMapTfm(hmp, desc, font_id * 3.1f, gt_lbl * 1.05f);
#endif
        hmpv.ShowHeightMap();
    }
    public float Run(float _L0, float _L1, float _L2)
    {
        if (hmp == null)
            return 1e10f;
        hmp.NextStep(_L0, _L1, _L2);
        hmpv.ShowHeightMap();
        return hmp.GetMatchVal();
    }
    public void Refresh()
    {
        if (hmp == null) return;
        hmp.Refresh();
        hmpv.ShowHeightMap();
    }
    public void SetSkipShow(bool _skip_show)
    {
        if (hmpv != null)
            hmpv.skip_show_hm = _skip_show;
    }
}
public class HeightMapBallMB : MonoBehaviour
{
    List<Pix28x28> mnist_items = new List<Pix28x28>();
    List<Pix28x28>[] fonts_items = new List<Pix28x28>[10];
    List<MNIST_HMP>[] hmp_arr = new List<MNIST_HMP>[10];
    HeightMapPair hmp;
    ViewHeightMapPair hmpv;
    //HeightMapPair hmp_gt;
    //ViewHeightMapPair hmpv_gt;

    //=========================================================
    int item_id = 0, prev_item_id = -1, curr_gt_lbl = -1;
    bool b_auto_play = false;
    bool b_only_show_default = false;
    float fall_speed = 0.35f;
    float sort_grid_speed = 0.5f;
    float drag_unmatch_lbd = 0.4f;
    #region debug-gui
    int debug_font_id = 0, debug_font_sub_id = 0, debug_gt_sub_id = 0;
    int debug_mnist_id = 1014;//5724;//6449,8622,7445,1014,3333,2989
    int debug_match_row = 15, debug_match_col = 7, debug_color_policy = 1;
    string[] color_policy_desc = new string[6] {
        "0-keep-shape", "1-covered", "2-match-back", "3-break-sca", "4-sub-hm-id", "5-line-info"
        };
    #endregion
    //=========================================================
    void Start()
    {
        int seed = (int)System.DateTime.Now.Ticks;
        Random.InitState(seed);

        //MathUtils.CheckCalcABCD();
        //HBall.TestGetAB();

        fonts_items = MNISTLoader.LoadStdFonts();
        mnist_items = MNISTLoader.ParseMNISTDataByBytes();
        for (int i = 0; i < fonts_items.Length; i++)
        {
            hmp_arr[i] = new List<MNIST_HMP>();
            for (int j = 0; j < fonts_items[i].Count; j++)
            {
                MNIST_HMP mh = new MNIST_HMP(i, j);
                hmp_arr[i].Add(mh);
                mh.Create();
                mh.SetBGFont(fonts_items[i][j]);
            }
        }
        hmp = new HeightMapPair();
        hmpv = new ViewHeightMapPair();
        //hmp.SetDefaultByte();
        hmpv.CreateHeightMapTfm(hmp, "Default", -3, -2);

        //hmp_gt = new HeightMapPair();
        //hmpv_gt = new ViewHeightMapPair();
        //hmpv_gt.CreateHeightMapTfm(hmp_gt, "Debug-GT", -1, -3);
        LoadFgFromMNIST();

        item_id = PlayerPrefs.GetInt("default_item_id", 0);
        debug_font_id = PlayerPrefs.GetInt("debug_font_id", -1);
        debug_font_sub_id = PlayerPrefs.GetInt("debug_font_sub_id", 0);
        debug_match_row = PlayerPrefs.GetInt("debug_match_row", 15);
        debug_match_col = PlayerPrefs.GetInt("debug_match_col", 15);
        debug_color_policy = PlayerPrefs.GetInt("debug_color_policy", 3);
        b_only_show_default = PlayerPrefs.GetInt("b_only_show_default", 0) == 1;
        ApplyOnlyShowDefault();
    }
    void ApplyOnlyShowDefault()
    {
        for (int i = 0; i < hmp_arr.Length; i++)
            for (int j = 0; j < hmp_arr[i].Count; j++)
            {
                GameObject go = hmp_arr[i][j].hmpv.hmv1.view_root;
                go.SetActive(!b_only_show_default);// || i == curr_gt_lbl);
                if (i == curr_gt_lbl)
                {
                    var _init_pos = hmp_arr[i][j].hmpv.root_pos;
                    if (b_only_show_default)
                        _init_pos.y = 0;
                    go.transform.position = _init_pos;
                }
            }
    }

    private void OnGUI()
    {
        GUILayout.BeginHorizontal("Box");
        if (mnist_items.Count > 0)
        {
            GUILayout.Box("Item:" + item_id + "/" + mnist_items.Count + ":" + curr_gt_lbl, GUILayout.Width(150));
            item_id = int.Parse(GUILayout.TextArea(item_id.ToString(), GUILayout.Width(40)));
            item_id = Mathf.RoundToInt(GUILayout.HorizontalSlider((float)item_id, 0, mnist_items.Count - 1, GUILayout.Width(150)));

            b_only_show_default = GUILayout.Toggle(b_only_show_default, "only-default", GUILayout.Width(50));
        }
        if (GUI.changed)
            LoadFgFromMNIST();
        GUILayout.EndHorizontal();

        GUILayout.BeginHorizontal("Box");
        {
            //if (GUILayout.Button("(D)ebugCase" + debug_mnist_id))
            //    DebugCase(debug_mnist_id);
            if (GUILayout.Button("Sort(G)rid"))
            {
                hmp.BallSortGrid();
                hmpv.ShowHeightMap();
            }
            if (GUILayout.Button("Sort(G)rid-Thin"))
            {
                hmp.BallSortGrid();
                hmpv.ShowHeightMap();

            }
            if (GUILayout.Button("(R)efresh"))
                Refresh();
            if (GUILayout.Button("(N)extStep"))
                NextStep();
            b_auto_play = GUILayout.Toggle(b_auto_play, "(A)utPlay");
        }
        GUILayout.EndHorizontal();

        GUILayoutOption optw = GUILayout.Width(55);
        GUILayout.BeginHorizontal("Box");
        {
            GUILayout.Box("Font:", optw);
            if (GUILayout.Button("C"))
                debug_font_id = -1;
            debug_font_id = int.Parse(GUILayout.TextField(debug_font_id.ToString(), GUILayout.Width(25)));
            debug_font_sub_id = int.Parse(GUILayout.TextField(debug_font_sub_id.ToString(), GUILayout.Width(25)));
            debug_gt_sub_id = int.Parse(GUILayout.TextField(debug_gt_sub_id.ToString(), GUILayout.Width(25)));
            if (GUILayout.Button("DebugAll:" + debug_all_fr + "-->" + debug_all_to))
                StartCoroutine(DebugAll());
            if (GUILayout.Button("StopDebug"))
                StopDebugAll();
        }
        GUILayout.EndHorizontal();

        GUILayout.BeginHorizontal("Box");
        {
            if (GUILayout.Button("d-" + item_id))
                DebugCase(item_id);
            int[] dbg_arr = new int[10] {  7, 8, 18, 48, 56, 59, 61, 65, 69, 80 };
            for (int kk = 0; kk < dbg_arr.Length; kk++)
                if (GUILayout.Button("" + dbg_arr[kk])) DebugCase(dbg_arr[kk]);
        }
        GUILayout.EndHorizontal();

        GUILayout.BeginHorizontal("Box");
        {
            //if (GUILayout.Button("(F)allDown"))
            //    hmp.BallFallDown();

            //if (GUILayout.Button("Update(H)"))
            //{
            //    hmp.UpdateHbyBall();
            //    hmpv.ShowHeightMap();
            //}
            float best_match_val = 10000;
            float best_gt_val = 10000;
            int best_font_id = -1;
            int best_gt_font_id = -1;
            int best_match_lbl = BestMatchFont(ref best_match_val, ref best_font_id, ref best_gt_val, ref best_gt_font_id, curr_gt_lbl);
            GUILayout.Box("itr[" + hmp.i_run_itr.ToString()
                + "] : " + hmp.GetMatchVal().ToString("F1")
                + "[" + hmp.hm1.hm_match_loss.GetMatchStr("F2", true)
#if def_hm2
                + "\n\t+" + hmp.hm2.fg_match_loss.GetMatchStr()
#endif
                + "]", GUILayout.Width(300));
            Color uic = GUI.color;
            if (curr_gt_lbl == best_match_lbl)
                GUI.color = Color.green;
            if (GUILayout.Button("item[" + item_id
                + "].gt[" + curr_gt_lbl + "/" + best_gt_font_id + "], val=" + best_gt_val.ToString("F1")
                + ", \n best match[" + best_match_lbl + "/" + best_font_id
                + "], val=" + best_match_val.ToString("F1"), GUILayout.Width(192)))
                DebugMatchValue();
            GUI.color = uic;
        }
        GUILayout.EndHorizontal();

        GUILayout.BeginHorizontal("Box");
        {
            GUILayout.Box("fall-speed=" + fall_speed.ToString("F2"), GUILayout.Width(115));
            fall_speed = GUILayout.HorizontalSlider(fall_speed, 0.1f, 1.5f, GUILayout.Width(65));
            GUILayout.Box("grid-speed=" + sort_grid_speed.ToString("F2"), GUILayout.Width(115));
            sort_grid_speed = GUILayout.HorizontalSlider(sort_grid_speed, 0.001f, 1.0f, GUILayout.Width(65));
            GUILayout.Box("drag-unmatch=" + drag_unmatch_lbd.ToString("F2"), GUILayout.Width(115));
            drag_unmatch_lbd = GUILayout.HorizontalSlider(drag_unmatch_lbd, 0.1f, 1.0f, GUILayout.Width(65));
        }
        GUILayout.EndHorizontal();

        GUILayout.BeginHorizontal("Box");
        {
            if (GUILayout.Button("dbg-neibor"))
                hmpv.hmv1.DebugNeiborRC(debug_match_row, debug_match_col);
            GUILayout.Box("debug-row-col");
            debug_match_row = int.Parse(GUILayout.TextArea(debug_match_row.ToString(), GUILayout.Width(35)));
            debug_match_col = int.Parse(GUILayout.TextArea(debug_match_col.ToString(), GUILayout.Width(35)));
        }
        GUILayout.EndHorizontal();
        GUILayout.BeginHorizontal("Box");
        {
            GUILayout.Box("debug-color:" + color_policy_desc[debug_color_policy]);
            debug_color_policy = Mathf.RoundToInt(GUILayout.HorizontalSlider(debug_color_policy, 0, color_policy_desc.Length - 0.51f, GUILayout.Width(100)));
            if (GUILayout.Button("DebugFg:"))
                DebugFg(debug_match_row, debug_match_col, debug_color_policy);
            if (GUILayout.Button("DebugBg:"))
                DebugBg(debug_match_row, debug_match_col, debug_color_policy);
        }
        GUILayout.EndHorizontal();
        GUILayout.BeginHorizontal("Box");
        {
            if (GUILayout.Button(new GUIContent("Itr0", "1.fgb to nearest bg\n 2.sort-grid X 3\n 3.sub-hm search in range-3")))
            {
                hmp.hm1.RunItr0();// 
                hmpv.ShowHeightMap();
                //hmpv.hmv1.FlatTfmZ();
            }
            if (GUILayout.Button(new GUIContent("Itr1", "sub-hm NextStep, sort-unactive, det-sub-id")))
            {
                hmp.hm1.RunItr1();
                hmpv.ShowHeightMap();
            }
            if (GUILayout.Button(new GUIContent("Itr1-LT", "sub-hm find Linear Transform")))
            {
                hmp.hm1.RunItr1_LT();// sub hm fit fg2bg Linear-Transform
                hmpv.ShowHeightMap();
            }

            if (GUILayout.Button(new GUIContent("Itr2", "exclude-sub-id by depth/reverse-depth")))
            {
                hmp.hm1.RunItr2();
                hmpv.ShowHeightMap();
                foreach(var sub_hmv in hmpv.hmv1.sub_hmv_list)
                    sub_hmv.ShowTempID();
            }
            if (GUILayout.Button(new GUIContent("Itr3:merge-pos", "Merge tfm pos")))
            {
                hmp.hm1.RunItr3(3, fall_speed, sort_grid_speed, drag_unmatch_lbd);
                hmpv.ShowHeightMap();
            }
            if (GUILayout.Button("All"))
            {
                hmp.hm1.RunItr0123();
                hmpv.ShowHeightMap();
                DebugFg(debug_match_row, debug_match_col, debug_color_policy);
            }
            if (GUILayout.Button("thin/Init"))
            {
                hmp.hm1.RunThin(true);
                hmpv.ShowHeightMap();
                hmpv.hmv1.ShowLineInfo();
            }
            if (GUILayout.Button("thin/RT"))
            {
                hmp.hm1.RunThin(false);
                hmpv.ShowHeightMap();
                hmpv.hmv1.ShowLineInfo();
            }
        }
        GUILayout.EndHorizontal();

        //GUILayout.BeginHorizontal("Box");
        //{
            
        //}
        //GUILayout.EndHorizontal();
        
        if (GUI.changed)
        {
            PlayerPrefs.SetInt("default_item_id", item_id);
            PlayerPrefs.SetInt("debug_font_id", debug_font_id);
            PlayerPrefs.SetInt("debug_font_sub_id", debug_font_sub_id);
            PlayerPrefs.SetInt("debug_match_row", debug_match_row);
            PlayerPrefs.SetInt("debug_match_col", debug_match_col);
            PlayerPrefs.SetInt("debug_color_policy", debug_color_policy);
            PlayerPrefs.SetInt("b_only_show_default", b_only_show_default ? 1 : 0);
            ApplyOnlyShowDefault();
        }
        if (GUI.tooltip.Length > 1)
        {
            GUI.Box(new Rect(Input.mousePosition.x - 150, Screen.height - Input.mousePosition.y - 50, 200, 50), GUI.tooltip);
        }
    }
    void Update()
    {
        if (Input.GetKeyDown(KeyCode.C))
            MathUtils.CheckCalcABCD();
        if (Input.GetKeyDown(KeyCode.A) || Input.GetKeyDown(KeyCode.P))
            b_auto_play = !b_auto_play;
        if (Input.GetKeyDown(KeyCode.N) || Input.GetKeyDown(KeyCode.Space) || b_auto_play)
            NextStep();
        if (Input.GetKeyDown(KeyCode.D))
            DebugCase(item_id);
        if (Input.GetKeyDown(KeyCode.LeftControl))
            for (int i = 0; i < 10; i++)
                NextStep();
        if (Input.GetKeyDown(KeyCode.LeftArrow))
        {
            item_id--;
            LoadFgFromMNIST();
        }
        if (Input.GetKeyDown(KeyCode.RightArrow))
        {
            item_id++;
            LoadFgFromMNIST();
        }

        if (Input.GetKeyDown(KeyCode.F))
            hmp.BallFallDown();
        if (Input.GetKeyDown(KeyCode.G))
            hmp.BallSortGrid();
        if (Input.GetKeyDown(KeyCode.R)) Refresh();
    }
    public int debug_all_fr = 0;
    public int debug_all_to = 100;
    StringBuilder debug_all_sbd = new StringBuilder();
    IEnumerator DebugAll(int each_itr = 10)
    {
        debug_all_sbd = new StringBuilder();
        SetSkipShowHM(true);
        StringBuilder sbd = debug_all_sbd;// new StringBuilder();
        //for (int i = 0; i < mnist_items.Count; i++)
        int bad_ans_cnt = 0;
        float t0 = Time.realtimeSinceStartup;
        for (int i = debug_all_fr; i < debug_all_to; i++)
        {
            item_id = i;
            LoadFgFromMNIST();
            for (int itr = 0; itr < each_itr; itr++)
                NextStep();
            float best_match_val = 10000;
            float best_gt_val = 10000;
            int best_font_id = -1;
            int best_gt_font_id = -1;
            int best_match_lbl = BestMatchFont(ref best_match_val, ref best_font_id, ref best_gt_val, ref best_gt_font_id, mnist_items[i].gt_lbl);
            if (best_match_lbl != mnist_items[i].gt_lbl)
            {
                bad_ans_cnt++;
                var gt_hmp = hmp_arr[mnist_items[i].gt_lbl][best_gt_font_id].hmp;
                var best_hmp = hmp_arr[best_match_lbl][best_font_id].hmp;
                sbd.Append("\n Item[" + item_id
                    + "]:\n\tgt-ans:[" + mnist_items[i].gt_lbl + "/" + best_gt_font_id
                    + "]\t(" + best_gt_val.ToString("F2") + "[" + gt_hmp.hm1.hm_match_loss.GetMatchStr() + "]"
                    + ")  \n\tbest-m:[" + best_match_lbl + "/" + best_font_id
                    + "]\t(" + best_match_val.ToString("F2") + "[" + best_hmp.hm1.hm_match_loss.GetMatchStr() + "]"
                    + ")");
            }
            yield return 0;
            if ((i % 10) == 0 && i > 0)
            {
                float t2 = Time.realtimeSinceStartup;
                Debug.Log("debug " + i + " / [" + debug_all_fr + "--->" + debug_all_to
                    + "] time cost :" + (t2 - t0).ToString("F3")
                    + ", : bad ans cnt =" + bad_ans_cnt);
                yield return new WaitForSeconds(0.1f);
            }
        }
        float t1 = Time.realtimeSinceStartup;
        Debug.Log("[" + debug_all_fr + "--->" + debug_all_to
            + "] time cost :" + (t1 - t0).ToString("F3") + ", bad-ans-cnt=" + bad_ans_cnt +
            sbd.ToString());
        SetSkipShowHM(false);
        yield return 0;
    }
    void StopDebugAll()
    {
        Debug.Log(debug_all_sbd.ToString());
        StopCoroutine(DebugAll());
        SetSkipShowHM(false);
    }
    void SetSkipShowHM(bool _skip_show)
    {
        hmpv.skip_show_hm = _skip_show;
        for (int i = 0; i < hmp_arr.Length; i++)
            for (int j = 0; j < hmp_arr[i].Count; j++)
                hmp_arr[i][j].SetSkipShow(_skip_show);
    }
    void DebugBg(int drow, int dcol, int color_p)
    {
        if (hmp == null || hmp.hm1 == null
           || hmp.hm1.fg_active_ball.Count == 0
            )
            return;
        color_p += 200;
        if (drow < 0) drow = 0;
        if (dcol < 0) dcol = 0;
        if (drow > HeightMap.ih1) drow = HeightMap.ih1;
        if (dcol > HeightMap.iw1) dcol = HeightMap.iw1;
        HBall bg_src = hmp.hm1.bg_ball[drow, dcol];
        hmpv.hmv1.DebugBgColor(color_p, drow, dcol);
        if (bg_src == null) return;


        var mat_fg = bg_src.match_self_fg;
        var root_pos = hmpv.hmv1.quat_bg.position;
        Debug.DrawLine(root_pos + new Vector3(0, drow, 0), root_pos + new Vector3(HeightMap.iw1, drow, 0), Color.white, 3);
        Debug.DrawLine(root_pos + new Vector3(dcol, 0, 0), root_pos + new Vector3(dcol, HeightMap.ih1, 0), Color.white, 3);

        if (mat_fg != null)
            Debug.DrawLine(root_pos + new Vector3(dcol, drow, 0), root_pos + mat_fg.tfm_pos, Color.green, 5);
        else
        {
            Debug.DrawLine(root_pos + new Vector3(dcol + 1, drow - 1, 0), root_pos + new Vector3(dcol - 1, drow + 1, 0), Color.green, 3);
            Debug.DrawLine(root_pos + new Vector3(dcol - 1, drow + 1, 0), root_pos + new Vector3(dcol + 1, drow - 1, 0), Color.green, 3);
        }
    }

    void DebugFg(int drow, int dcol, int color_p, bool reverse = false)
    {
        color_p += 100;
        if (hmp == null
            || hmp.hm1 == null
            || hmp.hm1.fg_active_ball.Count == 0
#if def_hm2
            || hmp.hm2 == null 
            || hmp.hm2.fg_active_ball.Count == 0
#endif
            ) return;
        if (drow < 0) drow = 0;
        if (dcol < 0) dcol = 0;
        if (drow > HeightMap.ih1) drow = HeightMap.ih1;
        if (dcol > HeightMap.iw1) dcol = HeightMap.iw1;
        HBall fg_src = hmp.hm1.fg_ball[drow, dcol];
        hmpv.hmv1.DebugFgColor(color_p, drow, dcol);
        if (fg_src == null) return;
        StringBuilder sbd = new StringBuilder();
        sbd.Append("debug fg: " + fg_src.GetDebugStr(1));
        foreach (var sub_hm in hmpv.hmv1.hm.sub_hm_list)
            sbd.Append("\n======== sub_hm["+ sub_hm.sub_hm_id+"].fg: " + sub_hm.fg_ball[drow, dcol].GetDebugStr(1));
        Debug.Log(sbd.ToString());

    }
    void DrawLineFromTo(Vector3 fr, Vector3 to, Color c, float lt = 15)
    {
        fr -= Vector3.forward * 0.1f;
        to -= Vector3.forward * 0.1f;
        Vector3 mid = (fr + to) * 0.5f;
        mid.z -= 5;
        Debug.DrawLine(fr, mid, c, lt);
        Debug.DrawLine(mid, to, c, lt);

        Debug.DrawLine(mid, to + Vector3.up * 0.05f, c, lt);
        Debug.DrawLine(mid, fr + Vector3.up * 0.05f, c, lt);
        Debug.DrawLine(mid, to - Vector3.up * 0.05f, c, lt);
        Debug.DrawLine(mid, fr - Vector3.up * 0.05f, c, lt);
        Debug.DrawLine(mid, to + Vector3.right * 0.05f, c, lt);
        Debug.DrawLine(mid, fr + Vector3.right * 0.05f, c, lt);
        Debug.DrawLine(mid, to - Vector3.right * 0.05f, c, lt);
        Debug.DrawLine(mid, fr - Vector3.right * 0.05f, c, lt);
    }
    void DebugMatchValue()
    {
        StringBuilder sbd = new StringBuilder("Debug Match value:");
        List<float> all_mv = new List<float>();
        List<int> all_mv_lbl = new List<int>();
        List<int> all_mv_fid = new List<int>();
        float min_mv = 1e10f;
        int min_mv_id = -1;
        for (int i = 0; i < 10; i++)
        {
            List<MNIST_HMP> hmlist = hmp_arr[i];
            sbd.Append("\n");
            for (int j = 0; j < hmlist.Count; j++)
            {
                float cmv = hmlist[j].hmp.GetMatchVal();
                all_mv.Add(cmv);
                all_mv_lbl.Add(i);
                all_mv_fid.Add(j);
                if (min_mv > cmv)
                {
                    min_mv = cmv;
                    min_mv_id = all_mv.Count - 1;
                }
                sbd.Append("\t lbl[" + i + "/" + j + "] : "
                    + (cmv > 10000 ? "inf" : cmv.ToString("F1"))
                    + "[" + hmlist[j].hmp.hm1.hm_match_loss.GetMatchStr()
#if def_hm2
                    + "+" + hmlist[j].hmp.hm2.fg_match_loss.GetMatchStr()
#endif
                    + "],\t ");
            }
        }
        if (min_mv_id >= 0)
        {
            int min_mv_lbl = all_mv_lbl[min_mv_id];
            int min_mv_fid = all_mv_fid[min_mv_id];
            HeightMapPair _minv_hmp = hmp_arr[min_mv_lbl][min_mv_fid].hmp;
            sbd.Append("\n\t item[" + item_id + "].gt_lbl[" + curr_gt_lbl + "].itr[" + _minv_hmp.i_run_itr
                + "] best match val = " + min_mv.ToString("F1")
                + "[" + min_mv_lbl + "/" + min_mv_fid + "], candidates:");
            for (int i = 0; i < all_mv.Count; i++)
            {
                if (all_mv[i] < (min_mv * 1.8f))
                {
                    int _lbl = all_mv_lbl[i];
                    int _fid = all_mv_fid[i];
                    HeightMapPair _hmp = hmp_arr[_lbl][_fid].hmp;
                    sbd.Append("\n\t\t" + all_mv[i].ToString("F1")
                        + "[" + _lbl + "/" + _fid + "]:\t "
                        + _hmp.hm1.hm_match_loss.GetMatchStr()
#if def_hm2
                        + ", " + _hmp.hm2.fg_match_loss.GetMatchStr()
#endif
                        );
                }
            }
        }
        Debug.Log(sbd.ToString());
    }
    int BestMatchFont(ref float mv, ref int font_id, ref float best_gt_val, ref int best_gt_font_id, int _gt_lbl)
    {
        mv = 100000;
        best_gt_val = 100000;
        int best_lbl = -1;
        font_id = -1;
        best_gt_font_id = -1;
        for (int i = 0; i < 10; i++)
        {
            List<MNIST_HMP> hmlist = hmp_arr[i];
            for (int j = 0; j < hmlist.Count; j++)
            {
                float cmv = hmlist[j].hmp.GetMatchVal();
                if (mv > cmv)
                {
                    mv = cmv;
                    best_lbl = i;
                    font_id = j;
                }
                if (i == _gt_lbl && best_gt_val > cmv)
                {
                    best_gt_val = cmv;
                    best_gt_font_id = j;
                }
            }
        }
        return best_lbl;
    }
    void DebugCase(int _tar_id)
    {
        if (!(_tar_id >= 0 && _tar_id < mnist_items.Count))
            return;
        item_id = _tar_id;
        LoadFgFromMNIST();

        Pix28x28 hw_item = mnist_items[item_id];
        if (debug_font_id < 0)
            debug_font_id = hw_item.gt_lbl;

        List<Pix28x28> candidates_item = fonts_items[debug_font_id];
        if (debug_font_sub_id >= candidates_item.Count)
            debug_font_sub_id = candidates_item.Count - 1;
        byte[,] ft_byt = candidates_item[debug_font_sub_id].raw_byte;
        List<byte[,]> sub_byt_list = candidates_item[debug_font_sub_id].sub_raw_byte;
        hmp.SetRawByte(ft_byt, hw_item.raw_byte, sub_byt_list);
        hmpv.CreateHeightMapTfm(hmp, "Default", -3, -2);

        List<Pix28x28> candidates_item_gt = fonts_items[curr_gt_lbl];
        if (debug_gt_sub_id >= candidates_item_gt.Count)
            debug_gt_sub_id = candidates_item_gt.Count - 1;
        byte[,] gt_byt = candidates_item_gt[debug_gt_sub_id].raw_byte;
        List<byte[,]> gt_sub_byt_list = candidates_item_gt[debug_gt_sub_id].sub_raw_byte;
        //hmp_gt.SetRawByte(gt_byt, hw_item.raw_byte, gt_sub_byt_list);
        //hmpv_gt.CreateHeightMapTfm(hmp_gt, "Debug-GT", 0, -3);

        hmpv.ShowHeightMap();
    }
    void NextStep()
    {
        hmp.NextStep(fall_speed, sort_grid_speed, drag_unmatch_lbd);
        hmpv.ShowHeightMap();

        //hmp_gt.NextStep(fall_speed, sort_grid_speed, drag_unmatch_lbd);
        //hmpv_gt.ShowHeightMap();

        float best_match_val = 1e10f;
        for (int i = 0; i < hmp_arr.Length; i++)
        {
            if (!b_only_show_default)// || i == curr_gt_lbl)
                for (int j = 0; j < hmp_arr[i].Count; j++)
                {
                    float cmv = hmp_arr[i][j].Run(fall_speed, sort_grid_speed, drag_unmatch_lbd);
                    if (best_match_val > cmv)
                        best_match_val = cmv;
                    //if (best_match_val < 15 && hmp.i_run_itr > 25)
                    //{
                    //    float plus_f = HBall.X01Y01(hmp.i_run_itr, 25, 50, 60, 30);
                    //    float mult_f = HBall.X01Y01(hmp.i_run_itr, 25, 50, 20, 5);
                    //    float tar_mv0 = best_match_val + plus_f;
                    //    float tar_mv1 = best_match_val * mult_f;
                    //    if (cmv > tar_mv0 || cmv > tar_mv1)
                    //        hmp_arr[i][j].hmp.failed_to_best = true;
                    //}
                }
        }
    }
    void LoadFgFromMNIST()
    {
        int _cnt = mnist_items.Count;
        if (_cnt == 0 || prev_item_id == item_id)
            return;
        Pix28x28 hw_item = mnist_items[(_cnt + item_id) % _cnt];
        prev_item_id = item_id;
        curr_gt_lbl = hw_item.gt_lbl;

        for (int i = 0; i < hmp_arr.Length; i++)
        {
            if (!b_only_show_default)// || i == curr_gt_lbl)
                for (int j = 0; j < hmp_arr[i].Count; j++)
                {
                    hmp_arr[i][j].SetFGHW(hw_item);
                    hmp_arr[i][j].Refresh();
                }
        }
    }
    void Refresh()
    {
        hmp.Refresh();
        hmpv.ShowHeightMap();
        //hmp_gt.Refresh();
        //hmpv_gt.ShowHeightMap();

        for (int i = 0; i < hmp_arr.Length; i++)
        {
            if (!b_only_show_default)// || i == curr_gt_lbl)
                for (int j = 0; j < hmp_arr[i].Count; j++)
                    hmp_arr[i][j].Refresh();
        }
    }

}
