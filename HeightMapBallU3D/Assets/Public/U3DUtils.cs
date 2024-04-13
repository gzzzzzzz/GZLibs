using UnityEngine;
using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using System.Threading;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;
using Newtonsoft.Json.Converters;
using System.Runtime.InteropServices;
public class FPSCounter
{
    int fr_cnt = 0;
    public float fps = 0;
    float init_time = 0, time_elaps = 0;
    float[] last20_time = new float[20];
    public void Update()
    {
        time_elaps = Time.realtimeSinceStartup - init_time;
        last20_time[fr_cnt % 20] = time_elaps;
        fr_cnt++;
        if (0 == (fr_cnt % 20))
            fps = 20 / (last20_time[(fr_cnt - 1) % 20] - last20_time[fr_cnt % 20]);
    }
}
public class PalmViewer
{
    public PalmViewer(string _desc)
    {
        desc = _desc;
        InitTfm();
    }
    string desc = "";
    float radius_w = 0.01f;
    public Transform[] palm = new Transform[21];
    public Transform[] palm_bone = new Transform[21];
    public Transform[] palm_ball = new Transform[21];
    public void SetPos(Vector3[] _rh_pos, Vector3 _wri_r)
    {
        if (_rh_pos != null && _rh_pos.Length == 21)
        {
            for (int i = 0; i < _rh_pos.Length; i++)
                palm[i].position = (_rh_pos[i] - _rh_pos[0]) + _wri_r;
        }
    }
    public void SetPos(Transform[] _rh_tfm, Vector3 _wri_r)
    {
        if (_rh_tfm != null && _rh_tfm.Length == 21)
        {
            for (int i = 0; i < _rh_tfm.Length; i++)
            {
                palm[i].position = (_rh_tfm[i].position - _rh_tfm[0].position) + _wri_r;
                palm[i].rotation = (_rh_tfm[i].rotation);
            }
        }
    }
    public void InitTfm()
    {
        for (int i = 0; i < 21; i++)
        {
            palm[i] = new GameObject(desc + "_finger_" + i).transform;
            GameObject _b = GameObject.CreatePrimitive(PrimitiveType.Sphere);
            U3DUtils.RemoveCollider(_b);
            palm_ball[i] = _b.transform;
            _b.GetComponent<MeshRenderer>().material = U3DUtils.CreateMaterial("Materials/blue-trans-mat");
            _b.name = palm[i].name + "_ball";
            _b.transform.localScale = Vector3.one * radius_w * 1f;
            _b.transform.parent = palm[i];
            if ((i % 4) == 1)
                palm[i].parent = palm[0];
            else if (i > 0)
                palm[i].parent = palm[i - 1];

            if (i > 0)
            {
                Transform _c = GameObject.CreatePrimitive(PrimitiveType.Cube).transform;
                U3DUtils.RemoveCollider(_c.gameObject);
                Material _mat = U3DUtils.CreateMaterial("Materials/blue-trans-mat");
                if (i < 5)
                    _mat.color *= 4;

                _c.GetComponent<MeshRenderer>().material = _mat;
                _c.name = desc + "_finger_bone_" + i;
                _c.localScale = Vector3.one * radius_w * 0.4f;
                _c.parent = palm[i].parent;
                palm_bone[i] = _c;
            }
        }
    }
    public void UpdateBone()
    {
        for (int i = 1; i < 21; i++)
        {
            int pfg = 0;
            if ((i % 4) != 1)
                pfg = i - 1;
            palm_bone[i].position = 0.5f * (palm[i].position + palm[pfg].position);
            Vector3 _dd = (palm[i].position - palm[pfg].position);
            palm_bone[i].localScale = new Vector3(_dd.magnitude, radius_w * 0.4f, radius_w * 0.4f);
            palm_bone[i].rotation = Quaternion.FromToRotation(Vector3.right, _dd);
        }
    }
    public void SetColor(Color _c)
    {
        for (int i = 1; i < palm.Length; i++)
        {
            Color _b = _c;
            if (i < 5)
                _b = Color.Lerp(_b, Color.black, 0.7f);
            Material _mat = palm_ball[i].GetComponent<MeshRenderer>().material;
            _mat.color = new Color(_b.r, _b.g, _b.b, _mat.color.a);
            _mat = palm_bone[i].GetComponent<MeshRenderer>().material;
            Color _d = _c * 0.8f;
            _mat.color = new Color(_d.r, _d.g, _d.b, _mat.color.a);
        }
    }
}
public class SkelShow
{
   
    public int tfm_cnt = 0;
    public int[] PID = new int[0];
    List<int> porder = new List<int>();
    public Transform root = null;
    public Transform[] tfms = null;
    public Transform[] balls = null;
    Transform[] cubes = null;
    Material mat_ball, mat_bone, mat_red;
    public PalmViewer l_palm, r_palm;
    public PalmViewer[] l_debug_palms = new PalmViewer[5];
    public PalmViewer[] r_debug_palms = new PalmViewer[5];
    int r_elb_id = 1;
    int r_wri_id = 2;
    int l_elb_id = 4;
    int l_wri_id = 5;
    public SkelShow()
    {
    }
    public SkelShow(int _tfm_cnt, int[,] bone_fr_to, string[] _name = null)
    {
        int[] _PID = new int[_tfm_cnt];
        for (int i = 0; i < _tfm_cnt; i++)
            _PID[i] = -1;
        int bone_cnt = bone_fr_to.GetLength(0);
        for (int i = 0; i < bone_cnt; i++)
            _PID[bone_fr_to[i, 0]] = bone_fr_to[i, 1];
        InitTfms(_tfm_cnt, _PID, _name);
    }
    public SkelShow(string root_name, int _tfm_cnt, int[] _PID, string[] _name = null)
    {
        InitTfms(_tfm_cnt, _PID, _name);
        root.name = root_name;
    }
    public void InitTfms(int _tfm_cnt, int[] _PID, string[] _name = null)
    {
        tfm_cnt = _tfm_cnt;
        PID = new int[_tfm_cnt];
        porder.Clear();
        bool[] pflag = new bool[_tfm_cnt];
        for (int i = 0; i < _tfm_cnt; i++)
        {
            PID[i] = (i < _PID.Length) ? _PID[i] : -1;
            pflag[i] = false;
        }
        for (int i = 0; i < _tfm_cnt; i++)
            AppendPOrder(i, pflag);
        if (porder.Count != tfm_cnt)
            Debug.Log("=========== bad p order:" + porder.Count);

        tfms = new Transform[tfm_cnt];
        balls = new Transform[tfm_cnt];
        cubes = new Transform[tfm_cnt];
        root = new GameObject("skel_root").transform;
        for (int i = 0; i < tfm_cnt; i++)
            tfms[i] = new GameObject().transform;
        for (int i = 0; i < tfm_cnt; i++)
        {
            string tfm_name = _name == null ? "tfm_" + i : _name[i];
            tfms[i].name = tfm_name;
            string b_name = _name == null ? "ball_" + i : _name[i];
            balls[i] = CreateGO(PrimitiveType.Sphere, Vector3.one * 0.03f, b_name);
            balls[i].parent = tfms[i];

            int cpid = PID[i];
            if (cpid >= 0)
            {
                tfms[i].parent = tfms[cpid];
                string c_name = _name == null ? "bone_" + i + "_" + cpid : _name[i] + _name[cpid];
                cubes[i] = CreateGO(PrimitiveType.Cube, Vector3.one * 0.02f, c_name);
                cubes[i].parent = tfms[cpid];
            }
            else
                tfms[i].parent = root;
        }
        for (int i = 0; i < tfms.Length; i++)
        {
            int _pid = PID[i];
            if (tfms[i] != null && _pid >= 0 && _pid < tfms.Length)
                tfms[i].parent = tfms[_pid];
        }

        l_palm = new PalmViewer("L_palm_" + root.name);
        r_palm = new PalmViewer("R_palm_" + root.name);
        l_palm.palm[0].parent = tfms[l_elb_id];
        r_palm.palm[0].parent = tfms[r_elb_id];
        l_palm.palm[0].gameObject.SetActive(false);
        r_palm.palm[0].gameObject.SetActive(false);
        for (int k = 0; k < 5; k++)
        {
            l_debug_palms[k] = new PalmViewer("L_debug_palm_" + k);
            r_debug_palms[k] = new PalmViewer("R_debug_palm_" + k);
            l_debug_palms[k].palm[0].parent = tfms[l_elb_id];
            r_debug_palms[k].palm[0].parent = tfms[r_elb_id];
            l_debug_palms[k].palm[0].gameObject.SetActive(false);
            r_debug_palms[k].palm[0].gameObject.SetActive(false);
        }
    }
    void AppendPOrder(int head, bool[] pflag)
    {
        if (head < 0 || head >= tfm_cnt || pflag[head])
            return;
        int pi = PID[head];
        AppendPOrder(pi, pflag);
        pflag[head] = true;
        porder.Add(head);
    }
    void UnParent()
    {
        for (int i = 0; i < tfms.Length; i++)
            if (tfms[i] != null)
                tfms[i].parent = root;
    }
    void SetParent()
    {
        for (int i = 0; i < tfms.Length; i++)
            if (tfms[i] != null && PID[i] >= 0)
                tfms[i].parent = tfms[PID[i]];
    }
    Transform CreateGO(PrimitiveType type, Vector3 local_scale, string name = "_gameObject")
    {
        Transform ret = GameObject.CreatePrimitive(type).transform;
        U3DUtils.RemoveCollider(ret.gameObject);
        ret.localScale = local_scale;
        ret.name = name;
        return ret;
    }
    public void UpdateLimb()
    {
        for (int i = 0; i < tfm_cnt; i++)
        {
            int _pid = PID[i];
            if (cubes[i] != null && _pid >= 0 && _pid < tfms.Length)
            {
                //if (!(tfms[i].gameObject.activeInHierarchy && tfms[_pid].gameObject.activeInHierarchy))
                if (!(tfms[i].gameObject.activeSelf && tfms[_pid].gameObject.activeSelf))
                {
                    cubes[i].gameObject.SetActive(false);
                    continue;
                }
                cubes[i].gameObject.SetActive(true);
                Vector3 dir = tfms[i].position - tfms[_pid].position;
                Vector3 mid = (tfms[i].position + tfms[_pid].position) * 0.5f;
                cubes[i].position = mid;
                cubes[i].rotation = Quaternion.FromToRotation(Vector3.right, dir);
                float blen = Mathf.Min(0.5f, dir.magnitude);
                cubes[i].localScale = new Vector3(dir.magnitude, blen * 0.05f, blen * 0.05f);
            }
        }
    }
    public void SetColor(Color _c)
    {
        if(mat_ball == null)            mat_ball = U3DUtils.CreateMaterial("");
        if(mat_bone == null)            mat_bone = U3DUtils.CreateMaterial("");
        mat_bone.color = _c;
        mat_ball.color = 0.5f * _c;
        if (balls != null)
            foreach (var _b in balls)
                if (_b != null)
                    _b.GetComponent<MeshRenderer>().material = mat_ball;
        if (cubes != null)
            foreach (var _b in cubes)
                if (_b != null)
                    _b.GetComponent<MeshRenderer>().material = mat_bone;
    }
    public void SetPalmColor()
    {
        Color[] palm_color = new Color[5] { Color.gray, Color.blue, Color.green, Color.yellow, Color.red };
        if (r_palm != null)
            r_palm.SetColor(Color.blue + Color.green);
        if (l_palm != null)
            l_palm.SetColor(Color.blue + Color.green);
        if (l_debug_palms != null)
            for (int i = 0; i < l_debug_palms.Length; i++)
                l_debug_palms[i].SetColor(palm_color[i % 5]);
        if (r_debug_palms != null)
            for (int i = 0; i < r_debug_palms.Length; i++)
                r_debug_palms[i].SetColor(palm_color[i % 5]);
    }
    public void SetColorVisible(int id, float visible_val, float thresh = 0.6f)
    {
        SetColorVisible(id, visible_val < thresh);
    }
    public void SetColorVisible(int id, bool red_or_orig)
    {
        if (mat_ball == null) mat_ball = U3DUtils.CreateMaterial("");
        if (mat_red == null)
        {
            mat_red = U3DUtils.CreateMaterial("");
            mat_red.color = Color.red;
        }
        if (balls != null && id >= 0 && id < balls.Length && balls[id] != null)
            balls[id].GetComponent<MeshRenderer>().material = red_or_orig ? mat_red : mat_ball;
    }
    public void UpdateTfm(Transform[] _tfms)
    {
        root.position = Vector3.zero;
        for(int i =0;  i < tfm_cnt && i < _tfms.Length; i++)
        {
            tfms[i].localPosition = _tfms[i].localPosition;
            tfms[i].localRotation = _tfms[i].localRotation;
        }
        root.position = _tfms[0].root.position;
        UpdateLimb();
    }
    public void UpdateTfm(Vector3[] _lcl_pos, Quaternion[] _lcl_rot, int [] cvt_order = null)
    {
        root.position = Vector3.zero;
        for(int j = 0; j < tfm_cnt; j++)
        {
            int i = porder[j];
            int k = i;
            if (cvt_order != null)
                k = cvt_order[i];
            tfms[i].localPosition = _lcl_pos[k];
            tfms[i].localRotation = _lcl_rot[k];

        }
        UpdateLimb();
    }
    public void UpdatePalmPos(int pj, Vector3[] _rh_wpos, Vector3[] _lh_wpos)
    {
        if (pj >= r_debug_palms.Length)
            return;
        PalmViewer rpv = r_debug_palms[pj];
        PalmViewer lpv = l_debug_palms[pj];
        for (int i = 0; i < 21; i++)
        {
            rpv.palm[i].position = tfms[r_wri_id].position + (_rh_wpos[i] - _rh_wpos[0]);
            lpv.palm[i].position = tfms[l_wri_id].position + (_lh_wpos[i] - _lh_wpos[0]);
        }
        lpv.UpdateBone();
        rpv.UpdateBone();
    }
    public void ShowPalmTfm(int pj, bool b)
    {
        r_palm.palm[0].gameObject.SetActive(b);
        l_palm.palm[0].gameObject.SetActive(b);
    }
    public void ShowPalm(int pj, bool b)
    {
        if (pj >= r_debug_palms.Length)
            return;
        PalmViewer rpv = r_debug_palms[pj];
        PalmViewer lpv = l_debug_palms[pj];
        rpv.palm[0].gameObject.SetActive(b);
        lpv.palm[0].gameObject.SetActive(b);
    }
    public void UpdatePalmTfm(Vector3[] _rh_lcl_pos, Quaternion[] _rh_lcl_rot, Vector3[] _lh_lcl_pos, Quaternion[] _lh_lcl_rot)
    {
        for (int i = 0; i < 21; i++)
        {
            l_palm.palm[i].localPosition = _lh_lcl_pos[i];
            l_palm.palm[i].localRotation = _lh_lcl_rot[i];

            r_palm.palm[i].localPosition = _rh_lcl_pos[i];
            r_palm.palm[i].localRotation = _rh_lcl_rot[i];
        }
        l_palm.UpdateBone();
        r_palm.UpdateBone();


    }
    public void UpdatePos(Vector3[] _pos, int[] cvt_order = null)
    {
        root.position = Vector3.zero;
        for(int j = 0; j < tfm_cnt; j++)
        {
            int i = porder[j];
            int k = i;
            if (cvt_order != null)
                k = cvt_order[i];
            tfms[i].position = _pos[k];
        }
        UpdateLimb();
    }
}
public class ViewFrameDataSln<T>
{
    public List<T> fd_list = new List<T>();
    public int frame_cnt = 0;
    public int curr_fid = 0;
    //string prev_folder = "";
    //string curr_folder = "";
    public void Start()
    {

    }
    public void Update()
    {
        if (frame_cnt == 0) return;
    }
    public void OnGUI()
    {
        if (frame_cnt == 0) return;
        GUILayout.BeginHorizontal("Box");
        GUILayout.Box(curr_fid.ToString("D03") + "/" + frame_cnt.ToString(), U3DUtils.WHOpt(70));
        int new_fid = (int)GUILayout.HorizontalSlider(curr_fid, 0, frame_cnt, U3DUtils.WHOpt(300));
        GUILayout.EndHorizontal();
        if (new_fid != curr_fid)
            SetFrameID(new_fid);
    }
    public void SetFrameID(int _fid)
    {
        if (frame_cnt == 0) return;
    }
}
public class DrawLineTool
{
    static Texture2D _aaLineTex;
    static Texture2D _lineTex;
    static void SetIfInTex(Texture2D tex, int x0, int y0, Color _color)
    {
        if (x0 >= 0 && x0 < tex.width && y0 >= 0 && y0 < tex.height)
            tex.SetPixel(x0, y0, _color);

    }
    public static void DrawLineThick(Texture2D tex, int thick, int x0, int y0, int x1, int y1, Color col, bool apply = true, bool wrap = false)
    {
        //DrawLine(tex, x0, y0, x1, y1, col, apply, wrap);
        int ht = thick >> 1;
        for (int i = -ht; i <= ht; i++)
        {
            DrawLine(tex, x0 + i, y0, x1 + i, y1, col, apply, wrap);
            DrawLine(tex, x0, y0 + i, x1, y1 + i, col, apply, wrap);
        }
    }
    public static void DrawLine(Texture2D tex, int x0, int y0, int x1, int y1, Color col, bool apply = true, bool wrap = false)
    {
        int dy = (int)(y1 - y0);
        int dx = (int)(x1 - x0);
        int stepx, stepy;

        if (dy < 0) { dy = -dy; stepy = -1; }
        else { stepy = 1; }
        if (dx < 0) { dx = -dx; stepx = -1; }
        else { stepx = 1; }
        dy <<= 1;
        dx <<= 1;

        float fraction = 0;

        tex.SetPixel(x0, y0, col);
        if (dx > dy)
        {
            fraction = dy - (dx >> 1);
            while (Mathf.Abs(x0 - x1) > 1)
            {
                if (fraction >= 0)
                {
                    y0 += stepy;
                    fraction -= dx;
                }
                x0 += stepx;
                fraction += dy;
                if (wrap)
                    tex.SetPixel(x0, y0, col);
                else
                    SetIfInTex(tex, x0, y0, col);
            }
        }
        else
        {
            fraction = dx - (dy >> 1);
            while (Mathf.Abs(y0 - y1) > 1)
            {
                if (fraction >= 0)
                {
                    x0 += stepx;
                    fraction -= dy;
                }
                y0 += stepy;
                fraction += dx;
                if (wrap)
                    tex.SetPixel(x0, y0, col);
                else
                    SetIfInTex(tex, x0, y0, col);
            }
        }
        if (apply)
            tex.Apply();
    }
    //Usethisforinitialization
    //static bool antiAlias = false;
    static Texture2D adLineTex
    {
        get
        {
            if (!_aaLineTex)
            {
                _aaLineTex = new Texture2D(1, 3, TextureFormat.ARGB32, true);
                _aaLineTex.SetPixel(0, 0, new Color(1, 1, 1, 0));
                _aaLineTex.SetPixel(0, 1, Color.white);
                _aaLineTex.SetPixel(0, 2, new Color(1, 1, 1, 0));
                _aaLineTex.Apply();
            }
            return _aaLineTex;
        }
    }

    static Texture2D lineTex
    {
        get
        {
            if (!_lineTex)
            {
                _lineTex = new Texture2D(1, 1, TextureFormat.ARGB32, true);
                _lineTex.SetPixel(0, 1, Color.white);
                _lineTex.Apply();
            }
            return _lineTex;
        }
    }

    public static void DrawLineMac(Vector2 pointA, Vector2 pointB, Color color, float width, bool antiAlias)
    {
        Color savedColor = GUI.color;
        Matrix4x4 savedMatrix = GUI.matrix;

        float oldWidth = width;

        if (antiAlias) width *= 3;
        float angle = Vector3.Angle(pointB - pointA, Vector2.right) * (pointA.y <= pointB.y ? 1 : -1);
        float m = (pointB - pointA).magnitude;

        if (m > 0.01f)
        {

            Vector3 dz = new Vector3(pointA.x, pointA.y, 0);
            Vector3 offset = new Vector3((pointB.x - pointA.x) * 0.5f, (pointB.y - pointA.y) * 0.5f, 0f);
            Vector3 tmp = Vector3.zero;

            if (antiAlias)
                tmp = new Vector3(-oldWidth * 1.5f * Mathf.Sin(angle * Mathf.Deg2Rad), oldWidth * 1.5f * Mathf.Cos(angle * Mathf.Deg2Rad));
            else
                tmp = new Vector3(-oldWidth * 0.5f * Mathf.Sin(angle * Mathf.Deg2Rad), oldWidth * 0.5f * Mathf.Cos(angle * Mathf.Deg2Rad));
            GUI.color = color;
            GUI.matrix = translationMatrix(dz) * GUI.matrix;
            GUIUtility.ScaleAroundPivot(new Vector2(m, width), new Vector2(-0.5f, 0));
            GUI.matrix = translationMatrix(dz) * GUI.matrix;
            GUIUtility.RotateAroundPivot(angle, Vector2.zero);
            GUI.matrix = translationMatrix(dz - tmp - offset) * GUI.matrix;

            GUI.DrawTexture(new Rect(0, 0, 1, 1), antiAlias ? adLineTex : lineTex);

        }
        GUI.matrix = savedMatrix;
        GUI.color = savedColor;

    }

    private static Matrix4x4 translationMatrix(Vector3 v)
    {
        return Matrix4x4.TRS(v, Quaternion.identity, Vector3.one);
    }
    public static bool SaveTexture2D(Texture2D _tex, string name)
    {
        if (_tex == null) return false;
        byte[] test = _tex.EncodeToJPG();//EncodeToJPG EncodeToPNG
        FileStream fs = new FileStream(name + ".jpg", FileMode.Create);
        fs.Write(test, 0, test.Length);
        fs.Close();
        return true;
    }

}
public class TextureScale
{
    public class ThreadData
    {
        public int start;
        public int end;
        public ThreadData(int s, int e)
        {
            start = s;
            end = e;
        }
    }

    private static Color[] texColors;
    private static Color[] newColors;
    private static int w;
    private static float ratioX;
    private static float ratioY;
    private static int w2;
    private static int finishCount;
    private static Mutex mutex;

    public static Texture2D Point(Texture2D tex, int newWidth, int newHeight)
    {
        return ThreadedScale(tex, newWidth, newHeight, false);
    }

    public static void Bilinear(Texture2D tex, int newWidth, int newHeight)
    {
        ThreadedScale(tex, newWidth, newHeight, true);
    }

    private static Texture2D ThreadedScale(Texture2D tex, int newWidth, int newHeight, bool useBilinear)
    {
        texColors = tex.GetPixels();
        newColors = new Color[newWidth * newHeight];
        if (useBilinear)
        {
            ratioX = 1.0f / ((float)newWidth / (tex.width - 1));
            ratioY = 1.0f / ((float)newHeight / (tex.height - 1));
            ratioX = ratioY;
        }
        else
        {
            ratioX = ((float)tex.width) / newWidth;
            ratioY = ((float)tex.height) / newHeight;
            ratioX = ratioY;
        }
        w = tex.width;
        w2 = newWidth;
        var cores = Mathf.Min(SystemInfo.processorCount, newHeight);
        var slice = newHeight / cores;

        finishCount = 0;
        if (mutex == null)
        {
            mutex = new Mutex(false);
        }
        if (cores > 1)
        {
            int i = 0;
            ThreadData threadData;
            for (i = 0; i < cores - 1; i++)
            {
                threadData = new ThreadData(slice * i, slice * (i + 1));
                ParameterizedThreadStart ts = useBilinear ? new ParameterizedThreadStart(BilinearScale) : new ParameterizedThreadStart(PointScale);
                Thread thread = new Thread(ts);
                thread.Start(threadData);
            }
            threadData = new ThreadData(slice * i, newHeight);
            if (useBilinear)
            {
                BilinearScale(threadData);
            }
            else
            {
                PointScale(threadData);
            }
            while (finishCount < cores)
            {
                Thread.Sleep(1);
            }
        }
        else
        {
            ThreadData threadData = new ThreadData(0, newHeight);
            if (useBilinear)
            {
                BilinearScale(threadData);
            }
            else
            {
                PointScale(threadData);
            }
        }
        //==
        Texture2D new_tex = new Texture2D(newWidth, newHeight, tex.format, false);
        new_tex.SetPixels(newColors);
        new_tex.Apply();
        //==
        //tex.Resize(newWidth, newHeight);
        //tex.SetPixels(newColors);
        //tex.Apply();
        //==

        texColors = null;
        newColors = null;
        return new_tex;
    }

    public static void BilinearScale(System.Object obj)
    {
        ThreadData threadData = (ThreadData)obj;
        for (var y = threadData.start; y < threadData.end; y++)
        {
            int yFloor = (int)Mathf.Floor(y * ratioY);
            var y1 = yFloor * w;
            var y2 = (yFloor + 1) * w;
            var yw = y * w2;

            for (var x = 0; x < w2; x++)
            {//w2//TODO, x < w2,x < 216,
                int xFloor = (int)Mathf.Floor(x * ratioX);
                var xLerp = x * ratioX - xFloor;//newColors[yw + x + 84]=
                newColors[yw + x] = ColorLerpUnclamped(ColorLerpUnclamped(texColors[y1 + xFloor], texColors[y1 + xFloor + 1], xLerp),
                                                       ColorLerpUnclamped(texColors[y2 + xFloor], texColors[y2 + xFloor + 1], xLerp),
                                                       y * ratioY - yFloor);
            }
        }

        mutex.WaitOne();
        finishCount++;
        mutex.ReleaseMutex();
    }

    public static void PointScale(System.Object obj)
    {
        ThreadData threadData = (ThreadData)obj;
        for (var y = threadData.start; y < threadData.end; y++)
        {
            var thisY = (int)(ratioY * y) * w;
            var yw = y * w2;//yw = y * w2
            for (var x = 0; x < w2; x++)
            {//TODO, (var x = 0;x < w2;x++), //for (var x = 0; x < 216; x++)384*1080/1920=216
                newColors[yw + x] = texColors[(int)(thisY + ratioX * x)];// newColors[yw + x + 84]
            }
        }

        mutex.WaitOne();
        finishCount++;
        mutex.ReleaseMutex();
    }

    private static Color ColorLerpUnclamped(Color c1, Color c2, float value)
    {
        return new Color(c1.r + (c2.r - c1.r) * value,
                          c1.g + (c2.g - c1.g) * value,
                          c1.b + (c2.b - c1.b) * value,
                          c1.a + (c2.a - c1.a) * value);
    }
}
//public class EM//convert YB <===> U3D
//{
//    public static Keyframe YBToU(YB.Keyframe kf)
//    {
//        return new Keyframe(kf.time, kf.value, kf.inTangent, kf.outTangent);
//    }
//    public static AnimationCurve YBToU(YB.AnimationCurve ac)
//    {
//        if (ac == null || ac.keys.Length == 0)
//            return null;
//        Keyframe[] kfs = new Keyframe[ac.keys.Length];
//        for (int i = 0; i < kfs.Length; i++)
//        {
//            kfs[i] = YBToU(ac.keys[i]);
//        }
//        return new AnimationCurve(kfs);
//    }
//    public static UnityEngine.Vector3 YBToU(YB.Vector3 v3)
//    {
//        return new UnityEngine.Vector3(v3.x, v3.y, v3.z);
//    }
//    public static YB.Vector3 UToYB(UnityEngine.Vector3 v3)
//    {
//        return new YB.Vector3(v3.x, v3.y, v3.z);
//    }
//    public static YB.Quaternion UToYB(UnityEngine.Quaternion q4)
//    {
//        return new YB.Quaternion(q4.x, q4.y, q4.z, q4.w);
//    }
//    public static UnityEngine.Quaternion YBToU(YB.Quaternion q)
//    {
//        return new Quaternion(q.x, q.y, q.z, q.w);
//    }
//    public static Matrix4x4 YBToU(YB.Matrix4x4 m)
//    {
//        Vector4 c0 = new Vector4(m.m00, m.m10, m.m20, m.m30);
//        Vector4 c1 = new Vector4(m.m01, m.m11, m.m21, m.m31);
//        Vector4 c2 = new Vector4(m.m02, m.m12, m.m22, m.m32);
//        Vector4 c3 = new Vector4(m.m03, m.m13, m.m23, m.m33);
//        return new Matrix4x4(c0, c1, c2, c3);
//    }

//}
public class U3DUtils
{
    public static float X01Y01(float val, float x0, float x1, float y0, float y1)
    {
        if (x0 == x1)
            x0 = x1 - 1e-6f;
        val = Mathf.Clamp(val, Mathf.Min(x0, x1), Mathf.Max(x0, x1));
        return ((y0 * x1 - y1 * x0) + val * (y1 - y0)) / (x1 - x0);
    }
    public static void MakeVert(Vector3 a, ref Vector3 b)
    {
        float vd = Vector3.Dot(a, b);
        float bl = b.magnitude;
        if (a.magnitude < 1e-6f || bl < 1e-6f || vd < 1e-6f)
            return;
        Vector3 c = Vector3.Cross(a, b);
        b = Vector3.Cross(c, a).normalized * bl;
    }
    public static float QuatAngle(Quaternion rot)
    {
        Quaternion r = rot.normalized;
        if (r.w < 0)
            r = new Quaternion(-r.x, -r.y, -r.z, -r.w);
        return Mathf.Acos(Mathf.Clamp(r.w, -1, 1)) * 2 * Mathf.Rad2Deg;

    }
    public static Vector3 QuatAxis(Quaternion rot)
    {
        Quaternion r = rot.normalized;
        if (r.w < 0)
            r = new Quaternion(-r.x, -r.y, -r.z, -r.w);
        return new Vector3(r.x, r.y, r.z).normalized;

    }
    public static Transform FindTfmByName(Transform[] tfms, string part_name)
    {
        for (int j = 0; j < tfms.Length; j++)
            if (tfms[j].name == part_name)
                return tfms[j];
        return null;
    }
    public static object LoadJson(string json_file)
    {
        if (!File.Exists(json_file))
        {
            Debug.Log("not exist: \n\t\t\t" + json_file);
            return null;
        }
        string json_text = File.ReadAllText(json_file);
        if (string.IsNullOrEmpty(json_text))
        {
            Debug.Log("null or empty: " + json_file);
            return null;
        }
        return JsonConvert.DeserializeObject(json_text);
    }
    public static GUILayoutOption[] WHOpt(float w = -1.0f, float h = -1.0f)
    {
        if (w < 0) return new GUILayoutOption[1] { GUILayout.Height(h) };
        if (h < 0) return new GUILayoutOption[1] { GUILayout.Width(w) };
        return new GUILayoutOption[2] { GUILayout.Width(w), GUILayout.Height(h) };
    }
    public static GameObject CreateGO(PrimitiveType ptype, float size = 0.1f, string name = "")
    {
        GameObject _go = GameObject.CreatePrimitive(ptype);
        RemoveCollider(_go);
        _go.transform.localScale = Vector3.one * size;
        _go.transform.name = name;
        return _go;
    }
    public static void RemoveCollider(GameObject _go)
    {
        BoxCollider bcd;
        if (_go.TryGetComponent<BoxCollider>(out bcd))
            GameObject.Destroy(bcd);
        SphereCollider scd;
        if (_go.TryGetComponent<SphereCollider>(out scd))
            GameObject.Destroy(scd);
    }
    public static Material CreateMaterial(string mat_name_in_resource_folder)
    {
        Material _mat = Resources.Load(mat_name_in_resource_folder) as Material;
        if (_mat == null)
            return new Material(Shader.Find("Diffuse"));
        else
            return GameObject.Instantiate(_mat) as Material;
    }
    public static string[] GetFilesInDir(string dir, bool remove_suffix = false)
    {
        DirectoryInfo info = new DirectoryInfo(dir);
        if (!info.Exists)
            return new string[0];
        FileInfo[] files = info.GetFiles();
        string[] ret = new string[files.Length];
        for (int i = 0; i < files.Length; i++)
        {
            ret[i] = files[i].Name;
            if (remove_suffix)
            {
                int dot_pos = ret[i].LastIndexOf('.');
                if (dot_pos != -1)
                    ret[i] = ret[i].Substring(0, dot_pos);
            }
            //Debug.Log("=======" + ret[i]);
        }
        return ret;
    }

    #region read/draw/save texture2D
    public static void SetTextureScale(Material mat, Vector2 val)
    {
        if (mat == null) return;
        string[] names = mat.GetTexturePropertyNames();
        foreach (var n in names)
        {
            mat.SetTextureScale(n, val);
            //Debug.Log("set " + n + " scale to :" + val.ToString());
        }
    }
    public static bool LoadTexture2D(Texture2D _tex, string _url)
    {
        FileInfo fileinfo = new FileInfo(_url);
        if (!File.Exists(fileinfo.FullName) || _tex == null)
            return false;
        byte[] fileData = File.ReadAllBytes(fileinfo.FullName);
        _tex.LoadImage(fileData);
        return true;
    }
    public static bool SaveTexture2D(Texture2D _tex, string name)
    {
        if (_tex == null) return false;
        //byte[] test = _tex.EncodeToPNG();
        byte[] test = _tex.EncodeToPNG();//EncodeToJPG EncodeToPNG
        FileStream fs = new FileStream(name + ".png", FileMode.Create);
        if (fs == null) return false;
        fs.Write(test, 0, test.Length);
        fs.Close();
        return true;
    }
    public static bool SaveTexture2D_JPG(Texture2D _tex, string name, bool with_postfix = true)
    {
        if (_tex == null || string.IsNullOrEmpty(name)) return false;
        if (!with_postfix)
            name += ".jpg";
        //byte[] test = _tex.EncodeToPNG();
        byte[] test = _tex.EncodeToJPG();//EncodeToJPG EncodeToPNG
        FileStream fs = new FileStream(name, FileMode.Create);
        if (fs == null) return false;
        fs.Write(test, 0, test.Length);
        fs.Close();
        return true;
    }
    public static void DrawRectOnTex(Texture2D _t, int x0, int y0, int x1, int y1, Color c, bool imme = true)
    {
        DrawVertLineOnTex(_t, x0, y0, x0, y1, c, imme);
        DrawVertLineOnTex(_t, x1, y0, x1, y1, c, imme);
        DrawHoriLineOnTex(_t, x0, y0, x1, y0, c, imme);
        DrawHoriLineOnTex(_t, x0, y1, x1, y1, c, imme);
    }
    public static void DrawVertLineOnTex(Texture2D _t, int x0, int y0, int x1, int y1, Color c, bool imme = true)
    {
        DrawLineTool.DrawLine(_t, x0, y0, x1, y1, c, imme);
        DrawLineTool.DrawLine(_t, x0 - 1, y0, x1 - 1, y1, c, imme);
        DrawLineTool.DrawLine(_t, x0 + 1, y0, x1 + 1, y1, c, imme);
    }
    public static void DrawHoriLineOnTex(Texture2D _t, int x0, int y0, int x1, int y1, Color c, bool imme = true)
    {
        DrawLineTool.DrawLine(_t, x0, y0, x1, y1, c, imme);
        DrawLineTool.DrawLine(_t, x0, y0 - 1, x1, y1 - 1, c, imme);
        DrawLineTool.DrawLine(_t, x0, y0 + 1, x1, y1 + 1, c, imme);
    }
    public static void DrawCrossOnTex(Texture2D _t, int cx, int cy, Color _c, int len = 30)
    {
        len = len >> 1;
        DrawLineTool.DrawLine(_t, cx - len, cy + 0, cx + len, cy + 0, _c);
        DrawLineTool.DrawLine(_t, cx - len, cy + 1, cx + len, cy + 1, _c);
        DrawLineTool.DrawLine(_t, cx - len, cy - 1, cx + len, cy - 1, _c);

        DrawLineTool.DrawLine(_t, cx + 0, cy - len, cx + 0, cy + len, _c);
        DrawLineTool.DrawLine(_t, cx + 1, cy - len, cx + 1, cy + len, _c);
        DrawLineTool.DrawLine(_t, cx - 1, cy - len, cx - 1, cy + len, _c);

    }
    public static GameObject DrawBallAt(Vector3 pos, string _name, float _size = 0.05f, string _parent = "debug_draw_root")
    {
        GameObject _ball = CreateDebugGO(PrimitiveType.Sphere, _size, _name, _parent);
        _ball.transform.position = pos;
        return _ball;
    }
    public static GameObject DrawBoxAt(Vector3 fr, Vector3 to, string _name, float _size = 0.05f, string _parent = "debug_draw_root")
    {
        Vector3 _dir = to - fr;
        GameObject _cube = CreateDebugGO(PrimitiveType.Cube, _size, _name, _parent);
        _cube.transform.position = 0.5f * (fr + to);
        _cube.transform.localScale = new Vector3(_dir.magnitude, _size, _size);
        _cube.transform.rotation = Quaternion.FromToRotation(Vector3.right, _dir);
        return _cube;
    }
    private static GameObject CreateDebugGO(PrimitiveType _type, float _size, string _name, string _parent = "debug_draw_root")
    {
        GameObject _root = GameObject.Find(_parent);
        if (_root == null)
            _root = new GameObject(_parent);
        GameObject rrr = GameObject.Find("debug_draw_root");
        if (rrr == null)
            rrr = new GameObject("debug_draw_root");
        if (rrr != _root)
            _root.transform.parent = rrr.transform;
        GameObject _ball = GameObject.Find(_name);
        if (_ball == null)
            _ball = CreateGO(_type, _size, _name);
        _ball.transform.parent = _root.transform;
        return _ball;
    }
    //*
    enum DebugColor { red, green, blue, black, white, yellow, orange };
    //[MonoPInvokeCallback(typeof(debugCallback))]
    public static void OnDebugCallback(IntPtr request, int color, int size)
    {
        //Ptr to string
        string debug_string = Marshal.PtrToStringAnsi(request, size);
        //Add Specified Color
        //string fmt_debug_string =
        //    String.Format("{0}{1}{2}{3}{4}",
        //    "<color=",
        //    ((DebugColor)color).ToString(),
        //    ">",
        //    debug_string,
        //    "</color>"
        //    );

        UnityEngine.Debug.Log(debug_string + "\n");
    }//*/
    #endregion

    public static float DistPos2NRay(UnityEngine.Vector3 pos, UnityEngine.Vector3[] _pos, UnityEngine.Vector3[] _d)
    {
        float ret = 0;
        for (int i = 0; i < _pos.Length; i++)
        {
            UnityEngine.Vector3 ab = pos - _pos[i];
            UnityEngine.Vector3 cr = UnityEngine.Vector3.Cross(ab, _d[i].normalized);
            ret += cr.magnitude;
        }
        return ret;
    }
    //public static UnityEngine.Vector3 MinDistPos2NRay(UnityEngine.Vector3[] _pos, UnityEngine.Vector3[] _d_input, float[] _w = null)
    //{
    //    if (_pos.Length != _d_input.Length || _d_input.Length < 2)
    //        UnityEngine.Debug.LogError("cacl min dist pos error, pos count != direction count, or less than 2");
    //    if (_w == null || _w.Length != _pos.Length)
    //    {
    //        _w = new float[_pos.Length];
    //        for (int i = 0; i < _w.Length; i++)
    //            _w[i] = 1.0f / _w.Length;
    //    }
    //    UnityEngine.Vector3[] _d = new UnityEngine.Vector3[_d_input.Length];
    //    for (int i = 0; i < _d.Length; i++)
    //        _d[i] = _d_input[i].normalized;

    //    UnityEngine.Vector3 ret = UnityEngine.Vector3.zero;
    //    YB.MyMat A = new YB.MyMat(3, 3, 0);
    //    YB.Vec B = new YB.Vec(3, 0);
    //    for (int i = 0; i < _pos.Length; i++)
    //    {
    //        UnityEngine.Vector3 p = _pos[i];
    //        UnityEngine.Vector3 d = _d[i];
    //        double wi = _w[i];
    //        double d01 = (double)(d.x * d.y);
    //        double d02 = (double)(d.x * d.z);
    //        double d12 = (double)(d.y * d.z);
    //        A.m[0, 0] += wi * (double)(d.y * d.y + d.z * d.z);
    //        A.m[0, 1] += wi * -d01;
    //        A.m[0, 2] += wi * -d02;
    //        A.m[1, 0] += wi * -d01;
    //        A.m[1, 1] += wi * (double)(d.x * d.x + d.z * d.z);
    //        A.m[1, 2] += wi * -d12;
    //        A.m[2, 0] += wi * -d02;
    //        A.m[2, 1] += wi * -d12;
    //        A.m[2, 2] += wi * (double)(d.x * d.x + d.y * d.y);

    //        B.m[0] += wi * (d.z * (d.z * p.x - d.x * p.z) - d.y * (d.x * p.y - d.y * p.x));
    //        B.m[1] += wi * (d.x * (d.x * p.y - d.y * p.x) - d.z * (d.y * p.z - d.z * p.y));
    //        B.m[2] += wi * (d.y * (d.y * p.z - d.z * p.y) - d.x * (d.z * p.x - d.x * p.z));
    //    }
    //    //UnityEngine.Debug.Log("A:(det="+A.Det3()+")" + A.ToString("F3"));
    //    //UnityEngine.Debug.Log("B:" + B.ToString("F3"));
    //    YB.MyMat Ainv = A.Inv33();
    //    YB.Vec ans = Ainv * B;
    //    //UnityEngine.Debug.Log("A*Ainv:" + (A*Ainv).ToString("F2"));
    //    ret.x = (float)ans.m[0];
    //    ret.y = (float)ans.m[1];
    //    ret.z = (float)ans.m[2];
    //    return ret;
    //}

    #region convert U3d.Matrix <===> U3d.Quaternion+U3d.Vector3
    //public static UnityEngine.Matrix4x4 Mat44FromRT(YB.MyMat R, YB.Vec T)
    //{
    //    UnityEngine.Matrix4x4 ret = new UnityEngine.Matrix4x4();
    //    ret.SetRow(0, new UnityEngine.Vector4((float)R.m[0, 0], (float)R.m[0, 1], (float)R.m[0, 2], (float)T.m[0]));
    //    ret.SetRow(1, new UnityEngine.Vector4((float)R.m[1, 0], (float)R.m[1, 1], (float)R.m[1, 2], (float)T.m[1]));
    //    ret.SetRow(2, new UnityEngine.Vector4((float)R.m[2, 0], (float)R.m[2, 1], (float)R.m[2, 2], (float)T.m[2]));
    //    ret.SetRow(3, new UnityEngine.Vector4(0, 0, 0, 1));
    //    return ret;
    //}
    public static Quaternion ReadCompressedRotation(BinaryReader reader)
    {
        const float FLOAT_PRECISION_MULT = 10000f;
        //if(false)
        //{
        //    float qx = reader.ReadSingle();
        //    float qy = reader.ReadSingle();
        //    float qz = reader.ReadSingle();
        //    float qw = reader.ReadSingle();
        //    return EM.QuatNorm(new Quaternion(qx, qy, qz, qw));
        //}
        // Read the index of the omitted field from the stream.
        byte maxIndex = reader.ReadByte();

        // Values between 4 and 7 indicate that only the index of the single field whose value is 1f was
        // sent, and (maxIndex - 4) is the correct index for that field.
        if (maxIndex >= 4 && maxIndex <= 7)
        {
            var x = (maxIndex == 4) ? 1f : 0f;
            var y = (maxIndex == 5) ? 1f : 0f;
            var z = (maxIndex == 6) ? 1f : 0f;
            var w = (maxIndex == 7) ? 1f : 0f;

            return new Quaternion(x, y, z, w);
        }

        // Read the other three fields and derive the value of the omitted field
        var a = (float)reader.ReadInt16() / FLOAT_PRECISION_MULT;
        var b = (float)reader.ReadInt16() / FLOAT_PRECISION_MULT;
        var c = (float)reader.ReadInt16() / FLOAT_PRECISION_MULT;
        float d = 0f;
        float abc_sqr = (a * a + b * b + c * c);
        if (abc_sqr > 1)
        {
            float abc_reci = 1.0f / Mathf.Sqrt(abc_sqr);
            a *= abc_reci;
            b *= abc_reci;
            c *= abc_reci;
        }
        else
            d = Mathf.Sqrt(1f - abc_sqr);

        if (maxIndex == 0)
            return new Quaternion(d, a, b, c);
        else if (maxIndex == 1)
            return new Quaternion(a, d, b, c);
        else if (maxIndex == 2)
            return new Quaternion(a, b, d, c);

        return new Quaternion(a, b, c, d);
    }
    public static UnityEngine.Quaternion QuatNorm(UnityEngine.Quaternion _q)
    {
        if (_q == null) return UnityEngine.Quaternion.identity;
        float _len = (_q.x * _q.x + _q.y * _q.y + _q.z * _q.z + _q.w * _q.w);
        if (Mathf.Abs(_len) < 1e-6f)
            return UnityEngine.Quaternion.identity;
        else
        {
            _len = 1.0f / UnityEngine.Mathf.Sqrt(_len);
            if (_q.w < 0)
                _len = -_len;
            return new UnityEngine.Quaternion(_q.x * _len, _q.y * _len, _q.z * _len, _q.w * _len);
        }
    }
    public static Matrix4x4 Q2M(Quaternion q)
    {
        Vector4 c0 = new Vector4(1 - 2 * q.y * q.y - 2 * q.z * q.z, 2 * q.x * q.y + 2 * q.z * q.w, 2 * q.x * q.z - 2 * q.y * q.w, 0);
        Vector4 c1 = new Vector4(2 * q.x * q.y - 2 * q.z * q.w, 1 - 2 * q.x * q.x - 2 * q.z * q.z, 2 * q.y * q.z + 2 * q.x * q.w, 0);
        Vector4 c2 = new Vector4(2 * q.x * q.z + 2 * q.y * q.w, 2 * q.y * q.z - 2 * q.x * q.w, 1 - 2 * q.x * q.x - 2 * q.y * q.y, 0);
        Vector4 c3 = new Vector4(0, 0, 0, 1);
        return new Matrix4x4(c0, c1, c2, c3);
    }
    public static Quaternion M2Q(Matrix4x4 m)
    {
        float m00 = m[0, 0];
        float m01 = m[0, 1];
        float m02 = m[0, 2];
        float m10 = m[1, 0];
        float m11 = m[1, 1];
        float m12 = m[1, 2];
        float m20 = m[2, 0];
        float m21 = m[2, 1];
        float m22 = m[2, 2];
        float qx = 0, qy = 0, qz = 0, qw = 1;
        float tr = m00 + m11 + m22;

        if (tr > 0)
        {
            float S = Mathf.Sqrt(tr + 1.0f) * 2; // S=4*qw 
            qw = 0.25f * S;
            qx = (m21 - m12) / S;
            qy = (m02 - m20) / S;
            qz = (m10 - m01) / S;
        }
        else if ((m00 > m11) & (m00 > m22))
        {
            float S = Mathf.Sqrt(1.0f + m00 - m11 - m22) * 2; // S=4*qx 
            qw = (m21 - m12) / S;
            qx = 0.25f * S;
            qy = (m01 + m10) / S;
            qz = (m02 + m20) / S;
        }
        else if (m11 > m22)
        {
            float S = Mathf.Sqrt(1.0f + m11 - m00 - m22) * 2; // S=4*qy
            qw = (m02 - m20) / S;
            qx = (m01 + m10) / S;
            qy = 0.25f * S;
            qz = (m12 + m21) / S;
        }
        else
        {
            float S = Mathf.Sqrt(1.0f + m22 - m00 - m11) * 2; // S=4*qz
            qw = (m10 - m01) / S;
            qx = (m02 + m20) / S;
            qy = (m12 + m21) / S;
            qz = 0.25f * S;
        }
        return new Quaternion(qx, qy, qz, qw);
    }

    public static bool r_z = true;
    public static Matrix4x4 swapYZ = new Matrix4x4(new Vector4(1, 0, 0, 0), new Vector4(0, 0, 1, 0), new Vector4(0, 1, 0, 0), new Vector4(0, 0, 0, 1));
    public static Matrix4x4 invertZM = Matrix4x4.TRS(Vector3.zero, Quaternion.identity, new Vector3(1, 1, -1));
    public static Matrix4x4 invertYM = Matrix4x4.TRS(Vector3.zero, Quaternion.identity, new Vector3(1, -1, 1));
    public static Matrix4x4 FromPQ(Vector3 p, Quaternion q)
    {
        return Matrix4x4.TRS(p, QuatNorm(q), Vector3.one);
    }

    public static Vector3 ExtractTranslationFromMatrix(Matrix4x4 matrix)
    {
        Vector3 translate;
        translate.x = matrix.m03;
        translate.y = matrix.m13;
        translate.z = matrix.m23;
        return translate;
    }
    public static Quaternion ExtractRotationFromMatrix(Matrix4x4 matrix)
    {
        Vector3 forward;
        forward.x = matrix.m02;
        forward.y = matrix.m12;
        forward.z = matrix.m22;

        Vector3 upwards;
        upwards.x = matrix.m01;
        upwards.y = matrix.m11;
        upwards.z = matrix.m21;

        return Quaternion.LookRotation(forward, upwards);
    }
    public static Vector3 ExtractScaleFromMatrix(Matrix4x4 matrix)
    {
        Vector3 scale;
        scale.x = new Vector4(matrix.m00, matrix.m10, matrix.m20, matrix.m30).magnitude;
        scale.y = new Vector4(matrix.m01, matrix.m11, matrix.m21, matrix.m31).magnitude;
        scale.z = new Vector4(matrix.m02, matrix.m12, matrix.m22, matrix.m32).magnitude;
        return scale;
    }

    public static void SetTransformFromMatrix(Transform transform, Matrix4x4 matrix)
    {
        transform.localPosition = ExtractTranslationFromMatrix(matrix);
        transform.localRotation = ExtractRotationFromMatrix(matrix);
        transform.localScale = ExtractScaleFromMatrix(matrix);
    }
    #endregion
}
