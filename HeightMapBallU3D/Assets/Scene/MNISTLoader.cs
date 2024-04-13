using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using System.Text;
public class Pix28x28
{
    public const int iw = 28, ih = 28;
    public const float inv_255 = 1.0f / byte.MaxValue;
    //public float[,] pix = new float[ih, iw];
    public byte[,] raw_byte = new byte[ih, iw];
    public List<byte[,]> sub_raw_byte = new List<byte[,]>();
    int iid = 0;
    public int gt_lbl = 0;
    //Texture2D tex;
    public Pix28x28(int _gt_lbl, Texture2D _tex, List<Texture2D> sub_tex = null)
    {
        iid = -1;
        gt_lbl = _gt_lbl;
        raw_byte = ReadTex2Bytes(_tex);
        if (sub_tex != null)
            foreach (var _st in sub_tex)
                sub_raw_byte.Add(ReadTex2Bytes(_st));
        //if (sub_raw_byte.Count > 0)
            //Debug.Log("gt lbl=" + _gt_lbl + ", sub tex:" + sub_raw_byte.Count);
    }
    static byte[,] ReadTex2Bytes(Texture2D _tex)
    {
        byte[,]  _byte = new byte[ih, iw];
        Color[] tex_pix = _tex.GetPixels();
        int bid = 0;
        for (int row = 0; row < ih; row++)
            for (int col = 0; col < iw; col++)
            {
                float _val = tex_pix[bid].r;
                //pix[row, col] = _val;
                byte _byt = (byte)(Mathf.Clamp(_val * 255, 0, 255));
                _byte[row, col] = _byt;
                bid++;
            }
        return _byte;
    }
    public Pix28x28(int id, byte[] raw_bytes, int _gt_lbl, int offset_byte16 = 16)
    {
        iid = id;
        gt_lbl = _gt_lbl;
        //tex = new Texture2D(iw, ih, TextureFormat.RGB24, false);

        int bid = id * iw * ih + offset_byte16;
        for (int row = 0; row < ih; row++)
            for (int col = 0; col < iw; col++)
            {
                byte val = raw_bytes[bid];
                raw_byte[ih - 1 - row, col] = val;
                //pix[ih - 1 - row, col] = val * inv_255;
                bid++;
            }
        //tex_pix[0] = Color.red;//left bottom
        //tex.SetPixels(tex_pix);
        //tex.Apply();
    }
}
class MNISTLoader
{
    public static string mnist_folder = "C:/gzWork/U3DProjects/CloudMatch/Docs/MNIST/";
    public const string _val_img_file = "t10k-images.idx3-ubyte";
    public const string _val_lbl_file = "t10k-labels.idx1-ubyte";
    public const string _train_img_file = "train-images.idx3-ubyte";
    public const string _train_lbl_file = "train-labels.idx1-ubyte";
    static int FourBytes2Int(byte[] a, int id)
    {
        return FourBytes2Int(a[id], a[id + 1], a[id + 2], a[id + 3]);
    }
    static int FourBytes2Int(byte a, byte b, byte c, byte d)
    {
        int ret = a;
        ret = (ret << 8) | b;
        ret = (ret << 8) | c;
        ret = (ret << 8) | d;
        return ret;
    }
    public static List<Pix28x28> ParseMNISTDataByBytes()//string _img_file, string _lbl_file)
    {
        List<Pix28x28> item_list = new List<Pix28x28>();
        byte[] _img_bytes = File.ReadAllBytes(mnist_folder + "/" + _val_img_file);
        byte[] _lbl_bytes = File.ReadAllBytes(mnist_folder + "/" + _val_lbl_file);

        int magic_num = FourBytes2Int(_img_bytes, 0);
        int img_cnt = FourBytes2Int(_img_bytes, 4);
        int img_wid = FourBytes2Int(_img_bytes, 8);
        int img_hei = FourBytes2Int(_img_bytes, 12);
        int magic_num2 = FourBytes2Int(_lbl_bytes, 0);
        int lbl_cnt = FourBytes2Int(_lbl_bytes, 4);
        Debug.Log("img bytes len:" + _img_bytes.Length
            + ", magic num=" + magic_num
            + ", img cnt=" + img_cnt
            + ", img size=" + img_wid + "x" + img_hei
            + "\nbytes len:" + _lbl_bytes.Length
            + ", magic num=" + magic_num2
            + ", lbl cnt=" + lbl_cnt);
        for (int i = 0; i < img_cnt; i++)
        {
            item_list.Add(new Pix28x28(i, _img_bytes, (int)_lbl_bytes[8 + i]));
        }
        return item_list;
    }

    public static List<Pix28x28>[] LoadStdFonts()
    {
        List<Pix28x28>[] std_fonts = new List<Pix28x28>[10];
        for (int i = 0; i < 10; i++)
            std_fonts[i] = LoadStdFonts(mnist_folder + "/fonts/", i);
        return std_fonts;
    }
    static List<Pix28x28> LoadStdFonts(string _folder, int num)
    {
        List<Pix28x28> ret = new List<Pix28x28>();
        _folder += "/num" + num + "/";
        string[] img_files = U3DUtils.GetFilesInDir(_folder, true);
        for (int i = 0; i < img_files.Length; i++)
        {
            Texture2D _tex = new Texture2D(0, 0);
            if (U3DUtils.LoadTexture2D(_tex, _folder + img_files[i] + ".png"))
            {
                string sub_folder = _folder + "/" + img_files[i] + "/";
                List<Texture2D> sub_tex_list = LoadSubTex(sub_folder);
                ret.Add(new Pix28x28(num, _tex, sub_tex_list));
            }
            else
                Debug.Log("load std fonts faild: " + img_files[i]);
        }
        return ret;
    }
    static List<Texture2D> LoadSubTex(string _sub_folder)
    {
        List<Texture2D> sub_tex_list = new List<Texture2D>();
        string[] sub_img_files = U3DUtils.GetFilesInDir(_sub_folder);
        foreach (var _sub_if in sub_img_files)
        {
            Texture2D _sub_tex = new Texture2D(0, 0);
            if (U3DUtils.LoadTexture2D(_sub_tex, _sub_folder + "/" + _sub_if))
                sub_tex_list.Add(_sub_tex);
        }
        return sub_tex_list;
    }
}
