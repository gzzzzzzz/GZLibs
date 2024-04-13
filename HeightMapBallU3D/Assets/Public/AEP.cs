using System.Collections;
using System.Collections.Generic;
using System.IO;
public static class BIP
{
    public enum EP
    {
        root = 0, pelvis = 1, spine = 2,
        spine1 = 3, spine2 = 4, neck = 5, head = 6,
        l_clavicle = 7, l_upper_arm = 8, l_fore_arm = 9, l_hand = 10,
        r_clavicle = 11, r_upper_arm = 12, r_fore_arm = 13, r_hand = 14,
        l_thigh = 15, l_calf = 16, l_foot = 17, l_toe0 = 18,
        r_thigh = 19, r_calf = 20, r_foot = 21, r_toe0 = 22,
        //e_spine3 = 23, e_spine4 = 24, e_spine5 = 25, e_spine6 = 26,
        l_heel = 23, r_heel = 24,
        bip_num = 25
    }
    public static int[] PID_normal = new int[(int)EP.bip_num]//for rotate structure
        {
            -1,                 //0 root->null
            (int)EP.root,       //1 pelvis-->root
            (int)EP.pelvis,     //2 e_spine-->pelvis
            (int)EP.spine,        //3 e_spine1--->e_spine
            (int)EP.spine1,        //4 e_spine2--->e_spine1
            (int)EP.spine2,        //5 e_neck--->e_spine2
            (int)EP.neck,       //6 e_head--->e_neck

            (int)EP.spine2,       //7 l_cla------------------------->spine2
            (int)EP.l_clavicle, //8 l_upper--->l_cla
            (int)EP.l_upper_arm,//9 l_fore--->l_upper
            (int)EP.l_fore_arm, //10 l_hand--->l_fore

            (int)EP.spine2,         //11 r_cla-------------------->spine2
            (int)EP.r_clavicle,   //12 r_upper--->r_cla
            (int)EP.r_upper_arm,  //13 r_fore--->r_upper
            (int)EP.r_fore_arm,   //14 r_hand--->r_fore

            (int)EP.pelvis,       //15 l_thigh------------------>pelvis
            (int)EP.l_thigh,      //16 l_calf--->l_thigh
            (int)EP.l_calf,       //17 l_foot--->l_calf
            (int)EP.l_foot,       //18 l_toe0--->l_foot
            
            (int)EP.pelvis,       //19 r_thigh------------------>pelvis
            (int)EP.r_thigh,      //20 r_calf--->r_thigh
            (int)EP.r_calf,       //21 r_foot--->r_calf
            (int)EP.r_foot,       //22 r_toe0--->r_foot

            (int)EP.l_foot,       //23 l_heel--->l_foot
            (int)EP.r_foot,       //23 r_heel--->r_foot
    };
    public static int[] PID = new int[(int)EP.bip_num]//3d-max bip structure
        {
            -1,                 //0 root->null
            (int)EP.root,       //1 pelvis-->root
            (int)EP.pelvis,     //2 e_spine-->pelvis
            (int)EP.spine,        //3 e_spine1--->e_spine
            (int)EP.spine1,        //4 e_spine2--->e_spine1
            (int)EP.spine2,        //5 e_neck--->e_spine2
            (int)EP.neck,       //6 e_head--->e_neck

            (int)EP.spine2,       //7 l_cla--->e_neck
            (int)EP.l_clavicle, //8 l_upper--->l_cla
            (int)EP.l_upper_arm,//9 l_fore--->l_upper
            (int)EP.l_fore_arm, //10 l_hand--->l_fore

            (int)EP.spine2,         //11 r_cla--->e_neck
            (int)EP.r_clavicle,   //12 r_upper--->r_cla
            (int)EP.r_upper_arm,  //13 r_fore--->r_upper
            (int)EP.r_fore_arm,   //14 r_hand--->r_fore

            (int)EP.pelvis,          //15 l_thigh--->e_spine
            (int)EP.l_thigh,      //16 l_calf--->l_thigh
            (int)EP.l_calf,       //17 l_foot--->l_calf
            (int)EP.l_foot,       //18 l_toe0--->l_foot
            
            (int)EP.pelvis,          //19 r_thigh--->e_spine
            (int)EP.r_thigh,      //20 r_calf--->r_thigh
            (int)EP.r_calf,       //21 r_foot--->r_calf
            (int)EP.r_foot,       //22 r_toe0--->r_foot

            (int)EP.l_foot,       //23 l_heel--->l_foot
            (int)EP.r_foot,       //23 r_heel--->r_foot
    };
    public static int[] B2A = new int[(int)EP.bip_num] {-1, (int)AIC.EP.sim_bip_pelvis, (int)AIC.EP.sim_bip_spine, 
        (int)AIC.EP.sp1, (int)AIC.EP.sp2, (int)AIC.EP.sim_bip_neck, (int)AIC.EP.sim_bip_head,
        (int)AIC.EP.sim_bip_cla_l, (int)AIC.EP.sho_l, (int)AIC.EP.elb_l, (int)AIC.EP.wri_l,
        (int)AIC.EP.sim_bip_cla_r, (int)AIC.EP.sho_r, (int)AIC.EP.elb_r, (int)AIC.EP.wri_r,
        (int)AIC.EP.hip_l, (int)AIC.EP.kne_l, (int)AIC.EP.ank_l, (int)AIC.EP.toe_l,
        (int)AIC.EP.hip_r, (int)AIC.EP.kne_r, (int)AIC.EP.ank_r, (int)AIC.EP.toe_r,

        (int)AIC.EP.heel_l, (int)AIC.EP.heel_r
    };

    public static string[] desc = new string[(int)EP.bip_num] {//For MaxScript Name
        "Biped Root", "Pelvis", "Spine",
        "Spine1", "Spine2", "Neck", "Head",                                 //[0,6]
        "Left Clavicle", "Left UpperArm", "Left Forearm", "Left Hand",      //[7,10]
        "Right Clavicle", "Right UpperArm", "Right Forearm", "Right Hand",  //[11,14]
        "Left Thigh", "Left Calf", "Left Foot", "Left Toe0",                //[15,18]
        "Right Thigh", "Right Calf", "Right Foot", "Right Toe0"             //[19,22]
        ,"Left Heel", "Right Heel"
    };
    public static float[] COM_w = new float[(int)EP.bip_num]
    {
        0, 0.07f, 0.05f, 0.05f, 0.05f,//root,plvs, sp012        22
        0.04f, 0.08f,                //neck, head               12
        0.04f, 0.04f, 0.03f, 0.02f,//clv,sho,elb,wri            
        0.04f, 0.04f, 0.03f, 0.02f,//                           26
        0.08f, 0.06f, 0.04f, 0.01f,//thigh,calf,foot,toe        
        0.08f, 0.06f, 0.04f, 0.01f,//thigh,calf,foot,toe        38
        0.01f, 0.01f,//heelx2
    };
    public static string E2S(EP ep) { return ep.ToString(); }
    public static string I2S(int i) { return ((EP)i).ToString(); }
    public const float spine_r = 0.4f;//plvs_sp / plvs_sp1
}
public static class AIC
{
    public const int part_14 = 14;
    public const int part_cnt = 18;
    public const int j3d_part_cnt = 18;
    public const int connects_cnt = 16;
    public const float ui_toe_length_bias = 0.0f;
    public const float ui_toe_vert_bias = 0.0f;
    public const float ui_toe_hori_bias = 0.0f;
    public const float ui_shoulder_width = 0.52f;
    public const float ui_clv_hori_bias = 0.05f;
    public const float ui_clv_vert_bias = 0.05f;
    public const float ui_neck_bias = 0.35f;
    public const float ui_head_bias = 0.12f;
    public enum EP
    {
        sho_r = 0, elb_r = 1, wri_r = 2,
        sho_l = 3, elb_l = 4, wri_l = 5,
        hip_r = 6, kne_r = 7, ank_r = 8,
        hip_l = 9, kne_l = 10, ank_l = 11,
        head = 12, neck = 13,
        //uv_cnt=14,
        sp2 = 14, sp1 = 15,
        toe_r = 16, toe_l = 17,
        //bp_cnt=18
        heel_l = 18, heel_r = 19, sim_bip_pelvis = 20,
        //on_ground_cnt=21
        sim_bip_neck = 21, sim_bip_head = 22,
        sim_bip_cla_l = 23, sim_bip_cla_r = 24, sim_bip_spine = 25,
        //sim_bip_cnt=26
        uv_cnt = 14,
        bp_cnt = 18,
        on_ground_cnt = 21,
        sim_bip_cnt = 26,
        total_cnt = 26,
    }

    public enum EC
    {
        neck_sho_r = 0, sho_elb_r = 1, elb_wri_r = 2,
        neck_sho_l = 3, sho_elb_l = 4, elb_wri_l = 5,

        pelvis_hip_r = 6, hip_kne_r = 7, kne_ank_r = 8,
        pelvis_hip_l = 9, hip_kne_l = 10, kne_ank_l = 11,

        head_neck = 12, neck_sp2 = 13, sp2_sp1 = 14, sp1_pelvis = 15,
        basic_cnt = 16,

        ank_toe_r = 16, ank_heel_r = 17, heel_toe_r = 18,
        ank_toe_l = 19, ank_heel_l = 20, heel_toe_l = 21,
        on_ground_cnt = 22,

        sim_head_neck = 22, sim_neck_sp2 = 23,
        sim_neck_cla_r = 24, sim_cla_sho_r = 25,
        sim_neck_cla_l = 26, sim_cla_sho_l = 27,

        sim_bip_cnt = 28,
        total_cnt = 28,
    }
    public static int[,] connects = new int[(int)EC.total_cnt, 2]
#region connects
    {
    { (int)EP.neck, (int)EP.sho_r }, //e_neck_sho_r = 0
    { (int)EP.sho_r, (int)EP.elb_r }, //e_sho_elb_r = 1
    { (int)EP.elb_r, (int)EP.wri_r }, //e_elb_wri_r = 2

    { (int)EP.neck, (int)EP.sho_l }, //e_neck_sho_l = 3
    { (int)EP.sho_l, (int)EP.elb_l }, //e_sho_elb_l = 4
    { (int)EP.elb_l, (int)EP.wri_l }, //e_elb_wri_l = 5

    { (int)EP.sim_bip_pelvis, (int)EP.hip_r }, //e_pelvis_hip_r = 6
    { (int)EP.hip_r, (int)EP.kne_r },       //e_hip_kne_r = 7
    { (int)EP.kne_r, (int)EP.ank_r },       //e_kne_ank_r = 8
    
    { (int)EP.sim_bip_pelvis, (int)EP.hip_l }, //e_pelvis_hip_l = 9
    { (int)EP.hip_l, (int)EP.kne_l },       //e_hip_kne_l = 10
    { (int)EP.kne_l, (int)EP.ank_l },       //e_kne_ank_l = 11

    { (int)EP.head, (int)EP.neck },         //e_head_neck = 12
    { (int)EP.neck, (int)EP.sp2 },          //e_neck_sp2 = 13
    { (int)EP.sp2, (int)EP.sp1 },           //e_sp2_sp1 = 14
    { (int)EP.sp1, (int)EP.sim_bip_pelvis },   //e_sp1_pelvis = 15
    //basic connects 16
    { (int)EP.ank_r, (int)EP.toe_r },       //e_ank_toe_r = 16
    { (int)EP.ank_r, (int)EP.heel_r },       //e_ank_heel_r = 17
    { (int)EP.heel_r, (int)EP.toe_r },       //e_heel_toe_r = 18

    { (int)EP.ank_l, (int)EP.toe_l },       //e_ank_toe_r = 19
    { (int)EP.ank_l, (int)EP.heel_l },       //e_ank_heel_r = 20
    { (int)EP.heel_l, (int)EP.toe_l },       //e_heel_toe_r = 21
    //adj_cnt 22
    { (int)EP.sim_bip_head, (int)EP.sim_bip_neck }, //sim_head_neck = 22
    { (int)EP.sim_bip_neck, (int)EP.sp2 },       //sim_neck_sp2 = 23
    { (int)EP.sim_bip_neck, (int)EP.sim_bip_cla_r },//sim_neck_cla_r = 24,
    { (int)EP.sim_bip_cla_r, (int)EP.sho_r},      //sim_cla_sho_r = 25,
    { (int)EP.sim_bip_neck, (int)EP.sim_bip_cla_l },//sim_neck_cla_l = 26,
    { (int)EP.sim_bip_cla_l, (int)EP.sho_l},      //sim_cla_sho_l = 27,
    };
#endregion
    public static int E2I(EP ep) { return (int)ep; }
    public static EP I2E(int i) { return (EP)i; }
    public static string E2S(EP ep) { return ep.ToString(); }
    public static string I2S(int i) { return I2E(i).ToString(); }

    //public static Color[] part_color = new Color[(int)EP.total_cnt] {
    //        new Color(1,0.6f,0), new Color(1,1,0), new Color(0.6f,1,0), 
    //        new Color(0.4f,0,1f), new Color(0,0,0.6f), new Color(0.4f,0,0.6f), 
    //        new Color(1,0.6f,0.2f), new Color(1,1,0.2f), new Color(0.6f,1,0.2f), 
    //        new Color(0.4f,0.2f,1f), new Color(0,0.2f,0.6f), new Color(0.4f,0.2f,0.6f), 
    //        new Color(1,0,1), new Color(0.7f,0,0.7f),               //head neck
    //        new Color(0,0,0), new Color(1,1,1),                     //sp2, sp1
    //        new Color(0.5f,0.5f,0.5f), new Color(0.8f,0.8f,0.8f),   //toe r l
    //        new Color(0.5f,0.5f,0.5f), new Color(0.8f,0.8f,0.8f),   //heel r l

    //        new Color(0,1,1),new Color(0,1,1),new Color(0,1,1),
    //        new Color(0,1,1),new Color(0,1,1),new Color(0,1,1)
    //};
    #region aic2bip
    public static readonly float[] SegRatio = new float[6]{
        0.4f,   //pelvis_sp1 / spine_len
        0.2f,   //sp1_sp2 / spine_len
        0.4f,   //sp2_neck / spine_len

        0.3f,   //BIP_neck_head_bottom / AIC_neck_head_top
        0.08f,  //neck_offset_len / AIC_spine_len
        0.2f    //clavicle_offset_len / AIC_neck_shoulder
    };
    public enum SegRatioEnum
    {
        plvs_sp1 = 0,       //0.4   	//pelvis_sp1 / spine_len
        sp1_sp2 = 1,        //0.2	    //sp1_sp2 / spine_len
        sp2_neck = 2,       //0.4	    //sp2_neck / spine_len
        neck_head = 3,      //0.3	    //BIP_neck_head_bottom / AIC_neck_head_top
        neck_offset = 4,    //0.08	//neck_offset_len / AIC_spine_len
        clavicle_offset = 5 //0.2	    //clavicle_offset_len / AIC_neck_shoulder
    }
    public static float GetBipSpineLen(float aic_spine_len)
    {
        return aic_spine_len * (1 - SegRatio[(int)SegRatioEnum.neck_offset]);
    }
    public static float GetAICSpineLen(float bip_spine_len)
    {
        return bip_spine_len / (1 - SegRatio[(int)SegRatioEnum.neck_offset]);
    }
    public static float GetBipNeckHeadLen(float aic_neck_head_len)
    {
        return aic_neck_head_len * SegRatio[(int)SegRatioEnum.neck_head];
    }
    public static float GetAICNeckHeadLen(float bip_neck_head_len)
    {
        return bip_neck_head_len / SegRatio[(int)SegRatioEnum.neck_head];
    }
    #endregion
}

public class AEP
{
    public int eid = 0;
    public int pid = -1;
    public string desc = "";
    public AEP(int _eid, string _desc, int _pid = -1)
    {
        eid = _eid;
        desc = _desc;
        pid = _pid;
    }
    public string ToStr()
    {
        return desc;
    }

    public static implicit operator int(AEP ep)
    {
        return ep.eid;
    }
    public static implicit operator string(AEP ep)
    {
        return ep.desc;
    }

    #region get parts PID[]
    public static int FindIndexInPartsToFitEid(int eid, AEP[] parts)
    {
        for (int i = 0; i < parts.Length; i++)
            if (parts[i].eid == eid)
                return i;
        return -1;
    }
    public static int FindParentInParts(int id_in_parts, AEP[] parts)
    {
        int ret = -1;
        int pid_in_all = parts[id_in_parts].pid;
        while (pid_in_all != -1)
        {
            ret = FindIndexInPartsToFitEid(pid_in_all, parts);
            if (ret != -1)
                break;
            pid_in_all = AEP.all_ep[pid_in_all].pid;
        }
        return ret;
    }
    public static int[] FindParentArr(AEP[] parts)
    {
        int[] to = new int[parts.Length];
        for (int i = 0; i < parts.Length; i++)
            to[i] = FindParentInParts(i, parts);
        return to;
    }
    public static int[] GetJ3dPID()
    {
        if (j3d_pid != null)
            return j3d_pid;
        j3d_pid = FindParentArr(j3d_order);
        return j3d_pid;
    }
    public static int[] GetJ2dPID()
    {
        if (j2d_pid != null)
            return j2d_pid;
        j2d_pid = FindParentArr(j2d_order);
        return j2d_pid;
    }
    public static int[] GetBipPID()
    {
        if (bip_pid != null)
            return bip_pid;
        bip_pid = FindParentArr(bip_order);
        return bip_pid;
    }
    public static int[] GetUHBPID()
    {
        if (uhb_pid != null)
            return uhb_pid;
        uhb_pid = FindParentArr(uhb_order);
        return uhb_pid;
    }
    public static int[] GetWBPID()
    {
        if (wb_pid != null)
            return wb_pid;
        wb_pid = FindParentArr(wb_order);
        return wb_pid;
    }
    #endregion

    #region all_ep
    public static AEP sho_r = new AEP(0, "sho_r", 18);
    public static AEP elb_r = new AEP(1, "elb_r", 0);
    public static AEP wri_r = new AEP(2, "wri_r", 1);

    public static AEP sho_l = new AEP(3, "sho_l", 19);
    public static AEP elb_l = new AEP(4, "elb_l", 3);
    public static AEP wri_l = new AEP(5, "wri_l", 4);

    public static AEP hip_r = new AEP(6, "hip_r", 14);
    public static AEP kne_r = new AEP(7, "kne_r", 6);
    public static AEP ank_r = new AEP(8, "ank_r", 7);

    public static AEP hip_l = new AEP(9, "hip_l", 14);
    public static AEP kne_l = new AEP(10, "kne_l", 9);
    public static AEP ank_l = new AEP(11, "ank_l", 10);

    public static AEP head = new AEP(12, "head", 13);
    public static AEP neck = new AEP(13, "neck", 17);
    public static AEP plvs = new AEP(14, "plvs", 27);
    public static AEP sp0 = new AEP(15, "sp0", 14);
    public static AEP sp1 = new AEP(16, "sp1", 15);
    public static AEP sp2 = new AEP(17, "sp2", 16);

    public static AEP clv_r = new AEP(18, "clv_r", 17);
    public static AEP clv_l = new AEP(19, "clv_l", 17);

    public static AEP toe_r = new AEP(20, "toe_r", 8);
    public static AEP toe_l = new AEP(21, "toe_l", 11);

    public static AEP heel_r = new AEP(22, "heel_r", 8);
    public static AEP heel_l = new AEP(23, "heel_l", 11);

    public static AEP palm_r = new AEP(24, "palm_r", 2);
    public static AEP palm_l = new AEP(25, "palm_l", 5);

    public static AEP nose = new AEP(26, "nose", 13);
    public static AEP root = new AEP(27, "root", -1);
    #endregion

    #region all ep arr
    public const int all_cnt = 28;
    public static AEP[] all_ep = new AEP[all_cnt] {
        sho_r, elb_r, wri_r,
        sho_l, elb_l, wri_l,
        hip_r, kne_r, ank_r,
        hip_l, kne_l, ank_l,

        head, neck,
        plvs, sp0, sp1, sp2,
        clv_r, clv_l,

        toe_r, toe_l,
        heel_r, heel_l,
        palm_r, palm_l,
        nose, root
        //new YBEP(28, "all")
    };
    #endregion
    #region gd order
    public const int gd_cnt = 4;
    public static AEP[] gd_order = new AEP[gd_cnt]{
        toe_r, heel_r, 
        toe_l, heel_l, 
    };
    public static int[] to_gd_order = new int[AEP.all_cnt]{
        -1,-1,-1,
        -1,-1,-1,
        -1,-1,-1,
        -1,-1,-1,
        -1,-1,
        -1,-1,  //plvs, sp0
        -1,-1,  //sp1,sp2
        -1,-1,  //clv r/l
        0,2,  //toe r/l
        1,3,  //heel r/l
        -1,-1,  //palm r/l
        -1,-1   //nose, root
    };
    #endregion
    #region j2d order
    public const int j2d_cnt = 20;
    public static int[] j2d_pid = null;
    public static AEP[] j2d_order = new AEP[j2d_cnt]{
        sho_r, elb_r, wri_r,
        sho_l, elb_l, wri_l,
        hip_r, kne_r, ank_r,
        hip_l, kne_l, ank_l,
        head, neck, nose,
        toe_r, toe_l,
        palm_r, palm_l, 
        sp0
    };
    public static int[] to_j2d_order = new int[all_cnt]
    {
        0,1,2,
        3,4,5,
        6,7,8,
        9,10,11,
        12,13,      //head neck
        -1,19,-1,-1,//plvs, sp0, sp1, sp2
        -1,-1,      //clv_r, clv_l,
        15,16,      //toe_r, toe_l
        -1,-1,      //heel_r, heel_l
        17,18,      //palm_r, palm_l
        14,-1       //nose, root
    };
    #endregion

    #region j3d order
    public const int j3d_cnt = 23;
    private static int[] j3d_pid = null;
    public static AEP[] j3d_order = new AEP[j3d_cnt]{//jid ===> ep
        sho_r, elb_r, wri_r, //0,1,2,
        sho_l, elb_l, wri_l, //3,4,5
        hip_r, kne_r, ank_r, //6,7,8
        hip_l, kne_l, ank_l, //9,10,11
        head, neck,         //12,13
        sp2, sp1,           //14,15
        toe_r, toe_l,       //16,17
        palm_r, palm_l,     //18,19
        nose,               //20
        heel_r, heel_l      //21,22
    };
    public static string[] j3d_desc = new string[j3d_cnt]
    {
        sho_r.desc, elb_r.desc, wri_r.desc, //0,1,2,
        sho_l.desc, elb_l.desc, wri_l.desc, //3,4,5
        hip_r.desc, kne_r.desc, ank_r.desc, //6,7,8
        hip_l.desc, kne_l.desc, ank_l.desc, //9,10,11
        head.desc, neck.desc,         //12,13
        sp2.desc, sp1.desc,           //14,15
        toe_r.desc, toe_l.desc,       //16,17
        palm_r.desc, palm_l.desc,     //18,19
        nose.desc,               //20
        heel_r.desc, heel_l.desc      //21,22
    };
    public static int[] to_j3d_order = new int[all_cnt] {//ep ===> jid
        0,1,2,
        3,4,5,
        6,7,8,
        9,10,11,
        12,13,
        -1,-1,  //plvs, sp0
        15,14,  //sp1,sp2
        -1,-1,  //clv r/l
        16,17,  //toe r/l
        21,22,  //heel r/l
        18,19,  //palm r/l
        20,-1   //nose, root
    };
    #endregion

    #region bip order
    public const int bip_cnt = 25;
    public static int[] bip_pid = null;
    public static AEP[] bip_order = new AEP[bip_cnt]{
        root, plvs, sp0,            //root = 0, pelvis = 1, spine = 2,
        sp1, sp2, neck, head,       //spine1 = 3, spine2 = 4, neck = 5, head = 6,
        clv_l, sho_l, elb_l, wri_l, //l_clavicle = 7, l_upper_arm = 8, l_fore_arm = 9, l_hand = 10,
        clv_r, sho_r, elb_r, wri_r, //r_clavicle = 11, r_upper_arm = 12, r_fore_arm = 13, r_hand = 14,
        hip_l, kne_l, ank_l, toe_l, //l_thigh = 15, l_calf = 16, l_foot = 17, l_toe0 = 18,
        hip_r, kne_r, ank_r, toe_r, //r_thigh = 19, r_calf = 20, r_foot = 21, r_toe0 = 22,
        heel_r, heel_l
    };
    #endregion

    #region up half body order
    public const int uhb_cnt = 15;
    public static int[] uhb_pid = null;
    public static AEP[] uhb_order = new AEP[uhb_cnt] {
        sho_r, elb_r, wri_r,    //0,1,2
        sho_l, elb_l, wri_l,    //3,4,5
        head, neck, sp2, nose,  //6,7,8,9
        hip_r, hip_l, plvs      //10,11,12
        , clv_r, clv_l
    };
    public static int[] to_uhb_order = new int[all_cnt]
    {
        0,1,2,
        3,4,5,
        10,-1,-1,//hip r
        11,-1,-1,//hip l
        6,7,    //head, neck
        12,-1,  //plvs, sp0
        -1,8,  //sp1,sp2
        13,14,  //clv r/l
        -1,-1,  //toe r/l
        -1,-1,  //heel r/l
        -1,-1,  //palm r/l
        9,-1   //nose, root
    };
    #endregion

    #region full body order
    public const int wb_cnt = 21;
    public static int[] wb_pid = null;
    public static AEP[] wb_order = new AEP[wb_cnt] {
        sho_r, elb_r, wri_r,    //0,1,2
        sho_l, elb_l, wri_l,    //3,4,5
        hip_r, kne_r, ank_r,      
        hip_l, kne_l, ank_l,     
        head, neck, nose,  
        toe_r, toe_l,
        sp2, plvs,
        clv_r, clv_l
    };
    #endregion
    /*
    public static Dictionary<string, YBEP> map_str2ep = new Dictionary<string, YBEP>{
         {all_ep[0].desc, all_ep[0]}, {all_ep[1].desc, all_ep[1]}, {all_ep[2].desc, all_ep[2]},
         {all_ep[3].desc, all_ep[3]}, {all_ep[4].desc, all_ep[4]}, {all_ep[5].desc, all_ep[5]},
         {all_ep[6].desc, all_ep[6]}, {all_ep[7].desc, all_ep[7]}, {all_ep[8].desc, all_ep[8]},
         {all_ep[9].desc, all_ep[9]}, {all_ep[10].desc, all_ep[10]}, {all_ep[11].desc, all_ep[11]},

         {all_ep[12].desc, all_ep[12]}, {all_ep[13].desc, all_ep[13]}, 
         {all_ep[14].desc, all_ep[14]}, {all_ep[15].desc, all_ep[15]}, 
         {all_ep[16].desc, all_ep[16]}, {all_ep[17].desc, all_ep[17]}, 
         {all_ep[18].desc, all_ep[18]}, {all_ep[19].desc, all_ep[19]}, 

         {all_ep[20].desc, all_ep[20]}, {all_ep[21].desc, all_ep[21]}, 
         {all_ep[22].desc, all_ep[22]}, {all_ep[23].desc, all_ep[23]}, 
         {all_ep[24].desc, all_ep[24]}, {all_ep[25].desc, all_ep[25]}, 
         {all_ep[26].desc, all_ep[26]}, {all_ep[27].desc, all_ep[27]}, 
    };
    //*/
}
