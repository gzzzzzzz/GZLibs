Shader "AG/Character Metallic R"
{
	Properties
	{
		[KeywordEnum(Opaque, Cutout, Transparent)] _Mode("Rendering Mode", Float) = 0
		[Space(15)]

		_MainTex("Albedo (RGBA)", 2D) = "white" {}
		_Color("Color (RGBA)", Color) = (1, 1, 1, 1)
		_Cutoff("Alpha Cutoff", Range(0.0, 1.0)) = 0
		[Space(15)]

		[Toggle(_SLIDE_BAR_METALLIC)] _EnableSlideBarMetallic("Use Slider Bar Value /Map R Channel As Metallic", Float) = 0
		_Metallic("Metallic", Range(0.0, 1.0)) = 0.5
		[Toggle(_SLIDE_BAR_SMOOTHNESS)] _EnableSlideBarSmoothness("Use Slider Bar Value /Map A Channel As Smoothness", Float) = 0
		_Glossiness("Smoothness", Range(0.0, 1.0)) = 0.5
		[Toggle(_MAPG_SPECULAR_INTENSITY)] _EnableSpecularIntensityMap("Use Specular Intensity From Map G Channel", Float) = 0
		[Toggle(_REFLECTION_ENVIRONMENT)] _EnableReflection("Enable /Disable Reflection Enviroment", Float) = 1
		[Toggle(_SLIDE_BAR_REFLECTIVITY)] _EnableSliderBarReflectivity("Use Slider Bar Value /Map B Channel As Reflectivity", Float) = 0
		_Reflectivity("Reflectivity", Range(0.0, 1.0)) = 0.5
		[NoScaleOffset] _MetallicGlossMap("Metallic, Specular Intensity, Reflectivity And Smoothness (RGBA)", 2D) = "grey" {}
		_SpecularScale("Specular Scale (Usually 1.0 For Non-metallic Object)", Float) = 1.0
		[Space(15)]

		[Toggle(_MAPB_FRESNEL_INTENSITY)] _EnableFresnelIntensityMap("Use Fresnel Intensity From Map", Float) = 0
		_FresnelIntensityMap("Fresnel Intensity (Only Alpha)", 2D) = "white" {}
		_FresnelColor("Fresnel Color", Color) = (1, 1, 1)
		[Space(15)]

		[Toggle(_NORMALMAP)] _EnableNormalMap("Enable Normal Map", Float) = 0
		[NoScaleOffset] _BumpMapTex("Normal Map", 2D) = "bump" {}
		[Space(15)]

		[Toggle(_EMISSION)] _EnableEmission("Enable Emission", Float) = 0
		_EmissionColor("Color", Color) = (0, 0, 0)
		_EmissionMap("Emission", 2D) = "white" {}
		[Space(15)]

		[Toggle(_DETAIL)] _EnableDetail("Enable Detail", Float) = 0
		[NoScaleOffset] _DetailMask("Detail Mask", 2D) = "white" {}
		_DetailAlbedoMap("Detail Albedo x2", 2D) = "grey" {}
		[NoScaleOffset] _DetailNormalMapTex("Normal Map", 2D) = "bump" {}
		[Enum(UV0, 0, UV1, 1)] _UVSec("UV Set for secondary textures", Float) = 0

		// Blending state
		[HideInInspector] _SrcBlend("__src", Float) = 1.0
		[HideInInspector] _DstBlend("__dst", Float) = 0.0
		[HideInInspector] _ZWrite("__zw", Float) = 1.0
		[HideInInspector] _ColorMask("__cm", Float) = 14.0
	}

	CGINCLUDE
		#undef UNITY_SPECCUBE_BOX_PROJECTION
		#define UNITY_SPECCUBE_BOX_PROJECTION 0
		#undef UNITY_OPTIMIZE_TEXCUBELOD
		#define UNITY_OPTIMIZE_TEXCUBELOD 0
		#undef UNITY_SPECCUBE_BLENDING
		#define UNITY_SPECCUBE_BLENDING 0

		#define _DYANMIC_OBJECT 1
	ENDCG

	SubShader
	{
		Tags{ "RenderType" = "Opaque" "PerformanceChecks" = "False" }
		LOD 600

		Pass
		{
			Name "FORWARD"
			Tags{ "LightMode" = "ForwardBase" }

			Blend[_SrcBlend][_DstBlend]
			ZWrite[_ZWrite]
			ColorMask[_ColorMask]

			CGPROGRAM
			#pragma target 3.0
			#pragma exclude_renderers gles

			#pragma shader_feature ___ _ALPHATEST_ON _ALPHAPREMULTIPLY_ON
			#pragma shader_feature ___ _SLIDE_BAR_SMOOTHNESS
			#pragma shader_feature ___ _REFLECTION_ENVIRONMENT
			#pragma shader_feature ___ _SLIDE_BAR_METALLIC
			#pragma shader_feature ___ _MAPG_SPECULAR_INTENSITY
			#pragma shader_feature ___ _MAPB_FRESNEL_INTENSITY
			#pragma shader_feature ___ _SLIDE_BAR_REFLECTIVITY
			#pragma shader_feature ___ _NORMALMAP
			#pragma shader_feature ___ _EMISSION
			#pragma shader_feature ___ _DETAIL
			#pragma shader_feature ___ _USE_BRDF1 _USE_BRDF2 _USE_BRDF3 _USE_BRDF4

			#pragma multi_compile_fwdbase
			#pragma multi_compile_fog

			#pragma vertex VertDynamicObjectPBS
			#pragma fragment FragCharacterPBSMetallicR

			#include "LK PBS.cginc"

			ENDCG
		}

		Pass
		{
			Name "FORWARD_DELTA"
			Tags{ "LightMode" = "ForwardAdd" }

			Blend One One
			ZWrite Off
			ZTest LEqual
			ColorMask[_ColorMask]

			Fog{ Color(0, 0, 0, 0) }

			CGPROGRAM
			#pragma target 3.0
			#pragma exclude_renderers gles

			#pragma shader_feature ___ _ALPHATEST_ON _ALPHAPREMULTIPLY_ON
			#pragma shader_feature ___ _SLIDE_BAR_SMOOTHNESS
			#pragma shader_feature ___ _SLIDE_BAR_METALLIC
			#pragma shader_feature ___ _MAPG_SPECULAR_INTENSITY
			#pragma shader_feature ___ _MAPB_FRESNEL_INTENSITY
			#pragma shader_feature ___ _SLIDE_BAR_REFLECTIVITY
			#pragma shader_feature ___ _NORMALMAP

			#pragma multi_compile_fwdadd_fullshadows
			#pragma multi_compile_fog

			#pragma vertex vertAddStandard
			#pragma fragment fragAddStandardMetallicR

			#include "LK PBS.cginc"

			ENDCG
		}

		Pass
		{
			Name "ShadowCaster"
			Tags{ "LightMode" = "ShadowCaster" }

			ZWrite On ZTest LEqual

			CGPROGRAM
			#pragma target 3.0
			#pragma exclude_renderers gles

			#pragma shader_feature ___ _ALPHATEST_ON _ALPHAPREMULTIPLY_ON
			#pragma multi_compile_shadowcaster

			#pragma vertex vertShadowCaster
			#pragma fragment fragShadowCaster

			#include "LK Shadow.cginc"

			ENDCG
		}
	}

	SubShader
	{
		Tags{ "RenderType" = "Opaque" "PerformanceChecks" = "False" "ForceNoShadowCasting" = "True" }
		LOD 150

		Pass
		{
			Name "FORWARD"
			Tags{ "LightMode" = "ForwardBase" }

			Blend[_SrcBlend][_DstBlend]
			ZWrite[_ZWrite]
			ColorMask[_ColorMask]

			CGPROGRAM
			#pragma target 2.0

			#pragma shader_feature ___ _ALPHATEST_ON _ALPHAPREMULTIPLY_ON
			#pragma shader_feature ___ _SLIDE_BAR_SMOOTHNESS
			#pragma shader_feature ___ _SLIDE_BAR_METALLIC
			#pragma shader_feature ___ _MAPG_SPECULAR_INTENSITY
			#pragma shader_feature ___ _MAPB_FRESNEL_INTENSITY
			#pragma shader_feature ___ _SLIDE_BAR_REFLECTIVITY
			#pragma shader_feature ___ _EMISSION
			#pragma skip_variants LIGHTMAP_ON SHADOWS_SCREEN

			#pragma multi_compile_fwdbase
			#pragma multi_compile_fog

			#pragma vertex VertDynamicObjectPBS
			#pragma fragment FragCharacterPBSMetallicRBrdf4

			#include "LK PBS.cginc"

			ENDCG
		}
	}

	FallBack "VertexLit"
	CustomEditor "PBSShaderScript"
}