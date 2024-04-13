// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

Shader "Cartoon Sihoutte"
{
	Properties
	{
		_Color("Color (RGBA)", Color) = (1, 1, 1, 1)
		_MainTex("Albedo (RGBA)", 2D) = "white" {}
		_ScaleSize("Scale Size Along Normals", Float) = 0.01
		_SillhoutteColor("Sillhoutte Color", Color) = (0, 0, 0, 1)
		_StruckColor("Struck Color", Color) = (1, 0, 0, 1)
		_TimeSinceApply("Time Since Apply", Float) = 0
		_TimePeriod("Time Period", Range(0.01, 2)) = 0.2

		_Angle ("Angle", Range(-3.14, 3.14)) = 0
		_StartAmount("StartAmount", Range(0.01, 1)) = 0.2
		_DissolveSpeed ("DissolveSpeed", float) = 0.3	//消融速度
		_ColorAnimate ("ColorAnimate", vector) = (0,1,0,0)
		_FireColor ("FireColor", Color) = (1, 1, 0, 1)
		_NoiseTex ("NoiseTex", 2D) = "white" {}
		
		_RedFadeSpeed ("RedSpeed", Range(0.1, 10.0)) = 0.5
		_GreenFadeSpeed ("GreenSpeed", Range(0.1, 10.0)) = 2
		_BlueFadeSpeed ("BlueSpeed", Range(0.1, 10.0)) = 2
	}
	SubShader
	{
		LOD 300
		Tags{ "RenderType" = "Opaque" "PerformanceChecks" = "False" }

		Pass
		{
			Name "SIMPLE"

			CGPROGRAM
			#pragma target 2.0

			#pragma vertex vertSimple
			#pragma fragment fragSimple
			#pragma multi_compile __ _GREY
			#pragma multi_compile __ _HURT
			#pragma multi_compile __ _DISSOLVE
			#include "UnityCG.cginc"

			struct VertexInputSimple
			{
				float4 vertex : POSITION;
				float2 uv0 : TEXCOORD0;
			};

			struct VertexOutputSimple
			{
				float4 pos : SV_POSITION;
				float2 tex : TEXCOORD0;

				float2 uv_NoiseTex : TEXCOORD1;
				//float2 uv_MatCap : TEXCOORD2;
			};

			sampler2D _MainTex;
			float4 _MainTex_ST;
			half4 _Color;
			half _TimeSinceApply;
			half _TimePeriod;

			float _DissolveSpeed;
			float _StartAmount;
			float _Angle;
			fixed4 _FireColor;
			half4 _ColorAnimate;
			static half3 static_Color = float3(1,1,1);
			sampler2D _NoiseTex;
			float4 _NoiseTex_ST;
			
			float _RedFadeSpeed;
			float _GreenFadeSpeed;
			float _BlueFadeSpeed;

			VertexOutputSimple vertSimple(VertexInputSimple v)
			{
				VertexOutputSimple o;
				UNITY_INITIALIZE_OUTPUT(VertexOutputSimple, o);

				o.pos = UnityObjectToClipPos(v.vertex);
				o.tex = TRANSFORM_TEX(v.uv0, _MainTex);
				#if _DISSOLVE
					o.uv_NoiseTex = TRANSFORM_TEX(v.uv0, _NoiseTex);
					/*
					float3 worldNorm = normalize(unity_WorldToObject[0].xyz * v.normal.x + unity_WorldToObject[1].xyz * v.normal.y + unity_WorldToObject[2].xyz * v.normal.z);
					worldNorm = mul((float3x3)UNITY_MATRIX_V, worldNorm);
					float2 direct;
					float s = sin(_Angle);
					float c = cos(_Angle);
					direct.x = dot(worldNorm.xy, float2(c, -s));
					direct.y = dot(worldNorm.xy, float2(s, c));
					o.uv_MatCap.xy = direct * 0.5 + 0.5;
					*/
				#endif
				return o;
			}

			half4 fragSimple(VertexOutputSimple i) : SV_Target
			{
				half4 albedo = tex2D(_MainTex, i.tex);
				#if _GREY
				    half c3 = 0.3333f * (albedo.r + albedo.g + albedo.b);
				    albedo.rgb = c3;
				#elif _HURT
					half progress  = (_Time.y-_TimeSinceApply) / _TimePeriod;//0~infinity
					albedo.g = lerp(1, albedo.g, saturate(progress*_GreenFadeSpeed));
					albedo.b = lerp(1, albedo.b, saturate(progress*_BlueFadeSpeed));
					albedo.r = lerp(1, albedo.r, saturate(progress*_RedFadeSpeed));
				#elif _DISSOLVE
					fixed4 noi =  tex2D(_NoiseTex, i.uv_NoiseTex);
					half timeVal = (_Time.y - _TimeSinceApply) * _DissolveSpeed;
					half clipAmount = noi.r - timeVal;
					clip(clipAmount);
					if (clipAmount < _StartAmount)
					{
						half lit = clipAmount / _StartAmount;
						static_Color.x = lerp(_FireColor.r, lit, _ColorAnimate.x);
						static_Color.y = lerp(_FireColor.g, lit, _ColorAnimate.y);
						static_Color.z = lerp(_FireColor.b, lit, _ColorAnimate.z);
						half mulFactor = static_Color.x + static_Color.y + static_Color.z;
						albedo.rgb *= static_Color.xyz * mulFactor * mulFactor;
					}
				#else
				    albedo *= _Color;
				#endif
				return half4(albedo.rgb, 0.9);
			}

			ENDCG
		}
		
		Pass
		{
			Name "SILHOUTTE"

			Cull Front

			CGPROGRAM
			#pragma target 2.0

			#pragma vertex vertSihoutte
			#pragma fragment fragSihoutte
			#pragma multi_compile __ _STRUCK
			#pragma multi_compile __ _DISSOLVE
			#include "UnityCG.cginc"

			struct VertexInputSihoutte
			{
				float4 vertex : POSITION;
				float3 normal : NORMAL;
			};

			struct VertexOutputSihoutte
			{
				float4 pos : SV_POSITION;
			};

			half _ScaleSize;
			half4 _SillhoutteColor;
			half _TimeSinceApply;
			half _TimePeriod;
			half4 _StruckColor;

			VertexOutputSihoutte vertSihoutte(VertexInputSihoutte v)
			{
				VertexOutputSihoutte o;
				UNITY_INITIALIZE_OUTPUT(VertexOutputSihoutte, o);
				half silhouetteSize = _ScaleSize;
				#if _STRUCK
					half progress  = (_Time.y-_TimeSinceApply) / _TimePeriod;//0~infinity
					half lamda = saturate(1 - progress);
					silhouetteSize = lerp(0.01f, _ScaleSize, lamda);
				#endif
				o.pos = UnityObjectToClipPos(half4(v.vertex.xyz + normalize(v.normal) * silhouetteSize, 1.0));
				return o;
			}

			half4 fragSihoutte(VertexOutputSihoutte i) : SV_Target
			{
			#if _STRUCK
				//half progress  = (_Time.y-_TimeSinceApply) / _TimePeriod;//0~infinity
				//half lamda = saturate(1 - progress);
				return _StruckColor;//lerp(_SillhoutteColor, _StruckColor, lamda);
			#endif
			#if _DISSOLVE
				clip(-1);
			#endif
				return _SillhoutteColor;
			}

			ENDCG
		}
	}
	SubShader
	{
		LOD 200
		Tags{ "RenderType" = "Opaque" "PerformanceChecks" = "False" }

		Pass
		{
			Name "SIMPLE"

			CGPROGRAM
			#pragma target 2.0

			#pragma vertex vertSimple
			#pragma fragment fragSimple
			#pragma multi_compile __ _GREY
			#pragma multi_compile __ _HURT
			#pragma multi_compile __ _DISSOLVE
			#include "UnityCG.cginc"

			struct VertexInputSimple
			{
				float4 vertex : POSITION;
				float2 uv0 : TEXCOORD0;
			};

			struct VertexOutputSimple
			{
				float4 pos : SV_POSITION;
				float2 tex : TEXCOORD0;

				float2 uv_NoiseTex : TEXCOORD1;
				//float2 uv_MatCap : TEXCOORD2;
			};

			sampler2D _MainTex;
			float4 _MainTex_ST;
			half4 _Color;
			half _TimeSinceApply;
			half _TimePeriod;

			float _DissolveSpeed;
			float _StartAmount;
			float _Angle;
			fixed4 _FireColor;
			half4 _ColorAnimate;
			static half3 static_Color = float3(1,1,1);
			sampler2D _NoiseTex;
			float4 _NoiseTex_ST;
			
			float _RedFadeSpeed;
			float _GreenFadeSpeed;
			float _BlueFadeSpeed;

			VertexOutputSimple vertSimple(VertexInputSimple v)
			{
				VertexOutputSimple o;
				UNITY_INITIALIZE_OUTPUT(VertexOutputSimple, o);

				o.pos = UnityObjectToClipPos(v.vertex);
				o.tex = TRANSFORM_TEX(v.uv0, _MainTex);
				#if _DISSOLVE
					o.uv_NoiseTex = TRANSFORM_TEX(v.uv0, _NoiseTex);
					/*
					float3 worldNorm = normalize(unity_WorldToObject[0].xyz * v.normal.x + unity_WorldToObject[1].xyz * v.normal.y + unity_WorldToObject[2].xyz * v.normal.z);
					worldNorm = mul((float3x3)UNITY_MATRIX_V, worldNorm);
					float2 direct;
					float s = sin(_Angle);
					float c = cos(_Angle);
					direct.x = dot(worldNorm.xy, float2(c, -s));
					direct.y = dot(worldNorm.xy, float2(s, c));
					o.uv_MatCap.xy = direct * 0.5 + 0.5;
					*/
				#endif
				return o;
			}

			half4 fragSimple(VertexOutputSimple i) : SV_Target
			{
				half4 albedo = tex2D(_MainTex, i.tex);
				#if _GREY
				    half c3 = 0.3333f * (albedo.r + albedo.g + albedo.b);
				    albedo.rgb = c3;
				#elif _HURT
					half progress  = (_Time.y-_TimeSinceApply) / _TimePeriod;//0~infinity
					albedo.g = lerp(1, albedo.g, saturate(progress*_GreenFadeSpeed));
					albedo.b = lerp(1, albedo.b, saturate(progress*_BlueFadeSpeed));
					albedo.r = lerp(1, albedo.r, saturate(progress*_RedFadeSpeed));
				#elif _DISSOLVE
					fixed4 noi =  tex2D(_NoiseTex, i.uv_NoiseTex);
					half timeVal = (_Time.y - _TimeSinceApply) * _DissolveSpeed;
					half clipAmount = noi.r - timeVal;
					clip(clipAmount);
					if (clipAmount < _StartAmount)
					{
						half lit = clipAmount / _StartAmount;
						static_Color.x = lerp(_FireColor.r, lit, _ColorAnimate.x);
						static_Color.y = lerp(_FireColor.g, lit, _ColorAnimate.y);
						static_Color.z = lerp(_FireColor.b, lit, _ColorAnimate.z);
						half mulFactor = static_Color.x + static_Color.y + static_Color.z;
						albedo.rgb *= static_Color.xyz * mulFactor * mulFactor;
					}
				#else
				    albedo *= _Color;
				#endif
				return half4(albedo.rgb, 0.9);
			}

			ENDCG
		}
	}
	FallBack "VertexLit"
}
