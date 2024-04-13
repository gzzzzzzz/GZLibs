#pragma once
#include "../Graphics/Transform.h"
#include "../Math/Rect.h"
#include "../Geometry/Ray.h"

namespace YB
{
	void RectfToViewport(const Rectf& r, int viewPort[4]);

	bool CameraUnProject(const Vector3f& p, const Matrix4x4f& cameraToWorld, const Matrix4x4f& clipToWorld, const int viewport[4], Vector3f& outP);

	bool CameraProject(const Vector3f& p, const Matrix4x4f& cameraToWorld, const Matrix4x4f& worldToClip, const int viewport[4], Vector3f& outP);


	class Camera
	{
	private:
		mutable Matrix4x4f  m_WorldToCameraMatrix;
		mutable Matrix4x4f  m_ProjectionMatrix;
		mutable Matrix4x4f	m_WorldToClipMatrix;

		Rectf m_NormalizedViewPortRect = Rectf(0, 0, 1, 1);

		mutable bool		m_DirtyWorldToCameraMatrix = true;
		mutable bool        m_DirtyProjectionMatrix = true;
		mutable bool        m_DirtyWorldToClipMatrix = true;

		bool                m_ImplicitWorldToCameraMatrix;
		bool                m_ImplicitProjectionMatrix;


		float                 m_Depth; ///< A camera with a larger depth is drawn on top of a camera with a smaller depth range {-100, 100}

	public:
		Transform transform;

		float target_width = 16;
		float target_height = 16;

		bool                  m_Orthographic = false;
		float                 m_OrthographicSize = 5;
		float                 m_FieldOfView = 60;	///< Field of view of the camera range { 0.00001, 179 }
		float                 m_NearClip = 0.3f;		///< Near clipping plane
		float                 m_FarClip = 1000.0f;		///< Far clipping plane
		float                 m_Aspect = 1.0f;


	public:

		Rectf GetCameraRect(bool zeroOrigin) const;

		Rectf GetScreenViewportRect() const { return GetCameraRect(true); }

		Vector3f WorldToScreenPoint(const Vector3f& worldSpacePoint, bool* canProject = NULL) const;

		Vector3f ScreenToWorldPoint(const Vector3f& v) const;

		Matrix4x4f& GetWorldToClipMatrix() const;

		void GetClipToWorldMatrix(Matrix4x4f& outMatrix) const;

		const Matrix4x4f& GetProjectionMatrix() const;

		const Matrix4x4f& GetWorldToCameraMatrix() const;

		Matrix4x4f GetCameraToWorldMatrix() const;

		// Converts a screen point into a world space ray
		Ray ScreenPointToRay(const Vector2f& screenPos) const;
		// Converts a viewport point into a world space ray
		Ray ViewportPointToRay(const Vector2f& viewportPos) const;

		Vector3f ScreenToViewportPoint(const Vector3f& screenPos) const;
		Vector3f ViewportToScreenPoint(const Vector3f& viewPortPos) const;

		void Reset();

		void SetTargetSize(float width, float height);
	};
};

