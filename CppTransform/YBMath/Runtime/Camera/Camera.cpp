#include "Camera.h"
#include "../Utilities/Word.h"
#include <iostream>

namespace YB
{
	void RectfToViewport(const Rectf& r, int viewPort[4])
	{
		// We have to take care that the viewport doesn't exceed the buffer size (case 569703).
		// Bad rounding to integer makes D3D11 crash and burn.
		viewPort[0] = RoundfToInt(r.x);
		viewPort[1] = RoundfToInt(r.y);
		viewPort[2] = RoundfToIntPos(r.GetRight()) - viewPort[0];
		viewPort[3] = RoundfToIntPos(r.GetBottom()) - viewPort[1];
	}


	bool CameraUnProject(const Vector3f& p, const Matrix4x4f& cameraToWorld, const Matrix4x4f& clipToWorld, const int viewport[4], Vector3f& outP)
	{
		// pixels to -1..1
		Vector3f in;
		in.x = (p.x - viewport[0]) * 2.0f / viewport[2] - 1.0f;
		in.y = (p.y - viewport[1]) * 2.0f / viewport[3] - 1.0f;
		// It does not matter where the point we unproject lies in depth; so we choose 0.95, which
		// is further than near plane and closer than far plane, for precision reasons.
		// In a perspective camera setup (near=0.1, far=1000), a point at 0.95 projected depth is about
		// 5 units from the camera.
		in.z = 0.95f;

#if UNITY_WP8
		RotatePointIfNeeded(in, true);
#endif

		Vector3f pointOnPlane;
		if (clipToWorld.PerspectiveMultiplyPoint3(in, pointOnPlane))
		{
			// Now we have a point on the plane perpendicular to the viewing direction. We need to return the one that is on the line
			// towards this point, and at p.z distance along camera's viewing axis.
			Vector3f cameraPos = cameraToWorld.GetPosition();
			Vector3f dir = pointOnPlane - cameraPos;

			// The camera/projection matrices follow OpenGL convention: positive Z is towards the viewer.
			// So negate it to get into Unity convention.
			Vector3f forward = -cameraToWorld.GetAxisZ();
			float distToPlane = Dot(dir, forward);
			if (Abs(distToPlane) >= 1.0e-6f)
			{
				bool isPerspective = (clipToWorld.m_Data[3] != 0.0f || clipToWorld.m_Data[7] != 0.0f || clipToWorld.m_Data[11] != 0.0f || clipToWorld.m_Data[15] != 1.0f);
				if (isPerspective)
				{
					dir *= p.z / distToPlane;
					outP = cameraPos + dir;
				}
				else
				{
					outP = pointOnPlane - forward * (distToPlane - p.z);
				}
				return true;
			}
		}

		outP.Set(0.0f, 0.0f, 0.0f);
		return false;
	}


	static inline Rectf GetCameraTargetRect(const Camera& camera, bool zeroOrigin)
	{
		return Rectf(0, 0, camera.target_width, camera.target_height);
	}

	bool CameraProject(const Vector3f& p, const Matrix4x4f& cameraToWorld, const Matrix4x4f& worldToClip, const int viewport[4], Vector3f& outP)
	{
		Vector3f clipPoint;
		if (worldToClip.PerspectiveMultiplyPoint3(p, clipPoint))
		{
			Vector3f cameraPos = cameraToWorld.GetPosition();
			Vector3f dir = p - cameraPos;
			// The camera/projection matrices follow OpenGL convention: positive Z is towards the viewer.
			// So negate it to get into Unity convention.
			Vector3f forward = -cameraToWorld.GetAxisZ();
			float dist = Dot(dir, forward);

#if UNITY_WP8
			RotatePointIfNeeded(clipPoint, false);
#endif

			outP.x = viewport[0] + (1.0f + clipPoint.x) * viewport[2] * 0.5f;
			outP.y = viewport[1] + (1.0f + clipPoint.y) * viewport[3] * 0.5f;
			//outP.z = (1.0f + clipPoint.z) * 0.5f;
			outP.z = dist;

			return true;
		}

		outP.Set(0.0f, 0.0f, 0.0f);
		return false;
	}


	Rectf Camera::GetCameraRect(bool zeroOrigin) const
	{
		// Get the screen rect from either the target texture or the viewport we're inside
		Rectf screenRect = GetCameraTargetRect(*this, zeroOrigin);

		// Now figure out how large this camera is depending on the normalized viewRect.
		Rectf viewRect = m_NormalizedViewPortRect;
		viewRect.Scale(screenRect.width, screenRect.height);
		viewRect.Move(screenRect.x, screenRect.y);
		viewRect.Clamp(screenRect);
		return viewRect;
	}

	const Matrix4x4f& Camera::GetProjectionMatrix() const
	{
		//if (m_DirtyProjectionMatrix && m_ImplicitProjectionMatrix)
		{
			if (!m_Orthographic)
				m_ProjectionMatrix.SetPerspective(m_FieldOfView, m_Aspect, m_NearClip, m_FarClip);
			else
				m_ProjectionMatrix.SetOrtho(-m_OrthographicSize * m_Aspect, m_OrthographicSize * m_Aspect, -m_OrthographicSize, m_OrthographicSize, m_NearClip, m_FarClip);


			//m_DirtyProjectionMatrix = false;
		}
		return m_ProjectionMatrix;
	}

	Matrix4x4f& Camera::GetWorldToClipMatrix() const
	{
		//if (m_DirtyWorldToClipMatrix)
		{
			MultiplyMatrices4x4(&GetProjectionMatrix(), &GetWorldToCameraMatrix(), &m_WorldToClipMatrix);
			//m_DirtyWorldToClipMatrix = false;
		}
		return m_WorldToClipMatrix;
	}

	void Camera::GetClipToWorldMatrix(Matrix4x4f& outMatrix) const
	{
		Matrix4x4f::Invert_Full(GetWorldToClipMatrix(), outMatrix);
	}

	Matrix4x4f Camera::GetCameraToWorldMatrix() const
	{
		Matrix4x4f m;
		Matrix4x4f::Invert_Full(GetWorldToCameraMatrix(), m);
		return m;
	}

	Vector3f Camera::WorldToScreenPoint(const Vector3f& v, bool* canProject) const
	{
		int viewPort[4];
		RectfToViewport(GetScreenViewportRect(), viewPort);

		Vector3f out;
		bool ok = CameraProject(v, GetCameraToWorldMatrix(), GetWorldToClipMatrix(), viewPort, out);
		if (canProject != NULL)
			*canProject = ok;
		return out;
	}

	Vector3f Camera::ScreenToWorldPoint(const Vector3f& v) const
	{
		int viewPort[4];
		RectfToViewport(GetScreenViewportRect(), viewPort);

		Vector3f out;
		Matrix4x4f clipToWorld;
		GetClipToWorldMatrix(clipToWorld);
		if (!CameraUnProject(v, GetCameraToWorldMatrix(), clipToWorld, viewPort, out))
		{
			std::cout << Format("Screen position out of view frustum (screen pos %f, %f, %f) (Camera rect %d %d %d %d)", v.x, v.y, v.z, viewPort[0], viewPort[1], viewPort[2], viewPort[3]) << std::endl;
		}
		return out;
	}

	void Camera::Reset()
	{
		m_NormalizedViewPortRect = Rectf(0, 0, 1, 1);

		m_Depth = 0.0F;
		m_NearClip = 0.3F;
		m_FarClip = 1000.0F;
		m_Aspect = 1.0F;
		m_Orthographic = false;

		m_OrthographicSize = 5.0F;
		m_FieldOfView = 60.0F;

		m_DirtyWorldToCameraMatrix = m_DirtyProjectionMatrix = m_DirtyWorldToClipMatrix = true;

	}


	const Matrix4x4f& Camera::GetWorldToCameraMatrix() const
	{
		/*if (m_DirtyWorldToCameraMatrix && m_ImplicitWorldToCameraMatrix)
		{
			m_WorldToCameraMatrix.SetScale(Vector3f(1.0F, 1.0F, -1.0F));
			m_WorldToCameraMatrix *= transform.GetWorldToLocalMatrixNoScale();
			m_DirtyWorldToCameraMatrix = false;
		}*/

		m_WorldToCameraMatrix.SetScale(Vector3f(1.0F, 1.0F, -1.0F));
		m_WorldToCameraMatrix *= transform.GetWorldToLocalMatrixNoScale();
		return m_WorldToCameraMatrix;
	}

	void Camera::SetTargetSize(float width, float height)
	{
		target_width = width;
		target_height = height;
		m_Aspect = width / height;
	}

	Ray Camera::ScreenPointToRay(const Vector2f& viewPortPos) const
	{
		int viewPort[4];
		RectfToViewport(GetScreenViewportRect(), viewPort);

		Ray ray;
		Vector3f out;
		Matrix4x4f clipToWorld;
		GetClipToWorldMatrix(clipToWorld);

		const Matrix4x4f& camToWorld = GetCameraToWorldMatrix();
		if (!CameraUnProject(Vector3f(viewPortPos.x, viewPortPos.y, m_NearClip), camToWorld, clipToWorld, viewPort, out))
		{
			if (viewPort[0] > 0 || viewPort[1] > 0 || viewPort[2] > 0 || viewPort[3] > 0)
			{
				std::cout << (Format("Screen position out of view frustum (screen pos %f, %f) (Camera rect %d %d %d %d)", viewPortPos.x, viewPortPos.y, viewPort[0], viewPort[1], viewPort[2], viewPort[3])) << endl;
			}
			return Ray(transform.GetPosition(), Vector3f(0, 0, 1));
		}
		ray.SetOrigin(out);
		if (!CameraUnProject(Vector3f(viewPortPos.x, viewPortPos.y, m_NearClip + 1.0f), camToWorld, clipToWorld, viewPort, out))
		{
			if (viewPort[0] > 0 || viewPort[1] > 0 || viewPort[2] > 0 || viewPort[3] > 0)
			{
				std::cout << (Format("Screen position out of view frustum (screen pos %f, %f) (Camera rect %d %d %d %d)", viewPortPos.x, viewPortPos.y, viewPort[0], viewPort[1], viewPort[2], viewPort[3])) << endl;
			}
			return Ray(transform.GetPosition(), Vector3f(0, 0, 1));
		}
		Vector3f dir = out - ray.GetOrigin();
		ray.SetDirection(Normalize(dir));

		return ray;
	}

	Ray Camera::ViewportPointToRay(const Vector2f& viewPortPos) const {
		Vector3f screenPos = ViewportToScreenPoint(Vector3f(viewPortPos.x, viewPortPos.y, 0.0F));
		return ScreenPointToRay(Vector2f(screenPos.x, screenPos.y));
	}

	Vector3f Camera::ScreenToViewportPoint(const Vector3f& screenPos) const
	{
		Rectf r = GetScreenViewportRect();
		float nx = (screenPos.x - r.x) / r.Width();
		float ny = (screenPos.y - r.y) / r.Height();
		return Vector3f(nx, ny, screenPos.z);
	}

	Vector3f Camera::ViewportToScreenPoint(const Vector3f& viewPos) const
	{
		Rectf r = GetScreenViewportRect();
		float nx = viewPos.x * r.Width() + r.x;
		float ny = viewPos.y * r.Height() + r.y;
		return Vector3f(nx, ny, viewPos.z);
	}
};
