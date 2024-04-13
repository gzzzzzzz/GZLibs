#ifndef BASEOBJECT_H
#define BASEOBJECT_H

#include <string>
using namespace std;

template<class T>
class PPtr
{
	SInt32	m_InstanceID;
#if !UNITY_RELEASE
	mutable T* m_DEBUGPtr;
#endif

protected:

	inline void AssignObject(const void* o);

private:
	static string s_TypeString;

public:

	static const char* GetTypeString();
	static bool IsAnimationChannel() { return false; }
	static bool MightContainPPtr() { return true; }
	static bool AllowTransferOptimization() { return false; }

	template<class TransferFunction>
	void Transfer(TransferFunction& transfer);

	// Assignment
	explicit PPtr(int instanceID)
	{
		m_InstanceID = instanceID;
#if !UNITY_RELEASE
		m_DEBUGPtr = NULL;
#endif
	}
	PPtr(const T* o) { AssignObject(o); }
	PPtr(const PPtr<T>& o)
	{
		m_InstanceID = o.m_InstanceID;
#if !UNITY_RELEASE
		m_DEBUGPtr = NULL;
#endif
	}

	PPtr()
	{
#if !UNITY_RELEASE
		m_DEBUGPtr = NULL;
#endif
		m_InstanceID = 0;
	}

	PPtr& operator = (const T* o) { AssignObject(o); return *this; }
	PPtr& operator = (const PPtr<T>& o)
	{
#if !UNITY_RELEASE
		m_DEBUGPtr = NULL;
#endif
		m_InstanceID = o.m_InstanceID; return *this;
	}
	
	void SetInstanceID(int instanceID) { m_InstanceID = instanceID; }
	int GetInstanceID()const { return m_InstanceID; }

	// Comparison
	bool operator <  (const PPtr& p)const { return m_InstanceID < p.m_InstanceID; }
	bool operator == (const PPtr& p)const { return m_InstanceID == p.m_InstanceID; }
	bool operator != (const PPtr& p)const { return m_InstanceID != p.m_InstanceID; }

	// MSVC gets confused whether it should use operator bool(), or operator T* with implicit
	// comparison to NULL. So we add explicit functions and use them instead.
	bool IsNull() const;
	bool IsValid() const;

	operator T* () const
	{
		Assert("PPtr");
		return nullptr;
	}
	T* operator -> () const;
	T& operator * () const;
};

template<class T>
class ImmediatePtr
{
	mutable intptr_t m_Ptr;
#if !UNITY_RELEASE
	mutable T* m_DEBUGPtr;
#endif

	void AssignInstanceID(int instanceID)
	{
		AssertIf(instanceID & 1); m_Ptr = instanceID | 1; AssertIf((m_Ptr & 1) == 0);
#if !UNITY_RELEASE
		m_DEBUGPtr = NULL;
#endif
	}
	void AssignObject(const T* o)
	{
		m_Ptr = (intptr_t)o; AssertIf(m_Ptr & 1);
#if !UNITY_RELEASE
		m_DEBUGPtr = const_cast<T*>(o);
#endif
	}
	void Load() const
	{
		AssertIf((m_Ptr & 1) == 0);
		T* loaded = PPtr<T>(m_Ptr & (~1));
		m_Ptr = (intptr_t)(loaded);
		AssertIf(m_Ptr & 1);
#if !UNITY_RELEASE
		m_DEBUGPtr = loaded;
#endif
	}

	inline T* GetPtr() const
	{
		if ((m_Ptr & 1) == 0)
		{
			return (T*)(m_Ptr);
		}
		else
		{
			Load();
			return (T*)(m_Ptr);
		}
	}

	static string s_TypeString;

public:

	bool IsLoaded() const;

	static const char* GetTypeString();
	static bool IsAnimationChannel() { return false; }
	static bool MightContainPPtr() { return true; }
	static bool AllowTransferOptimization() { return false; }

	template<class TransferFunction>
	void Transfer(TransferFunction& transfer);

	// Assignment
	ImmediatePtr(const T* o) { AssignObject(o); }
	ImmediatePtr(const ImmediatePtr<T>& o) { m_Ptr = o.m_Ptr; }
	ImmediatePtr() { m_Ptr = 0; }

	ImmediatePtr& operator = (const T* o) { AssignObject(o); return *this; }

	void SetInstanceID(int instanceID) { AssignInstanceID(instanceID); }
	int GetInstanceID()const
	{
		if ((m_Ptr & 1) == 0 && m_Ptr != 0)
		{
			T* o = (T*)(m_Ptr);
			SInt32 instanceID = o->GetInstanceID();
			AssertIf(instanceID & 1);
			return instanceID;
		}
		else
			return m_Ptr & (~1);
	}

	inline bool operator == (const T* p)const { return GetPtr() == p; }
	inline bool operator != (const T* p)const { return GetPtr() != p; }

	inline operator T* () const { return GetPtr(); }
	inline T* operator -> () const { T* o = GetPtr(); AssertIf(o == NULL); return o; }
	inline T& operator * () const { T* o = GetPtr(); AssertIf(o == NULL); /*ANALYSIS_ASSUME(o);*/ return *o; }
};


#endif
