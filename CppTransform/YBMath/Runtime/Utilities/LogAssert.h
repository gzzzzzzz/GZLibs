#pragma once
#include "Runtime/mecanim/types.h"

#if DEBUGMODE

#define AssertIf(x)												\
	do {															\
		if(x)														\
		{															\
			DebugStringToFile (#x, 0, __FILE_STRIPPED__, __LINE__, kAssert);	\
			ASSERT_BREAK;											\
			ANALYSIS_ASSUME(!(x));									\
		}															\
	} while(0)

#define AssertIfObject(x,o)																	\
	do {																						\
		if(x)																					\
		{																						\
			DebugStringToFile (#x, 0, __FILE_STRIPPED__, __LINE__, kAssert, (o)?(o)->GetInstanceID():0);	\
			ASSERT_BREAK;																		\
			ANALYSIS_ASSUME(!(x));																\
		}																						\
	} while(0)


#define Assert(x)												\
	do {															\
		if(!(x))													\
		{															\
			DebugStringToFile (#x, 0, __FILE_STRIPPED__, __LINE__, kAssert);	\
			ASSERT_BREAK;											\
			ANALYSIS_ASSUME(x);										\
		}															\
	} while(0)

#define AssertMsg(x,...)															\
	do {																				\
		if(!(x))																		\
		{																				\
			DebugStringToFile (Format(__VA_ARGS__), 0, __FILE_STRIPPED__, __LINE__, kAssert);	\
			ASSERT_BREAK;																\
			ANALYSIS_ASSUME(x);															\
		}																				\
	} while(0)

#define AssertMsgObject(x,o,...)													\
	do {																				\
		if(!(x))																		\
		{																				\
			DebugStringToFile (Format(__VA_ARGS__), 0, __FILE_STRIPPED__, __LINE__, kAssert, (o) ? (o)->GetInstanceID() : 0);	\
			ASSERT_BREAK;																\
			ANALYSIS_ASSUME(x);															\
		}																				\
	} while(0)


#else

#define AssertIf(x) 		do { (void)sizeof(x); } while(0)
#define AssertIfObject(x,o) do { (void)sizeof(x); (void)sizeof(o); } while(0)
#define Assert(x) 			do { (void)sizeof(x); } while(0)
#define AssertMsg(x,...)	do { (void)sizeof(x); } while(0)
#define AssertMsgObject(x,o,...)	do { (void)sizeof(x); } while(0)

#endif

#define DebugAssertIf(x) AssertIf(x)
#define DebugAssert(x) Assert(x)
#define DebugAssertMsg(x, ...) AssertMsg(x, __VA_ARGS__)

#define AssertBreak(x)											\
	do {															\
		if(!(x))													\
		{															\
			DEBUG_BREAK;											\
			Assert(x);												\
		}															\
	} while(0)

typedef unsigned int UInt32;
typedef int SInt32;
typedef unsigned long UInt64;
typedef unsigned short UInt16;
typedef unsigned char UInt8;

#define EXPORT_COREMODULE
