#include "Word.h"
#include <cstdarg>

std::string VFormat(const char* format, va_list ap)
{
	va_list zp;
	va_copy(zp, ap);
	char buffer[1024 * 10];
	vsnprintf(buffer, 1024 * 10, format, zp);
	va_end(zp);
	return buffer;
}

std::string Format(const char* format, ...)
{
	va_list va;
	va_start(va, format);
	std::string formatted = VFormat(format, va);
	va_end(va);
	return formatted;
}