#pragma once
#include <string>

std::string VFormat(const char* format, va_list ap);
std::string Format(const char* format, ...);

