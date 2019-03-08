/**
 * Copyright (C) 2013 Uni Research AS.
 * This code is licensed under The GNU General Public License v3.0 or later.
 */
#include <opm/grid/verteq/utility/visibility.h>
#include <opm/grid/verteq/utility/exc.hpp>

#include <boost/exception/error_info.hpp>
#include <boost/exception/get_error_info.hpp>
#include <boost/exception/info.hpp>
#include <boost/exception/diagnostic_information.hpp>
#include <cstdarg>	// va_list, va_start, va_end
#include <cstdio>		// vsnprintf
#include <string>		// std::string
#ifdef __MSVC__
#	include <malloc.h>											// alloca
#endif

namespace Opm {
namespace Exc {

// this symbol must be exported so that other DSO can get it
struct SYMBOL_IS_EXPORTED message;

typedef boost::error_info <message, std::string> message_t;

Base&
Base::operator () (char const* fmt, ...) {
	// number of bytes necessary to render this string
	va_list params;
	va_start (params, fmt);
	int count = vsnprintf (0, 0, fmt, params) + 1;
	va_end (params);

	// allocate these bytes; this is necessary since we aren't allowed
	// to touch the internal buffer of a std::string. we make up for it
	// by allocating on the stack instead of on the free store!
	size_t bytes = count * sizeof (char);
	char* buf = static_cast<char*> (alloca (bytes));

	// print the string to those bytes. note that it is important to
	// restart the argument list, it cannot be reused from above!
	va_start (params, fmt);
	vsnprintf (buf, bytes, fmt, params);
	va_end (params);

	// put a copy of the buffer inside the exception structure
	// see <http://marknelson.us/2007/11/13/no-exceptions/>
	*this << message_t (std::string (buf));

	// allow the exception object to be chained
	return *this;
}

// this message is used if nothing else is specified. it is in the
// anonymous namespace so that no-one else is able to use it.
namespace {
char const* UNSPECIFIED = "<unspecified>";
}

char const*
Base::what () const throw () {
	// retrieve the stored reason, or a generic message if there is none
	std::string const* str = boost::get_error_info <message_t> (*this);
	return str ? str->c_str () : UNSPECIFIED;
}

std::string
diag_what (std::exception const& ex) {
	// header
	std::stringstream buf;
	buf << "Error";
	// std::exception info; this contains the reason
	char const* what = ex.what();
	if (what != UNSPECIFIED) {
		buf << ": " << what;
	}
	// boost::exception info; this contains the location
	boost::exception const* bex = dynamic_cast <boost::exception const*> (&ex);
	if (bex) {
		char const* const* file = boost::get_error_info <boost::throw_file> (*bex);
		int const* line = boost::get_error_info <boost::throw_line> (*bex);
		char const* const* func = boost::get_error_info <boost::throw_function> (*bex);
		if (func) {
			buf << ", in function " << *func;
		}
		if (file && line) {
			buf << ", at " << *file << ":" << *line;
		}
	}
	return buf.str();
}

} /* namespace Exc */
} /* namespace Opm */
