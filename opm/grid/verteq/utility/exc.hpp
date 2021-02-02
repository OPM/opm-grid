/**
 * Copyright (C) 2013 Uni Research AS.
 * This code is licensed under The GNU General Public License v3.0 or later.
 */
#include <exception> // std::exception
#ifndef UUID_274DA366004E11DCB1DDFE2E56D89593
#include <boost/exception/exception.hpp> // boost::exception
#endif
#ifndef UUID_8D22C4CA9CC811DCAA9133D256D89593
#include <boost/exception/info.hpp> // boost::throw_{function,file,line}
#endif

#ifndef OPM_VERTEQ_VISIBILITY_INCLUDED
#include <opm/grid/verteq/utility/visibility.h>
#endif
#if defined (opmverteq_EXPORTS)
#  define OPM_VERTEQ_PUBLIC SYMBOL_IS_EXPORTED
#else
#  define OPM_VERTEQ_PUBLIC SYMBOL_IS_IMPORTED
#endif
#define OPM_VERTEQ_PRIVATE SYMBOL_IS_LOCALDEF

// enable printf-style attribute for MSVC
#ifdef _MSC_VER
#  if _MSC_VER >= 1400
#    define _USE_ATTRIBUTES_FOR_SAL 1
#    include <sal.h>
#  endif
#endif

namespace Opm {
namespace Exc {
/**
 * Base class for exceptions thrown from the OPM framework.
 *
 * Throw exceptions like this:
 *
 *	throw OPM_EXC ("Magic value = %d", 42);
 *
 * Catch exceptions like this:
 *
 *	catch (std::exception& ex) {
 *		std::cerr << OPM_WHAT (ex) << std::endl;
 *		exit (-1);
 *	}
 */
struct OPM_VERTEQ_PUBLIC Base
	: public virtual std::exception
	, public virtual boost::exception {

	/**
	 * @brief Add a message to the string
	 * @param fmt printf-style format string. The rest of the parameters
	 *	are formatted according to this string.
	 * @return Same exception object so it can be chained. (In particular,
	 *	more boost::error_info objects may be added).
	 */
	virtual Base& operator () (
#ifdef _MSC_VER
#  if _MSC_VER >= 1400
		__format_string
#  endif
#endif
			char const* fmt, ...)
	// there is an implicit this parameter at the start, so the format
	// string is the second parameter, and the variable argument list
	// starts with the third.
#ifdef __GNUC__
	__attribute__ ((format (printf, 2, 3)))
#endif
	;

	/**
	 * @brief Message created at the throw-site about the error.
	 * @return String that can be printed in the log about the exception.
	 */
	virtual char const* what () const throw ();
};

/**
 * @brief Retrieve information about the code that failed
 * @param ex Exception that was thrown
 * @return Text containing error information and location
 */
std::string OPM_VERTEQ_PUBLIC diag_what (std::exception const& ex);

/**
* Create a new exception object, possibly filled with location
* information if a debug build was done.
*/
#ifdef DEBUG
// const_cast is necessary because the operator<< returns a
// const, and we may want it mutable to add more information
	#define OPM_EXC const_cast <Opm::Exc::Base&>(Opm::Exc::Base ()\
		<< ::boost::throw_function (BOOST_CURRENT_FUNCTION)\
		<< ::boost::throw_file (__FILE__)\
		<< ::boost::throw_line (static_cast <int> (__LINE__))\
)
#else
	#define OPM_EXC Opm::Exc::Base ()
#endif

#define OPM_WHAT(ex) Opm::Exc::diag_what (ex)

} /* namespace Opm::Exc */
} /* namespace Opm */
