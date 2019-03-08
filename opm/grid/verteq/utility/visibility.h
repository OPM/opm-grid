#ifndef OPM_VERTEQ_VISIBILITY_INCLUDED
#define OPM_VERTEQ_VISIBILITY_INCLUDED

/**
 * Macros to encapsulate symbol visibility on various platforms.
 * You need to define a separate macro for your package; symbols
 * that are exported in one dynamic shared object, may of course
 * be imported by another one!
 *
 * Usage:
 *
 * #include <opm/verteq/utility/visibility.h>
 *
 * #if defined (foo_EXPORTS)
 * #  define FOO_PUBLIC SYMBOL_IS_EXPORTED
 * #else
 * #  define FOO_PUBLIC SYMBOL_IS_IMPORTED
 * #endif
 * #define FOO_PRIVATE SYMBOL_IS_LOCALDEF
 *
 * struct FOO_PUBLIC Bar {
 * };
 *
 * int FOO_PRIVATE mumble ();
 */
#if defined (_WIN32)
#  define SYMBOL_IS_EXPORTED __declspec (dllexport)
#  define SYMBOL_IS_IMPORTED __declspec (dllimport)
#  define SYMBOL_IS_LOCALDEF
#else
#  if __GNUC__ >= 4
#    define SYMBOL_IS_EXPORTED __attribute__ ((visibility ("default")))
#    define SYMBOL_IS_IMPORTED __attribute__ ((visibility ("default")))
#    define SYMBOL_IS_LOCALDEF __attribute__ ((visibility ("hidden")))
#  else
#    define SYMBOL_IS_EXPORTED
#    define SYMBOL_IS_IMPORTED
#    define SYMBOL_IS_LOCALDEF
#  endif
#endif

#endif /* OPM_VERTEQ_VISIBILITY_INCLUDED */
