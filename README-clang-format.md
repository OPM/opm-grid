# What is clang-format

The clang-format program is a tool to reformat source code files,
according to a given rule set. In particular, it is useful to reformat
C++, as it leverages the clang/llvm frontend to parse the code,
meaning that the tool does not understand the code as pure text, but
as code. This can be leveraged to make different rules for the
formatting of functions, classes, namespaces and more.

Since clang-format understands the structure of code, it should never
change the meaning of it, making it a safe tool to apply in an
automated setting.

# Installing clang-format

On most platforms installing clang-format should be fairly easy, using
pre-built packages.

## Linux

Tested on Ubuntu 18.10 (cosmic):

    sudo apt-get install clang-format

This installs clang-format version 7.0.0 at the time of writing
(September 2019). If you get an error message like "Unable to locate
package clang-format", you may need to update the package database
first:

    sudo apt-get update

and then run the install command above.

## macOS

Using Homebrew:

    brew install clang-format

or if you need to upgrade an old installation

    brew upgrade clang-format

This installs clang-format version 8.0.0 at the time of writing
(September 2019).

## Windows

Visual Studio 2017 includes and integrates clang-format, so it should
not be necessary to download it separately if you use this IDE, unless
you want to run it on the command line.

# Running clang-format on the command line

Running the following command modifies File.hpp in-place:

    clang-format -i File.hpp

You can the use git diff to see what was changed. If you want to
avoid changing the file, running

    clang-format File.hpp

will print the reformatted code to the terminal, you can redirect it
to a new file as follows:

    clang-format File.hpp > ReformattedFile.hpp

# Controlling clang-format

Clang-format is controlled by the rule-set selected. Typically such a
rule-set is read from a ".clang-format" file. If such a file is found
in the current directory it will be applied, otherwise it will look in
the parent directories for such a file until it is found. Since such a
file is present in the "opm-grid" directory, all code contained in
that directory will be formatted according to it (when running
clang-format on a file there). Putting a ".clang-format" file in your
home directory would make it apply to all your files, unless there are
more local files that apply, such as the one in opm-grid.

Rulesets are typically based on a broad style, the one in opm-grid is
based on WebKit. On top of that, a large number of style options are
available, see
https://clang.llvm.org/docs/ClangFormatStyleOptions.html
for complete documentation.


# Integration with your editor

## emacs

Download "clang-format.el" from
https://llvm.org/svn/llvm-project/cfe/trunk/tools/clang-format/
and put it in your emacs lisp folder, typically named ".emacs.d/lisp".

Then add the following to your ".emacs" file:

    ;; clang-format-region macro
    (require 'clang-format)

You can now invoke clang format by selecting a region and doing

    M-x clang-format-region

or to change the entire file

    M-x clang-format-buffer

## Visual Studio 2017

See the notes by Magne Sjaastad:

https://magnesj.github.io/resinsight-system-doc/system/clang-format-notes

# Integration with git

It is possible to have clang-format apply automatically via a git
commit hook. We do not at this point enforce this, as it can have some
annoying consequences, for example if developers have different
versions of clang-format.
