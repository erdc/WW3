# c++ >= 11 save a lot of work. So let's use it!
QMAKE_CXXFLAGS += -std=c++17

# Keine Ahnung warum man das manuell angeben muss
QMAKE_CXXFLAGS += -fPIC


#
# -Wextra                          This are some extra warning flags that are not enabled by -Wall
#
# Enable this tests from time to time
#
# -Wsuggest-attribute=const        Warn about functions that might be candidates for attributes const.
# g++ 5 only begin
# -Wsuggest-final-types            Warn about types with virtual methods where code quality would be improved if the type were declared with the C++11 final specifier
# -Wsuggest-final-methods          Warn about virtual methods where code quality would be improved if the method were declared with the C++11 final specifier
#                                  Use -Wsuggest-final-types. Use  link time optimization for both.
# g++ 5 only end
# -Weffc++                         Warn about violations of the following style guidelines from Scott Meyers' Effective C++ series of books
# -Wdouble-promotion               Give a warning when a value of type float is implicitly promoted to double.
# -Wconversion                     Warn for implicit conversions that may alter a value
# -Wzero-as-null-pointer-constant  Warn when a literal '0' is used as null pointer constant.
# -Wlogical-op                     Warn about suspicious uses of logical operators in expressions.
# -Wpadded                         Warn if padding is included in a structure
# -Wredundant-decls                Warn if anything is declared more than once in the same scope
# -Wvector-operation-performance   Warn if vector operation is not implemented via SIMD capabilities

{
    GCC_VERSION = $$system($$QMAKE_CXX " -dumpversion")
    VERSIONS = $$split(GCC_VERSION, .)
    VERSION_MAJ = $$member(VERSIONS)

    contains(VERSION_MAJ, 9) {
        message( "g++ version 9.x found" )
        CONFIG += g++9
    }
    contains(VERSION_MAJ, 8) {
        message( "g++ version 8.x found" )
        CONFIG += g++8
    }
    contains(VERSION_MAJ, 7) {
        message( "g++ version 7.x found" )
        CONFIG += g++7
    }
    contains(VERSION_MAJ, 6) {
        message( "g++ version 6.x found" )
        CONFIG += g++6
    }
    contains(VERSION_MAJ, 5) {
        message( "g++ version 5.x found" )
        CONFIG += g++5
    }
    contains(VERSION_MAJ, 4) {
        message( "g++ version 4.x found" )
        CONFIG += g++4
    }
}

g++8: QMAKE_CXXFLAGS += -Werror=class-memaccess
g++9: QMAKE_CXXFLAGS += -Werror=init-list-lifetime -Werror=redundant-move -Werror=pessimizing-move -Werror=class-conversion

QMAKE_CXXFLAGS += -Wpedantic -Wextra -Wvector-operation-performance -Werror=multichar -Werror=switch -Werror=switch-enum -Werror=switch-bool -Werror=switch-unreachable
QMAKE_CXXFLAGS += -Werror=pedantic -pedantic-errors -Werror=unused-value -Werror=overflow -Werror=address
QMAKE_CXXFLAGS += -Werror=narrowing -Werror=return-type -Werror=unused-value -Werror=misleading-indentation -Werror=strict-aliasing -Werror=sign-promo -Werror=dangling-else
QMAKE_CXXFLAGS += -Werror=parentheses -Werror=logical-not-parentheses -Werror=logical-op -Werror=bool-compare -Werror=uninitialized -Werror=maybe-uninitialized
QMAKE_CXXFLAGS += -Werror=unused-result -Werror=reorder -Werror=duplicated-cond -Werror=duplicated-branches -Werror=type-limits -Werror=missing-field-initializers
# lot of work
# QMAKE_CXXFLAGS += -Werror=conversion
# we love pain
# QMAKE_CXXFLAGS += -Werror=shadow
# GCC 9
# QMAKE_CXXFLAGS += -Werror=init-list-lifetime -Werror=redundant-move -Werror=pessimizing-move -Werror=class-conversion


# Disable asserts for release build
QMAKE_CXXFLAGS_RELEASE += -DNDEBUG

# enable OpenMP for release build
QMAKE_CXXFLAGS_RELEASE +=  -fopenmp
QMAKE_LFLAGS_RELEASE +=  -fopenmp

#QMAKE_CXXFLAGS_DEBUG +=  -fsignaling-nans -fsanitize-undefined-trap-on-error -fsanitize=undefined -fsanitize=float-divide-by-zero

# more opti
# QMAKE_CXXFLAGS_RELEASE -= -O2
# QMAKE_CXXFLAGS_RELEASE += -O3

#profile
# QMAKE_CXXFLAGS += -pg
# QMAKE_LFLAGS += -pg

# link time optimazion
#QMAKE_CXXFLAGS += -flto
#QMAKE_LFLAGS += -flto

# GNU C++ library special modes
# https://gcc.gnu.org/onlinedocs/libstdc++/manual/debug_mode.html
# -D_GLIBCXX_DEBUG                Safe iterators, pre-conditions in the algorithms

# https://gcc.gnu.org/onlinedocs/libstdc++/manual/profile_mode.html
# -D_GLIBCXX_PROFILE              collect statistics, select optimal std::containser
QMAKE_CXXFLAGS_DEBUG +=  -D_GLIBCXX_DEBUG

# The QT library is not free of warnings. If we enable e.g. -Wconversion we see a lot
# of warnings from Qt we cant fix nor care. Here is a simple way to disable warnings
# for external framwork. Just  include library headers using -isystem instead of -I.
# This will make them "system headers" and GCC won't report warnings for them.
# This still trigger some warnings due to inline of code.
# If a standard system include directory, or a directory specified with -isystem, is also specified with -I, the -I option is ignored.
QMAKE_CXXFLAGS += -isystem$$[QT_INSTALL_HEADERS]/QtCore
QMAKE_CXXFLAGS += -isystem$$[QT_INSTALL_HEADERS]/QtGui
QMAKE_CXXFLAGS += -isystem$$[QT_INSTALL_HEADERS]/QtNetwork
QMAKE_CXXFLAGS += -isystem$$[QT_INSTALL_HEADERS]/QtOpenGL
QMAKE_CXXFLAGS += -isystem$$[QT_INSTALL_HEADERS]/QtScript
QMAKE_CXXFLAGS += -isystem$$[QT_INSTALL_HEADERS]/QtTest


