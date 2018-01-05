#
# spec file for package opm-grid
#

%define tag rc1

Name:           opm-grid
Version:        2017.10
Release:        0
Summary:        Cornerpoint grid management module for OPM
License:        GPL-3.0
Group:          Development/Libraries/C and C++
Url:            http://www.opm-project.org/
Source0:        https://github.com/OPM/%{name}/archive/release/%{version}/%{tag}.tar.gz#/%{name}-%{version}.tar.gz
BuildRequires: blas-devel lapack-devel dune-common-devel devtoolset-6-toolchain
BuildRequires:  git suitesparse-devel doxygen bc ecl-devel opm-common-devel
BuildRequires:  tinyxml-devel dune-istl-devel dune-grid-devel
%{?el6:BuildRequires: cmake3 boost148-devel}
%{!?el6:BuildRequires: cmake boost-devel}
BuildRequires:  opm-parser-devel opm-material-devel opm-output-devel
BuildRoot:      %{_tmppath}/%{name}-%{version}-build
Requires:       libopm-grid1 = %{version}

%description
This module enables working with corner-point or, more
generally, pillar grids.  A standard grid type in the petroleum
industry, corner-point grids fill space with a relatively low number
of cells while still providing sufficient flexibility to model faults,
fractures and erosion.  The grid format was originally designed with
an eye towards geological modelling rather than numerical simulation
and this design choice does limit the number of feasible numerical
methods.

%package -n libopm-grid1
Summary:        Cornerpoint grid management module for OPM
Group:          System/Libraries

%description -n libopm-grid1
This module enables working with corner-point or, more
generally, pillar grids.  A standard grid type in the petroleum
industry, corner-point grids fill space with a relatively low number
of cells while still providing sufficient flexibility to model faults,
fractures and erosion.  The grid format was originally designed with
an eye towards geological modelling rather than numerical simulation
and this design choice does limit the number of feasible numerical
methods.

%package devel
Summary:        Development and header files for opm-grid
Group:          Development/Libraries/C and C++
Requires:       %{name} = %{version}
Requires:       blas-devel
Requires:       lapack-devel
Requires:       suitesparse-devel
Requires:       libopm-grid1 = %{version}

%description devel
This package contains the development and header files for opm-grid

%package doc
Summary:        Documentation files for opm-grid
Group:          Documentation
BuildArch:	noarch

%description doc
This package contains the documentation files for opm-grid

%package bin
Summary:        Applications in opm-grid
Group:          Scientific
Requires:       %{name} = %{version}
Requires:       libopm-grid1 = %{version}

%description bin
This package contains the applications for opm-grid

%prep
%setup -q -n %{name}-release-%{version}-%{tag}

# consider using -DUSE_VERSIONED_DIR=ON if backporting
%build
scl enable devtoolset-6 bash
%{?el6:cmake3} %{!?el6:cmake} -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gcc -DCMAKE_Fortran_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gfortran %{?el6:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=%{_includedir}/boost148}
make

%install
make install DESTDIR=${RPM_BUILD_ROOT}
make install-html DESTDIR=${RPM_BUILD_ROOT}

%clean
rm -rf %{buildroot}

%post -n libopm-grid1 -p /sbin/ldconfig

%postun -n libopm-grid1 -p /sbin/ldconfig

%files
%doc COPYING README.md

%files doc
%{_docdir}/*

%files -n libopm-grid1
%defattr(-,root,root,-)
%{_libdir}/*.so.*

%files devel
%defattr(-,root,root,-)
%{_libdir}/*.so
%{_libdir}/dunecontrol/*
%{_libdir}/pkgconfig/*
%{_includedir}/*
%{_datadir}/cmake/*
%{_datadir}/opm/cmake/Modules/*

%files bin
%{_bindir}/*
