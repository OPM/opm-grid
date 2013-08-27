#
# spec file for package opm-core
#

Name:           dune-cornerpoint
Version:        2013.03
Release:        0
Summary:        Cornerpoint grid management module for DUNE
License:        GPL-3.0
Group:          Development/Libraries/C and C++
Url:            http://www.opm-project.org/
Source0:        %{name}-%{version}.tar.gz
BuildRequires:  blas-devel gcc-c++ gcc-gfortran lapack-devel dune-common-devel
BuildRequires:  boost-devel git suitesparse-devel cmake28 doxygen bc
BuildRequires:  tinyxml-devel dune-istl-devel opm-core-devel dune-grid-devel
BuildRoot:      %{_tmppath}/%{name}-%{version}-build
Requires:       libdune-cornerpoint1 = %{version}

%description
This module enables working with corner-point or, more
generally, pillar grids.  A standard grid type in the petroleum
industry, corner-point grids fill space with a relatively low number
of cells while still providing sufficient flexibility to model faults,
fractures and erosion.  The grid format was originally designed with
an eye towards geological modelling rather than numerical simulation
and this design choice does limit the number of feasible numerical
methods.

%package -n libdune-cornerpoint1
Summary:        Cornerpoint grid management module for DUNE
Group:          System/Libraries

%description -n libdune-cornerpoint1
This module enables working with corner-point or, more
generally, pillar grids.  A standard grid type in the petroleum
industry, corner-point grids fill space with a relatively low number
of cells while still providing sufficient flexibility to model faults,
fractures and erosion.  The grid format was originally designed with
an eye towards geological modelling rather than numerical simulation
and this design choice does limit the number of feasible numerical
methods.

%package devel
Summary:        Development and header files for dune-cornerpoint
Group:          Development/Libraries/C and C++
Requires:       %{name} = %{version}
Requires:       blas-devel
Requires:       lapack-devel
Requires:       suitesparse-devel
Requires:       libdune-cornerpoint1 = %{version}

%description devel
This package contains the development and header files for dune-cornerpoint

%package doc
Summary:        Documentation files for dune-cornerpoint
Group:          Documentation
BuildArch:	noarch

%description doc
This package contains the documentation files for dune-cornerpoint

%package bin
Summary:        Applications in dune-cornerpoint
Group:          Scientific
Requires:       %{name} = %{version}
Requires:       libdune-cornerpoint1 = %{version}

%description bin
This package contains the applications for dune-cornerpoint

%prep
%setup -q

# consider using -DUSE_VERSIONED_DIR=ON if backporting
%build
cmake28 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DWHOLE_PROG_OPTIM=ON -DUSE_RUNPATH=OFF
make

%install
make install DESTDIR=${RPM_BUILD_ROOT}
make install-html DESTDIR=${RPM_BUILD_ROOT}

%clean
rm -rf %{buildroot}

%post -n libdune-cornerpoint1 -p /sbin/ldconfig

%postun -n libdune-cornerpoint1 -p /sbin/ldconfig

%files
%doc COPYING README

%files doc
%{_docdir}/*

%files -n libdune-cornerpoint1
%defattr(-,root,root,-)
%{_libdir}/*.so.*

%files devel
%defattr(-,root,root,-)
%{_libdir}/*.so
%{_libdir}/dunecontrol/*
%{_libdir}/pkgconfig/*
%{_includedir}/*
%{_datadir}/cmake/*

%files bin
%{_bindir}/*
