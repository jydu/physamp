%define _basename physamp
%define _version 1.0.0
%define _release 1
%define _prefix /usr

URL: http://home.gna.org/comap/

Name: %{_basename}
Version: %{_version}
Release: %{_release}
License: GPL-3.0
Vendor: Julien Dutheil
Source: %{_basename}-%{_version}.tar.gz
Summary: The PhySamp package
Group: Productivity/Scientific/Other

Requires: libbpp-phyl9 = 2.3.0
Requires: libbpp-seq9 = 2.3.0
Requires: libbpp-core2 = 2.3.0

BuildRoot: %{_builddir}/%{_basename}-root
BuildRequires: cmake >= 2.8.11
BuildRequires: gcc-c++ >= 4.7.0
BuildRequires: groff
BuildRequires: texinfo >= 4.0.0
BuildRequires: libbpp-core2 = 2.3.0
BuildRequires: libbpp-core-devel = 2.3.0
BuildRequires: libbpp-seq9 = 2.3.0
BuildRequires: libbpp-seq-devel = 2.3.0
BuildRequires: libbpp-phyl9 = 2.3.0
BuildRequires: libbpp-phyl-devel = 2.3.0


AutoReq: yes
AutoProv: yes
%if 0%{?mdkversion}
%if 0%{?mdkversion} >= 201100
BuildRequires: xz
%define zipext xz
%else
BuildRequires: lzma
%define zipext lzma
%endif
%else
BuildRequires: gzip
%define zipext gz
%endif

%description
Includes programs:
 - bppalnoptim, optimize the amount of missing data in a sequence alignment,
 - bppphysamp, sample species according phylogeny

%prep
%setup -q

%build
CFLAGS="-I%{_prefix}/include $RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix}"
if [ %{_lib} == 'lib64' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DLIB_SUFFIX=64"
fi
if [ %{zipext} == 'lzma' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DDOC_COMPRESS=lzma -DDOC_COMPRESS_EXT=lzma"
fi
if [ %{zipext} == 'xz' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DDOC_COMPRESS=xz -DDOC_COMPRESS_EXT=xz"
fi

cmake $CMAKE_FLAGS .
make
make info

%install
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root)
%doc AUTHORS LICENSE INSTALL ChangeLog
%{_prefix}/bin/bppalnoptim
%{_prefix}/bin/bppphysamp
%{_prefix}/share/info/physamp.info.%{zipext}
%{_prefix}/share/man/man1/bppalnoptim.1.%{zipext}
%{_prefix}/share/man/man1/bppphysamp.1.%{zipext}

%changelog
* Fri May 26 2017 Julien Dutheil <dutheil@evolbio.mpg.de> 1.0.0-1
- Initial release

