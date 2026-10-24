# RPM recipe for RCV (used by packaging/linux/build_linux.sh, which passes rcv_version and rcv_stage).
# The package holds the ready-made PyInstaller build: nothing is compiled, stripped or scanned for libraries
# (the bundled libraries are private to /opt/rcv and must not become package dependencies or provides).
%global debug_package %{nil}
%global __os_install_post %{nil}
%global _build_id_links none
%define _binary_payload w6.xzdio

Name:           rcv
Version:        %{rcv_version}
Release:        1
Summary:        RNA Conformation Viewer: predict, view and compare RNA conformational ensembles
License:        GPL-3.0-only
URL:            https://www.rnaconfviewer.com/
BuildArch:      x86_64
AutoReqProv:    no
Requires:       glibc >= %{rcv_glibc}
Requires:       xdg-utils

%description
RCV predicts RNA conformational ensembles from a sequence (ViennaRNA secondary structures, 3D assembly and
refinement) and shows them in 3D and 2D: morphs between conformations, plots, palettes, surfaces,
annotations, figure and video export, plugins. Includes the demo ensembles and the user guide.

%install
mkdir -p %{buildroot}
cp -a %{rcv_stage}/. %{buildroot}/

%files
%defattr(-,root,root,-)
/opt/rcv
/usr/bin/rcv
/usr/share/applications/rcv.desktop
/usr/share/icons/hicolor/*/apps/rcv.png
/usr/share/mime/packages/rcv.xml
/usr/share/doc/rcv
