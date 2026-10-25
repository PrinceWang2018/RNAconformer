; Inno Setup script for the RCV Windows installer (built by packaging\build_windows.ps1).
;   Inputs: dist\RCV (PyInstaller one-folder build) and build\docs_html (python packaging\make_docs.py build\docs_html).
#ifndef AppVersion
  #define AppVersion "0.7.0"
#endif

[Setup]
AppId={{A7D1F3C2-5B8E-4F0A-9C61-2E4D7B9A0C35}
AppName=RCV - RNA Conformation Viewer
AppVersion={#AppVersion}
AppVerName=RCV {#AppVersion}
AppPublisher=Zixiang Wang
AppPublisherURL=https://www.rnaconfviewer.com/
AppSupportURL=https://github.com/PrinceWang2018/RNAconformer/issues
AppUpdatesURL=https://www.rnaconfviewer.com/#download
AppContact=princewang2018@gmail.com
AppCopyright=Copyright (C) 2026 Zixiang Wang
VersionInfoVersion={#AppVersion}
VersionInfoDescription=RCV - RNA Conformation Viewer setup
DefaultDirName={autopf}\RCV
DefaultGroupName=RCV
OutputDir=..\..\dist\installer
OutputBaseFilename=RCV-{#AppVersion}-Windows-x64-Setup
SetupIconFile=..\icons\rcv.ico
UninstallDisplayIcon={app}\RCV.exe
UninstallDisplayName=RCV {#AppVersion} - RNA Conformation Viewer
LicenseFile=..\..\LICENSE
Compression=lzma2/ultra64
SolidCompression=yes
LZMAUseSeparateProcess=yes
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
PrivilegesRequiredOverridesAllowed=dialog
WizardStyle=modern
ChangesAssociations=yes
ShowLanguageDialog=auto

; The installer follows the Windows display language when a translation exists (RCV itself has its own
; language setting with more languages, including Chinese).
[Languages]
Name: "en"; MessagesFile: "compiler:Default.isl"
Name: "de"; MessagesFile: "compiler:Languages\German.isl"
Name: "es"; MessagesFile: "compiler:Languages\Spanish.isl"
Name: "fr"; MessagesFile: "compiler:Languages\French.isl"
Name: "it"; MessagesFile: "compiler:Languages\Italian.isl"
Name: "ja"; MessagesFile: "compiler:Languages\Japanese.isl"
Name: "ko"; MessagesFile: "compiler:Languages\Korean.isl"
Name: "ptbr"; MessagesFile: "compiler:Languages\BrazilianPortuguese.isl"
Name: "ru"; MessagesFile: "compiler:Languages\Russian.isl"

[CustomMessages]
en.FileTypes=File types:
en.AssocRcv=Open .rcv files with RCV
en.UserGuide=RCV User Guide
en.Demo=RCV with demo ensembles
en.WhatsNew=What is new in RCV
en.PluginGuide=RCV Plugin Developer Guide
de.FileTypes=Dateitypen:
de.AssocRcv=.rcv-Dateien mit RCV öffnen
de.UserGuide=RCV-Benutzerhandbuch
de.Demo=RCV mit Demo-Ensembles
de.WhatsNew=Neu in RCV
de.PluginGuide=RCV-Plugin-Entwicklerhandbuch
es.FileTypes=Tipos de archivo:
es.AssocRcv=Abrir los archivos .rcv con RCV
es.UserGuide=Guía del usuario de RCV
es.Demo=RCV con conjuntos de demostración
es.WhatsNew=Novedades de RCV
es.PluginGuide=Guía para desarrollar complementos de RCV
fr.FileTypes=Types de fichiers :
fr.AssocRcv=Ouvrir les fichiers .rcv avec RCV
fr.UserGuide=Guide de l'utilisateur de RCV
fr.Demo=RCV avec les ensembles de démonstration
fr.WhatsNew=Nouveautés de RCV
fr.PluginGuide=Guide du développeur d'extensions RCV
it.FileTypes=Tipi di file:
it.AssocRcv=Apri i file .rcv con RCV
it.UserGuide=Guida utente di RCV
it.Demo=RCV con gli ensemble dimostrativi
it.WhatsNew=Novità di RCV
it.PluginGuide=Guida per sviluppatori di plugin RCV
ja.FileTypes=ファイルの種類:
ja.AssocRcv=.rcv ファイルを RCV で開く
ja.UserGuide=RCV ユーザーガイド
ja.Demo=RCV（デモアンサンブル付き）
ja.WhatsNew=RCV の新機能
ja.PluginGuide=RCV プラグイン開発ガイド
ko.FileTypes=파일 형식:
ko.AssocRcv=.rcv 파일을 RCV로 열기
ko.UserGuide=RCV 사용자 가이드
ko.Demo=RCV (데모 앙상블 포함)
ko.WhatsNew=RCV의 새로운 기능
ko.PluginGuide=RCV 플러그인 개발자 가이드
ptbr.FileTypes=Tipos de arquivo:
ptbr.AssocRcv=Abrir arquivos .rcv com o RCV
ptbr.UserGuide=Guia do usuário do RCV
ptbr.Demo=RCV com conjuntos de demonstração
ptbr.WhatsNew=Novidades do RCV
ptbr.PluginGuide=Guia do desenvolvedor de plug-ins do RCV
ru.FileTypes=Типы файлов:
ru.AssocRcv=Открывать файлы .rcv в RCV
ru.UserGuide=Руководство пользователя RCV
ru.Demo=RCV с демонстрационными ансамблями
ru.WhatsNew=Что нового в RCV
ru.PluginGuide=Руководство разработчика плагинов RCV

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"
Name: "assoc"; Description: "{cm:AssocRcv}"; GroupDescription: "{cm:FileTypes}"

[Files]
Source: "..\..\dist\RCV\*"; DestDir: "{app}"; Flags: recursesubdirs createallsubdirs ignoreversion
; user guide, plugin guide and release notes as web pages (also inside RCV: Help menu)
Source: "..\..\build\docs_html\*"; DestDir: "{app}\docs"; Flags: recursesubdirs createallsubdirs ignoreversion
Source: "..\..\LICENSE"; DestDir: "{app}"; DestName: "LICENSE.txt"; Flags: ignoreversion
Source: "..\THIRD-PARTY-NOTICES.txt"; DestDir: "{app}"; Flags: ignoreversion

[InstallDelete]
Type: filesandordirs; Name: "{app}\_internal"
Type: filesandordirs; Name: "{app}\docs"

[Registry]
Root: HKA; Subkey: "Software\Classes\.rcv"; ValueType: string; ValueName: ""; ValueData: "RCV.Ensemble"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\RCV.Ensemble"; ValueType: string; ValueName: ""; ValueData: "RCV conformational ensemble"; Flags: uninsdeletekey; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\RCV.Ensemble\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\RCV.exe,0"; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\RCV.Ensemble\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\RCV.exe"" ""%1"""; Tasks: assoc

[Icons]
Name: "{group}\RCV"; Filename: "{app}\RCV.exe"
Name: "{group}\{cm:Demo}"; Filename: "{app}\RCV.exe"; Parameters: "--demo"
Name: "{group}\{cm:UserGuide}"; Filename: "{app}\docs\index.html"
Name: "{group}\{cm:WhatsNew}"; Filename: "{app}\docs\changelog.html"
Name: "{group}\{cm:PluginGuide}"; Filename: "{app}\docs\plugins.html"
Name: "{group}\{cm:UninstallProgram,RCV}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\RCV"; Filename: "{app}\RCV.exe"; Tasks: desktopicon

[Run]
Filename: "{app}\RCV.exe"; Description: "{cm:LaunchProgram,RCV}"; Flags: nowait postinstall skipifsilent
Filename: "{app}\docs\index.html"; Description: "{cm:UserGuide}"; Flags: shellexec nowait postinstall skipifsilent unchecked
