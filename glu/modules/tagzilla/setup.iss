[Setup]
AppName=TagZilla
AppVerName=TagZilla version 1.0
DefaultDirName = {pf}\TagZilla
DefaultGroupName = TagZilla
OutputBaseFilename = TagZilla100_Win32_Installer

[Files]
Source: "dist\*"; DestDir: {app}

[Registry]
Root: HKLM; Subkey: "SOFTWARE\Microsoft\Windows\CurrentVersion\App Paths\TagZilla.exe"; ValueType: string; ValueName: ""; ValueData: "{app}\TagZilla.exe"; Flags: uninsdeletekey

[Icons]
Name: "{group}\TagZilla"; Filename: "{app}\TagZilla.exe"
Name: "{group}\UnInstall"; Filename: "{app}\{uninstallexe}"
