; Script generated by the Inno Setup Script Wizard.
; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!

; Set this path to where you've gathered the Deployable's files.
#define workdirectory "H:\NGSE"

#define MyAppName "My Program"
#define InstallerName "distributable Shiny app setup"
#define MyAppVersion "1.0"
; #define MyAppPublisher "My Company, Inc."
; #define MyAppURL "http://www.example.com/"
#define MyAppExeName "run.vbs"

[Setup]
; NOTE: The value of AppId uniquely identifies this application.
; Do not use the same AppId value in installers for other applications.
; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)
AppId={{CC150627-0A91-47D4-AD78-500C1E3FB13D}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
; AppVerName={#MyAppName} {#MyAppVersion}
; AppPublisher={#MyAppPublisher}
; AppPublisherURL={#MyAppURL}
; AppSupportURL={#MyAppURL}
; AppUpdatesURL={#MyAppURL}
PrivilegesRequired=lowest
UsePreviousAppDir=no
DefaultDirName={localappdata}\{#MyAppName}
DisableProgramGroupPage=yes
OutputDir={#workdirectory}
OutputBaseFilename={#InstallerName}
Compression=lzma
SolidCompression=yes

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"
; Name: "dutch"; MessagesFile: "compiler:Languages\Dutch.isl"

; [Tasks]
; Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: {#workdirectory}\*; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs; Excludes: {#InstallerName}
;Source: "C:\Users\pathreg\Desktop\sample dist complete\run.vbs"; DestDir: "{app}"; Flags: ignoreversion
; Source: "C:\Users\pathreg\Desktop\sample dist complete\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
; NOTE: Don't use "Flags: ignoreversion" on any shared system files

[Icons]
Name: "{app}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"
Name: "{userdesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"
; Tasks: desktopicon
; NOTE: To use a custom file icon, use  IconFilename: "{app}\[icon].ico"; (where [icon] is the name of your custom icon).