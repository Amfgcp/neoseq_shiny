Rexe           = "R-Portable\App\R-Portable\bin\Rscript.exe"
Ropts          = "--no-save --no-environ --no-init-file --no-restore --no-Rconsole "
RScriptFile    = "runShinyApp.R"
Outfile        = "ShinyApp.log"
startChrome    = "GoogleChromePortable\GoogleChromePortable.exe --app=http://127.0.0.1:7777"
strCommand     = Rexe & " " & Ropts & " " & RScriptFile '& " 1> " & Outfile & " 2>&1"'

intWindowStyle = 0     ' Hide the window and activate another window.'
bWaitOnReturn  = False ' continue running script after launching R   '

' Terminate any currently running R scripts. There can be only one!
'CreateObject("WScript.Shell").Run "taskkill /im Rscript.exe", intWindowStyle, True'

' the following is a Sub call, so no parentheses around arguments'
CreateObject("Wscript.Shell").Run strCommand, intWindowStyle, bWaitOnReturn
CreateObject("Wscript.Shell").Run startChrome, intWindowStyle, bWaitOnReturn