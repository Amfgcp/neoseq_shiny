require(shiny)
#The app found in the /shiny/ folder will be launched on this port.
#The .vbs script launches a ChromePortable browser with this port for the App to run in.
shiny::runApp('./shiny/', port=7777)