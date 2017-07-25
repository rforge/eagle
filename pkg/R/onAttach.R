.onAttach <- function(libname, pkgname) {
  # Runs when attached to search() path such as by library() or require()
  if (interactive()) {


##-------------------------------
## Opening message
##--------------------------------
packageStartupMessage(' ')
packageStartupMessage('   +---------------------------------------------------+ ')
packageStartupMessage('   |              Eagle has been loaded.               | ')               
packageStartupMessage('   |  Documentation: vignette("QuickStart", "Eagle")   | ')
packageStartupMessage('   |  Run OpenUI() to open web-based user interface.   | ')
packageStartupMessage('   +---------------------------------------------------+ ')
packageStartupMessage(' ')
}

}

