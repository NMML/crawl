#' @title Start a shiny app to check data stored in a .csv file for model fitting with \code{crwMLE} function.
#' @export
#' @importFrom shiny runApp
check_csv = function(){
  appDir <- system.file("shiny_apps", "check_project_data", package = "crawl")
  if (appDir == "") stop("Could not find shiny app directory. Try re-installing `crawl`.", call. = FALSE)
  shiny::runApp(appDir, display.mode = "normal")
}