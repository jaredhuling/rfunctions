
  
KnitPost <- function(input, base.url = "/") {
  #stolen from here:
  #http://jfisher-usgs.github.io/r/2012/07/03/knitr-jekyll/
  #used just like knit() function
  require(knitr)
  opts_knit$set(base.url = base.url)
  fig.path <- paste0("figs/", sub(".Rmd$", "", basename(input)), "/")
  opts_chunk$set(fig.path = fig.path)
  opts_chunk$set(fig.cap = "center")
  render_jekyll()
  knit(input, envir = parent.frame())
}