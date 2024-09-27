library(connectapi)
client <- connect(server = Sys.getenv("CONNECT_SERVER"),
                  api_key = Sys.getenv("CONNECT_API_KEY"))

install.packages("toxEval", repos = "https://rpkg.chs.usgs.gov/prod-cran/latest")

rsconnect::writeManifest(appDir = "./inst/shiny")
bundle <- bundle_dir("./inst/shiny")

content <- client %>% 
  deploy(bundle, name = "toxEval") %>% 
  poll_task()
