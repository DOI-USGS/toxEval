library(connectapi)
client <- connect(server = Sys.getenv("CONNECT_SERVER"),
                  api_key = Sys.getenv("CONNECT_API_KEY"))

rsconnect::writeManifest(appDir = "./inst/shiny")
bundle <- bundle_dir("./inst/shiny")

content <- client %>% 
  deploy(bundle, name = "toxEval") %>% 
  poll_task()
