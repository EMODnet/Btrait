
#MWTLext <- merge(MWTL$density, Taxonomy, by="taxon")
#TaxMWTL <- subset(Taxonomy, subset=taxon %in% unique(MWTLext$taxon))
#Traits_niozext <- extendTrait(Traits_nioz, taxonomy=TaxMWTL[,1:3])


# ----------------------------------------------
# ----------------------------------------------
# SHINY interface for the northsea benthic data
# ----------------------------------------------
# ----------------------------------------------

shinyUI(dashboardPage(
  dashboardHeader(title = "Northsea Benthos MWTL data",
      titleWidth = 400),
  
#===================================
# Sidebar with interactive inputs
#===================================
  dashboardSidebar(disable=TRUE),  # no sidebar

  dashboardBody(
    # Boxes put in two rows
    fluidRow(
      box(status = "primary", title = "map", width=4, 
         plotOutput ("spMap", height="500px"),
        splitLayout(cellWidths = c("30%", "70%"),
        checkboxInput("clog", "log scale", value=FALSE),
        radioButtons(inputId="color", label="Color scheme", inline=TRUE,
             choices=c("blue/yellow/red", "green/blue")))
        ),
      
      box(width=4, status="warning", title = "xy plots", 
         plotOutput ("xyPlot", height="550px")
        ),
      
      box(status = "primary", title = "select", width=4, 
        splitLayout(cellWidths = c("25%", "25%","25%", "25%"),
      div(img(src="Small_tube_dwellers.jpg", height=65, width=65)),
      div(img(src="Shallow_shell.jpg", height=65, width=65)),
      div(img(src="Echinocardium_zeeklit.png", height=65, width=65)),
      div(img(src="Bathyporeia_small.jpg", height=65, width=65))), 
      wellPanel(  
         
       radioButtons(inputId="settaxon", label="What",
           choices=c("taxa", "abiotics", "traits", "groups", "summary"), selected="taxa", inline=TRUE),
       uiOutput("SettaxonList"),  # list of taxa generated in server
       conditionalPanel(condition="input.settaxon == 'taxa'",
         checkboxInput("all", "list all taxa", value=FALSE)),
  
       conditionalPanel(condition="input.settaxon != 'abiotics'",
        radioButtons(inputId="what", label="Variable",
           choices=c("density", "biomass"), selected="density", inline=TRUE)),
       uiOutput("SetabioticsList"),  # list of abiotic data generated in server
       radioButtons(inputId="logaxes", label="Logarithmic axes",
           choices=c("none", "x-axis", "y-axis", "both"),
           selected="none", inline=TRUE) 
      )
    )
    ),
    fluidRow(
      box(status = "primary", title = "Time series", width=5, 
#      wellPanel(  
       plotOutput ("tsPlot", height="300px"))
#      )
    ,
      box(status="warning", title = "Info", width=5, 
         verbatimTextOutput ("printTraits"),
         verbatimTextOutput ("printTaxonomy")
        ),
      box(status = "info", title = "Acknowledgements", width=2, 
      wellPanel(  
       box(title="Created by:", width=NULL, background="light-blue",
           "Karline Soetaert, 2022", 
            h4("illustrations:"), h5("Anna van der Kaaden")),
        splitLayout(cellWidths = c("40%", "40%", "20%"),
            div(img(src="Rlogo.png", width=60)),
            div(img(src="EMODnet.png", width=60)),
          div(img(src="RWSlogo.png", width=40)))
        )
      )
    )
  )
))
    