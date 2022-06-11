# Se cargan las librerias
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(igraph)
library(tidyverse)
library(tidyr)
library(CINNA)
library(bibliometrix)
library(rebus)
library(tm)
library(wordcloud)
library(shinyWidgets)
library(tidygraph)
source('verbs2.R')


# ------ CREACION DE  FUNCIONES AUXILIARES PARA ANALISIS DE CITACIONES---------

#-------------------------------------------------------------------------------
# ------ Se crea la funcion para generar los wordclouds-------------------------
#-------------------------------------------------------------------------------

wordclouds <- function (subarea_1){

  # Se crea el corpus
  corp <- Corpus(VectorSource(subarea_1$TI %>% na.omit()))
  paperCorp <- corp

  # Se realiza la limpieza el corpus
  paperCorp <- tm_map(paperCorp, removePunctuation)
  paperCorp <- tm_map(paperCorp, removeNumbers)
  paperCorp <- tm_map(paperCorp, content_transformer(tolower))
  paperCorp <- tm_map(paperCorp, removeWords, stopwords("english"))
  paperCorp <- tm_map(paperCorp, stripWhitespace)
  #paperCorp <- tm_map(paperCorp, stemDocument)
  paperCorp <- tm_map(paperCorp, removeWords, c("viral", "market"))
}


get_data <- function(scopus_dataframe){

  references_df        <- get_references(scopus_dataframe)
  citation_network     <- get_citation_network(scopus_dataframe,
                                               references_df = references_df)
  citation_network_tos <- get_citation_network_tos(citation_network)
  return(citation_network_tos)
}

# ------------------ CREACION DE LA INTERFAZ DE USUARIO ------------------------


title1 <- "Core of science"

ui <- fluidPage(

  dashboardPage(skin = "green",
                dashboardHeader(title = title1),

                dashboardSidebar(
                  fileInput("file1", "Choose .bib File",
                            accept = c(
                              "text/csv",
                              "text/comma-separated-values,text/plain",
                              ".bib")
                  ),


                  sidebarMenu(
                    menuItem("Introduction", tabName = "introduction", icon = icon("th")),
                    menuItem("Importance"  , tabName = "importance"  , icon = icon("th")),
                    menuItem("Evolution - ToS"   , tabName = "evolution"   , icon = icon("th")),
                    menuItem("Subfields"   , tabName = "subfields"   , icon = icon("th"))
                  )
                ),
                dashboardBody(

                  tabItems(
                    tabItem(tabName = "introduction",h2("Introduction"),
                            box(width = 8,h3("Tree of Science ToS"),
                                "Tree of Science (ToS) is a Web based tool for scientific articles selection. ToS was created in order to solve the difficulties to find relevant articles in a research topics and to make easier the process on writing the theoretical framework. ToS has three advantages: decreases the time interval bias in the search, decreases bias of the databases indexed in the search, diminishes the rigor of keywords. ToS is directed to all academic community and researchers, and students completing a short-term research project.
                                 ToS philosophy has been based on three pillars: simplicity, effectiveness and innovation. Simplicity is based on organic concepts to help the user understanding about the structure of science from the tree metaphor. The effectiveness is based on the accuracy of the results of scientific articles. Finally, innovation part of a continuous improvement of the services provided so ToS can surprise the users." ),
                            #actionButton("do", "Realizar calculos")
                            ),
                    tabItem(tabName = "importance",h2("Articles importance"),

                            fluidRow( box(title = "History publication",
                                          collapsible = TRUE, closable = F,
                                          downloadButton("downloadData2", strong("Download History publication")),
                                          actionButton(inputId='ab1', label="Learn More",
                                                       icon = icon("th"),
                                                       onclick ="window.open('https://www.sciencedirect.com/science/article/abs/pii/S1751157717300500?via%3Dihub')"),
                                          plotOutput("grafico1")),
                                      box(title = "Most productive authors",
                                          collapsible = TRUE, closable = F,
                                          downloadButton("downloadData3", strong("Download Most productive authors")),
                                          tableOutput("tabla1"))),

                            fluidRow(box(title = "Most popular journals",
                                         collapsible = TRUE, closable = F,
                                         downloadButton("downloadData4", strong("Download Most popular journals")),
                                         tableOutput("tabla2")))

                    ),
                    tabItem(tabName = "evolution",h2("Evolution ToS"),
                            fluidRow(
                              box(title = "Evolution ToS",collapsible = TRUE, closable = F,
                                  selectInput("tosid",
                                              label = "Choose a TOS category",
                                              choices = list("All",
                                                             "Root",
                                                             "Trunk",
                                                             "Leaves"),
                                              selected = "All"),
                                  box(title = "Tree of Science ToS",
                                      DT::dataTableOutput("contents")
                                  )
                              ),

                              box(title = "Download ToS",collapsible = TRUE, closable = F,
                                  downloadButton("downloadData", strong("Download Tree of science")),
                                  actionButton(inputId='ab2', label="Learn More",
                                               icon = icon("th"),
                                               onclick ="window.open('https://revistas.unal.edu.co/index.php/ingeinv/article/view/77718/74279')")
                              )

                            )
                    ),
                    tabItem(tabName = "subfields",h2("Subfields"),
                            fluidRow(
                              box(title = "Evolution ToS",collapsible = TRUE, width = 7,closable = F,
                                  selectInput("grupo",
                                              label = "Choose a subfield group",
                                              choices = list("Group 1",
                                                             "Group 2",
                                                             "Group 3"),
                                              selected = "Group 1"),
                                  downloadButton("downloadData1", strong("Download Subfields")),
                                  actionButton(inputId='ab3', label="Learn More",
                                               icon = icon("th"),
                                               onclick ="window.open('https://revistas.usb.edu.co/index.php/Psychologia/article/download/4230/3520/')"),
                                  box(title = "Subfields clasification",
                                      DT::dataTableOutput("subareas1"))
                              ),
                              box(title = "Worcloud",collapsible = TRUE, width = 4, closable = F,
                                  sliderInput("freq",
                                              "Minimum Frequency:",
                                              min = 10,  max = 100, value = 50),
                                  sliderInput("max",
                                              "Maximum Number of Words:",
                                              min = 1,  max = 300,  value = 50),
                                  plotOutput("wordcloud1"),
                                  selectInput("dinamyc", "Select words to remove",
                                              choices  = NULL,
                                              selected = NULL,
                                              multiple = TRUE),
                                  actionButton("reset", "Reset"),
                                  downloadButton("down", strong("Download wordcloud"))
                              ),

                            )
                    )

                  )
                )
  )
)


# ------------------ CREACION DE LA FUNCION PRINCIPAL  ------------------------

server <- function(input, output,session) {



  # Se carga el archivo
  data <- reactive({
    req(input$file1)
    scopus_dataframe <- convert2df(file = input$file1$datapath,
                                   dbsource = "scopus",
                                   format   = "bibtex")
    citation_network_tos <- get_data(scopus_dataframe)})



  citation_tos1 <- reactive({
    TOS <- get_sap(data())
    return(TOS)
  })

  observeEvent(input$file1, {
    showNotification("Computing!")
    withProgress(message = 'Processing data', value = 100,{
    data()
    citation_tos1()
    subareas.t()
    })
  })

  output$tabla11 <- renderTable({
    as.tibble(data())
  })

  # Se crea el arbol TOS
  output$contents <- DT::renderDataTable({
    TOS <- citation_tos1()
    TOS$link <- paste0("https://www.google.com/search?q=",TOS$Title)
    TOS$link <- paste0("<a href='",TOS$link,"'>",TOS$link,"</a>")

    if (input$tosid == "All"){
      return(DT::datatable(TOS,escape = FALSE))
    }
    if (input$tosid == "Root"){
      return(DT::datatable(TOS[TOS$TOS == "Root",],escape = FALSE))
    }
    if (input$tosid == "Leaves"){
      return(DT::datatable(TOS[TOS$TOS == "Leaves",],escape = FALSE))
    }
    if (input$tosid == "Trunk"){
      return(DT::datatable(TOS[TOS$TOS == "Trunk",],escape = FALSE))
    }
  })


  output$downloadData <- downloadHandler(
    filename = function() {
      paste("TreeOfScience", ".csv", sep = "")
    },
    content = function(file) {
      #write.csv(TOS, file, row.names = FALSE)
      write.table(citation_tos1(),file,sep=";")
    }
  )

  # Se crea el grafico de produccion cientifica por año
  output$grafico1 <- renderPlot({
    scopus_dataframe  <- as.tibble(data())
    fecha.publicacion <- scopus_dataframe$PY         # Año de publicacion
    N                 <- length(fecha.publicacion)   # Numero de articulos
    frec.publicacion  <- table(fecha.publicacion)
    plot(x    = frec.publicacion,  # Datos
         main = "Scientific anual production",
         ylab = "Papers",
         xlab = "Year",
         col  = "Red",
         type = "o")
  })

  output$tabla1 <- renderTable({
    scopus_dataframe  <- as.tibble(data()) %>% separate(col = 'name', into = c("AU", "year", "SO"), sep = ',')
    autores <- scopus_dataframe$AU
    Author  <- c()

    for (w in autores) {
      for (w1 in w){
        Author <- c(Author,w1)
      }
    }
    frec.autores <- as.data.frame(table(Author)) %>%
      arrange(desc(Freq))


    return(frec.autores[1:10,])

  })

  output$tabla2 <- renderTable({

    scopus_dataframe  <- as.tibble(data()) %>% separate(col = 'name', into = c("AU", "year", "SO"), sep = ',')
    Journals <- as.data.frame(table(scopus_dataframe$SO)) %>%
      arrange(desc(Freq)) %>%
      rename(Journal = "Var1")

    return(Journals[1:10,])

  })

  # Determinacion de las subareas

  subareas.t <- reactive({
    subareas <- extract_subareas(data())
  })

  output$subareas1 <- DT::renderDataTable({
    subareas      <- subareas.t()
    subareai      <- subareas[subareas$group == input$grupo, ]
    subareai$link <- paste0("https://www.google.com/search?q=",subareai$TI)
    subareai$link <- paste0("<a href='",subareai$link,"'>",subareai$link,"</a>")

    return(DT::datatable(subareai[,-c(5)],escape = FALSE))
  })

  wordcloud.con <- reactive({
    subareas    <- subareas.t()
    wordcloud1  <- wordclouds(subareas[subareas$group == input$grupo,])
    dtm <- TermDocumentMatrix(wordcloud1)
    m   <- as.matrix(dtm)
    v   <- sort(rowSums(m),decreasing=TRUE)
    d   <- data.frame(word = names(v),freq=v)
  })

  output$wordcloud1 <- renderPlot({

    wordcloud1 <- wordcloud.con()

    if (!is.null(input$dinamyc)){
      wordcloud1 <- wordcloud1[-c(which(wordcloud1$word %in% input$dinamyc)),]
    }


    wordcloud(words = wordcloud1$word,
              freq  = wordcloud1$freq,
              min.freq  = input$freq,
              max.words = input$max,
              random.order = FALSE,
              rot.per=0.35,
              colors=brewer.pal(8,"Dark2"))
  })

  observeEvent(input$reset, {
    wordcloud1 <- wordcloud.con()
    palabras   <- wordcloud1$word[1:10]
    updateSelectInput(session,"dinamyc",choices = palabras)
  })

  output$downloadData1 <- downloadHandler(
    filename = function() {
      paste("subfields", ".csv", sep = "")
    },
    content = function(file) {
      #write.csv(TOS, file, row.names = FALSE)
      write.table(subareas.t(),file,sep=";")
    }
  )

  output$down <- downloadHandler(
    filename =  function() {
      paste("wordcloud", ".png", sep="")
    },

    content = function(file) {
      png(file)
      wordcloud1 <- wordcloud.con()

      if (!is.null(input$dinamyc)){
        wordcloud1 <- wordcloud1[-c(which(wordcloud1$word %in% input$dinamyc)),]
      }

      print( wordcloud(words = wordcloud1$word,
                       freq  = wordcloud1$freq,
                       min.freq  = input$freq,
                       max.words = input$max,
                       random.order = FALSE,
                       rot.per=0.35,
                       colors=brewer.pal(8,"Dark2")))
      dev.off()
    }
  )

  output$downloadData2 <- downloadHandler(
    filename =  function() {
      paste("historyPublication", ".png", sep="")
    },

    content = function(file) {
      png(file)
      scopus_dataframe  <- data()
      fecha.publicacion <- scopus_dataframe$PY         # Año de publicacion
      N                 <- length(fecha.publicacion)   # Numero de articulos
      frec.publicacion  <- table(fecha.publicacion)
      plot(x    = frec.publicacion,  # Datos
           main = "Scientific anual production",
           ylab = "Papers",
           xlab = "Year",
           col  = "Red",
           type = "o")
      dev.off()
    }
  )

  output$downloadData3 <- downloadHandler(
    filename = function() {
      paste("Most_productive_authors", ".csv", sep = "")
    },
    content = function(file) {
      scopus_dataframe  <- data()
      autores <- strsplit( x = scopus_dataframe$AU, split = ";", fixed = FALSE )
      Author  <- c()

      for (w in autores) {
        for (w1 in w){
          Author <- c(Author,w1)
        }
      }
      frec.autores <- as.data.frame(table(Author)) %>%
        arrange(desc(Freq))

      write.table(frec.autores,file,sep=";")
    }
  )

  output$downloadData4 <- downloadHandler(
    filename = function() {
      paste("Most_popular_journals", ".csv", sep = "")
    },
    content = function(file) {
      scopus_dataframe  <- data()
      Journals <- as.data.frame(table(scopus_dataframe$SO)) %>%
        arrange(desc(Freq)) %>%
        rename(Journal = "Var1")

      write.table(Journals,file,sep=";")
    }
  )

}

# ------------------ CREACION DE LA APLICACION--------  ------------------------

shinyApp(ui, server)
