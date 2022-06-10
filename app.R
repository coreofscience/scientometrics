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
library(shinythemes)
library(tidygraph)


# ------ CREACION DE  FUNCIONES AUXILIARES PARA ANALISIS DE CITACIONES---------

#-------------------------------------------------------------------------------
# ------ Se crea la funcion para la generacion de la edgelist ------------------
#-------------------------------------------------------------------------------
edgelista.bib <- function(scopus_dataframe){

  pattern_authors <-
    SPC %R%
    one_or_more(WRD) %R%
    SPC %R%
    one_or_more(or(WRD, ANY_CHAR))

  pattern_titles <-
    OPEN_PAREN %R%
    repeated(DGT, 4) %R%
    CLOSE_PAREN %R%
    one_or_more(or(WRD,ANY_CHAR))

  pattern_year <-
    OPEN_PAREN %R%
    repeated(DGT, 4) %R%
    CLOSE_PAREN

  pattern_journal <-
    one_or_more(or(WRD,SPC))

  pattern_volume <-
    one_or_more(or(WRD, SPC))

  pattern_pages <-
    "PP. " %R%
    one_or_more(or(DGT, ANY_CHAR))

  cited_references <-
    scopus_dataframe %>%
    separate_rows(CR, sep = "; ") %>%
    select(SR_TOS,
           CR) %>%
    mutate(CR_AUTHOR = str_remove(CR, pattern_authors),
           CR_TITLE_1 = str_extract(CR, pattern_authors),
           CR_TITLE = str_remove(CR_TITLE_1, pattern_titles),
           CR_TITLE = str_trim(CR_TITLE),
           CR_YEAR_1 <- str_extract(CR_TITLE_1, pattern_titles),
           CR_YEAR = str_extract(CR_YEAR_1, repeated(DGT, 4)),
           CR_JOURNAL_1 = str_remove(CR_YEAR_1, pattern_year),
           CR_JOURNAL = str_extract(CR_JOURNAL_1, pattern_journal),
           CR_JOURNAL = str_trim(CR_JOURNAL),
           CR_VOLUME_1 = str_remove(CR_JOURNAL_1, pattern_journal),
           CR_VOLUME = str_extract(CR_VOLUME_1, pattern_volume),
           CR_PAGES = str_extract(CR_VOLUME_1, pattern_pages),
           CR_PAGES = str_remove(CR_PAGES, "PP. ")) %>%
    select(SR_TOS,
           CR,
           CR_AUTHOR,
           CR_TITLE,
           CR_YEAR,
           CR_JOURNAL,
           CR_VOLUME,
           CR_PAGES) %>%
    mutate(lastname = sub("\\., .*", "", CR),
           lastname = sub(",", "", lastname),
           lastname = sub("\\.", "", lastname),
           CR_SO = str_c(lastname,
                         ", ",
                         CR_YEAR,
                         ", ",
                         CR_JOURNAL)) %>%
    select(-lastname)


  edge_list <-
    cited_references %>%
    select(SR_TOS,
           CR_SO) %>%
    na.omit() %>%
    unique()

  return(edge_list)

}


#-------------------------------------------------------------------------------
# ------ Se crea la funcion para depurar el grafo ------------------------------
#-------------------------------------------------------------------------------

crear.grafo <- function(edgelist){

  graph <- graph.data.frame(edgelist) %>%
    simplify()

  # Se eliminan los vertices con indegree = 1 y con outdegree = 0
  graph_1 <- delete.vertices(graph,
                             which(degree(graph, mode = "in") == 1 &
                                     degree(graph, mode = "out") == 0))

  # Se escoge el componente mas grande conectado
  graph_2 <- giant_component_extract(graph_1, directed = TRUE)
  graph_2 <- graph_2[[1]]

  subareas <-
    as.undirected(graph_2,
                  mode = "each") %>%
    cluster_louvain()

  graph_2 <-
    graph_2 %>%
    set_vertex_attr(name = "sub_area",
                    value = membership(subareas))


}

#-------------------------------------------------------------------------------
# ------ Se crea la funcion algoritmo SAP para analisis de citaciones ----------
#-------------------------------------------------------------------------------
algoritmoSAP <- function(graph_2){

  # Se crean la metricas de la red

  metricas.red <- tibble(
    id        = V(graph_2)$name,
    indegree  = degree(graph_2, mode = "in"),
    outdegree = degree(graph_2, mode = "out"),
    bet       = betweenness(graph_2))

  metricas.red <- metricas.red %>%
    mutate(year = as.numeric(str_extract(id, "[0-9]{4}")))

  # Clasificacion de las raices

  Raices <- metricas.red[metricas.red$outdegree == 0, c("id","indegree")] %>%
    arrange(desc(indegree))
  Raices <- Raices[1:10,]


  # Clasificacion de las hojas
  Hojas.ext <- metricas.red[metricas.red$indegree == 0, c("id","outdegree","year")]
  act.year  <- as.numeric(format(Sys.Date(),'%Y'))
  Hojas.ext <- Hojas.ext %>%
    mutate(antiguedad = act.year - year) %>%
    arrange(antiguedad)
  Hojas     <- filter(Hojas.ext, antiguedad <= 5)

  # Se determina el numero del vertice de las Hojas
  num.vertices.hojas <- c()
  for (vertice in Hojas$id){
    num.vertices.hojas <- c(num.vertices.hojas,which(metricas.red$id == vertice))
  }

  # Se determina el numero del vertice de las raices
  num.vertices.raices <- c()
  for (vertice in Raices$id){
    num.vertices.raices <- c(num.vertices.raices,which(metricas.red$id == vertice))
  }

  # Calculo del SAP de las Hojas
  SAP_hojas <- c()
  for (vert in Hojas$id){
    h <- get.all.shortest.paths(graph_2,
                                from = vert,
                                to   = Raices$id,
                                mode = "out")

    SAP_hojas   <- c(SAP_hojas, length(h[[1]]))
  }
  Hojas <- Hojas %>%
    mutate(SAP = SAP_hojas) %>%
    arrange(desc(SAP))

  Hojas <- Hojas[1:60,] %>%
    filter(SAP > 0)

  Caminos   <- c()
  for (vert in Hojas$id){
    h <- get.all.shortest.paths(graph_2,
                                from = vert,
                                to   = Raices$id,
                                mode = "out")
    lista.nodos <- unique(unlist(h[1]))
    lista.nodos <- lista.nodos[!(lista.nodos %in% num.vertices.raices)]
    lista.nodos <- lista.nodos[!(lista.nodos %in% num.vertices.hojas)]
    Caminos     <- c(Caminos,lista.nodos)
  }



  # Seleccion del tronco

  Tronco     <- metricas.red[unique(Caminos), c("id","indegree","year")]
  mas.nuevo  <- max(Tronco$year, na.rm = TRUE)
  Tronco     <- Tronco %>%
    mutate(antiguedad = mas.nuevo - year)

  # Tree of science
  Raices$TOS <- "Root"
  Hojas$TOS  <- "Leaves"
  Tronco$TOS <- "Trunk"

  TOS   <- rbind(Raices[,c(1,3)], Tronco[,c(1,5)], Hojas[,c(1,6)])
  return(TOS)
}


#-------------------------------------------------------------------------------
# ------ Se crea la funcion para generar subareas ------------------------------
#-------------------------------------------------------------------------------
subarea.grafo <- function(citation_tos){

  citation_tos <- as_tibble(citation_tos)
  subfields <- citation_tos %>% group_by(subfield) %>% count() %>%
    arrange(desc(n)) %>%
    head(3) %>%
    select(subfield)

  group1             <- citation_tos[citation_tos$subfield == subfields$subfield[1],] %>%
    mutate(group = 'Group 1')
  group2             <- citation_tos[citation_tos$subfield == subfields$subfield[2],] %>%
    mutate(group = 'Group 2')

  group3             <- citation_tos[citation_tos$subfield == subfields$subfield[3],] %>%
    mutate(group = 'Group 3')

  referencias.grupo  <- rbind(group1,group2, group3)

  return(referencias.grupo)
}

#-------------------------------------------------------------------------------
# ------ Se crea la funcion para generar referencias citadas -------------------
#-------------------------------------------------------------------------------

cited_references <- function(scopus_dataframe) {

  pattern_authors <-
    SPC %R%
    one_or_more(WRD) %R%
    SPC %R%
    one_or_more(or(WRD, ANY_CHAR))

  pattern_titles <-
    OPEN_PAREN %R%
    repeated(DGT, 4) %R%
    CLOSE_PAREN %R%
    one_or_more(or(WRD,ANY_CHAR))

  pattern_year <-
    OPEN_PAREN %R%
    repeated(DGT, 4) %R%
    CLOSE_PAREN

  pattern_journal <-
    one_or_more(or(WRD,SPC))

  pattern_volume <-
    one_or_more(or(WRD, SPC))

  pattern_pages <-
    "PP. " %R%
    one_or_more(or(DGT, ANY_CHAR))

  cited_references <-
    scopus_dataframe %>%
    separate_rows(CR, sep = "; ") %>%
    select(SR_TOS,
           CR) %>%
    mutate(CR_AUTHOR = str_remove(CR, pattern_authors),
           CR_TITLE_1 = str_extract(CR, pattern_authors),
           CR_TITLE = str_remove(CR_TITLE_1, pattern_titles),
           CR_TITLE = str_trim(CR_TITLE),
           CR_YEAR_1 <- str_extract(CR_TITLE_1, pattern_titles),
           CR_YEAR = str_extract(CR_YEAR_1, repeated(DGT, 4)),
           CR_JOURNAL_1 = str_remove(CR_YEAR_1, pattern_year),
           CR_JOURNAL = str_extract(CR_JOURNAL_1, pattern_journal),
           CR_JOURNAL = str_trim(CR_JOURNAL),
           CR_VOLUME_1 = str_remove(CR_JOURNAL_1, pattern_journal),
           CR_VOLUME = str_extract(CR_VOLUME_1, pattern_volume),
           CR_PAGES = str_extract(CR_VOLUME_1, pattern_pages),
           CR_PAGES = str_remove(CR_PAGES, "PP. ")) %>%
    select(SR_TOS,
           CR,
           CR_AUTHOR,
           CR_TITLE,
           CR_YEAR,
           CR_JOURNAL,
           CR_VOLUME,
           CR_PAGES) %>%
    mutate(lastname = sub("\\., .*", "", CR),
           lastname = sub(",", "", lastname),
           lastname = sub("\\.", "", lastname),
           CR_SO = str_c(lastname,
                         ", ",
                         CR_YEAR,
                         ", ",
                         CR_JOURNAL)) %>%
    select(-lastname)

  return(cited_references)

}

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



# ------------------ CREACION DE LA INTERFAZ DE USUARIO ------------------------


title1 <- "Core of Science"

ui <- fluidPage(

  dashboardPage(skin = "green",
                dashboardHeader(title = title1),

                dashboardSidebar(
                  fileInput("file1", "Choose .bib File",
                            accept = c(
                              "text/csv",
                              "text/comma-separated-values,text/plain",
                              ".bib")
                  ), actionButton("go", "Process your data"),


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
                                         "Leaves",
                                         "Trunk"),
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
                          min = 10,  max = 100, value = 20),
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




shiny_busy <- function() {
  HTML("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;", paste0(
    '<span data-display-if="',
    '$(&#39;html&#39;).attr(&#39;class&#39;)==&#39;shiny-busy&#39;',
    '">',
    '<i class="fa fa-spinner fa-pulse fa-fw" style="color:orange"></i>',
    '</span>'
  ))
}

get_data <- function(scopus_data){

  scopus_dataframe <- convert2df(file = scopus_data,
                                 dbsource = "scopus",
                                 format   = "bibtex")
  references_df        <- get_references(scopus_dataframe)
  citation_network     <- get_citation_network(scopus_dataframe,
                                               references_df = references_df)
  citation_network_tos <- get_citation_network_tos(citation_network)
  tos                  <- get_sap(citation_network_tos)

  return(list(scopus_dataframe = scopus_dataframe,
              references_df = references_df,
              citation_network_tos = citation_network_tos,
              tos = tos))
}
# ------------------ CREACION DE LA FUNCION PRINCIPAL  ------------------------

server <- function(input, output, session) {

  # Se carga el archivo
  data <- reactive({
    req(input$file1)

    showModal(modalDialog("Processing data...wait a minute",
                          shiny_busy(), easyClose = FALSE))

    # scopus_dataframe <- convert2df(file = input$file1$datapath,
    #                                dbsource = "scopus",
    #                                format   = "bibtex")

    new_data <- get_data(input$file1$datapath)

    if(is.null(new_data)) {
      showModal(modalDialog("there was a problem"))
    }else{
      showModal(modalDialog("Your Tree of Science is done"))
    }
    # TOS <- get_sap(new_data)
    return(new_data)

  })

  observeEvent(input$go, {
    data()
  })

  ## ------------------ Evolution ToS ------------------------

  # citation_tos <- reactive({
  #   scopus_dataframe     <- new_data$tos
  #   new_data <- get_data(scopus_dataframe)
  #   return(new_data)
  # })
  #
  # citation_tos1 <- reactive({
  #   scopus_dataframe     <- data()
  #   new_data <- get_data(scopus_dataframe)
  #   TOS <- get_sap(new_data)
  #   return(TOS)
  # })

  output$contents <- DT::renderDataTable({

    TOS <- new_data$tos |>
      select(-id, -name)

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

  # citation_tos <- reactive({
  #   scopus_dataframe     <- data()
  #   references_df        <- get_references(scopus_dataframe)
  #   citation_network     <- get_citation_network(scopus_dataframe, references_df = references_df)
  #   citation_network_tos <- get_citation_network_tos(citation_network)
  # })
  #
  # datatos <- reactive({
  #   TOS                  <- get_sap(citation_tos())
  # })

  ## ------------------ Download ToS  ------------------------

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("TreeOfScience", ".csv", sep = "")
    },
    content = function(file) {
      #write.csv(TOS, file, row.names = FALSE)
      write.table(new_data$tos,file,sep=";")
    }
  )

  ## ------------------ Importance  ------------------------

  # Se crea el grafico de produccion cientifica por año
  output$grafico1 <- renderPlot({
    scopus_dataframe  <- new_data$scopus_dataframe
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
    scopus_dataframe  <- new_data$scopus_dataframe
    autores <- strsplit( x = scopus_dataframe$AU, split = ";", fixed = FALSE )
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

    scopus_dataframe  <- new_data$scopus_dataframe
    Journals <- as.data.frame(table(scopus_dataframe$SO)) %>%
      arrange(desc(Freq)) %>%
      rename(Journal = "Var1")

    return(Journals[1:10,])

  })
  #
  # # Determinacion de las subareas
  #
  # # subareas.t <- reactive({
  # #   subareas <- subarea.grafo(citation_tos())
  # # })

  ## ------------------ Subfields  ------------------------

  output$subareas1 <- DT::renderDataTable({
    subareas      <- new_data$citation_network_tos |>
      activate(nodes) |>
      as_tibble() |>
      select(group = subfield, TI, PY) |>
      mutate(group = str_c("Group", group, sep = " "))

    subareai      <- subareas[subareas$group == input$grupo, ]
    #subareai$link <- paste0("https://www.google.com/search?q=",subareai$CR_TITLE)
    #subareai$link <- paste0("<a href='",subareai$link,"'>",subareai$link,"</a>")

    return(DT::datatable(subareai,escape = FALSE))
  })

  ## ------------------ Wordcloud  ------------------------

  wordcloud.con <- reactive({
    subareas      <- new_data$citation_network_tos |>
      activate(nodes) |>
      as_tibble() |>
      select(group = subfield, TI, PY) |>
      mutate(group = str_c("Group", group, sep = " "))

    wordcloud1  <- wordclouds(subareas[subareas$group == input$grupo,])
    dtm <- TermDocumentMatrix(wordcloud1)
    m   <- as.matrix(dtm)
    v   <- sort(rowSums(m),decreasing=TRUE)
    d   <- data.frame(word = names(v),freq=v)
  })
  #
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
  #
  observeEvent(input$reset, {
    wordcloud1 <- wordcloud.con()
    palabras   <- wordcloud1$word[1:10]
    updateSelectInput(session,"dinamyc",choices = palabras)
  })


  ## ------------------ Download Subfields  ------------------------

  output$downloadData1 <- downloadHandler(
    filename = function() {
      paste("subfields", ".csv", sep = "")
    },
    content = function(file) {
      #write.csv(TOS, file, row.names = FALSE)
      write.table(subareas      <- new_data$citation_network_tos |>
                    activate(nodes) |>
                    as_tibble() |>
                    select(group = subfield, TI, PY) |>
                    mutate(group = str_c("Group", group, sep = " ")),
                  file,
                  sep=";")
    }
  )
  #
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
  ## ------------------ Download importance data ------------------

  output$downloadData2 <- downloadHandler(
    filename =  function() {
      paste("historyPublication", ".png", sep="")
    },

    content = function(file) {
      png(file)
      scopus_dataframe  <- new_data$scopus_dataframe
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
      scopus_dataframe  <- new_data$scopus_dataframe
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
      scopus_dataframe  <- new_data$scopus_dataframe
      Journals <- as.data.frame(table(scopus_dataframe$SO)) %>%
        arrange(desc(Freq)) %>%
        rename(Journal = "Var1")

      write.table(Journals,file,sep=";")
    }
  )

}

# ------------------ CREACION DE LA APLICACION--------  ------------------------

shinyApp(ui, server)
