get_references <- function(data) {

  references_df <-
    data |>
    select(SR, CR) |>
    na.omit()  |>
    separate_rows(CR, sep = "; ") |>
    mutate(PY = str_extract(CR, "\\([0-9]{4}\\)"),
           PY = str_remove_all(PY, "[\\(\\)]"),
           PY = as.numeric(PY)) |>
    na.omit() |>
    mutate(AU = str_extract(CR, ".*\\([0-9]{4}\\)"),
           AU = str_extract(AU, ".*\\.,"),
           AU = gsub("([^,]+,[^,]+),", "\\1;", AU),
           AU = str_sub(AU, 1, nchar(AU)-1),
           AU = str_replace_all(AU,
                                pattern = "; ",
                                replacement = ";"),
           AU = str_remove_all(AU, pattern = "\\."),
           AU = str_remove_all(AU, pattern = ",")) |>
    mutate(type_ref = if_else(str_detect(CR,
                                         "\\., \\("), 2, # books
                              if_else(str_detect(CR,
                                                 "^\\([0-9]{4}\\)"), 3,
                                      if_else(str_detect(CR,
                                                         " \\([0-9]{4}\\), "), 4,
                                              1)))) |> # papers
    mutate(TI = if_else(type_ref == 1,
                        str_extract(CR,
                                    ".*\\([0-9]{4}\\)"),
                        CR)) |>
    mutate(TI = if_else(type_ref == 1,
                        str_remove(TI, "\\([0-9]{4}\\)"),
                        TI)) |>
    mutate(TI = if_else(type_ref == 1,
                        str_remove(TI, ".*\\., "),
                        TI)) |>
    mutate(TI = if_else(type_ref == 1,
                        str_trim(TI),
                        TI)) |>
    mutate(TI = if_else(type_ref == 1,
                        str_remove(TI, "\""),
                        TI)) |>
    mutate(TI = if_else(type_ref == 1,
                        str_remove(TI, "\"$"),
                        TI)) |>
    mutate(TI = if_else(type_ref == 2,
                        str_extract(CR,
                                    "\\([0-9]{4}\\).*"),
                        TI)) |>
    mutate(TI = if_else(type_ref == 2,
                        str_remove(TI, "\\([0-9]{4}\\)"),
                        TI)) |>
    mutate(TI = if_else(type_ref == 2,
                        str_remove(TI, ", [0-9].*"),
                        TI)) |>
    mutate(TI = if_else(type_ref == 2,
                        str_trim(TI),
                        TI)) |>
    mutate(TI = if_else(type_ref == 3,
                        str_remove(CR, "\\([0-9]{4}\\)"),
                        TI)) |>
    mutate(TI = if_else(type_ref == 3,
                        str_remove(TI, ", ,.*"),
                        TI)) |>
    mutate(TI = if_else(type_ref == 3,
                        str_trim(TI),
                        TI)) |>
    mutate(TI = if_else(type_ref == 4,
                        str_extract(CR,
                                    ".* \\([0-9]{4}\\) "),
                        TI)) |>
    mutate(TI = if_else(type_ref == 4,
                        str_remove(TI, "\\([0-9]{4}\\)"),
                        TI)) |>
    mutate(TI = if_else(type_ref == 4,
                        str_remove(TI, ".*\\., "),
                        TI)) |>
    mutate(TI = if_else(type_ref == 4,
                        str_trim(TI),
                        TI)) |>
    mutate(JI = if_else(type_ref == 1,
                        str_remove(CR, ".*\\([0-9]{4}\\)"),
                        CR)) |>
    mutate(JI = if_else(type_ref == 1,
                        str_remove(JI, ", .*"),
                        JI)) |>
    mutate(JI = if_else(type_ref == 1,
                        str_trim(JI),
                        JI)) |>
    filter(JI != "") |>
    mutate(JI = str_remove_all(JI, "\\.")) |>
    mutate(SR_ref = gsub("^(.*?);.*", "\\1", AU),
           SR_ref = str_c(SR_ref, ", ", PY, ", ", JI, sep = "")) |>
    select(SR, SR_ref, TI, AU, JI, PY,
           ref_type = type_ref, CR_ref = CR) |>
    filter(!is.na(SR_ref))
  # Finding right SR_ref in main dataset
  references_df_1 <- tibble()

  list_ref_TI <-
    references_df |>
    filter(ref_type == 1) |>
    select(TI) |>
    dplyr::distinct()

  # for (i in list_ref_TI$TI) {
  #
  #   df_1 <-
  #     scopus_df |>
  #     filter(TI == i)
  #
  #   if (length(df_1$TI) != 0) {
  #
  #     ref_1 <- references_df %>% filter(!is.na(TI))
  #     ref_1$SR_ref[ref_1$TI == df_1$TI] <- df_1$SR
  #
  #     references_df$SR_ref[references_df$TI == df_1$TI] <- df_1$SR
  #
  #     references_df_1 <-
  #       references_df |>
  #       mutate(SR_ref = replace(SR_ref,
  #                               TI == i,
  #                               df_1 |>
  #                                 filter(TI == i) |>
  #                                 select(SR) |>
  #                                 pull()))
  #   }
  #
  # }
  #
  # Finding duplicate titles in references with different SR_ref

  TI_ref_duplicates <-
    references_df |>
    filter(ref_type == 1) |>
    group_by(TI, SR_ref) |>
    count(TI, sort = TRUE) |>
    ungroup() |>
    select(TI) |>
    group_by(TI) |>
    count(TI) |>
    filter(n > 1) |>
    select(TI)

  # Create the

  df_1 <- references_df

  for (i in TI_ref_duplicates$TI) {

    SR_ref_1 <-
      references_df |>
      filter(TI == i) |>
      select(SR_ref) |>
      group_by(SR_ref) |>
      count(SR_ref, sort = TRUE) |>
      ungroup() |>
      slice(1) |>
      select(SR_ref)

    df_1 <-
      df_1 |>
      mutate(SR_ref = if_else(TI == i, SR_ref_1$SR_ref,
                              SR_ref))

  }

  return(references_df = df_1)
}


get_asn <- function(data) {

  data_tidied_scopus_ref <-
    data |>
    dplyr::select(SR, CR) |>
    separate_rows(CR, sep = "; ") |>
    na.omit() |>
    mutate(PY = str_extract(CR, "\\([0-9]{4}\\)"),
           PY = str_remove_all(PY, "[\\(\\)]")) |>
    na.omit() |>
    mutate(AU = str_extract(CR, ".*\\([0-9]{4}\\)"),
           AU = str_extract(AU, ".*\\.,"),
           AU = gsub("([^,]+,[^,]+),", "\\1;", AU),
           AU = str_sub(AU, 1, nchar(AU)-1),
           AU = str_replace_all(AU,
                                pattern = "; ",
                                replacement = ";"),
           AU = str_remove_all(AU, pattern = "\\."),
           AU = str_remove_all(AU, pattern = ",")) |>
    na.omit()|>
    rename(main_ref = SR,
           id_ref = CR)

  dummy <-
    data_tidied_scopus_ref |>
    dplyr::select(AU, PY) |>
    mutate(PY = as.numeric(PY)) |>
    bind_rows(data |>
                dplyr::select(AU, PY)) |>
    unique()

  asn <-
    biblioNetwork(M = data.frame(dummy),
                  analysis = "collaboration",
                  network = "authors") |>
    graph_from_adjacency_matrix(mode = "undirected",
                                weighted = TRUE) |>
    simplify() |>
    as_tbl_graph() |>
    activate(nodes) |>
    mutate(communities = group_components(type = "weak")) |>
    filter(communities == 1)

  return(asn)
}

get_citation_network <- function(scopus_df, references_df) {

  sca_citation_network <-
    references_df |>
    filter(ref_type == 1) |>
    select(SR, SR_ref) |>
    na.omit() |>
    dplyr::distinct() |>
    as_tbl_graph()

  sca_data_cleaned_1 <-
    scopus_df |>
    select(SR, TI, PY) |>
    bind_rows(references_df |>
                select(SR = SR_ref,
                       TI,
                       PY)) |>
    dplyr::distinct() |>
    dplyr::rename(name = SR) |>
    filter(!duplicated(name))

  sca_citation_network_1 <-
    sca_citation_network |>
    activate(nodes) |>
    left_join(sca_data_cleaned_1, by = "name") |>
    activate(edges)
}

get_citation_network_tos <- function(citation_network) {

  citation_network_gc <-
    citation_network |>
    activate(nodes) |>
    mutate(component = group_components(type = "weak")) |>
    filter(component == 1) %>%
    mutate(in_degree = centrality_degree(mode = "in"),
           out_degree = centrality_degree(mode = "out")) |>
    filter(!(in_degree == 1 & out_degree == 0)) |>
    select(-in_degree, -out_degree)

  subfields <-
    citation_network_gc |>
    tidygraph::to_undirected() |>
    activate(nodes) |>
    mutate(subfield = tidygraph::group_louvain()) |>
    as_tibble() |>
    select(name, subfield)

  citation_network_subfield <-
    citation_network_gc |>
    activate(nodes) |>
    left_join(subfields, by = "name")

  df_tos = tibble()

  for (i in unique(subfields$subfield) ) {

    df_1 <-
      citation_network_subfield |>
      activate(nodes) |>
      filter(subfield == i) |>
      mutate(in_degree = centrality_degree(mode = "out"),
             out_degree = centrality_degree(mode = "in"),
             bet = centrality_degree(mode = "in")) |>
      as_tibble() |>
      select(name,
             in_degree,
             out_degree
      )

    df_tos <-
      df_tos |>
      bind_rows(df_1)

  }

  citation_network_tos <-
    citation_network_subfield |>
    activate(nodes) |>
    left_join(df_tos, by = "name") |>
    select(-component)

}

get_sap <- function(citation_network) {
  nodes <-
    citation_network %>%
    activate(nodes) %>%
    data.frame() %>%
    rownames_to_column("rowid") %>%
    mutate(rowid = as.integer(rowid))

  edges <-
    citation_network %>%
    activate(edges) %>%
    data.frame()

  for (i in 1 : nrow(edges)){
    from = edges[i, 1]
    to   = edges[i, 2]
    edges[i, 1] = nodes[from, 'name']
    edges[i, 2] = nodes[to, 'name']
  }

  g <- graph_from_data_frame(edges, directed = TRUE) %>%
    simplify()

  # Se eliminan los vertices con indegree = 1 y con outdegree = 0
  g1 <- delete.vertices(g,
                        which(degree(g, mode = "in") == 1 &
                                degree(g, mode = "out") == 0))

  # Se escoge el componente mas grande conectado
  g2 <- giant_component_extract(g1, directed = TRUE)
  g2 <- g2[[1]]


  metricas.red <- tibble(
    id        = V(g2)$name,
    indegree  = degree(g2, mode = "in"),
    outdegree = degree(g2, mode = "out"),
    bet       = betweenness(g2))


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
    h <- get.all.shortest.paths(g2,
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
    h <- get.all.shortest.paths(g2,
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

  tos_articles <- c(TOS['id'])
  titles       <- c()
  PY           <- c()
  for (id in tos_articles$id){
    mask <- nodes['name'] == id
    titles <- c(titles, nodes[mask,'TI'])
    PY     <- c(PY, nodes[mask,'PY'])
  }
  TOS <- TOS %>%
    mutate(Title = titles, PY = PY)

  return(TOS)

}
