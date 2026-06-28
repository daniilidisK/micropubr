# Interactive exploration --------------------------------------------------

#' Interactively explore a phyloseq object
#'
#' Launches a Shiny gadget that wraps [mp_alpha()], [mp_ordination()] and
#' [mp_taxa()] behind interactive controls, with a live preview and the
#' equivalent reproducible `micropubr` code. Closing the gadget with **Done**
#' returns the currently displayed `ggplot` object so it can be saved or
#' further edited.
#'
#' @param object A [phyloseq::phyloseq] object.
#'
#' @return Invisibly, the last `ggplot` shown (when run interactively).
#' @export
#' @examples
#' \dontrun{
#' data(GlobalPatterns, package = "phyloseq")
#' p <- mp_explore(GlobalPatterns)
#' }
mp_explore <- function(object) {
  .mp_check_phyloseq(object)
  .mp_require(c("shiny", "miniUI"), "Interactive exploration")

  vars <- c("(none)", phyloseq::sample_variables(object))
  ranks <- phyloseq::rank_names(object)
  has_tree <- !is.null(phyloseq::phy_tree(object, errorIfNULL = FALSE))

  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("micropubr explorer"),
    miniUI::miniTabstripPanel(
      miniUI::miniTabPanel(
        "Plot", icon = shiny::icon("chart-area"),
        miniUI::miniContentPanel(
          shiny::fillRow(
            flex = c(1, 3),
            shiny::tagList(
              shiny::selectInput("kind", "Plot",
                                 c("Alpha diversity" = "alpha",
                                   "Ordination" = "ordination",
                                   "Composition" = "taxa")),
              shiny::selectInput("group", "Group by", vars),
              shiny::conditionalPanel(
                "input.kind == 'alpha'",
                shiny::selectInput("measure", "Measure",
                                   c("Shannon", "Observed", "Simpson",
                                     "InvSimpson", "Chao1"))),
              shiny::conditionalPanel(
                "input.kind == 'ordination'",
                shiny::selectInput("method", "Method",
                                   c("PCoA", "NMDS", "PCA")),
                shiny::selectInput("distance", "Distance",
                                   c("bray", "jaccard", "aitchison",
                                     if (has_tree) c("unifrac", "wunifrac")))),
              shiny::conditionalPanel(
                "input.kind == 'taxa'",
                shiny::selectInput("rank", "Rank", ranks,
                                   selected = utils::tail(ranks, 1)),
                shiny::sliderInput("topn", "Top N taxa", 3, 30, 10))
            ),
            shiny::plotOutput("plot", height = "100%")
          )
        )
      ),
      miniUI::miniTabPanel(
        "Code", icon = shiny::icon("code"),
        miniUI::miniContentPanel(shiny::verbatimTextOutput("code"))
      )
    )
  )

  server <- function(input, output, session) {
    grp <- shiny::reactive(if (input$group == "(none)") NULL else input$group)

    current <- shiny::reactive({
      switch(
        input$kind,
        alpha = mp_alpha(object, measure = input$measure, group.by = grp()),
        ordination = mp_ordination(object, method = input$method,
                                   distance = input$distance,
                                   group.by = grp()),
        taxa = mp_taxa(object, tax.level = input$rank, top.n = input$topn,
                       group.by = grp())
      )
    })

    output$plot <- shiny::renderPlot(print(current()))
    output$code <- shiny::renderText(.mp_explore_code(input, grp()))

    shiny::observeEvent(input$done, {
      shiny::stopApp(invisible(current()))
    })
    shiny::observeEvent(input$cancel, shiny::stopApp(invisible(NULL)))
  }

  shiny::runGadget(ui, server,
                   viewer = shiny::dialogViewer("micropubr", width = 1000,
                                                height = 680))
}

#' Generate the reproducible call shown in the explorer's Code tab
#' @noRd
.mp_explore_code <- function(input, group) {
  g <- if (is.null(group)) "" else sprintf(", group.by = \"%s\"", group)
  switch(
    input$kind,
    alpha = sprintf("mp_alpha(ps, measure = \"%s\"%s)", input$measure, g),
    ordination = sprintf(
      "mp_ordination(ps, method = \"%s\", distance = \"%s\"%s)",
      input$method, input$distance, g),
    taxa = sprintf(
      "mp_taxa(ps, tax.level = \"%s\", top.n = %d%s)",
      input$rank, input$topn, g)
  )
}
