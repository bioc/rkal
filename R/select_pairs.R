#' Select pairs and replicates for an experiment.
#'
#' Used to select a pairs for paired-end experiments and replicate samples. Please follow prompt, ensuring correct file name matching and
#' end-type experiment identification.
#'
#' @param data_dir Directory with raw fastq.gz RNA-Seq files.
#'
#' @return data.frame with character column \code{'File Name'} and integer columns \code{'Pair'} (if pair-ended) and \code{'Replicate'} (if added).
#' @export
#'
#' @examples
#'
select_pairs <- function(data_dir) {

  # setup ----

  group_colors <- c("#C7E9C0", "#C6DBEF", "#FCBBA1", "#FDD0A2", "#BCBDDC", "#D9D9D9", "#F6E8C3", "#DC143C",
                    "#A1D99B", "#9ECAE1", "#FC9272", "#FDAE6B", "#9E9AC8", "#BDBDBD", "#DFC27D", "#F0EAD6",
                    "#C3B091", "#007FFF", "#00FFFF", "#7FFFD4", "#228B22", "#808000", "#7FFF00", "#BFFF00",
                    "#FFD700", "#DAA520", "#FF7F50", "#FA8072","#FC0FC0", "#CC8899", "#E0B0FF", "#B57EDC", "#843179")

  ncolors <- length(group_colors)

  background <- 'url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAPklEQVQoU43Myw0AIAgEUbdAq7VADCQaPyww55dBKyQiHZkzBIwQLqQzCk9E4Ytc6KEPMnTBCG2YIYMVpHAC84EnVbOkv3wAAAAASUVORK5CYII=) repeat'


  fastqs <- list.files(data_dir, '.fastq.gz')
  pdata <- tibble::tibble('File Name' = fastqs)

  pdata <- tibble::add_column(pdata, Pair = NA_character_, Replicate = NA_character_, .before = 1)

  # auto-detect if paired
  fastq_id1s <- get_fastq_id1s(file.path(data_dir, fastqs))
  paired <- detect_paired(fastq_id1s)

  # select and mark auto-detected pair type
  end_types <- c('single-ended', 'pair-ended')
  if (paired) end_types <- end_types[c(2, 1)]
  end_types[1] <- paste(end_types[1], '(detected)')

  message("For paired-end experiments:
            - select and add paired rows (usually contain _1 and _2 in filename)\n")

  message("For paired-end and single-ended experiments:
            - confirm that file names (in column 'File Name') are assigned to the correct sample rows
            - confirm that the experiment was correctly identified as single-ended or paired-end
            - select samples to treat as a single library (if any - e.g. same sample sequenced in replicate)\n")

  message("Click 'Done' to continue with kallisto quantification.")


  # things we will update/return to user
  pairs <- reps <- rep(NA, nrow(pdata))


  #  user interface ----

  ui <- miniUI::miniPage(
    shinyjs::useShinyjs(),
    shiny::tags$head(
      shiny::tags$style("#pdata {white-space: nowrap;}"), # table text on 1 line
      shiny::tags$style(".dt-fake-height {height: 1px;}"), # to make 100% height div work
      shiny::tags$style("td.dt-nopad {padding: 0px !important; height: 100%;}"), # td for bg color group column
      shiny::tags$style("td.dt-nopad div {height: 100%; width: 100%; text-align: center;}"), # div inside td for bg color group column
      shiny::tags$style("td.dt-nopad span {display: inline-block; padding: 8px 10px; color: white;}") # span inside div inside td for bg color group column
    ),
    # title bar
    miniUI::gadgetTitleBar(shiny::textOutput('title', inline=TRUE), left = miniUI::miniTitleBarButton("reset", "Reset")),
    miniUI::miniContentPanel(
      shiny::fillCol(flex = c(NA, NA, 1),
                     shiny::selectizeInput('end_type',
                                           'Confirm end-type:',
                                           choices = end_types),
                     shiny::hr(),
                     DT::dataTableOutput("pdata")
      )
    ),
    miniUI::miniButtonBlock(
      shiny::actionButton("pair", "Pair Samples"),
      shiny::actionButton("rep", "Mark Replicates")
    )
  )

  # server ----

  server <- function(input, output, session) {

    # make reactive state value to keep track of ctrl vs test group
    state <- shiny::reactiveValues(pdata = 0)



    # show phenotype data
    output$pdata <- DT::renderDataTable({

      DT::datatable(
        isolate(pdata_r()),
        class = 'cell-border dt-fake-height',
        rownames = FALSE,
        escape = FALSE, # to allow HTML in table
        options = list(
          columnDefs = list(list(className = 'dt-nopad', targets = c(0, 1))),
          scrollY = FALSE,
          paging = FALSE,
          bInfo = 0
        )
      )
    })

    # pdata reactive so that will update Pair/Replicate column
    pdata_r <- shiny::eventReactive(state$pdata, {

      # update pdata Replicate column
      rep_nums <- sort(unique(setdiff(reps, NA)))
      for (rep_num in rep_nums) {
        color <- group_colors[ncolors - rep_num]
        rows <- which(reps == rep_num)
        pdata[rows, 'Replicate'] <<- paste('<div style="background-color:', color, ';"></div>')
      }

      # update pdata Pair column
      if (paired) {
        pair_nums <- sort(unique(setdiff(pairs, NA)))
        for (pair_num in pair_nums) {
          color <- group_colors[pair_num]
          rows <- which(pairs == pair_num)
          pdata[rows, 'Pair'] <<- paste('<div style="background:', color, background, ';"></div>')
        }
      } else {
        pdata[1:nrow(pdata), 'Pair'] <<- NA
      }

      return(pdata)
    })

    # proxy used to replace data
    pdata_proxy <- DT::dataTableProxy("pdata")
    shiny::observe({
      DT::replaceData(pdata_proxy, pdata_r(), rownames = FALSE)
    })


    # click 'Pair Samples' ----

    shiny::observeEvent(input$pair, {

      # get rows
      rows  <- input$pdata_rows_selected

      # check for incomplete/wrong input
      if (is.null(validate_pairs(pairs, rows, reps))) {

        # add rows as a pair
        pair_num <- length(unique(setdiff(pairs, NA))) + 1
        pairs[rows] <<- pair_num

        # update states to trigger updates
        state$pdata <- state$pdata + 1
      }
    })

    # click 'Mark Replicates' ----

    shiny::observeEvent(input$rep, {

      # get rows
      rows  <- input$pdata_rows_selected
      if (validate_reps(pairs, rows, reps)) {
        # add rows as replicates
        rep_num <- length(unique(setdiff(reps, NA))) + 1
        reps[rows] <<- rep_num

        # update states to trigger updates
        state$pdata <- state$pdata + 1
      }
    })

    # click 'Done' ----

    shiny::observeEvent(input$done, {
      if (paired & sum(is.na(pdata$Pair)) != 0) {
        message("All samples must be in a pair for paired-end experiments.")
      } else {

        # replace html with integers identifying pairs/replicates
        if (paired) {
          pdata$Pair <<- pairs
        } else {
          pdata$Pair <<- NULL
        }

        if (all(is.na(reps))) {
          pdata$Replicate <<- NULL
        } else {
          pdata$Replicate <<- reps
        }

        shiny::stopApp(pdata)
      }
    })

    # select 'end-type' ----
    shiny::observeEvent(input$end_type, {

      # single-ended experiments are not paired, so disable
      # hide pair column labels if select single-ended

      if (grepl('^single', input$end_type)) {
        shinyjs::disable("pair")
        paired <<- FALSE
        state$pdata <- state$pdata + 1

      } else {
        shinyjs::enable("pair")
        paired <<- TRUE
        state$pdata <- state$pdata + 1
      }
    })


    # click 'Reset' ----

    shiny::observeEvent(input$reset, {
      pairs <<- rep(NA, nrow(pdata))
      reps <<- rep(NA, nrow(pdata))
      pdata$Pair <<- NA
      pdata$Replicate <<- NA

      # remove groups from table and reset control
      state$pdata <- state$pdata + 1
    })
  }

  shiny::runGadget(shiny::shinyApp(ui, server), viewer = shiny::paneViewer())
}
