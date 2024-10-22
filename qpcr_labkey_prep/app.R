library(shiny)
library(tidyverse)
library(bslib)

# Define UI for slider demo app ----
ui <- page_sidebar(

    # App title ----
    title = "Prepare Instrument qPCR Results for LabKey Upload",

    # Sidebar panel for inputs ----
    sidebar = sidebar(

        # Input: Select a file ----
        fileInput(
            "file1",
            "Choose Instrument qPCR Results",
            multiple = TRUE,
            accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv"
            )
        ),

        # Input: Checkbox if file has header ----
        checkboxInput("header", "Header", TRUE),

        # Input: Select separator ----
        radioButtons(
            "sep",
            "Separator",
            choices = c(
                Comma = ",",
                Semicolon = ";",
                Tab = "\t"
            ),
            selected = "\t"
        ),

        # Input: Define rows to skip, if any
        numericInput(
            inputId = "skip",
            label = "Rows at the top to skip:",
            value = 0,
            min = 0,
        ),

        # Input: Select quotes ----
        radioButtons(
            "quote",
            "Quote",
            choices = c(
                None = "",
                "Double Quote" = '"',
                "Single Quote" = "'"
            ),
            selected = '"'
        ),

        # Horizontal line ----
        tags$hr(),

        # Input: Choose an instrument
        selectInput(
            inputId = "instrument",
            label = "Choose an instrument:",
            choices = c("", "LC480", "LC96"),
            selected = "",
        ),

        # Input: the current experiment number
        numericInput(
            inputId = "experiment",
            label = "Experiment Number:",
            value = NULL,
            min = 1,
        ),

        # Input: the current date
        dateInput(
            inputId = "date",
            label = "qPCR Run Date:",
            value = NULL,
        ),

        # Input: Choose a pathogen
        selectInput(
            inputId = "pathogen",
            label = "Choose a pathogen:",
            choices = c("", "IAV", "SC2"),
            selected = "",
        ),

        # Input: Choose the lab
        textInput(
            inputId = "lab",
            label = "PCR Lab:",
            value = "AVRL"
        ),

        # Input: Choose the cq_conf amp_score
        textInput(
            inputId = "cq_conf",
            label = "Cq configuration:",
            value = NULL
        ),

        # Input: Choose the amp_score
        numericInput(
            inputId = "amp_score",
            label = "Amplification score:",
            value = NULL,
            min = 0,
        ),

        # Input: Add comments:
        textInput(
            inputId = "comments",
            label = "Additional comments:",
            value = NULL,
        ),

        # Input: the name of the assay
        selectInput(
            inputId = "assay",
            label = "Name of the assay run:",
            choices = c(
                "",
                "IAV M gene",
                "IAV M gene AVRL-20240316-multiplex",
                "CDC N1",
                "CDC N2",
                "CDC N1 AVRL-20240316-multiplex",
                "CDC N1 AVRL-20240710-multiplex",
                "IAV M gene-20240710-multiplex"
            ),
            selected = "",
        ),

        # Input: the dilution factor
        numericInput(
            inputId = "dilution",
            label = "Dilution factor",
            value = NULL,
            min = 0,
        ),

        # Horizontal line ----
        tags$hr(),

        # Button
        downloadButton("downloadData", "Download Dataset"),
    ),

    # Output: Data file ----
    tableOutput("table"),
)

processLC480 <- function(plate_results, input) {

    labkey_format <- plate_results |>
        # Add new columns and values for experiment, pcr_lab, pcr_date, comments
        mutate(
            experiment = input$experiment,
            pcr_lab = input$lab,
            pcr_date = as.character(input$date),
            cq_conf = input$cq_conf,
            amp_score = input$amp_score,
            comments = input$comments,
            pcr_assay = input$assay
        ) |>
        # Change column names for Sample, Target, Cq, Cq.Conf, and Amp.Score
        # to match labkey input columns
        rename(
            cartridge_id = Name,
            pcr_ct = Cp,
            pcr_copies = Concentration
        ) |>
        # Replace pcr_ct values with values accepted by labkey
        mutate(
            pcr_ct = if_else(is.na(pcr_copies), 99, pcr_ct)
        ) |>
        # add dilution factor
        mutate(
            dilution_factor = if_else(grepl("AE", cartridge_id), input$dilution, NA)
        ) |>
        # reorder columns required for labkey
        select(
            experiment, cartridge_id, pcr_lab, pcr_date, pcr_assay, pcr_ct, pcr_copies,
            dilution_factor, cq_conf, amp_score, comments
        )

    return(labkey_format)
}

processLC96 <- function(plate_results, input) {

    labkey_format <- plate_results |>
        # Add new columns and values for experiment, pcr_lab, pcr_date, comments
        mutate(
            experiment = input$experiment,
            pcr_lab = input$lab,
            pcr_date = as.character(input$date),
            cq_conf = input$cq_conf,
            amp_score = input$amp_score,
            comments = input$comments,
        ) |>
        # Change column names for Sample, Target, Cq, Cq.Conf, and Amp.Score to
        # match labkey input columns
        rename(
            cartridge_id = Sample.Name,
            pcr_assay = Gene.Name,
            pcr_ct = Cq.Mean,
            pcr_copies = Concentration.Mean,
        ) |>
        # Replace pcr_ct, pcr_copies, and pcr_assays with values accepted by labkey
        mutate(
            pcr_ct = str_replace(pcr_ct, "-", "99") |> as.numeric(),
            pcr_copies = str_remove(pcr_copies, "-"),
            pcr_assay = if_else(
                pcr_assay %in% c("N1", "N2"),
                paste("CDC", pcr_assay, sep = " "),
                pcr_assay
            )
        ) |>
        # add dilution factor
        mutate(
            dilution_factor = if_else(grepl("AE", cartridge_id), input$dilution, NA)
        ) |>
        select(
            experiment, cartridge_id, pcr_lab, pcr_date, pcr_assay, pcr_ct,
            pcr_copies, dilution_factor, cq_conf, amp_score, comments
        )

    return(labkey_format)

}


# Define server logic to read selected file ----
server <- function(input, output) {

    datasetInput <- reactive({

        req(input$file1)
        plate_results <- read.csv(
            input$file1$datapath,
            header = input$header,
            sep = input$sep,
            quote = input$quote,
            skip = input$skip,
        )

        if (input$instrument == "LC480") {

            labkey_format <- processLC480(plate_results, input)

        } else if (input$instrument == "LC96") {

            labkey_format <- processLC96(plate_results, input)

            return(labkey_format)

        } else {

            return(plate_results)

        }

    })

    output$table <- renderTable({
        datasetInput()
    })

    output$downloadData <- downloadHandler(
        filename = function() {

            # Use the selected dataset as the suggested file name
            paste(input$experiment, "_labkey_", input$pathogen, "_results.csv", sep = "")
        },
        content = function(file) {

            # Write the dataset to the `file` that will be downloaded
            write_csv(datasetInput(), file, na = "")
        }
    )
}

# Create Shiny app ----
shinyApp(ui, server)
