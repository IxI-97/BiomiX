library(shiny)
library(dplyr)
library(ggplot2)
library(missForest)
library(mixOmics)
library(tidyr)
library(mice)
library(data.table)



get_default_browser <- function() {
        print(getwd())
        directory <- getwd()
        browser_analysis <- readLines(paste(directory,'/_INSTALL/CHOISE_BROWSER_pre-view',sep=""), n = 1)
        print(browser_analysis)
        
        if (.Platform$OS.type == "windows") {
                if (browser_analysis == "chrome") {
                        return("C:/Program Files/Google/Chrome/Application/chrome.exe")
                } else if (browser_name == "firefox") {
                        return("C:/Program Files/Mozilla Firefox/firefox.exe")
                } else if (browser_name == "chromium") {
                        return("C:/Program Files/Chromium/chromium.exe")
                } else {
                        return(NULL)
                }
        } else if (.Platform$OS.type == "unix") {
                if (Sys.info()[["sysname"]] == "Darwin") {
                        if (browser_analysis == "firefox") {
                                return("/Applications/Firefox.app/Contents/MacOS/firefox")
                        } else if (browser_analysis == "chrome" ) {
                                return("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")  # or "chromium-browser" based on your installation
                        } else if (browser_analysis == "chromium-browser" ) {
                                return("/Applications/Chromium.app/Contents/MacOS/Chromium")  # or "chromium-browser" based on your installation
                        } else {
                                return(NULL)
                        }} else {        
                                if (browser_analysis == "firefox") {
                                        return("/usr/bin/firefox")
                                } else if (browser_analysis == "chrome" ) {
                                        return("/usr/bin/chrome")  # or "chromium-browser" based on your installation
                                } else if (browser_analysis == "chromium-browser" ) {
                                        return("/usr/bin/chromium-browser")  # or "chromium-browser" based on your installation
                                } else {
                                        return(NULL)
                                }}
        } else {
                stop("Unsupported OS")
        }
}

        
impute_data <- function(data, method) {
        df <- data  # Start with the original data    
                
        if (method == "Replace with 0") {
                df[is.na(df)] <- 0
        } else if (method == "Mean") {
                df <- df %>%
                        mutate_all(~ifelse(is.na(.), mean(., na.rm = TRUE), .))
        } else if (method == "Median") {
                df <- df %>%
                        mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .))
        } else if (method == "Lasso") {
                # Implement Lasso-based imputation using mice
                imputed_data <- mice(df, method = "lasso.norm", m = 1, maxit = 5, printFlag = FALSE)
                df <- complete(imputed_data)  # Extract the completed dataset
        } else if (method == "Random Forest") {
                # Implement Random Forest Imputation
                df <- missForest(df)$ximp
        } else if (method == "NIPALS") {
                # Check if the matrix has missing values
                if (any(is.na(df))) {
                        # Impute using NIPALS from mixOmics
                        # num_comp is the number of principal components to compute (you can adjust it)
                        nipals_result <- nipals(df, ncomp = 2)
                        
                        # Get the matrix with missing values imputed
                        df <- nipals_result$rec  # Reconstructed matrix with imputed values
                } else {
                        message("No missing values detected, skipping NIPALS imputation.")
                }
        } else if (method == "Replace NA with 0") {
                df[is.na(df)] <- 0  # Replace only NA with 0
        }
        
        return(df)
}


library(shiny)
library(dplyr)
library(tidyr)  # For pivot_longer
library(ggplot2)
library(randomForest)
library(glmnet)
library(mixOmics)

        
        # Define UI for app
        ui <- fluidPage(
                titlePanel("Missing Data Imputation"),
                
                sidebarLayout(
                        sidebarPanel(
                                # Imputation Methods
                                selectInput("imputation_method", "Choose Imputation Method:",
                                            choices = c("Replace with 0", "Mean", "Median", "Lasso", "Random Forest", "NIPALS", "Replace NA with 0")),
                                
                                
                                # Thresholds for filtering variables and samples
                                numericInput("var_threshold", "Variable Missingness Threshold (%)", 50, min = 0, max = 100),
                                numericInput("sample_threshold", "Sample Missingness Threshold (%)", 50, min = 0, max = 100),
                                actionButton("Impute", "Impute the matrix"),
                                actionButton("Filter", "Filter Variables and samples"),
                                br(),
                                br(),
                                actionButton("resetBtn", "Reset to Original Data"),  # Reset button
                                br(),
                                br(),
                                # Close App button
                                actionButton("closeApp", "Close App and Return Data")
                        ),
                        
                        mainPanel(
                                tabsetPanel(
                                        tabPanel("Upload Data",
                                                 fileInput("file", "Upload a File", accept = c(".csv",".tsv",".xls",".xlsx")),
                                                 checkboxInput("transpose", "Transpose Data", value = FALSE),
                                                 checkboxInput("remove_first_row", "Remove First Row", value = FALSE),
                                                 checkboxInput("remove_first_col", "Remove First Column", value = FALSE),
                                                 actionButton("load_data", "Load Data"),
                                                 uiOutput("warningMessage")
                                                 
                                        ),
                                        tabPanel("Boxplot by Variables", plotOutput("boxplot")),
                                        tabPanel("Barplot by Samples", plotOutput("barplot")),
                                        tabPanel("Summary Table", DT::dataTableOutput("summary_table")),
                                        tabPanel("Download",
                                                 downloadButton("downloadData", "Download Normalized Data")
                                        )
                                )
                        )
                ),
                
                # JavaScript handler for closing the app
                tags$script(HTML("
        Shiny.addCustomMessageHandler('closeApp', function(message) {
            window.close();
        });
    "))
        )
        
        server <- function(input, output, session) {
                
                # Make combined_data a reactiveValues object to track changes
                combined_data <- reactiveValues(
                        numeric = NULL)
                
                # Store the original dataset (numeric and metadata)
                original_data <- reactiveValues(numeric = NULL)
                
                
                # Upload data when the "Load Data" button is clicked
                observeEvent(input$load_data, {

                        file_ext <- tools::file_ext(input$file$name)
                        
                        # Read the uploaded file based on its extension
                        if (file_ext == "csv") {
                                df <- fread(input$file$datapath, header = TRUE)  # Using data.table's fread for CSV
                        } else if (file_ext == "tsv") {
                                df <- fread(input$file$datapath, header = TRUE, sep = "\t")  # TSV
                        } else if (file_ext == "xlsx" | file_ext == "xls") {
                                df <- read_excel(input$file$datapath, col_names = TRUE)  # Excel files
                        } else {
                                showNotification("Unsupported file format! Please upload CSV, TSV, or Excel files.", type = "error")
                                return()
                        }
                        
                        
                        # Identify the rows with NA values
                        rows_with_na <- which(is.na(colnames(df)))  # Change 'A' to the relevant column
                        # Generate new row names for the rows with NA values
                        if (length(rows_with_na) > 1){
                        colnames(df)[rows_with_na] <- paste("missing", seq_along(rows_with_na), sep="_")}
                        colnames(df) <- make.unique(colnames(df), sep = "_")
                        
                        # Optionally transpose the data
                        if (input$transpose) {
                                df <- t(df)
                                colnames(df) <- as.character(unlist(df[1, ]))
                                df <- df[-1, ]  # Remove the row used for column names
                                df <- as.data.frame(df, stringsAsFactors = FALSE)  # Convert back to data frame
                                
                        }
                        
                        # Remove the first row if selected
                        if (input$remove_first_row) {
                                df <- df[-1, ]
                        }
                        
                        # Remove the first column if selected
                        if (input$remove_first_col) {
                                sav<- df[, 1]
                                #print(dim(sav))
                                df <- df[, -1]
                                df <- as.data.frame(df)
                                #print(dim(df))
                                #print("OK")
                                sav <- as.vector(sav[[1]])
                                # Identify the rows with NA values
                                rows_with_na <- which(is.na(sav))  # Change 'A' to the relevant column
                                #print(rows_with_na)
                                # Generate new row names for the rows with NA values
                                if (length(rows_with_na) > 1){
                                        sav[rows_with_na] <- paste("missing", seq_along(rows_with_na), sep="_")}
                                #print(sav)
                                sav <- make.unique(sav, sep = "_")
                                
                                row.names(df) <- as.vector(sav)
                                
                        }
                        
                        # # Use the first row as column names

                        
                        # Update the combined_data reactive value
                        #combined_data$numeric <- df
                        combined_data$numeric <- df # Update original data with uploaded data
                        #print(combined_data$numeric)
                        #print(original_data)
                        original_data$numeric <- df # Update original data with uploaded data
                        #print(original_data$numeric)
                        
                        
                        
                })
                
                
                
                # Check if any column contains non-numeric values
                output$warningMessage <- renderUI({

                        # Get matrix
                        plot_data <- combined_data$numeric
                        
                        # Check if any column is non-numeric
                        is_numeric <- sapply(plot_data, is.numeric)
                        if (any(!is_numeric)) {
                                warning_cols <- paste(names(plot_data)[!is_numeric], collapse = ", ")
                                div(
                                        style = "color: red; font-weight: bold;",
                                        paste("Warning: Non-numeric values found in columns:", warning_cols)
                                )
                        } else {
                                return(NULL)  # No warning if all columns are numeric
                        }
                })
                
                
                
                # Barplot: Missing/zero values by variable
                output$boxplot <- renderPlot({
                        plot_data <- combined_data$numeric
                        #print(plot_data)
                        df <- plot_data
                        total_samples <- nrow(df)
                        
                        # Check for NA or 0 values in each variable
                        stats <- df %>%
                                summarise_all(~sum(is.na(.) | . == 0)) %>%
                                pivot_longer(everything(), names_to = "variable", values_to = "count_missing") %>%
                                mutate(percent_missing = (count_missing / total_samples) * 100)
                        variable_stats <- stats
                        
                        ggplot(variable_stats, aes(x = variable, y = percent_missing)) +
                                geom_bar(stat = "identity", fill = "steelblue") +  # Bar plot
                                labs(title = "Percentage of Missing/Zero Values by Variable", y = "% Missing or Zero", x = "Variable") +
                                theme(axis.text.x = element_text(angle = 90, hjust = 1))
                })
                
                # Barplot: Missing/zero values by sample
                output$barplot <- renderPlot({
                        plot_data <- combined_data$numeric
                        df <- plot_data
                        total_vars <- ncol(df)
                        
                        # Calculate number of missing/zero values for each sample (row)
                        sample<-rownames(df)
                        stats <- df %>%
                                rowwise() %>%
                                mutate(count_missing = sum(is.na(c_across()) | c_across() == 0)) %>%
                                ungroup() %>%
                                mutate(percent_missing = (count_missing / total_vars) * 100)
                        
                        stats <- stats[, c(ncol(stats) - 1, ncol(stats))]
                        rownames(stats) <- sample
                        sample_stats <- stats
                        
                        ggplot(sample_stats, aes(x = factor(rownames(stats)), y = percent_missing)) +
                                geom_bar(stat = "identity", fill = "darkorange") +  # Bar plot
                                labs(title = "Percentage of Missing/Zero Values by Sample", y = "% Missing or Zero", x = "Sample") +
                                theme(axis.text.x = element_text(angle = 90, hjust = 1))
                })
                
                # Reset to original data when reset button is clicked
                observeEvent(input$resetBtn, {
                        combined_data$numeric <- original_data$numeric # Reset the numeric data

                        # Clear any warning messages or notifications
                        output$warningMessage <- renderText({ NULL })
                        
                        # Show a notification to confirm data reset
                        showNotification("Data has been reset to its original state", type = "message", duration = 3)
                })
                
                # Filter data based on thresholds
                observeEvent(input$Filter, {
                        df <- combined_data$numeric
                        
                        # Filter variables based on missing percentage
                        var_stats <- variable_stats
                        vars_to_keep <- var_stats %>%
                                filter(percent_missing <= input$var_threshold) %>%
                                pull(variable)
                        
                        df_filtered <- df[, c(vars_to_keep)]
                        
                        # Filter samples based on missing percentage
                        samples_to_keep <- which(sample_stats$percent_missing <= input$sample_threshold)
                        df_filtered <- df_filtered[samples_to_keep, ]
                        
                        # Update reactiveValues object
                        combined_data$numeric <- df_filtered
                        })
                
                # Summary table of the filtered and imputed data
                output$summary_table <- DT::renderDataTable({
                        #print(combined_data$numeric)
                        plot_data <- combined_data$numeric
                        if (ncol(plot_data) > 50 | nrow(plot_data) > 50) {
                                if (ncol(plot_data) > 50) {
                                        showNotification("More than 50 variables detected. Table displaying only the top 50 variable", type = "warning", duration = 5)
                                        plot_data[,1:50] }
                                if (nrow(plot_data) > 50) {
                                        showNotification("More than 50 samples detected. Table displaying only the top 50 samples", type = "warning", duration = 5)
                                        plot_data[1:50,] }
                        } else {
                                plot_data
                        }
                })
                
                # Show data in an interactive table
                output$dataTable <- DT::renderDataTable({
                        combined_data$numeric  # Explicitly trigger updates with reactiveValues
                })
                
                # Imputation methods
                observeEvent(input$Impute, {
                        df <- combined_data$numeric
                        method <- input$imputation_method
                        combined_data$numeric <- impute_data(df, method)
                })
                
                # Download handler to download normalized data as CSV
                output$downloadData <- downloadHandler(
                        filename = function() { paste("imputated_data", Sys.Date(), ".tsv", sep = "") },
                        content = function(file) {
                                saved<-cbind(rownames(combined_data$numeric), combined_data$numeric)
                                colnames(saved)[1] <- "ID"
                                write.table(saved, file, row.names = FALSE, quote = FALSE, sep = "\t")
                        }
                )
                
                # Closing the app and returning the data
                observeEvent(input$closeApp, {
                        session$sendCustomMessage("closeApp", "Closing app")
                        stopApp()
                })
        }
        
        # Get the browser path and set it
        browser_path <- get_default_browser()
        if (is.null(browser_path)) {
                stop("Browser path could not be determined.")
        }
        
        # Launch the Shiny app
        options(browser = get_default_browser())
        shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))  # Launch browser