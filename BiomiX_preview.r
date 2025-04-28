library(shiny)
library(ggplot2)
library(DT)
library(preprocessCore)  # Quantile normalization
library(stats)           # Loess normalization
library(pheatmap)        # Heatmap generation
library(umap)            # UMAP
library(MASS)            # For Mahalanobis distance




get_default_browser <- function() {
        if (.Platform$OS.type == "windows") {
                
                if (browser_analysis == "chrome") {
                        return("C:/Program Files/Google/Chrome/Application/chrome.exe")
                } else if (browser_analysis == "firefox") {
                        return("C:/Program Files/Mozilla Firefox/firefox.exe")
                } else if (browser_analysis == "chromium") {
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




# Function to get the top N most variable features
get_top_variable_genes <- function(data, top_n = 500) {
        # Compute variance of each feature
        feature_variances <- apply(data, 2, var)
        
        # Get indices of top N most variable genes
        top_genes <- order(feature_variances, decreasing = TRUE)[1:top_n]
        
        # Return the column names of the top N genes
        return(colnames(data)[top_genes])
}



get_top_variable_features <- function(data, top_n = 1000) {
        # Calculate variance for each feature (gene)
        variances <- apply(data, 2, var, na.rm = TRUE)
        
        # Get indices of the top N most variable features
        top_features <- order(variances, decreasing = TRUE)[1:top_n]
        
        # Subset the data to keep only the top N features
        return(data[, top_features])
}


# Function to normalize data
normalize_data <- function(data, methods) {
        result <- data  # Start with the original data
        #print(names(result))
        
        # Check for negative values for logarithmic transformation
        if (any(result < 0, na.rm = TRUE ) && "Logarithmic" %in% methods) {
                warning("Cannot apply logarithmic transformation due to negative values in the dataset.")
        } else {
                for (method in methods) {
                        if (method == "Logarithmic") {
                                sample_names <- rownames(result)
                                result <- as.data.frame(lapply(result, function(x) log(x + 1)))
                                rownames(result) <- sample_names
                        } else if (method == "Median Centralization") {
                                sample_names <- rownames(result)
                                result <- as.data.frame(lapply(result, function(x) x - median(x, na.rm = TRUE)))
                                rownames(result) <- sample_names
                        } else if (method == "Mean Centralization") {
                                sample_names <- rownames(result)
                                result <- as.data.frame(lapply(result, function(x) x - mean(x, na.rm = TRUE)))
                                rownames(result) <- sample_names
                        } else if (method == "Z-Score") {
                                sample_names <- rownames(result)
                                result <- as.data.frame(lapply(result, scale))
                                rownames(result) <- sample_names
                        } else if (method == "MAD") {
                                sample_names <- rownames(result)
                                result <- as.data.frame(lapply(result, function(x) (x - median(x, na.rm = TRUE)) / mad(x, na.rm = TRUE)))
                                rownames(result) <- sample_names
                        } else if (method == "Quantile Normalization") {
                                sample_names <- rownames(result)
                                result_values <- preprocessCore::normalize.quantiles(as.matrix(result))
                                result <- as.data.frame(result_values)
                                colnames(result) <- colnames(data)[1:ncol(result)]  # Preserve column names
                                rownames(result) <- sample_names
                                # Check for NaN values after normalization
                                if (any(is.na(result))) {
                                        warning("Quantile normalization resulted in NaN values.")
                                }
                        } else if (method == "Loess Normalization") {
                                for (col_name in names(result)) {
                                        #print(col_name)
                                        if (any(is.na(result[[col_name]]))){  #Keep the NA after the Loess normalization
                                                # Define the positions where you want to insert NA (in ascending order)
                                                positions <- which(is.na(result[[col_name]]))
                                                # Insert NA values at specified positions, adjusting for the shift after each insertion
                                                filled_by_na = loess(result[[col_name]] ~ seq_along(result[[col_name]]), span = 0.5)$fitted
                                                for (i in positions ) {
                                                        filled_by_na <- append(filled_by_na, NA, after = as.numeric(i) - 1)
                                                }
                                                result[[col_name]] <- filled_by_na
                                        }else{
                                        result[[col_name]] <- loess(result[[col_name]] ~ seq_along(result[[col_name]]), span = 0.5)$fitted
                                        }
                                }
                        } else if (method == "VST") {
                                library(DESeq2)  # Make sure DESeq2 is loaded
                                
                                result_matrix <- as.matrix(result)
                                #print(result)
                                result_matrix<-t(result_matrix)
                                #print(head(result_matrix, n=1))
                                sample_names <- rownames(result_matrix)
                                
                                #VST normalization
                                result_matrix<- apply(result_matrix,2,round)
                                colData <- data.frame(condition = rep(1, ncol(result_matrix)))  # Create a data frame for colData
                                rownames(colData) <- colnames(result_matrix)  # Ensure the row names match the column names of countData
                                dds <- DESeqDataSetFromMatrix(countData = result_matrix,
                                                              colData = colData,
                                                              design = ~ 1)
                                dds
                                
                                dds <- estimateSizeFactors(dds)
                                sizeFactors(dds)
                                
                                results <- vst(assay(dds)) #Normalization VST
                                result <- t(results)
                                
                        }
                }
        }
        
        return(result)
}

# Wrap your shiny app in the shinyApp() function
runShinyApp <- function(numeric_data, metadata) {
        print("TIME Data upload")
        #print(Sys.time())
        # Combine into a single reactive value for the app
        original_data <- list(numeric = numeric_data, metadata = metadata)
        combined_data <- reactiveVal(original_data)
        
        # Shiny UI layout with tabs
        ui <- fluidPage(
                titlePanel("Advanced Data Normalization and Visualization"),
                
                tabsetPanel(
                        tabPanel("Plot",
                                 sidebarLayout(
                                         sidebarPanel(
                                                 selectInput("x_var", "Select X Variable:", choices = NULL),
                                                 selectInput("y_var", "Select Y Variable:", choices = NULL),
                                                 actionButton("applyScatter", "Update Scatter Plot"),
                                                 br(),
                                                 br(),
                                                 tags$h3("Transformation methods"),
                                                 br(),
                                                 checkboxGroupInput("transform_methods", NULL,
                                                                    choices = list(
                                                                            "Logarithmic" = "Logarithmic", 
                                                                            "Median Centralization" = "Median Centralization", 
                                                                            "Mean Centralization" = "Mean Centralization"
                                                                    ),
                                                                    selected = "None"
                                                 ),
                                                 # Title for the normalization section
                                                 tags$h3("Normalization methods"),
                                                 br(),
                                                 checkboxGroupInput("norm_methods", NULL,
                                                                    choices = list(
                                                                            "Z-Score" = "Z-Score", 
                                                                            "MAD" = "MAD", 
                                                                            "Quantile Normalization" = "Quantile Normalization", 
                                                                            "Loess Normalization" = "Loess Normalization",
                                                                            "VST (Variance Stabilizing Transformation)" = "VST"
                                                                    ),
                                                                    selected = "None"
                                                 ),
                                                 selectInput("color_col", "Choose Column for Coloring:", choices = names(metadata), selected = "condition"),
                                                 actionButton("applyBtn", "Apply Normalization"),
                                                 actionButton("resetBtn", "Reset to Original Data"),  # Reset button
                                                 downloadButton("downloadPlot", "Download Current Plot"),
                                                 # Existing inputs...
                                                 actionButton("closeAppBtn", "Close App"),  # Add this line
                                                 tags$script(HTML("
        Shiny.addCustomMessageHandler('closeApp', function(message) {
            window.close();
        });
    "))
                                                 
                                         ),
                                         mainPanel(
                                                 plotOutput("dataPlot"),
                                                 textOutput("selectedNorms"),
                                                 verbatimTextOutput("warningMessage"),  # Output for warnings
                                                 plotOutput("pcaPlot"),
                                                 plotOutput("corrHeatmap"),
                                                 plotOutput("umapPlot"),
                                                 plotOutput("sumBoxplot")  # Boxplot output
                                         )
                                 )
                        ),
                        tabPanel("Data Table",
                                 DT::dataTableOutput("dataTable")
                        ),
                        
                        # Updated Outliers Tab
                        tabPanel("Outliers",
                                 sidebarLayout(
                                         sidebarPanel(
                                                 h3("Coefficient of Variation"),
                                                 sliderInput("cvThreshold", "CV Threshold (percentile)", min = 0, max = 100, value = 100),
                                                 actionButton("filterCVBtn", "Filter Variables with High CV"),
                                                 plotOutput("cvPlot"),
                                                 h3(" PCA-Based Outlier Detection"),
                                                 numericInput("pvalueThreshold", "Outlier Threshold (p-value)", value = 0.0001, min = 0, max = 1, step = 0.01),
                                                 actionButton("calcOutliersBtn", "Calculate Outliers"),
                                                 verbatimTextOutput("outliers"),
                                                 selectInput("removeSamples", "Select Samples to Remove:", choices = NULL, multiple = TRUE),
                                                 actionButton("removeBtn", "Remove Selected Samples")
                                         ),
                                         mainPanel(
                                                 DT::dataTableOutput("cvTable"),
                                                 plotOutput("mahalanobisPlot")
                                         )
                                 )
                        ),
                        
                        tabPanel("Download",
                                 downloadButton("downloadData", "Download Normalized Data")
                        )
                )
                ,div(
                        style = "background-color: #fffacd; padding: 15px; border-radius: 5px; margin-bottom: 20px;",
                        p(HTML("<b>Suggestions for the user</b>")),
                        p("Normalization is a critical step in the multiomics analysis. The methods provided here can guarantee a statistical normalization but cannot manage the instrument normalization (such as pooled samples in metabolomics). We strongly suggest performing the instrument normalization before using BiomiX to avoid biased results."),
                        p("Here is a brief explanation of the available normalization methods:"),
                        tags$ul(
                                tags$li(HTML("<b>Z-Score Normalization:</b> Standardizes the data by subtracting the <b>mean</b> and dividing by the <b>standard deviation</b>, resulting in a dataset with a mean of 0 and a standard deviation of 1.")),
                                tags$li(HTML("<b>MAD (Median Absolute Deviation):</b> Scales the data based on the <b>median absolute deviation</b>, making it robust to outliers.")),
                                tags$li(HTML("<b>Quantile Normalization:</b> Transforms the data to match a reference distribution, ensuring comparable statistical properties across samples.")),
                                tags$li(HTML("<b>Loess Normalization:</b> Applies local polynomial regression fitting (LOESS) to smooth the data and correct systematic biases. (Suggested for Linear data)")),
                                tags$li(HTML("<b>VST (Variance Stabilizing Transformation):</b> Stabilizes variance across different levels of intensity, reducing heteroscedasticity for more accurate analysis. <b>(Suggested for Transcriptomics)</b>"))
                        ),
                        p("Here is a brief explanation of the available transformation methods:"),
                        tags$ul(
                                tags$li(HTML("<b>Logarithmic:</b> Reduces skewness and compresses the scale of the data by taking the logarithm of each data point. This transformation stabilizes variance and makes relationships more linear, especially for data with exponential growth. <b>(Suggested for Metabolomics)</b>")),
                                tags$li(HTML("<b>Mean Centralization:</b> Centers the data by subtracting the mean value from each data point. This method is effective for normally distributed data but is sensitive to outliers.")),
                                tags$li(HTML("<b>Median Centralization:</b> Centers the data by subtracting the median value from each data point, reducing skewness and making the data more symmetric, especially in the presence of outliers."))
                        )
                ),
                )
        
        
        # Shiny server logic
        server <- function(input, output, session) {
                
                
                # Watch for the combined data
                observe({
                        plot_data <- combined_data()$numeric  # Assume combined_data() provides the numeric data
                        
                        # Check the number of columns (variables) in the numeric data
                        if (ncol(plot_data) > 1000) {
                                print("TIME variables")
                                #print(Sys.time())
                                # If more than 1000 columns, select the top 500 most variable genes
                                top_genes <- get_top_variable_genes(plot_data, top_n = 100)
                                #print(Sys.time())
                        } else {
                                # Otherwise, use all columns
                                top_genes <- colnames(plot_data)
                        }
                        
                        # Update the choices for x_var and y_var select inputs dynamically
                        updateSelectInput(session, "x_var", choices = top_genes)
                        updateSelectInput(session, "y_var", choices = top_genes)
                })
                
                
                # Apply normalization when button is clicked
                observeEvent(input$applyBtn, {
                        numeric_matrix <- combined_data()$numeric
                        output$warningMessage <- renderText({ NULL })  # Clear previous warnings
                        # Check for negative values before normalization
                        if (any(numeric_matrix < 0, na.rm = TRUE) && "Logarithmic" %in% input$norm_methods) {
                                output$warningMessage <- renderText("Warning: Cannot apply logarithmic transformation due to negative values in the dataset.")
                        } else {
                                set.seed(123)  # Add seed here for reproducibility
                                methods_join <- c(input$transform_methods, input$norm_methods)
                                normalized_data <- normalize_data(numeric_matrix, methods_join)
                                #print(head(normalized_data,n=5))
                                combined_data(list(numeric = normalized_data, metadata = combined_data()$metadata))
                        }
                })
                
                # Reset to original data when reset button is clicked
                observeEvent(input$resetBtn, {
                        combined_data(original_data)  # Reset to original unnormalized data
                        output$warningMessage <- renderText({ NULL })  # Clear warnings on reset
                })
                
                # Display the selected normalization methods
                output$selectedNorms <- renderText({
                        if(is.null(input$transform_methods) | is.null(input$norm_methods)){
                                paste("Selected Transformation: None"," ","/"," ", "Selected Normalizations: None", sep="\t")
                        }else{
                                paste(paste("Selected Transformation:", paste(input$transform_methods, collapse = ", "))," ","/"," ", paste("Selected Normalizations:", paste(input$norm_methods, collapse = ", ")), sep="\t")
                                }
                })
                
                # Close app when button is clicked
                observeEvent(input$closeAppBtn, {
                        current_data <- combined_data()
                        matrix_preview <- current_data$numeric
                        metadata_preview <- current_data$metadata
                        long_data = NULL
                        corr_matrix = NULL
                        umap_data  = NULL
                        umap_results = NULL
                        pca_data = NULL
                        pca = NULL
                        combined_data = NULL
                        updated_metadata  = NULL
                        pca_scores  = NULL
                        current_data = NULL
                        gc()
                        session$sendCustomMessage(type = 'closeApp', message = 'close')
                        # Stop the Shiny app gracefully
                        stopApp(list(matrix = matrix_preview, metadata = metadata_preview))
                })
                
                # Render scatter plot with conditional coloring
                output$dataPlot <- renderPlot({
                        set.seed(123)  # Reproducibility
                        #print("TIME Scatter")
                        #print(Sys.time())
                        plot_data <- combined_data()

                        req(input$x_var, input$y_var)

                        
                        # Use the selected variables for x and y axes
                        ggplot(data = cbind(plot_data$numeric, plot_data$metadata), aes_string(x = input$x_var, y = input$y_var, color = input$color_col)) +
                                geom_point() +
                                ggtitle("Scatter Plot of Selected Variables") +
                                theme_minimal()
                        #print(Sys.time())
                })
                
                # Render PCA plot
                output$pcaPlot <- renderPlot({
                        set.seed(123)  # Reproducibility
                        #print("TIME PCA")
                        #print(Sys.time())
                        plot_data <- combined_data()$numeric
                        # If there are more than 1000 genes, select the top 1000 based on variance
                        if (ncol(plot_data) > 1000) {
                                plot_data <- get_top_variable_features(plot_data, top_n = 1000)
                        }
                        
                        if (any(is.na(plot_data))) {
                                # Replace NA with the mean of the respective column
                                plot_data <- apply(plot_data, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
                        }
                        pca <- prcomp(plot_data, scale. = FALSE)
                        pca_data <- data.frame(pca$x, condition = combined_data()$metadata[[input$color_col]])
                        ggplot(pca_data, aes(PC1, PC2, color = condition)) +
                                geom_point() +
                                ggtitle("PCA Plot") +
                                theme_minimal()
                        #print(Sys.time())
                })
                
                # Render Correlation Heatmap
                output$corrHeatmap <- renderPlot({
                        #print("TIME HEAT")
                        #print(Sys.time())
                        plot_data <- combined_data()$numeric  # Assume combined_data() provides the numeric data
                        # Ensure the heatmap is only rendered if the number of variables (columns) is less than 1000
                        # Check if there are too many variables and subset to the top 500 most variable features
                        if (ncol(plot_data) > 100) {
                                # Show a notification to the user
                                showNotification("More than 100 variables detected. Heatmap displaying only the top 100 most variable features.", type = "warning", duration = 10)
                                #Sys.time()
                                plot_data <- get_top_variable_features(plot_data, top_n = 100)
                                #Sys.time()
                        }
                        set.seed(123)  # Reproducibility
                        
                        if (any(is.na(plot_data))) {
                                # Replace NA with the mean of the respective column
                                plot_data <- apply(plot_data, 2, function(x) ifelse(is.na(x), 0, x))
                        }
                        
                        corr_matrix <- cor(plot_data)
                        t<-pheatmap(corr_matrix, main = "Correlation Heatmap", cluster_rows = TRUE, cluster_cols = TRUE)
                        print(t)
                        #print(Sys.time())
                })
                
                # Render UMAP plot
                output$umapPlot <- renderPlot({
                        #print("TIME Umap")
                        #print(Sys.time())
                        set.seed(123)  # Reproducibility
                        plot_data <- combined_data()$numeric
                        
                        if (any(is.na(plot_data))) {
                                # Replace NA with the mean of the respective column
                                plot_data <- apply(plot_data, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
                        }
                        
                        umap_result <- umap(plot_data)
                        umap_data <- data.frame(UMAP1 = umap_result$layout[, 1], UMAP2 = umap_result$layout[, 2], condition = combined_data()$metadata[[input$color_col]])
                        ggplot(umap_data, aes(UMAP1, UMAP2, color = condition)) +
                                geom_point() +
                                ggtitle("UMAP Plot") +
                                theme_minimal()
                        #print(Sys.time())
                })
                
                
                # Render Boxplot of Sums for Each Sample
                # Render Boxplot of Values for Each Sample with individual points
                output$sumBoxplot <- renderPlot({
                        #print("TIME Boxplot")
                        #print(Sys.time())
                        plot_data <- as.data.frame(combined_data()$numeric)
                        #print(plot_data)
                        metadata <- combined_data()$metadata  # Fetch metadata
                        
                        # Convert numeric data to long format manually using base R
                        long_data <- stack(plot_data)
                        long_data$Sample <- rep(row.names(plot_data), each = ncol(plot_data))  # Add sample names to the data
                        
                        # Add grouping information from the metadata (e.g., condition, gender, etc.)
                        long_data$Group <- rep(metadata[[input$color_col]], each = ncol(plot_data))
                        
                        
                        # Create a boxplot where color depends on the grouping variable (Group)
                        ggplot(long_data, aes(x = Sample, y = values, fill = Group)) +
                                geom_boxplot(outlier.shape = NA) +  # Disable default outlier points
                                geom_jitter(width = 0.2, alpha = 0.5, aes(color = Group)) +  # Add jittered points with color based on group
                                ggtitle("Boxplot of Individual Variable Values per Sample, Colored by Group") +
                                theme_minimal() +
                                xlab("Samples") +
                                ylab("Values of Variables") +
                                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate sample labels for clarity
                                scale_fill_brewer(palette = "Set2") +  # Color palette for fill
                                scale_color_brewer(palette = "Set2")  # Same color palette for jitter points
                        #print(Sys.time())
                })
                
                
                # # Show data in an interactive table
                # output$dataTable <- DT::renderDataTable({
                #         combined_data()$numeric
                # })
                
                output$dataTable <- DT::renderDataTable({
                        #print(combined_data$numeric)
                        plot_data <- combined_data()$numeric
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
                
                # 1. Coefficient of Variation (CV) Filtering Logic
                observeEvent(input$filterCVBtn, {
                        plot_data <- combined_data()$numeric
                        
                        # Step 1: Calculate CV for each variable (column)
                        cv_values <- apply(plot_data, 2, function(x) {
                                sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100
                                
                        })
                        
                        # Step 2: Plot CV distribution with user-selected threshold
                        output$cvPlot <- renderPlot({
                                ggplot(data.frame(CV = cv_values), aes(CV)) +
                                        geom_density(fill = "lightblue") +
                                        geom_vline(xintercept = quantile(cv_values, probs = input$cvThreshold / 100), linetype = "dashed", color = "red") +
                                        labs(title = "Distribution of Coefficient of Variation (CV)", x = "CV (%)", y = "Density") +
                                        theme_minimal()
                        })
                        
                        # Step 3: Filter variables based on the CV threshold
                        
                        #the NaN generated including variables with the same value (0)
                        #will be assigned with a value 0.001 to avoid the NaN
                        cv_values[is.nan(cv_values)] <- 1
                        cv_threshold <- quantile(cv_values, probs = input$cvThreshold / 100, na.rm = TRUE)
                        filtered_data <- plot_data[, cv_values <= cv_threshold]
                        
                        # Update combined data with filtered variables
                        combined_data(list(numeric = filtered_data, metadata = combined_data()$metadata))
                        
                        # Display the updated CV values in a table
                        output$cvTable <- DT::renderDataTable({
                                data.frame(Variable = names(cv_values), CV = cv_values, Retained = cv_values <= cv_threshold)
                        })
                })
                
                # PCA-Based Outlier Detection Logic
                # This method calculates outliers based on distances in the PCA-transformed space.
                
                observeEvent(input$calcOutliersBtn, {
                        plot_data <- combined_data()$numeric
                        
                        if (ncol(plot_data) > 1000) {
                                plot_data <- get_top_variable_features(plot_data, top_n = 1000)
                        }
                        
                        #If missing data are present are computated
                        if (any(is.na(plot_data))) {
                                # Replace NA with the mean of the respective column
                                plot_data <- apply(plot_data, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
                        }
                        
                        # Step 1: Perform PCA on the data with scaling (each variable gets the same weight).
                        pca_result <- prcomp(plot_data, scale. = FALSE)  # PCA with scaling
                        pca_scores <- as.data.frame(pca_result$x)  # Extract scores in PCA space (principal component values)
                        
                        # Step 2: Calculate Euclidean distances from the centroid (origin) in PCA space.
                        # Each sample's distance from the origin is computed as sqrt(sum of squares of the PC scores).
                        pca_distances <- sqrt(rowSums(pca_scores^2))  # Calculate distance for each sample
                        
                        # Step 3: Identify outliers based on an empirical threshold.
                        # We calculate the threshold distance as a quantile of the empirical distance distribution.
                        # For example, if pvalueThreshold = 0.05, the threshold corresponds to the 95th percentile of distances.
                        threshold_distance <- quantile(pca_distances, probs = 1 - input$pvalueThreshold)
                        
                        # Identify outliers as samples with distances greater than the threshold.
                        outliers <- rownames(plot_data)[pca_distances > threshold_distance]
                        
                        # Step 4: Visualize the distances and highlight outliers.
                        output$mahalanobisPlot <- renderPlot({
                                plot(pca_distances, pch = 19, col = ifelse(pca_distances > threshold_distance, "red", "black"),
                                     main = "PCA Distance Plot", xlab = "Sample Index", ylab = "PCA Distance")
                                abline(h = threshold_distance, col = "blue", lty = 2)  # Threshold line
                        })
                        
                        # Output the identified outliers (display sample names).
                        output$outlierText <- renderText({
                                paste("Outliers identified:", paste(outliers, collapse = ", "))
                        })
                        
                        
                        # Update select input for removing outliers
                        updateSelectInput(session, "removeSamples", choices = rownames(plot_data), selected = outliers)
                })
                
                # Remove selected samples from the dataset
                observeEvent(input$removeBtn, {
                        plot_data <- combined_data()$numeric
                        metadata <- combined_data()$metadata
                        
                        # Remove the selected samples
                        samples_to_remove <- input$removeSamples
                        updated_data <- plot_data[!rownames(plot_data) %in% samples_to_remove, ]
                        updated_metadata <- metadata[!metadata$ID %in% samples_to_remove, ]
                        print(dim(updated_metadata))
                        print(dim(updated_data))
                        
                        # Update the combined data
                        combined_data(list(numeric = updated_data, metadata = updated_metadata))
                        
                        # Clear outliers and reset select input
                        output$outliers <- renderText({ NULL })
                        updateSelectInput(session, "removeSamples", choices = rownames(updated_data), selected = NULL)
                })
                
                
                
                # Download handler to download normalized data as CSV
                output$downloadData <- downloadHandler(
                        filename = function() { paste("normalized_data", Sys.Date(), ".tsv", sep = "") },
                        content = function(file) {
                                plot_data <- combined_data()$numeric
                                Samples <- colnames(plot_data)
                                print(head(plot_data))
                                variables <- rownames(plot_data)
                                plot_data <- as.data.frame(t(plot_data))
                                plot_data <- cbind(Samples, plot_data)
                                colnames(plot_data) <- c("ID", variables)
                                write.table(plot_data, file, row.names = FALSE, quote = FALSE, sep = "\t")
                        }
                )
                
                # Download handler for the current plot
                output$downloadPlot <- downloadHandler(
                        filename = function() { paste("plot", Sys.Date(), ".pdf", sep = "") },
                        content = function(file) {
                                pdf(file, width = 18, height = 14)  # Open a PDF device
                                
                                # Scatter Plot
                                set.seed(123)  # Reproducibility
                                plot_data <- combined_data()
                                
                                # Use the selected variables for x and y axes
                                ggplot(data = cbind(plot_data$numeric, plot_data$metadata), aes_string(x = input$x_var, y = input$y_var, color = input$color_col)) +
                                        geom_point() +
                                        ggtitle("Scatter Plot of Selected Variables") +
                                        theme_minimal()
                                print(last_plot())  # Print to the PDF
                                
                                # PCA Plot
                                set.seed(123)  # Reproducibility
                                plot_data <- combined_data()$numeric
                                # If there are more than 1000 genes, select the top 1000 based on variance
                                if (ncol(plot_data) > 1000) {
                                        plot_data <- get_top_variable_features(plot_data, top_n = 1000)
                                }
                                pca <- prcomp(plot_data, scale. = FALSE)
                                pca_data <- data.frame(pca$x, condition = combined_data()$metadata[[input$color_col]])
                                ggplot(pca_data, aes(PC1, PC2, color = condition)) +
                                        geom_point() +
                                        ggtitle("PCA Plot") +
                                        theme_minimal()
                                print(last_plot())  # Print to the PDF
                                
                                # Correlation Heatmap
                                plot_data <- combined_data()$numeric  # Assume combined_data() provides the numeric data
                                # Ensure the heatmap is only rendered if the number of variables (columns) is less than 1000
                                # Check if there are too many variables and subset to the top 500 most variable features
                                if (ncol(plot_data) > 100) {
                                        # Show a notification to the user
                                        showNotification("More than 100 variables detected. Heatmap displaying only the top 100 most variable features.", type = "warning", duration = 10)
                                        #Sys.time()
                                        plot_data <- get_top_variable_features(plot_data, top_n = 100)
                                        #Sys.time()
                                }
                                set.seed(123)  # Reproducibility
                                corr_matrix <- cor(plot_data)
                                pheatmap(corr_matrix, main = "Correlation Heatmap", cluster_rows = TRUE, cluster_cols = TRUE)
                                print(last_plot())  # Print to the PDF
                                
                                # UMAP Plot
                                set.seed(123)  # Reproducibility
                                plot_data <- combined_data()$numeric
                                umap_result <- umap(plot_data)
                                umap_data <- data.frame(UMAP1 = umap_result$layout[, 1], UMAP2 = umap_result$layout[, 2], condition = combined_data()$metadata[[input$color_col]])
                                ggplot(umap_data, aes(UMAP1, UMAP2, color = condition)) +
                                        geom_point() +
                                        ggtitle("UMAP Plot") +
                                        theme_minimal()
                                print(last_plot())  # Print to the PDF
                                
                                
                                #Boxplot values for each sample
                                plot_data <- as.data.frame(combined_data()$numeric)
                                #print(plot_data)
                                metadata <- combined_data()$metadata  # Fetch metadata
                                
                                # Convert numeric data to long format manually using base R
                                long_data <- stack(plot_data)
                                long_data$Sample <- rep(row.names(plot_data), each = ncol(plot_data))  # Add sample names to the data
                                
                                # Add grouping information from the metadata (e.g., condition, gender, etc.)
                                long_data$Group <- rep(metadata[[input$color_col]], each = ncol(plot_data))
                                
                                
                                # Create a boxplot where color depends on the grouping variable (Group)
                                ggplot(long_data, aes(x = Sample, y = values, fill = Group)) +
                                        geom_boxplot(outlier.shape = NA) +  # Disable default outlier points
                                        geom_jitter(width = 0.2, alpha = 0.5, aes(color = Group)) +  # Add jittered points with color based on group
                                        ggtitle("Boxplot of Individual Variable Values per Sample, Colored by Group") +
                                        theme_minimal() +
                                        xlab("Samples") +
                                        ylab("Values of Variables") +
                                        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate sample labels for clarity
                                        scale_fill_brewer(palette = "Set2") +  # Color palette for fill
                                        scale_color_brewer(palette = "Set2")  # Same color palette for jitter pointsscale_color_brewer(palette = "Set2")  # Same color palette for jitter points
                                print(last_plot())  # Print to the PDF
                                
                                dev.off()  # Close the PDF device
                        }
                )
                
        }
        
        # Return the shiny app object
        shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
}



# Now directly call runShinyApp() to run the app and return the matrix and metadata
#result <- runShinyApp(numeric_data, metadata)

# Extract the matrix and metadata after the app closes
#matrix_data <- result$matrix
#metadata_data <- result$metadata
