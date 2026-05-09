# BF 591 Final Project - Huntington's Disease RNA-Seq Explorer
# Dataset: GSE64810 (20 HD vs 49 neurologically normal controls)
# BA9 prefrontal cortex, bulk RNA-Seq

library(shiny)
library(tidyverse)
library(DT)
library(pheatmap)
library(patchwork)
library(colourpicker)

# raise Shiny's default 5MB file upload limit to 100MB for the counts CSV
options(shiny.maxRequestSize = 100 * 1024^2)

# helper: read a CSV or TSV and check that required columns are present
validate_and_read <- function(file_input, required_cols) {
  req(file_input)
  
  # check that the extension is csv or tsv
  ext <- tolower(sub(".*\\.", "", file_input$name))
  validate(need(ext %in% c("csv", "tsv"),
                "File must be a .csv or .tsv"))
  
  # read using the matching parser
  df <- if (ext == "csv") {
    read_csv(file_input$datapath, show_col_types = FALSE, progress = FALSE)
  } else {
    read_tsv(file_input$datapath, show_col_types = FALSE, progress = FALSE)
  }
  
  # check that required columns exist
  validate(need(all(required_cols %in% colnames(df)),
                paste("File must contain columns:",
                      paste(required_cols, collapse = ", "))))
  df
}

ui <- fluidPage(
  titlePanel("Huntington's Disease RNA-Seq Explorer (GSE64810)"),
  
  tabsetPanel(
    # --- Tab 1: Sample Information ---
    tabPanel("Samples",
             sidebarLayout(
               sidebarPanel(
                 helpText("Explore the sample metadata for this experiment: ",
                          "column summaries, a sortable data table, ",
                          "and distribution plots of continuous variables."),
                 fileInput("sample_file", "Upload sample info CSV",
                           accept = c(".csv", ".tsv"))
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary", tableOutput("samples_summary")),
                   tabPanel("Table",   DTOutput("samples_table")),
                   tabPanel("Plots",
                            uiOutput("samples_plot_col_ui"),
                            plotOutput("samples_plot"))
                 )
               )
             )
    ),
    
    # --- Tab 2: Counts Matrix ---
    tabPanel("Counts",
             sidebarLayout(
               sidebarPanel(
                 helpText("Explore the normalized counts matrix. Use the sliders ",
                          "to filter genes by variance and by how many samples ",
                          "they are detected in. Both filters must pass."),
                 fileInput("counts_file", "Upload normalized counts CSV",
                           accept = c(".csv", ".tsv")),
                 sliderInput("var_pct", "Minimum variance percentile:",
                             min = 0, max = 100, value = 50, step = 1),
                 sliderInput("nonzero_n", "Minimum non-zero samples:",
                             min = 0, max = 69, value = 10, step = 1)
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Filter Summary", verbatimTextOutput("counts_summary")),
                   tabPanel("Diagnostic Plots", plotOutput("counts_diag", height = "600px")),
                   tabPanel("Heatmap",
                            checkboxInput("heatmap_log", "Log2-transform counts", value = TRUE),
                            plotOutput("counts_heatmap", height = "700px")),
                   tabPanel("PCA",
                            radioButtons("pca_mode", "Plot type:",
                                         choices = c("Scatter (pick 2 PCs)" = "scatter",
                                                     "Beeswarm (top N PCs)" = "beeswarm"),
                                         inline = TRUE),
                            uiOutput("pca_controls_ui"),
                            plotOutput("counts_pca", height = "600px"))
                 )
               )
             )
    ),
    
    # --- Tab 3: Differential Expression ---
    tabPanel("DE",
             sidebarLayout(
               sidebarPanel(
                 helpText("Explore differential expression results. Pick axis ",
                          "variables, adjust the significance threshold, and ",
                          "customize the colors. The table is searchable."),
                 fileInput("de_file", "Upload DE results CSV",
                           accept = c(".csv", ".tsv")),
                 radioButtons("x_axis", "X-axis variable",
                              choices = c("baseMean", "HD.mean", "Control.mean",
                                          "log2FoldChange", "lfcSE", "stat",
                                          "pvalue", "padj"),
                              selected = "log2FoldChange"),
                 radioButtons("y_axis", "Y-axis variable",
                              choices = c("baseMean", "HD.mean", "Control.mean",
                                          "log2FoldChange", "lfcSE", "stat",
                                          "pvalue", "padj"),
                              selected = "padj"),
                 colourInput("base_color", "Base point color", value = "grey"),
                 colourInput("hl_color",   "Highlight point color", value = "red"),
                 sliderInput("padj_slider", "Magnitude of adjusted p-value (10^X):",
                             min = -50, max = 0, value = -5, step = 1)
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Volcano Plot", plotOutput("volcano_plot", height = "600px")),
                   tabPanel("Table", DTOutput("de_table"))
                 )
               )
             )
    ),
    
    # --- Tab 4: Individual Gene Visualization ---
    tabPanel("Gene Viz",
             sidebarLayout(
               sidebarPanel(
                 helpText("Visualize a single gene's expression across groups. ",
                          "Upload both files, choose a grouping variable, pick a ",
                          "gene, select a plot type, and click Plot."),
                 fileInput("gv_counts", "Upload normalized counts CSV",
                           accept = c(".csv", ".tsv")),
                 fileInput("gv_meta",   "Upload sample info CSV",
                           accept = c(".csv", ".tsv")),
                 uiOutput("gv_cat_ui"),
                 uiOutput("gv_gene_ui"),
                 radioButtons("gv_plot_type", "Plot type:",
                              choices = c("Bar" = "bar",
                                          "Boxplot" = "box",
                                          "Violin" = "violin",
                                          "Beeswarm" = "beeswarm"),
                              selected = "box",
                              inline = TRUE),
                 actionButton("gv_go", "Plot", icon = icon("chart-simple"))
               ),
               mainPanel(plotOutput("gv_plot", height = "600px"))
             )
    )
  )
)

server <- function(input, output, session) {
  
  # --- Tab 1: Sample Information ---
  
  # read the sample info file; must have a sample_id column
  sample_data <- reactive({
    validate_and_read(input$sample_file, required_cols = "sample_id")
  })
  
  # summary table shows type + mean/sd for numeric, distinct values for categorical
  output$samples_summary <- renderTable({
    df <- sample_data()
    tibble(
      `Column Name` = colnames(df),
      Type          = sapply(df, function(x) class(x)[1]),
      `Mean (sd) or Distinct Values` = sapply(df, function(x) {
        if (is.numeric(x)) {
          paste0(round(mean(x, na.rm = TRUE), 2),
                 " (+/- ", round(sd(x, na.rm = TRUE), 2), ")")
        } else {
          paste(unique(x), collapse = ", ")
        }
      })
    )
  })
  
  output$samples_table <- renderDT({
    datatable(sample_data(), options = list(pageLength = 10))
  })
  
  # dropdown is populated with only the numeric columns of the uploaded file
  output$samples_plot_col_ui <- renderUI({
    df <- sample_data()
    numeric_cols <- colnames(df)[sapply(df, is.numeric)]
    selectInput("samples_plot_col", "Choose a continuous variable:",
                choices = numeric_cols)
  })
  
  output$samples_plot <- renderPlot({
    req(input$samples_plot_col)
    df <- sample_data()
    ggplot(df, aes(x = .data[[input$samples_plot_col]])) +
      geom_histogram(bins = 30, fill = "steelblue", color = "white") +
      labs(title = paste("Distribution of", input$samples_plot_col),
           x = input$samples_plot_col, y = "Count") +
      theme_minimal()
  })
  
  # --- Tab 2: Counts Matrix ---
  
  # must have a 'gene' column; the rest are sample columns
  counts_data <- reactive({
    validate_and_read(input$counts_file, required_cols = "gene")
  })
  
  # compute per-gene stats once and cache
  gene_stats <- reactive({
    df <- counts_data()
    mat <- as.matrix(df[, -1])
    rownames(mat) <- df$gene
    tibble(
      gene       = df$gene,
      median_cnt = apply(mat, 1, median),
      variance   = apply(mat, 1, var),
      n_nonzero  = apply(mat, 1, function(x) sum(x != 0)),
      n_zero     = apply(mat, 1, function(x) sum(x == 0))
    )
  })
  
  # TRUE/FALSE for whether each gene passes both slider filters
  passing <- reactive({
    stats <- gene_stats()
    var_cutoff <- quantile(stats$variance, probs = input$var_pct / 100, na.rm = TRUE)
    (stats$variance >= var_cutoff) & (stats$n_nonzero >= input$nonzero_n)
  })
  
  # filtered counts matrix used for heatmap and PCA
  filtered_matrix <- reactive({
    df <- counts_data()
    mat <- as.matrix(df[, -1])
    rownames(mat) <- df$gene
    mat[passing(), , drop = FALSE]
  })
  
  output$counts_summary <- renderPrint({
    df <- counts_data()
    n_samples <- ncol(df) - 1
    n_genes   <- nrow(df)
    n_pass    <- sum(passing())
    n_fail    <- n_genes - n_pass
    cat("Number of samples:", n_samples, "\n")
    cat("Total number of genes:", n_genes, "\n")
    cat("Genes passing current filter:", n_pass,
        paste0("(", round(100 * n_pass / n_genes, 2), "%)"), "\n")
    cat("Genes not passing current filter:", n_fail,
        paste0("(", round(100 * n_fail / n_genes, 2), "%)"), "\n")
  })
  
  output$counts_diag <- renderPlot({
    stats <- gene_stats()
    stats$passing <- passing()
    
    p1 <- ggplot(stats, aes(x = median_cnt, y = variance, color = passing)) +
      geom_point(alpha = 0.6, size = 0.7) +
      scale_x_log10() + scale_y_log10() +
      scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgray")) +
      labs(title = "Median count vs Variance",
           x = "Median count (log10)", y = "Variance (log10)") +
      theme_minimal()
    
    p2 <- ggplot(stats, aes(x = median_cnt, y = n_zero, color = passing)) +
      geom_point(alpha = 0.6, size = 0.7) +
      scale_x_log10() +
      scale_color_manual(values = c("TRUE" = "black", "FALSE" = "lightgray")) +
      labs(title = "Median count vs Number of zeros",
           x = "Median count (log10)", y = "Number of zeros") +
      theme_minimal()
    
    p1 + p2 + plot_layout(ncol = 2)
  })
  
  output$counts_heatmap <- renderPlot({
    mat <- filtered_matrix()
    req(nrow(mat) > 1)
    if (input$heatmap_log) {
      mat <- log2(mat + 1)
    }
    pheatmap(mat,
             scale         = "row",
             show_rownames = FALSE,
             show_colnames = TRUE,
             color         = colorRampPalette(c("navy", "white", "firebrick3"))(100))
  })
  
  # prcomp expects samples as rows and genes as columns, so we transpose
  pca_result <- reactive({
    mat <- filtered_matrix()
    req(nrow(mat) > 1)
    prcomp(t(mat), scale. = TRUE)
  })
  
  # PCA sub-tab controls switch between two modes
  output$pca_controls_ui <- renderUI({
    req(pca_result())
    n_pcs <- ncol(pca_result()$x)
    if (input$pca_mode == "scatter") {
      tagList(
        selectInput("pc_x", "X-axis PC:", choices = 1:n_pcs, selected = 1),
        selectInput("pc_y", "Y-axis PC:", choices = 1:n_pcs, selected = 2)
      )
    } else {
      sliderInput("top_n_pcs", "Number of top PCs:",
                  min = 2, max = min(10, n_pcs), value = 5, step = 1)
    }
  })
  
  output$counts_pca <- renderPlot({
    pca <- pca_result()
    var_explained <- summary(pca)$importance["Proportion of Variance", ]
    
    if (input$pca_mode == "scatter") {
      req(input$pc_x, input$pc_y)
      i <- as.numeric(input$pc_x)
      j <- as.numeric(input$pc_y)
      scores <- as.data.frame(pca$x)
      ggplot(scores, aes(x = .data[[paste0("PC", i)]],
                         y = .data[[paste0("PC", j)]])) +
        geom_point(size = 2) +
        labs(x = paste0("PC", i, " (", round(100 * var_explained[i], 1), "%)"),
             y = paste0("PC", j, " (", round(100 * var_explained[j], 1), "%)"),
             title = paste0("PC", i, " vs PC", j)) +
        theme_minimal()
    } else {
      req(input$top_n_pcs)
      n <- input$top_n_pcs
      scores <- as.data.frame(pca$x[, 1:n, drop = FALSE])
      scores$sample <- rownames(scores)
      long <- pivot_longer(scores, cols = -sample, names_to = "PC", values_to = "score")
      long$PC <- factor(long$PC, levels = paste0("PC", 1:n))
      ggplot(long, aes(x = PC, y = score)) +
        geom_jitter(position = position_jitter(width = 0.15, seed = 42),
                    size = 2, alpha = 0.7) +
        labs(title = paste("Top", n, "PCs"), x = "", y = "Score") +
        theme_minimal()
    }
  })
  
  # --- Tab 3: Differential Expression ---
  
  # must have padj (used for significance threshold)
  de_data <- reactive({
    validate_and_read(input$de_file, required_cols = "padj")
  })
  
  # points are colored 'highlight' if padj below 10^slider, else 'base'
  output$volcano_plot <- renderPlot({
    df <- de_data()
    threshold <- 10 ^ input$padj_slider
    df$sig <- ifelse(df$padj < threshold, "highlight", "base")
    
    ggplot(df, aes(x = .data[[input$x_axis]],
                   y = -log10(.data[[input$y_axis]]),
                   color = sig)) +
      geom_point(alpha = 0.7, size = 1) +
      scale_color_manual(values = c("base"      = input$base_color,
                                    "highlight" = input$hl_color)) +
      labs(x = input$x_axis,
           y = paste0("-log10(", input$y_axis, ")"),
           title = "Volcano Plot",
           color = "Significance") +
      theme_minimal()
  })
  
  # table also filters by padj threshold; DT provides a built-in search box
  output$de_table <- renderDT({
    df <- de_data()
    threshold <- 10 ^ input$padj_slider
    df_filtered <- df %>% filter(padj < threshold)
    df_filtered$pvalue <- formatC(df_filtered$pvalue, format = "e", digits = 5)
    df_filtered$padj   <- formatC(df_filtered$padj,   format = "e", digits = 5)
    datatable(df_filtered, options = list(pageLength = 10))
  })
  
  # --- Tab 4: Individual Gene Visualization ---
  
  gv_counts_data <- reactive({
    validate_and_read(input$gv_counts, required_cols = "gene")
  })
  
  gv_meta_data <- reactive({
    validate_and_read(input$gv_meta, required_cols = "sample_id")
  })
  
  # only non-numeric columns are shown as grouping options
  output$gv_cat_ui <- renderUI({
    df <- gv_meta_data()
    cat_cols <- colnames(df)[sapply(df, function(x) !is.numeric(x))]
    selectInput("gv_cat", "Group by categorical variable:",
                choices = cat_cols,
                selected = if ("diagnosis" %in% cat_cols) "diagnosis" else cat_cols[1])
  })
  
  # selectInput with many choices auto-enables search on the client side
  output$gv_gene_ui <- renderUI({
    df <- gv_counts_data()
    selectInput("gv_gene", "Choose a gene (type to search):",
                choices = df$gene,
                selected = df$gene[1])
  })
  
  # eventReactive: capture ALL plot parameters at button click.
  # this prevents the plot from updating until Plot is clicked.
  gv_plot_spec <- eventReactive(input$gv_go, {
    counts <- gv_counts_data()
    meta   <- gv_meta_data()
    req(input$gv_gene, input$gv_cat)
    
    gene_row <- counts %>% filter(gene == input$gv_gene)
    req(nrow(gene_row) == 1)
    
    # reshape the one gene's row into long format and join with sample metadata
    gene_long <- gene_row %>%
      pivot_longer(cols = -gene, names_to = "sample_id", values_to = "count")
    plot_data <- gene_long %>% left_join(meta, by = "sample_id")
    
    list(
      data = plot_data,
      gene = input$gv_gene,
      cat  = input$gv_cat,
      type = input$gv_plot_type
    )
  })
  
  output$gv_plot <- renderPlot({
    spec      <- gv_plot_spec()
    df        <- spec$data
    cat_var   <- spec$cat
    plot_type <- spec$type
    
    base_labs <- labs(title = paste("Expression of", spec$gene),
                      x = cat_var, y = "Normalized count")
    
    if (plot_type == "bar") {
      # bar = mean +/- standard error per group
      summary_df <- df %>%
        group_by(.data[[cat_var]]) %>%
        summarise(mean_count = mean(count, na.rm = TRUE),
                  se = sd(count, na.rm = TRUE) / sqrt(n()),
                  .groups = "drop")
      ggplot(summary_df, aes(x = .data[[cat_var]], y = mean_count)) +
        geom_col(fill = "steelblue") +
        geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se),
                      width = 0.2) +
        base_labs + theme_minimal()
    } else {
      p <- ggplot(df, aes(x = .data[[cat_var]], y = count))
      if (plot_type == "box") {
        p <- p + geom_boxplot(fill = "lightblue")
      } else if (plot_type == "violin") {
        p <- p + geom_violin(fill = "lightgreen")
      } else if (plot_type == "beeswarm") {
        # seed makes jitter positions reproducible across re-renders
        p <- p + geom_jitter(position = position_jitter(width = 0.2, seed = 42),
                             size = 2, alpha = 0.7)
      }
      p + base_labs + theme_minimal()
    }
  })
}

shinyApp(ui, server)