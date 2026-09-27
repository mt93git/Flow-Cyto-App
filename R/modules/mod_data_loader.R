# R/modules/mod_data_loader.R

# --- Parsing Helper Function ---
parse_legacy_facs <- function(file_path, population_mapping = NULL) {
  # Read all data
  raw_df <- read_excel(file_path)
  
  # Ensure columns exist
  if (!all(c("Name", "Statistic", "#Cells") %in% names(raw_df))) {
    stop("Invalid file format. Expected columns: Name, Statistic, #Cells")
  }
  
  # Initialize storage
  results_mfi <- list()
  results_counts <- list()
  
  current_sample <- NA_character_
  current_sample_id <- NA_character_
  current_population <- NA_character_
  
  # Regex patterns
  sample_pattern <- "^(.*_(\\d+)\\.fcs)$" 
  mfi_pattern <- "Geometric Mean\\s*:\\s*(.*)"
  
  cat("Starting parse of", nrow(raw_df), "rows...\n")
  
  for (i in 1:nrow(raw_df)) {
    row_name <- raw_df$Name[i]
    row_stat <- raw_df$Statistic[i]
    row_count <- raw_df$`#Cells`[i]
    
    if (is.na(row_name)) next
    
    # Check for New Sample Block
    if (grepl("\\.fcs$", row_name) && !grepl("/", row_name)) {
      current_sample <- row_name
      
      neu_match <- str_match(row_name, "BM_Neu_(\\d+)")
      if (!is.na(neu_match[1,2])) {
          current_sample_id <- as.numeric(neu_match[1,2])
      } else {
          match <- str_match(row_name, sample_pattern)
          if (!is.na(match[1,3])) {
            current_sample_id <- as.numeric(match[1,3])
          } else {
            fallback_match <- str_extract(row_name, "\\d+(?=\\.fcs)")
            if (!is.na(fallback_match)) {
                current_sample_id <- as.numeric(fallback_match)
            } else {
                current_sample_id <- NA
            }
          }
      }
      next
    }
    
    if (is.na(current_sample)) next
    
    # 1. Check for MFI Line
    if (grepl("Geometric Mean", row_name)) {
      clean_name <- sub("Geometric Mean\\s*:\\s*", "", row_name)
      clean_name <- sub("\\s*=.*", "", clean_name)
      
      results_mfi[[length(results_mfi) + 1]] <- data.frame(
        Sample = current_sample,
        SampleID = current_sample_id,
        Population = current_population,
        Parameter = clean_name,
        Value = row_stat,
        Type = "MFI",
        stringsAsFactors = FALSE
      )
    } 
    # 2. Check for Population Count Line
    else if (!is.na(row_count)) {
        parts <- strsplit(row_name, "/")[[1]]
        pop_name <- tail(parts, 1)
        
        if (!is.null(population_mapping) && pop_name %in% names(population_mapping)) {
            pop_name <- population_mapping[[pop_name]]
        }
        
        if (pop_name == "immune") pop_name <- "CD45+"
        
        current_population <- pop_name
        
        results_counts[[length(results_counts) + 1]] <- data.frame(
          Sample = current_sample,
          SampleID = current_sample_id,
          Population = pop_name,
          Count = row_count,
          Type = "Count",
          stringsAsFactors = FALSE
        )
    }
  }
  
  df_mfi <- do.call(rbind, results_mfi)
  df_counts <- do.call(rbind, results_counts)
  
  return(list(mfi = df_mfi, counts = df_counts))
}

# --- Core Normalization & Integration Pipeline ---
execute_flow_pipeline <- function(facs_path, meta_path, marker_path, bead_pattern = "CountBright|Bead", add_log = message) {
  add_log("--- Starting Analytical Pipeline ---")
  
  # 1. Parse FACS
  add_log("Parsing FACS file...")
  poly_map <- list(
      "immune" = "CD45+",
      "All cells" = "Cells"
  )
  parsed <- parse_legacy_facs(facs_path, population_mapping = poly_map)
  
  unique_pops <- unique(parsed$counts$Population)
  add_log(paste("Topology Report: Found", length(unique_pops), "unique populations:"))
  for(p in unique_pops) {
       add_log(paste("  [POP] ->", p))
  }
  add_log(paste("Parsed", nrow(parsed$mfi), "MFI rows and", nrow(parsed$counts), "count rows."))
  
  # 2. Load Metadata
  add_log("Loading Metadata...")
  meta <- read_excel(meta_path)
  if("sample_no" %in% names(meta)) meta$sample_no <- as.numeric(meta$sample_no)
  
  # 3. Load Marker Map (with robust BOM handling)
  add_log("Loading Marker Map...")
  markers <- read.csv(marker_path, fileEncoding = "UTF-8-BOM", check.names = FALSE)
  colnames(markers) <- gsub("^[^a-zA-Z0-9]+", "", colnames(markers))
  
  # 4. Join MFI with Metadata
  df_mfi <- parsed$mfi
  df_combined <- df_mfi %>%
    inner_join(meta, by = c("SampleID" = "sample_no"))
  
  # 5. Join with Marker Map
  df_combined$CleanFluor <- gsub("Comp-", "", df_combined$Parameter)
  df_combined$CleanFluor <- trimws(df_combined$CleanFluor)
  df_final <- df_combined %>%
    left_join(markers, by = c("CleanFluor" = "fluor"))
  
  # 6. Absolute Counts Normalization
  counts_df <- parsed$counts
  bead_raw_subset <- counts_df %>%
    filter(grepl(bead_pattern, Population, ignore.case=TRUE))
  
  detected_names <- unique(bead_raw_subset$Population)
  add_log(paste("Normalization: Found", length(detected_names), "unique Reference Standard labels"))
  
  bead_df <- bead_raw_subset %>%
    select(SampleID, BeadCount = Count) %>%
    group_by(SampleID) %>%
    summarise(BeadCount = mean(BeadCount, na.rm=TRUE), .groups="drop")
  
  counts_normalized <- counts_df %>%
    left_join(bead_df, by="SampleID") %>%
    left_join(meta, by=c("SampleID" = "sample_no"))
  
  if (all(c("volume_total", "volume_sample_Neu", "beads_input") %in% names(counts_normalized))) {
       counts_normalized <- counts_normalized %>%
         mutate(
           AbsCount = (Count / BeadCount) * beads_input * (volume_total / volume_sample_Neu)
         )
       add_log("Absolute Counts Successfully Normalized using Reference Standards.")
  } else {
       counts_normalized$AbsCount <- NA
       add_log("Notice: Missing bead spike volume metadata; skipping absolute cell counts.")
  }
  
  # 7. Transformation & Z-Score Standardization
  df_final$Value_Transformed <- asinh(df_final$Value / 150)
  
  df_final <- df_final %>%
    group_by(marker, Population) %>%
    mutate(Z_Score = (Value_Transformed - mean(Value_Transformed, na.rm=TRUE)) / sd(Value_Transformed, na.rm=TRUE)) %>%
    ungroup()
  
  add_log("Processing Complete! Data ready for interactive Heatmap & Clustering.")
  return(list(data = df_final, count_data = counts_normalized))
}

mod_data_loader_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(
      class = "card p-3 mb-3 bg-light border-primary shadow-sm",
      h6(class="text-primary fw-bold mb-2", "⚡ Quick Start Demo"),
      p(class="small text-muted mb-2", "Click below to evaluate the complete pipeline out-of-the-box using the built-in LPS vs SteadyState dataset."),
      actionButton(ns("load_demo_btn"), "📊 Load Built-in Demo Dataset", class = "btn-success w-100 fw-bold")
    ),
    
    hr(),
    h6(class="fw-bold", "📁 Custom Data Upload"),
    fileInput(ns("facs_file"), "Upload FACS Data (.xls)", accept = c(".xls")),
    fileInput(ns("meta_file"), "Upload Metadata (.xlsx)", accept = c(".xlsx")),
    fileInput(ns("marker_file"), "Upload Marker Map (.csv)", accept = c(".csv")),
    
    hr(),
    textInput(ns("bead_pattern"), "Bead Population Regex", value = "CountBright|Bead", width="100%"),
    actionButton(ns("process_btn"), "Process & Normalize Custom Data", class = "btn-primary w-100"),
    div(class="mt-2", downloadButton(ns("download_counts"), "Download Abs Counts (CSV)", class = "btn-sm w-100 btn-info")),
    hr(),
    verbatimTextOutput(ns("status_log"))
  )
}

mod_data_loader_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    processed_data <- reactiveVal(NULL)
    log_text <- reactiveVal("Waiting for input or click 'Load Built-in Demo Dataset'...")
    
    add_log <- function(msg) {
      log_text(paste(log_text(), msg, sep = "\n"))
    }
    
    output$status_log <- renderText({ log_text() })
    
    # --- Handler 1: Built-in Demo Dataset ---
    observeEvent(input$load_demo_btn, {
      tryCatch({
        base_dir <- getwd()
        demo_facs <- file.path(base_dir, "demo_data", "demo_facs_data.xls")
        demo_meta <- file.path(base_dir, "demo_data", "demo_metadata.xlsx")
        demo_marker <- file.path(base_dir, "demo_data", "demo_marker_map.csv")
        
        if(!file.exists(demo_facs)) {
          demo_facs <- file.path("..", "demo_data", "demo_facs_data.xls")
          demo_meta <- file.path("..", "demo_data", "demo_metadata.xlsx")
          demo_marker <- file.path("..", "demo_data", "demo_marker_map.csv")
        }
        
        if(!file.exists(demo_facs)) {
          stop("Demo files not found in demo_data/")
        }
        
        add_log("--- Loading Built-in Demo Dataset ---")
        res <- execute_flow_pipeline(
          facs_path = demo_facs,
          meta_path = demo_meta,
          marker_path = demo_marker,
          bead_pattern = input$bead_pattern,
          add_log = add_log
        )
        processed_data(res)
        showNotification("Demo Dataset loaded successfully!", type = "message")
      }, error = function(e) {
        add_log(paste("ERROR:", e$message))
        showNotification(paste("Error loading demo data:", e$message), type = "error")
      })
    })
    
    # --- Handler 2: Custom Uploaded Data ---
    observeEvent(input$process_btn, {
      req(input$facs_file, input$meta_file, input$marker_file)
      
      tryCatch({
        res <- execute_flow_pipeline(
          facs_path = input$facs_file$datapath,
          meta_path = input$meta_file$datapath,
          marker_path = input$marker_file$datapath,
          bead_pattern = input$bead_pattern,
          add_log = add_log
        )
        processed_data(res)
        showNotification("Custom Dataset processed successfully!", type = "message")
      }, error = function(e) {
        add_log(paste("ERROR:", e$message))
        showNotification(paste("Error:", e$message), type = "error")
      })
    })
    
    output$download_counts <- downloadHandler(
      filename = function() { paste("AbsCounts_", Sys.Date(), ".csv", sep="") },
      content = function(file) {
        data <- processed_data()
        req(data$count_data)
        write.csv(data$count_data, file, row.names = FALSE)
      }
    )
    
    return(list(
        data = processed_data,
        log = log_text
    ))
  })
}
