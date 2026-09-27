# R/modules/mod_heatmap.R

mod_heatmap_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # CSS for compact controls
    tags$style(HTML(paste0("#", ns("controls_card"), " .card-body { font-size: 0.85em; }"))),
    
    # 1. Resizable Plot Container
    div(class="card mb-3",
        div(class="card-body p-0",
            # Resizable wrapper
            jqui_resizable(
                plotOutput(ns("heatmap_plot"), width = "100%", height = "600px"),
                options = list(minHeight = 400, minWidth = 400)
            )
        )
    ),
    
    # 2. Control Ribbon (Tabbed Interface)
    card(
      id = ns("controls_card"), 
      card_header("Visualization Controls"),
      style = "min-height: 300px;",
      card_body(
        navset_card_tab(
            
            # --- TAB 1: HIERARCHY ---
            nav_panel("Hierarchy",
                div(class="d-flex gap-2 mb-2",
                    actionButton(ns("select_all"), "All", size="xs", class="btn-xs btn-light"),
                    actionButton(ns("deselect_all"), "None", size="xs", class="btn-xs btn-light")
                ),
                uiOutput(ns("anno_selector_ui"))
            ),
            
            # --- TAB 2: SORTING & ORDER ---
            nav_panel("Sorting",
                layout_columns(
                    col_widths = c(6, 6),
                    
                    # Left: Main Clustering Logic
                    div(
                        strong("Matrix Organization"),
                        radioButtons(ns("col_sort_method"), "Column Order:", 
                                     choices = c("By Hierarchy (Manual)"="hierarchy", "Cluster (Auto)"="cluster"), 
                                     selected = "hierarchy"),
                        checkboxInput(ns("cluster_rows"), "Cluster Rows", value = TRUE),
                        p(class="text-muted", style="font-size:0.8em; line-height:1.1;",
                          "Note: 'By Hierarchy' uses the order defined in the Hierarchy tab.")
                    ),
                    
                    # Right: Intra-Group Level Sequencing
                    div(
                        strong("Column Ordering"),
                        p(class="text-muted", style="font-size:0.8em; margin-bottom:5px;", 
                          "Reorder column for the selected variable."),
                        selectInput(ns("level_target_var"), NULL, choices=NULL, width="100%"),
                        uiOutput(ns("level_order_ui"))
                    )
                )
            ),
            
            # --- TAB 3: AESTHETICS ---
            nav_panel("Aesthetics",
                 layout_columns(
                     col_widths = c(6, 6),
                     
                     # Left: Color Editor
                     div(
                        strong("Annotation Colors"),
                        selectInput(ns("color_target_var"), "Target Variable:", choices=NULL, width="100%"),
                        uiOutput(ns("dynamic_colors_ui"))
                     ),
                     
                     # Right: Heatmap Settings
                     div(
                        strong("Heatmap & Grid"),
                        div(class="d-flex gap-1 mb-2",
                            colourInput(ns("col_low"), "Low", value = "blue", showColour = "both", width="32%"),
                            colourInput(ns("col_mid"), "Mid", value = "white", showColour = "both", width="32%"),
                            colourInput(ns("col_high"), "High", value = "red", showColour = "both", width="32%")
                        ),
                        div(class="d-flex gap-2",
                            colourInput(ns("grid_col"), "Grid", value = "white", showColour="both", width="80px"),
                            radioButtons(ns("cell_label_type"), "Labels:",
                                     choices = c("None"="none", "Raw"="raw", "Z"="z"),
                                     selected = "none", inline = TRUE)
                        )
                     )
                 )
            ),
            
            # --- TAB 4: EXPORT ---
            nav_panel("Export",
                div(class="d-flex flex-column align-items-center",
                    numericInput(ns("height"), "Height (px)", value=600, step=50, width="50%"),
                    numericInput(ns("width"), "Width (px)", value=800, step=50, width="50%"),
                    downloadButton(ns("download_pdf"), "Download PDF", class="btn-primary w-50 mt-3")
                )
            )
        )
      )
    )
  )
}

mod_heatmap_server <- function(id, data_reactive) {
  moduleServer(id, function(input, output, session) {
    
    # helper for matrix conversion
    make_matrix <- function(df, val_var="Z_Score") {
      # Handle duplicates: Combine marker and Population for highly resolved heatmap rows
      mat_df <- df %>%
        mutate(Feature = paste0(marker, " [", Population, "]")) %>%
        select(Feature, SampleID, !!sym(val_var)) %>%
        group_by(Feature, SampleID) %>%
        summarise(val = mean(!!sym(val_var), na.rm=TRUE), .groups="drop") %>%
        pivot_wider(names_from = SampleID, values_from = val) %>%
        tibble::column_to_rownames("Feature")
      
      return(as.matrix(mat_df))
    }
    
    # --- 0. Metadata Detection Engine ---
    detect_metadata <- reactive({
      df <- data_reactive()
      req(df)
      
      # Get Sample-level metadata
      meta <- df %>% 
        select(-marker, -Value, -Value_Transformed, -Z_Score, -CleanFluor, -Parameter) %>%
        distinct(SampleID, .keep_all=TRUE)
      
      # Identify columns
      candidates <- names(meta)[sapply(meta, function(x) {
         (is.character(x) || is.factor(x)) && length(unique(x)) < 24
      })]
      
      candidates <- setdiff(candidates, c("SampleID", "Sample", "remark", "machine"))
      
      return(list(meta=meta, cols=candidates))
    })
    
    # --- 1A. Hierarchy & Selection ---
    active_annos <- reactiveVal(NULL)
    
    observeEvent(detect_metadata(), {
        info <- detect_metadata()
        defaults <- intersect(c("condition", "bone_type"), info$cols)
        if(is.null(active_annos())) active_annos(defaults)
    }, once = TRUE)
    
    observeEvent(input$select_all, { active_annos(detect_metadata()$cols) })
    observeEvent(input$deselect_all, { active_annos(character(0)) })
    observeEvent(input$selected_annos, {
        if(!identical(input$selected_annos, active_annos())) active_annos(input$selected_annos)
    }, ignoreNULL = FALSE)
    
    output$anno_selector_ui <- renderUI({
       info <- detect_metadata()
       ns <- session$ns
       current_active <- intersect(active_annos(), info$cols)
       available <- setdiff(info$cols, current_active)
       
       bucket_list(
          header = NULL, 
          group_name = ns("bucket_group"),
          orientation = "vertical",
          add_rank_list("Active (Ordered)", labels = current_active, input_id = ns("selected_annos"), options = sortable_options(multiDrag=TRUE)),
          add_rank_list("Available / Hidden", labels = available, input_id = ns("available_annos"), options = sortable_options(multiDrag=TRUE))
       )
    })
    
    # --- 1B. Level Sequencing (Manual Order) ---
    # Store user preferences for level order per variable
    level_prefs <- reactiveValues()
    
    # Update choices for Level Sequencer
    observe({
        req(active_annos())
        updateSelectInput(session, "level_target_var", choices = active_annos())
    })
    
    # Render Rank List for Current Target
    output$level_order_ui <- renderUI({
        target <- input$level_target_var
        req(target)
        info <- detect_metadata()
        if(!(target %in% info$cols)) return(NULL)
        
        # Get levels: check if we have a stored pref, otherwise alphanumeric
        default_levels <- sort(unique(info$meta[[target]]))
        
        # If stored pref exists, use it validly
        if(!is.null(level_prefs[[target]])) {
            # Intersect to ensure validity (e.g. if data changed)
            current_pref <- level_prefs[[target]]
            # If all current levels are in pref, respect pref order
            if(all(default_levels %in% current_pref)) {
                # Filter pref to only existing, append any new
                existing_ordered <- intersect(current_pref, default_levels)
                new_items <- setdiff(default_levels, current_pref)
                final_levels <- c(existing_ordered, new_items)
            } else {
                final_levels <- default_levels
            }
        } else {
            final_levels <- default_levels
        }
        
        ns <- session$ns
        rank_list(
            text = NULL,
            labels = final_levels,
            input_id = ns("manual_level_order")
        )
    })
    
    # Save preferences when user reorders
    observeEvent(input$manual_level_order, {
        target <- input$level_target_var
        req(target)
        # Only update if meaningful
        level_prefs[[target]] <- input$manual_level_order
    })
    
    
    # --- 1C. Aesthetics (Colors) ---
    observe({
        req(active_annos())
        updateSelectInput(session, "color_target_var", choices = active_annos())
    })
    
    output$dynamic_colors_ui <- renderUI({
      info <- detect_metadata()
      target <- input$color_target_var
      req(target)
      if(!(target %in% info$cols)) return(NULL)
      ns <- session$ns
      
      # Base Palettes
      known_colors <- list(
          "SteadyState"="#808080", "LPS72h"="#FF0000", "PDAC"="#0000FF",
          "WTnaive24" = "#3498DB",   "WTzymozan24" = "#E74C3C",
          "Femur"="#1B9E77", "Lumbar"="#D95F02", "Sternum"="#7570B3", "Calvaria"="#E7298A"
      )
      get_color <- function(lvl) {
        if(lvl %in% names(known_colors)) return(known_colors[[lvl]])
        chk <- sum(utf8ToInt(as.character(lvl)))
        grDevices::rainbow(30, s=0.8, v=0.9)[(chk %% 30) + 1]
      }
      
      levels <- sort(unique(info$meta[[target]]))
      
      lapply(levels, function(lvl) {
         safe_col <- gsub("[^A-Za-z0-9]", "", target)
         safe_lvl <- gsub("[^A-Za-z0-9]", "", lvl)
         id <- paste0("cp_", safe_col, "_", safe_lvl)
         
         cur_val <- input[[id]]
         def_val <- if(!is.null(cur_val)) cur_val else get_color(lvl)
         
         div(style="display:inline-block; margin-right:5px; margin-bottom: 5px;",
             colourInput(ns(id), lvl, value = def_val, showColour="both", width="100px")
         )
      }) %>% div(class="d-flex flex-wrap", .)
    })

    # 1D. Render Heatmap (Generic + Secure)
    plot_obj <- reactive({
      req(data_reactive())
      # Use input directly to avoid sync issues
      # but verify it exists
      selected_annos <- input$selected_annos
      
      df <- data_reactive()
      validate(need(nrow(df) > 0, "No data available"))
      
      # Filter NAs and Inf
      df <- df %>% 
        filter(!is.na(marker), !is.na(Z_Score)) %>%
        filter(is.finite(Z_Score))
      
      if(nrow(df) == 0) return(NULL)

      # NOTE: Removed tryCatch to expose errors and fix syntax
      
          # Create Matrices
          mat <- make_matrix(df, "Z_Score")
          mat_raw <- make_matrix(df, "Value")
          
          # Get Metadata Map
          info <- detect_metadata()
          meta_map <- info$meta
          
          # Align
          common_samples <- intersect(colnames(mat), as.character(meta_map$SampleID))
          if(length(common_samples) == 0) return(NULL)
          
          mat <- mat[, common_samples, drop=FALSE]
          mat_raw <- mat_raw[, common_samples, drop=FALSE]
          meta_map <- meta_map %>% filter(as.character(SampleID) %in% common_samples)
          
          # --- DYNAMIC SORTING ---
          cluster_cols_arg <- TRUE
          cluster_rows_arg <- input$cluster_rows 
          if(is.null(cluster_rows_arg)) cluster_rows_arg <- TRUE
          
          # Only perform hierarchical sorting if we have selected annotations
          if (input$col_sort_method == "hierarchy" && length(selected_annos) > 0) {
              # Sort by the ORDER of selected_annos
              for(col in selected_annos) {
                  # Default alphanumeric levels
                  default_levels <- sort(unique(meta_map[[col]]))
                  
                  # Check for manual preference
                  if(!is.null(level_prefs[[col]])) {
                      pref <- level_prefs[[col]]
                      # Validate: only present levels
                      valid_pref <- intersect(pref, default_levels)
                      remainder <- setdiff(default_levels, valid_pref)
                      final_levs <- c(valid_pref, remainder)
                  } else {
                      final_levs <- default_levels
                  }
                  
                  meta_map[[col]] <- factor(meta_map[[col]], levels = final_levs)
              }
              
              meta_map <- meta_map %>% arrange(across(all_of(selected_annos)))
              
              sorted_samples <- as.character(meta_map$SampleID)
              mat <- mat[, sorted_samples, drop=FALSE]
              mat_raw <- mat_raw[, sorted_samples, drop=FALSE]
              cluster_cols_arg <- FALSE
              
              # Row Staircase Logic
              if(isTRUE(input$cluster_rows)) {
                   # If hierarchy is active, we might want staircase or cluster
                   # Original logic: staircase sort (max index)
                   max_indices <- apply(mat, 1, which.max)
                   sort_df <- data.frame(Marker = rownames(mat), MaxIndex = max_indices) %>% arrange(MaxIndex, Marker)
                   mat <- mat[sort_df$Marker, , drop=FALSE]
                   mat_raw <- mat_raw[sort_df$Marker, , drop=FALSE]
                   cluster_rows_arg <- FALSE
              }
          }
          
          # --- DYNAMIC ANNOTATION CONSTRUCTION ---
          ha <- NULL
          if(length(selected_annos) > 0) {
              # Build Annotation DF
              anno_df <- meta_map %>% select(all_of(selected_annos)) %>% as.data.frame()
              
              # Build Color List
              col_list <- list()
              
              known_colors_internal <- list(
                  "SteadyState"="#808080", "LPS72h"="#FF0000", "PDAC"="#0000FF",
                  "WTnaive24" = "#3498DB",   "WTzymozan24" = "#E74C3C",
                  "Femur"="#1B9E77", "Lumbar"="#D95F02", "Sternum"="#7570B3", "Calvaria"="#E7298A"
              )
              get_color_internal <- function(lvl) {
                if(lvl %in% names(known_colors_internal)) return(known_colors_internal[[lvl]])
                chk <- sum(utf8ToInt(as.character(lvl)))
                idx <- (chk %% 30) + 1
                grDevices::rainbow(30, s=0.8, v=0.9)[idx]
              }
              
              for(col in selected_annos) {
                   levels_vec <- as.character(unique(anno_df[[col]]))
                   
                   cols_vec <- sapply(levels_vec, function(lvl) {
                       safe_col <- gsub("[^A-Za-z0-9]", "", col)
                       safe_lvl <- gsub("[^A-Za-z0-9]", "", lvl)
                       id <- paste0("cp_", safe_col, "_", safe_lvl)
                       
                       val <- input[[id]]
                       if(is.null(val)) return(get_color_internal(lvl)) else return(val)
                   })
                   col_list[[col]] <- cols_vec
              }
              
              ha <- HeatmapAnnotation(df = anno_df, col = col_list)
          }
          
          # Dynamic Color Ramp
          col_fun <- colorRamp2(c(-2, 0, 2), c(input$col_low, input$col_mid, input$col_high))
          
          # Define Cell Function
          cell_fun_arg <- NULL
          if(input$cell_label_type == "raw") {
              cell_fun_arg <- function(j, i, x, y, width, height, fill) {
                  val <- mat_raw[i, j]
                  z_val <- mat[i, j]
                  txt_col <- if(abs(z_val) > 1.5) "white" else "black"
                  grid.text(round(val, 0), x, y, gp = gpar(fontsize = 8, col = txt_col))
              }
          } else if(input$cell_label_type == "z") {
              cell_fun_arg <- function(j, i, x, y, width, height, fill) {
                  val <- mat[i, j]
                  txt_col <- if(abs(val) > 1.5) "white" else "black"
                  grid.text(round(val, 2), x, y, gp = gpar(fontsize = 8, col = txt_col))
              }
          }

          Heatmap(mat,
                  name = "Z-Score",
                  col = col_fun,
                  top_annotation = ha,
                  cluster_rows = cluster_rows_arg,
                  cluster_columns = cluster_cols_arg, 
                  rect_gp = gpar(col = input$grid_col, lwd = 1),
                  cell_fun = cell_fun_arg,
                  column_title = paste("Samples (n=", ncol(mat), ")", sep=""),
                  row_title = "Markers",
                  show_column_names = TRUE)
    })
    
    # 3. Render Plot
    output$heatmap_plot <- renderPlot({
      p <- plot_obj()
      req(p)
      draw(p)
    })
    
    # 3b. Sync Resize to PDF Inputs (WYSIWYG)
    observeEvent(input$heatmap_plot_size, {
      sz <- input$heatmap_plot_size
      req(sz$height, sz$width)
      
      updateNumericInput(session, "height", value = round(sz$height))
      updateNumericInput(session, "width", value = round(sz$width))
    })
    
    # 4. Download Handler
    output$download_pdf <- downloadHandler(
      filename = function() {
        paste("Heatmap_ZScore_", Sys.Date(), ".pdf", sep="")
      },
      content = function(file) {
        req(plot_obj())
        w_in <- input$width / 72
        h_in <- input$height / 72
        
        pdf(file, width = w_in, height = h_in) 
        draw(plot_obj())
        dev.off()
      }
    )
    
  })
}
