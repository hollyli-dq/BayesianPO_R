# sushi_functions.R
# Utility functions for working with the Sushi preference dataset

#' Load sushi item features and names
#'
#' @param idata_path Path to the sushi3.idata file (100 rows × 7 numeric attributes)
#' @param item_name_path Optional path to item_mapping.txt (0-99 ⇒ Japanese-Roman names)
#' @return A data.frame with sushi item features
#' @details 
#' Returns a data.frame with:
#' - item_id: int (0-99)
#' - item_name: str (if mapping given)
#' - style: 0 = maki roll | 1 = otherwise
#' - major_group: 0 = seafood | 1 = other
#' - minor_group: 0-11 detailed class
#' - heaviness: 0-4 (0 = heavy/oily)
#' - freq_eat: 0-3 (3 = eaten frequently)
#' - price_norm: float [0,1]
#' - freq_sell: 0-1 (1 = most frequently sold)
load_sushi3_item_mapping <- function(idata_path, item_name_path = NULL) {
  # Check if main data file exists
  if (!file.exists(idata_path)) {
    stop(paste("File not found:", idata_path))
  }
  
  # 1) Read the raw data with flexible whitespace separation
  df <- read.table(
    idata_path,
    header = FALSE,
    sep = "",
    stringsAsFactors = FALSE,
    strip.white = TRUE
  )
  
  # 2) Extract item_id and name from first two columns
  item_ids <- as.integer(df[[1]])
  item_names <- df[[2]]
  
  # 3) Convert remaining numeric columns
  numeric_data <- df[, 3:ncol(df)]
  numeric_data <- as.data.frame(lapply(numeric_data, as.numeric))
  
  # 4) Combine with proper column names
  col_names <- c("style", "major_group", "minor_group",
                 "heaviness", "freq_eat", "price_norm", "freq_sell")
  
  result_df <- data.frame(
    item_id = item_ids,
    item_name = item_names,
    numeric_data,
    stringsAsFactors = FALSE
  )
  names(result_df)[3:ncol(result_df)] <- col_names
  
  # 5) Convert specific columns to integers
  int_cols <- c("style", "major_group", "minor_group", "heaviness")
  result_df[int_cols] <- lapply(result_df[int_cols], as.integer)
  
  # 6) Handle item name mapping if provided
  if (!is.null(item_name_path)) {
    if (!file.exists(item_name_path)) {
      stop(paste("File not found:", item_name_path))
    }
    
    # Read mapping file, skipping header lines
    lines <- readLines(item_name_path)
    name_rows <- list()
    
    for (line in lines) {
      line <- trimws(line)
      if (line == "" || startsWith(line, "*") || startsWith(line, "3.")) {
        next
      }
      
      if (grepl(":", line)) {
        parts <- strsplit(line, ":", fixed = TRUE)[[1]]
        if (length(parts) >= 2) {
          item_id <- as.integer(trimws(parts[1]))
          item_name <- trimws(strsplit(parts[2], "(", fixed = TRUE)[[1]][1])
          name_rows[[length(name_rows) + 1]] <- list(item_id = item_id, item_name = item_name)
        }
      }
    }
    
    if (length(name_rows) > 0) {
      names_df <- do.call(rbind, lapply(name_rows, as.data.frame))
      
      # Merge with existing data
      result_df <- merge(
        result_df, 
        names_df, 
        by = "item_id", 
        all.x = TRUE,
        suffixes = c("_data", "")
      )
      
      # Prefer names from mapping file
      if ("item_name_data" %in% names(result_df)) {
        result_df$item_name <- ifelse(
          is.na(result_df$item_name), 
          result_df$item_name_data, 
          result_df$item_name
        )
        result_df$item_name_data <- NULL
      }
      
      # Reorder columns
      result_df <- result_df[, c("item_id", "item_name", col_names)]
    }
  }
  
  return(result_df)
}

#' Visualize sushi features through exploratory plots
#'
#' @param df The sushi features dataframe from load_sushi3_item_mapping
#' @return Generates three ggplot2 visualizations
show_sushi_visualizations <- function(df) {
  # 1. Correlation Heatmap
  numeric_cols <- c('heaviness', 'freq_eat', 'price_norm', 'freq_sell')
  corr <- cor(df[numeric_cols], use = "complete.obs")
  
  print(
    ggplot(melt(corr), aes(x = Var1, y = Var2, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 0, limit = c(-1,1)) +
      geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      coord_fixed() +
      ggtitle("Feature Correlations")
  )
  
  # 2. Price Distribution
  print(
    ggplot(df, aes(x = price_norm)) +
      geom_histogram(bins = 20, fill = "steelblue", color = "black") +
      ggtitle("Distribution of Sushi Prices") +
      xlab("Normalized Price") +
      theme_minimal()
  )
  
  # 3. Heaviness vs Price
  print(
    ggplot(df, aes(x = heaviness, y = price_norm, color = factor(major_group))) +
      geom_point(size = 3) +
      ggtitle("Heaviness vs Price by Major Group") +
      xlab("Heaviness") +
      ylab("Normalized Price") +
      scale_color_discrete(name = "Major Group") +
      theme_minimal()
  )
}

#' Load sushi rankings from data file
#'
#' @param order_path Path to the sushi3a.5000.10.order file
#' @return List of rankings, where each ranking is a vector of 10 item IDs
#' @details 
#' Each ranking is a vector of 10 item IDs, where the first item is most preferred
load_sushi_rankings <- function(order_path) {
  # Check if file exists
  if (!file.exists(order_path)) {
    stop(paste("File not found:", order_path))
  }
  
  # Read all lines from file
  lines <- readLines(order_path)
  
  # Skip header line (first line)
  lines <- lines[-1]
  
  # Process each line
  rankings <- lapply(lines, function(line) {
    # Split line into components
    items <- strsplit(trimws(line), "\\s+")[[1]]
    
    # Extract item IDs (skip first two columns)
    item_ids <- items[3:length(items)]
    
    # Convert to integers
    as.integer(item_ids)
  })
  
  return(rankings)
}



process_sushi_numerical <- function(features_df) {
  # Select only numerical columns
  numerical_cols <- c("freq_eat", "price_norm", "freq_sell", "heaviness")
  
  # Create a copy with just numerical features
  numerical_features <- features_df[numerical_cols]
  
  # Normalize each column (skip if no variance)
  for(col in numerical_cols) {
    if(sd(numerical_features[[col]], na.rm = TRUE) > 0) {
      numerical_features[[col]] <- scale(numerical_features[[col]])
    }
  }
  
  return(t(as.matrix(numerical_features)))  # Transpose to match Python structure
}

filter_selected_items <- function(rankings, item_features, selected_items) {
  # Input validation
  if (!is.list(rankings)) stop("rankings must be a list")
  if (!is.data.frame(item_features)) stop("item_features must be a dataframe")
  if (!"item_id" %in% names(item_features)) stop("item_features must contain 'item_id' column")
  if (length(selected_items) < 2) stop("Must select at least 2 items")
  
  # Filter features
  filtered_features <- item_features[item_features$item_id %in% selected_items, ]
  
  # Filter rankings (keep only selected items in original order)
  filtered_rankings <- lapply(rankings, function(r) {
    r[r %in% selected_items]
  })
  
  # Remove empty rankings
  non_empty <- sapply(filtered_rankings, length) > 0
  filtered_rankings <- filtered_rankings[non_empty]
  
  # Return results
  list(
    filtered_rankings = filtered_rankings,
    filtered_features = filtered_features
  )
}
# Example usage:
# selected <- c(3, 7, 15, 23, 42)  # Example item IDs to keep
# filtered <- filter_selected_items(rankings, item_features, selected)
# 
# filtered$filtered_rankings  # Rankings with only selected items
# filtered$filtered_features  # Features for selected items only
# filtered$item_mapping       # Shows original and new item IDs
