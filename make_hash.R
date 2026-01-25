library(digest)

# Function to compute SHA-256 hash of a file
compute_file_hash <- function(filepath) {
  if (!file.exists(filepath)) {
    stop("File does not exist: ", filepath)
  }
  
  hash <- digest(file = filepath, algo = "sha256")
  return(hash)
}

# Function to save hash to CHECKSUMS file
save_hash_to_file <- function(filepath, output_file = "CHECKSUMS.md") {
  hash <- compute_file_hash(filepath)
  filename <- basename(filepath)
  
  # Append hash and filename to output file
  line <- sprintf("%s  %s\n", hash, filename)
  cat(line, file = output_file, append = TRUE)
  
  invisible(hash)
}

compute_and_save_hash <- function(filepath, output_file = "CHECKSUMS.md"){
  hash <- compute_file_hash(filepath)
  filename <- basename(filepath)
  
  # Append hash and filename to output file
  line <- sprintf("%s  %s\n", hash, filename)
  cat(line, file = output_file, append = TRUE)
  
  invisible(hash)
}