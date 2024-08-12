# Check if "data" directory exists
DATA_DIR <- "data"

if (!file.exists(DATA_DIR)) {
    # Create "data" directory
    dir.create(DATA_DIR)
}

# TODO: Set logs level and logs destination

