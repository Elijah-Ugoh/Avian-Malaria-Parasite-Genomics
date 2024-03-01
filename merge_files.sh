#!/bin/bash

# Check if correct number of arguments provided
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <output_file> <header_file> <data_file>"
    exit 1
fi

# Output file
output_file=$1
header_file=$2
data_file=$3

# Check if output file already exists
if [ -e "$output_file" ]; then
    echo "Output file already exists. Please specify a different name."
    exit 1
fi

# Concatenate files with a newline between them
cat "$header_file" "$data_file" > "$output_file"

echo "Files merged successfully into $output_file."