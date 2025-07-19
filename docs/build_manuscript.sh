#!/bin/bash
# Universal Distance Dilation Theory - Manuscript Build Script
# Converts markdown to PDF using pandoc with citations

echo "Building UDT manuscript..."

# Check if pandoc is installed
if ! command -v pandoc &> /dev/null; then
    echo "Error: pandoc is not installed. Please install pandoc first."
    exit 1
fi

# Input and output files
INPUT_FILE="Universal_Distance_Dilation_Theory_Manuscript.md"
OUTPUT_FILE="docs/udt.pdf"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: $INPUT_FILE not found in current directory"
    exit 1
fi

# Create docs directory if it doesn't exist
mkdir -p docs

# Convert markdown to PDF with citations
echo "Converting $INPUT_FILE to $OUTPUT_FILE..."

pandoc \
    --citeproc \
    --pdf-engine=xelatex \
    --variable geometry:margin=1in \
    --variable fontsize=11pt \
    --variable documentclass=article \
    --variable colorlinks=true \
    --variable linkcolor=blue \
    --variable urlcolor=blue \
    --variable citecolor=blue \
    --number-sections \
    --table-of-contents \
    --standalone \
    -o "$OUTPUT_FILE" \
    "$INPUT_FILE"

# Check if conversion was successful
if [ $? -eq 0 ]; then
    echo "✅ Successfully built manuscript: $OUTPUT_FILE"
    echo "📄 File size: $(ls -lh $OUTPUT_FILE | awk '{print $5}')"
else
    echo "❌ Error: Failed to build manuscript"
    exit 1
fi

echo "Build complete!"