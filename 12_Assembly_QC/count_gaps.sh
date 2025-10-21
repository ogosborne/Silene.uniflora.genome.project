# For each fasta, count the number of gaps (runs of Ns >= 10 bp)
for f in *.fasta; do
  if [ -f "$f" ]; then
    # Step 1: Convert to single-line FASTA using awk, store in a temp file
    temp_file=$(mktemp)
    awk '/^>/ { if (NR > 1) { print ""; } printf "%s\n", tolower($0); next; } { printf "%s", tolower($0); } END { print ""; }' "$f" > "$temp_file"

    # Step 2: Process the temp file to count gaps
    num_gaps=$(grep -v '^>' "$temp_file" | grep -o 'n\{10,\}' | wc -l)
    
    # Clean up the temporary file
    rm "$temp_file"
    
    # Print the filename and the gap count
    echo -e "$f\t$num_gaps"
  fi
done > gap_report.txt
