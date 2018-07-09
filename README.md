*scantop* scans sequences for instances of motifs.

Usage:
`scantop <ace_file> <motif_number> <sequence_file> <background_GC> <site_count> <output_file>`
- ace_file: AlignACE output file with motifs
- motif_number:	Motif number within ace_file to use
- sequence_file: File with sequences to scan (FASTA format)
- background_GC: Background GC content (should be between 0 and 1)
- site_count: Number of top scoring instances to output per sequence
- output_file: Output file with motif instances


