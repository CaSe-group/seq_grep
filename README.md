grep for reads
===

If you need grep with mutliple colors and searchterms to visual inspect reads for problems

## What does it do
* quick grep tool to highlight sequence patterns in fastq files
* e.g. searching and vis for multi barcodes or primers or remaining adapters
* supports reverse complement search `--rc`
* colors are red, blue, green, yellow and up to 9 search terms per color
    * `-r1` would be a search term colored in red
* can multithread if .fastq was used as an input; .fastq.gz is slow

example code:
```
python3 searcher.py all.fastq -r1 "GGTTACACAAACCCTGGACA" -b2 "GGATTCATTCCCACGGTAAC" -rc --limit 20
```

## terminal output

![alt text](image.png)