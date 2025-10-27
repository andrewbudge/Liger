package main

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"strings"
)

type Gene struct {
	Name      string
	Length    int
	Sequences map[string]string // key is the FULL original header
}

func parseFasta(filename string) (map[string]string, int, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, 0, err
	}
	defer file.Close()

	sequences := make(map[string]string)
	var currentHeader string
	var currentSeq strings.Builder
	maxLen := 0

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		
		if strings.HasPrefix(line, ">") {
			// Save previous sequence
			if currentHeader != "" {
				seq := strings.ToUpper(currentSeq.String())
				sequences[currentHeader] = seq
				if len(seq) > maxLen {
					maxLen = len(seq)
				}
			}
			
			// Store the FULL header (everything after >)
			currentHeader = line[1:]
			currentSeq.Reset()
		} else if line != "" {
			currentSeq.WriteString(line)
		}
	}

	// Save last sequence
	if currentHeader != "" {
		seq := strings.ToUpper(currentSeq.String())
		sequences[currentHeader] = seq
		if len(seq) > maxLen {
			maxLen = len(seq)
		}
	}

	return sequences, maxLen, scanner.Err()
}

// findSequenceByPattern searches for a taxon name pattern in headers (like grep -p)
func findSequenceByPattern(sequences map[string]string, pattern string) (string, bool) {
	for header, seq := range sequences {
		// Case-insensitive substring match, like grep
		if strings.Contains(strings.ToLower(header), strings.ToLower(pattern)) {
			return seq, true
		}
	}
	return "", false
}

func loadTaxaList(filename string) ([]string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var taxa []string
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		taxon := strings.TrimSpace(scanner.Text())
		if taxon != "" {
			taxa = append(taxa, taxon)
		}
	}

	return taxa, scanner.Err()
}

func main() {
	if len(os.Args) < 3 {
		fmt.Fprintln(os.Stderr, "Usage: liger [INPUT FILES] [TAXA LIST]")
		fmt.Fprintln(os.Stderr, "Output: Supermatrix to stdout, partition file to stderr")
		os.Exit(1)
	}

	inputFiles := os.Args[1 : len(os.Args)-1]
	taxaListFile := os.Args[len(os.Args)-1]

	// Load taxa list
	taxa, err := loadTaxaList(taxaListFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}

	// Load all genes
	genes := make([]Gene, 0, len(inputFiles))
	for _, fastaFile := range inputFiles {
		sequences, length, err := parseFasta(fastaFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error reading %s: %v\n", fastaFile, err)
			os.Exit(1)
		}

		geneName := strings.TrimSuffix(filepath.Base(fastaFile), filepath.Ext(fastaFile))
		genes = append(genes, Gene{
			Name:      geneName,
			Length:    length,
			Sequences: sequences,
		})
	}

	// Pre-generate missing data strings
	missingData := make(map[int]string)
	for i, gene := range genes {
		missingData[i] = strings.Repeat("N", gene.Length)
	}

	// Generate supermatrix
	for _, taxon := range taxa {
		var concat strings.Builder
		
		for i, gene := range genes {
			// Search for taxon pattern in headers (like grep -p)
			if seq, found := findSequenceByPattern(gene.Sequences, taxon); found {
				concat.WriteString(seq)
			} else {
				concat.WriteString(missingData[i])
			}
		}
		
		fmt.Printf(">%s\n%s\n", taxon, concat.String())
	}

	// Print partition file to stderr
	fmt.Fprintln(os.Stderr, "#NEXUS")
	fmt.Fprintln(os.Stderr, "begin sets;")
	
	currentPos := 1
	for _, gene := range genes {
		endPos := currentPos + gene.Length - 1
		fmt.Fprintf(os.Stderr, "  charset %s = %d-%d;\n", gene.Name, currentPos, endPos)
		currentPos = endPos + 1
	}
	
	fmt.Fprintln(os.Stderr, "end;")
}
