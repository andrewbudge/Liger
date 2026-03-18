// Liger: A smart concatination tool
// version: 0.2.0
// written and developed by Andrew Budge

package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"sort"
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
	expectedLength := -1

	scanner := bufio.NewScanner(file)

	// Future proof reading large seqs on a single line
	const maxCap = 1024 * 1024 * 1024
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, maxCap)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		if strings.HasPrefix(line, ">") {
			// Save previous sequence
			if currentHeader != "" {
				seq := strings.ToUpper(currentSeq.String())
				if expectedLength == -1 {
					expectedLength = len(seq)
				} else {
					if len(seq) != expectedLength {
						return nil, 0, fmt.Errorf("Error: unequal sequence length found in %s", filename)
					}
				}
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
		if expectedLength == -1 {
			expectedLength = len(seq)
		} else {
			if len(seq) != expectedLength {
				return nil, 0, fmt.Errorf("Error: unequal sequence length found in %s", filename)
			}
		}
		sequences[currentHeader] = seq
		if len(seq) > maxLen {
			maxLen = len(seq)
		}
	}

	return sequences, maxLen, scanner.Err()
}

// findSequenceByPattern searches for a taxon name pattern in headers (like grep -p)
func findSequenceByPattern(sequences map[string]string, pattern string, claimed map[string]bool) (string, string, bool) {
	normalizedPattern := strings.ToLower(pattern)
	for header, seq := range sequences {
		if claimed[header] {
			continue
		}
		normalizedHeader := strings.ToLower(header)
		if strings.Contains(normalizedHeader, normalizedPattern) {
			return seq, header, true
		}
	}
	return "", "", false
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

func printNexus(taxa []string, supermatrix map[string]string, genes []Gene, fillChar string) {
	totalLen := 0
	for _, gene := range genes {
		totalLen += gene.Length
	}

	maxTaxonLen := 0
	for _, taxon := range taxa {
		if len(taxon) > maxTaxonLen {
			maxTaxonLen = len(taxon)
		}
	}

	fmt.Println("#NEXUS")
	fmt.Println()

	fmt.Println("BEGIN TAXA;")
	fmt.Printf("    DIMENSIONS NTAX=%d;\n", len(taxa))
	fmt.Println("    TAXLABELS")
	for _, taxon := range taxa {
		fmt.Printf("        %s\n", taxon)
	}
	fmt.Println("    ;")
	fmt.Println("END;")
	fmt.Println()

	fmt.Println("BEGIN CHARACTERS;")
	fmt.Printf("    DIMENSIONS NCHAR=%d;\n", totalLen)
	fmt.Printf("    FORMAT DATATYPE=DNA MISSING=%s GAP=-;\n", fillChar)
	fmt.Println("    MATRIX")
	for _, taxon := range taxa {
		fmt.Printf("        %-*s  %s\n", maxTaxonLen, taxon, supermatrix[taxon])
	}
	fmt.Println("    ;")
	fmt.Println("END;")
	fmt.Println()

	fmt.Println("BEGIN SETS;")
	currentPos := 1
	for _, gene := range genes {
		endPos := currentPos + gene.Length - 1
		fmt.Printf("    CHARSET %s = %d-%d;\n", gene.Name, currentPos, endPos)
		currentPos = endPos + 1
	}
	fmt.Println("END;")
}

func main() {

	formatFlag := flag.String("f", "fasta", "Output format: fasta or nexus")
	missingCharFlag := flag.String("m", "N", "Character for missing data")

	flag.Usage = func() {
		fmt.Fprintln(os.Stderr, "Usage: liger [FLAGS] [TAXA LIST] [INPUT FILES]")
		fmt.Fprintln(os.Stderr, "\nFlags:")
		fmt.Fprintln(os.Stderr, "  -f | choose output format fasta or nexus (default fasta)")
		fmt.Fprintln(os.Stderr, "  -m | missing symbol (default N)")
		fmt.Fprintln(os.Stderr, "\nOutput:")
		fmt.Fprintln(os.Stderr, "  FASTA mode: Supermatrix to stdout, partitions to stderr")
		fmt.Fprintln(os.Stderr, "  NEXUS mode: Complete NEXUS file to stdout")
	}

	flag.Parse()

	args := flag.Args()
	if len(args) < 2 {
		flag.Usage()
		os.Exit(1)
	}

	taxaListFile := args[0]
	inputFiles := args[1:]

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
		missingData[i] = strings.Repeat(*missingCharFlag, gene.Length)
	}

	// Sort taxa list by length, for hierarchical matching
	sort.SliceStable(taxa, func(i, j int) bool {
		return len(taxa[i]) > len(taxa[j])
	})

	// claim seqs to prevent double assignment of seqs to taxa (very bad!)
	claimed := make([]map[string]bool, len(genes))
	for i := range claimed {
		claimed[i] = make(map[string]bool)
	}

	// Generate supermatrix
	supermatrix := make(map[string]string)

	for _, taxon := range taxa {
		var concat strings.Builder

		for i, gene := range genes {
			if seq, header, found := findSequenceByPattern(gene.Sequences, taxon, claimed[i]); found {
				concat.WriteString(seq)
				claimed[i][header] = true
			} else {
				concat.WriteString(missingData[i])
			}
		}

		supermatrix[taxon] = concat.String()
	}

	sort.Strings(taxa)

	if *formatFlag == "nexus" {
		printNexus(taxa, supermatrix, genes, *missingCharFlag)
	} else if *formatFlag == "fasta" {
		for _, taxon := range taxa {
			fmt.Printf(">%s\n%s\n", taxon, supermatrix[taxon])
		}

		fmt.Fprintln(os.Stderr, "#NEXUS")
		fmt.Fprintln(os.Stderr, "begin sets;")

		currentPos := 1
		for _, gene := range genes {
			endPos := currentPos + gene.Length - 1
			fmt.Fprintf(os.Stderr, "  charset %s = %d-%d;\n", gene.Name, currentPos, endPos)
			currentPos = endPos + 1
		}

		fmt.Fprintln(os.Stderr, "end;")
	} else {
		fmt.Fprintf(os.Stderr, "Error: invalid format '%s'. Use 'fasta' or 'nexus'\n", *formatFlag)
		os.Exit(1)
	}
}
