#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Download and read in PANTHER file
# ----------------------------------------------------------------------------------------

# PANTHER's macaque found at:
panther.macaque.ftp.path = paste0(	"ftp://ftp.pantherdb.org//sequence_classifications/",
									"current_release/",
									"PANTHER_Sequence_Classification_files/",
									"PTHR9.0_macaque")
# Local path to PANTHER database
panther.macaque.path = "data/macaque.panther.db"

if (file.exists(panther.macaque.path) == FALSE) {
	download.file(panther.macaque.ftp.path, panther.macaque.path)
}

panther = read.delim(panther.macaque.path, sep="\t", header=FALSE, stringsAsFactors=FALSE)
 
names(panther) = c(	"gene.id", "protein.id",
					"panther.sf.id", "panther.family", "panther.subfamily", 
					"go.mf", "go.bp", "go.cc",
					"panther.proteinclass", "panther.pathway")

# ----------------------------------------------------------------------------------------
# --- Parse ID string
# ----------------------------------------------------------------------------------------

# Separate gene.id into its parts: Ensembl ID, Uniprot accession
panther$ensembl.id = gsub(".*Ensembl=([^\\|]*)\\|*.*", "\\1", panther$gene.id)
panther$uniprot.id = gsub(".*UniProtKB=([^\\|]*)\\|*.*", "\\1", panther$gene.id)

# ----------------------------------------------------------------------------------------
# --- Parse lists of GO ID strings
# ----------------------------------------------------------------------------------------

# Function to make list of GO IDs from each ontology
parse.go = function(go.column) {
	go.list = regmatches(go.column, gregexpr("GO:[0-9]*", go.column))
	names(go.list) = panther$gene.id
 	return(go.list)
}

go.mf = parse.go(panther$go.mf)
go.bp = parse.go(panther$go.bp)
go.cc = parse.go(panther$go.cc)

# ----------------------------------------------------------------------------------------
# --- Parse pathway data
# ----------------------------------------------------------------------------------------

# From README:
#	Example Pathway:
#	Inflammation mediated by [...]#Inflammation mediated by [...]#P00031 [no break]
#	>Integrin#Integrin#P00853;Integrin signal [...]#Integrin signal [...]#P00034
#	>Integrin  alpha#Integrin alpha#P00941
#
#	The format of the pathway information is: 
#	pathway_long_name#pathway_short_name#pathway_accession [no break]
#	>component_long_name#component_short_name#component_accession
#
#	Explanation of pathway accessions:
#	Gxxxxx   Gene or RNA
#	Pxxxxx   Protein
#	Sxxxxx   small molecules
#	Uxxxxx   others, such as "unknown", etc.

# Odd values are pathway IDs, even values are component IDs
pathway.list = regmatches(	panther$panther.pathway, 
							gregexpr("#[G,P,S,U][0-9]*>*", panther$panther.pathway))

names(pathway.list) = panther$gene.id

# ----------------------------------------------------------------------------------------
# --- Parse protein class data
# ----------------------------------------------------------------------------------------

protein.list = regmatches(	panther$panther.proteinclass, 
							gregexpr("#[^;]+", panther$panther.proteinclass))

names(protein.list) = panther$gene.id

# ----------------------------------------------------------------------------------------
# --- Package parsed PANTHER data up as an object
# ----------------------------------------------------------------------------------------

panther.classification = list(  data=panther,
								go.mf=go.mf,
								go.bp=go.bp,
								go.cc=go.cc,
								panther.pathway=pathway.list,
								panther.proteinclass=protein.list)

class(panther.classification) = "panther.classification"
