#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------------
# --- Convert gene names to Ensembl IDs by parsing Uniprot reference proteome DB
# ------------------------------------------------------------------------------------

# Uniprot reference proteomes data found at:
uniprot.ftp.path = paste0(	"ftp://ftp.ebi.ac.uk/pub/databases/reference_proteomes/",
							"previous_releases/qfo_release-2013_04/",
							"QfO_release_2013_04.tar.gz")
# Local path to Uniprot proteome database
uniprot.path = "data/QfO_release_2013_04.tar.gz"

if (file.exists(uniprot.path) == FALSE) {
	download.file(uniprot.ftp.path, uniprot.path)
}

if (file.exists("data/9544.fasta") == FALSE) {
	# Just decompress macaque (taxon ID = 9544) files
	untar(uniprot.path, exdir="data", compress="gzip", 
			extras="--wildcards --no-anchored '9544*'")
}

uni.raw = readLines(file("data/9544.fasta"))
uni.head = unique(grep("^>", uni.raw, value = TRUE))

# Grab Uniprot ID and gene name
uni.ids    = gsub(".*\\|([\\w\\d]+)\\|.*", "\\1", uni.head, perl=TRUE)
gene.names = gsub(".*GN=([\\w\\d-\\(\\)\\.\\/]+)\\s.*", "\\1", uni.head, perl=TRUE)

# Combine and remove rows that lack gene names
assoc.ids = cbind(uni.ids, gene.names)[-grep(">", gene.names),]

# Write these associated identifiers to a file
write.table(assoc.ids, 
			file="data/9544_macaca_mulatta.uniprot2genename", 
			quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
