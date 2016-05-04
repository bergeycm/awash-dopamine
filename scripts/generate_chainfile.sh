#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Steps to generate a baboon to rhesus LiftOver chain file
# ----------------------------------------------------------------------------------------

# After:
#  http://iamphioxus.org/2013/06/25/using-liftover-to-convert-genome-assembly-coordinates/
# See also:
#  http://genomewiki.ucsc.edu/index.php/Minimal_Steps_For_LiftOver

module load kent/intel/06042014 

# --- Make directory to do chain file generation work ------------------------------------

mkdir /scratch/cmb433/liftOver_chain_generation
cd /scratch/cmb433/liftOver_chain_generation

# --- Download genomes, get rid of unassembled stuff, and sort chromosomes numerically. --

SOURCE_FA=./papAnu2.RM.fa
TARGET_FA=./rheMac3.RM.fa

SOURCE_NAME=papAnu2
TARGET_NAME=rheMac3

if [ ! -e $TARGET_FA ]; then

    wget \
        'ftp://hgdownload.cse.ucsc.edu/goldenPath/rheMac3/bigZips/rheMac3.fa.masked.gz' \
        -O rheMac3.RM.fa.gz

    gunzip rheMac3.RM.fa.gz
    
    echo "Getting rid of unassembled stuff..." >&2

    LAST_OK_LINE_TARGET=$((`grep -n "^>[^c]" $TARGET_FA | head -n 1 | cut -d":" -f 1` - 1))
    if [ $LAST_OK_LINE_TARGET -gt 0 ]; then
        mv $TARGET_FA ${TARGET_FA}.backup
        head -n $LAST_OK_LINE_TARGET ${TARGET_FA}.backup > ${TARGET_NAME}.RM.fa
        rm ${TARGET_FA}.backup
    fi
    
    mkdir tmp_for_sort
    faSplit byname ${TARGET_FA} tmp_for_sort/
    cd tmp_for_sort/;
    ls -v | xargs cat > ../${TARGET_FA}
    cd ..
    rm -r tmp_for_sort
fi

if [ ! -e $SOURCE_FA ]; then
        
    wget \
        'ftp://hgdownload.cse.ucsc.edu/goldenPath/papAnu2/bigZips/papAnu2.fa.masked.gz' \
        -O papAnu2.RM.fa.gz

    gunzip papAnu2.RM.fa.gz
    
    echo "Getting rid of unassembled stuff..." >&2
    
    LAST_OK_LINE_SOURCE=$((`grep -n "^>[^c]" $SOURCE_FA | head -n 1 | cut -d":" -f 1` - 1))
    if [ $LAST_OK_LINE_SOURCE -gt 0 ]; then
        mv $SOURCE_FA ${SOURCE_FA}.backup
        head -n $LAST_OK_LINE_SOURCE ${SOURCE_FA}.backup > ${SOURCE_NAME}.RM.fa
        rm ${SOURCE_FA}.backup
    fi
    
    mkdir tmp_for_sort
    faSplit byname ${SOURCE_FA} tmp_for_sort/
    cd tmp_for_sort/;
    ls -v | xargs cat > ../${SOURCE_FA}
    cd ..
    rm -r tmp_for_sort
fi

TARGET_N=`grep "^>" $TARGET_FA | wc -l`

# --- Download earlier version of target genome too. -------------------------------------

TARGET_2_FA=./rheMac2.RM.fa
TARGET_2_NAME=rheMac2

if [ ! -e $TARGET_2_FA ]; then

    wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/rheMac2/bigZips/chromFaMasked.tar.gz'

    tar -zxvf chromFaMasked.tar.gz
    rm hardMask/chrUr.fa.masked

    ls -v hardMask/* | xargs cat > ${TARGET_2_FA}
    rm -r hardMask
fi

# --- Define function to run chainfile generating steps ----------------------------------

function make_chainfile {

    # --- Split genome into BLASTable hunks and make lift files --------------------------
    
    echo "Splitting genome..." >&2
    
    mkdir -p split_genomes/${TARGET_NAME}
    
    # Split up genome by chromosome
    faSplit sequence $TARGET_FA $TARGET_N split_genomes/${TARGET_NAME}/${TARGET_NAME}_
    
    # Split further into 3kb hunks and generate lift files
    for (( c=0; c<$TARGET_N; c++ )); do
        # Add leading zero, if needed, to match chromosome name
        chr=$(printf "%02d" $c)
        faSplit size split_genomes/${TARGET_NAME}/${TARGET_NAME}_$chr.fa 3000 \
            split_genomes/${TARGET_NAME}/${TARGET_NAME}_$chr.split \
            -lift=split_genomes/${TARGET_NAME}/${TARGET_NAME}_$chr.lft -oneFile
    done;
    
    # --- BLAT target genomic hunks against source genome --------------------------------
    
    echo "BLATing genome..." >&2
    
    mkdir blat_psl
    for (( c=0; c<$TARGET_N; c++ )); do
        # Add leading zero, if needed, to match chromosome name
        chr=$(printf "%02d" $c)
        blat $SOURCE_FA split_genomes/${TARGET_NAME}/${TARGET_NAME}_$chr.split.fa \
            -t=dna -q=dna -tileSize=12 -fastMap -minIdentity=90 -noHead -minScore=100 \
            ./blat_psl/${TARGET_NAME}_$chr.psl
    done
    
    # --- Change coordinates of PSL files ------------------------------------------------
    
    echo "Changing coordinates of PSL files..." >&2
    
    mkdir liftup
    for (( c=0; c<$TARGET_N; c++ )); do
        # Add leading zero, if needed, to match chromosome name
        chr=$(printf "%02d" $c)
        liftUp -pslQ liftup/${TARGET_NAME}_$chr.liftup.psl \
            split_genomes/${TARGET_NAME}/${TARGET_NAME}_$chr.lft \
            warn blat_psl/${TARGET_NAME}_$chr.psl
    done
    
    # --- Make chain files ---------------------------------------------------------------
    
    echo "Making chain files..." >&2
    
    mkdir chain_raw
    for (( c=0; c<$TARGET_N; c++ )); do
        # Add leading zero, if needed, to match chromosome name
        chr=$(printf "%02d" $c)
        axtChain -linearGap=medium  -faQ -faT -psl liftup/${TARGET_NAME}_$chr.liftup.psl \
            $SOURCE_FA $TARGET_FA chain_raw/$chr.chain
    done
     
    # --- Merge and sort chain files -----------------------------------------------------
    
    echo "Merging and sorting chain files..." >&2
    
    chainMergeSort chain_raw/*.chain | chainSplit chain_split stdin
    
    # --- Figure out sizes of genomes ----------------------------------------------------
    
    echo "Figuring out genome size..." >&2
    
    faSize $SOURCE_FA -detailed > $SOURCE_NAME.chr_length.txt
    faSize $TARGET_FA -detailed > $TARGET_NAME.chr_length.txt
     
    # --- Make alignment nets from chain files -------------------------------------------
    
    echo "Making alignment nets from chain files..." >&2
    
    mkdir nets
    for i in chain_split/*.chain; do
        tag=${i/chain_split\//}
        touch nets/$tag.net
        chainNet $i $SOURCE_NAME.chr_length.txt \
            $TARGET_NAME.chr_length.txt nets/$tag.net /dev/null
    done
     
    # --- Create liftOver chain file -----------------------------------------------------
    
    echo "Creating liftOver chain file..." >&2
    
    mkdir over
    for i in chain_split/*.chain; do
        tag=${i/chain_split\//}
        netChainSubset nets/$tag.net $i over/$tag.chain
    done
     
    cat over/*.chain > ${SOURCE_NAME}_to_${TARGET_NAME}.over.chain
    
    # --- Clean up -----------------------------------------------------------------------

    rm -r split_genomes
    rm -r blat_psl
    rm -r liftup
    rm -r chain_raw
    rm -r chain_split
    rm -r nets
    rm -r over
 
}

# --- Call function to make chainfile source-to-target -----------------------------------

make_chainfile

# Test:
# liftOver baboon.bed  ${SOURCE_NAME}_to_${TARGET_NAME}.over.chain rhesus.bed unMapped

# --- Call function again to make reverse chainfile target-to-source ---------------------

# Switch target and source

ORIG_SOURCE_FA=$SOURCE_FA		# papAnu2.RM.fa
ORIG_TARGET_FA=$TARGET_FA		# rheMac3.RM.fa

ORIG_SOURCE_NAME=$SOURCE_NAME	# papAnu2
ORIG_TARGET_NAME=$TARGET_NAME	# rheMac3

SOURCE_FA=$ORIG_TARGET_FA		# rheMac3.RM.fa
TARGET_FA=$ORIG_SOURCE_FA		# papAnu2.RM.fa

SOURCE_NAME=$ORIG_TARGET_NAME	# rheMac3
TARGET_NAME=$ORIG_SOURCE_NAME	# papAnu2

TARGET_N=`grep "^>" $TARGET_FA | wc -l`

make_chainfile

# --- Call function again to use other version of target (rheMac2) -----------------------

SOURCE_FA=$ORIG_SOURCE_FA		# papAnu2.RM.fa
TARGET_FA=$TARGET_2_FA			# rheMac2.RM.fa

SOURCE_NAME=$ORIG_SOURCE_NAME	# papAnu2
TARGET_NAME=$TARGET_2_NAME		# rheMac2

TARGET_N=`grep "^>" $TARGET_FA | wc -l`

make_chainfile

# --- Call function again to make reverse chainfile [other version] of target-to-source --

# Switch target and source

SOURCE_FA=$TARGET_2_FA			# rheMac2.RM.fa
TARGET_FA=$ORIG_SOURCE_FA		# papAnu2.RM.fa

SOURCE_NAME=$TARGET_2_NAME		# rheMac2
TARGET_NAME=$ORIG_SOURCE_NAME	# papAnu2

TARGET_N=`grep "^>" $TARGET_FA | wc -l`

make_chainfile
