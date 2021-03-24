# nuc_progressive_aln.R
# Progressive alignment of nucleotide data using Needleman-Wunsch pairwise algorithm
#   with nodes (both leafs and internal) represented as profile matrices.
# Author: Juan J. Garcia Mesa

suppressMessages(library(Biostrings)) # default Needleman-Wunsch to construct distance matrix
suppressMessages(library(seqinr))     # read, write (handle) fasta sequences
suppressMessages(library(ape))        # create guide tree

progressive_aln = function(seqs_file) {
  # read FASTA in an XStringSet object
  fasta = readDNAStringSet(seqs_file, format="fasta", nrec=-1L, skip=0L,
                           seek.first.rec=FALSE, use.names=TRUE)

  # fixed substitution matrix
  mat = nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
  
  # names of sequences
  fasta_names = names(fasta)
  nseqs = length(fasta_names)
  
  # define sequence order by pairwise score (highest to lowest)
  pairwise_distance = c()
  for(j in 1:nseqs) {
    for(i in 1:nseqs) {
      if(i > j) {
        similarity = pairwiseAlignment(fasta[i],fasta[j], substitutionMatrix = mat,
                                  gapOpening = 0, gapExtension = 1, scoreOnly=TRUE)
        # convert similarity values to distance values (https://doi.org/10.1016/j.tcs.2009.02.023)
        distance = (fasta[i]@ranges@width*1+fasta[j]@ranges@width*1)/2-similarity
        pairwise_distance = c(pairwise_distance, distance)
      }
    }
  }
  
  # construct neighbor-joining tree
  M = matrix(0,nseqs,nseqs)
  M[lower.tri(M)] = pairwise_distance
  dimnames(M) = list(fasta_names)
  tree = nj(M)
  
  # newick = write.tree(tree)
  sub_trees = subtrees(tree)
  nclades = length(sub_trees)
  nodes_order = rep(list(NA),nclades)
  
  # get order in which sequences need to be (progressively) aligned
  for(i in 1:(nclades-1)) {
    nodes_order[[nclades-i+1]] = sub_trees[[i]]$tip.label[!sub_trees[[i]]$tip.label %in% sub_trees[[i+1]]$tip.label]
  }
  nodes_order[[1]] = sub_trees[[nclades]]$tip.label
  
  # core of progressive alignment using Needleman-Wunsch algorithm with profile matrices
  aln = NULL
  interm_aln = NULL
  seqs = read.fasta(seqs_file, forceDNAtolower = FALSE)
  
  for(i in 1:nclades) {
    j = length(nodes_order[[i]])
    if(j == 1) {
      # if there is only 1 sequence in clade, align it to main alignment (aln)
      aln = NW(aln, getSequence(seqs[[which(names(seqs) == nodes_order[[i]][j])]]))
    } else {
      while(j > 1) {
        if(is.null(aln)) {
          # if no alignment has been done yet, align first two sequences
          aln = NW(getSequence(seqs[[which(names(seqs) == nodes_order[[i]][j-1])]]),
                   getSequence(seqs[[which(names(seqs) == nodes_order[[i]][j])]]))
          j = j - 2
        } else if(is.null(interm_aln)) {
          # if alignment has been initiated but intermediate matrix is null, align next 2 seqs
          interm_aln = NW(getSequence(seqs[[which(names(seqs) == nodes_order[[i]][j-1])]]),
                          getSequence(seqs[[which(names(seqs) == nodes_order[[i]][j])]]))
          j = j - 2
        } else {
          # if both aln matrices have been used, align next seq in clade to intermediate (aux) aln
          interm_aln = NW(interm_aln, getSequence(seqs[[which(names(seqs) == nodes_order[[i]][j])]]))
          j = j - 1
        }
      }
      if(j==1) {
        # if 1 more seq (node) remaining in clade, align it with interm matrix and assign to aln
        interm_aln = NW(interm_aln, getSequence(seqs[[which(names(seqs) == nodes_order[[i]][j])]]))
        aln = NW(aln, interm_aln)
      } else if(!is.null(interm_aln)) {
        # if no more sequences in clade, align interm to aln
        aln = NW(aln, interm_aln)
      }
    }
    interm_aln = NULL
  }
  
  # correct sequence order
  alignment = aln
  for(i in 1:length(names(seqs))) {
    alignment[i,] = aln[which(unlist(nodes_order) == names(seqs)[i]),]
  }
  
  # save alignment into fasta file
	write.fasta(sequences = as.list(as.data.frame(t(alignment))), names = names(seqs),
	            file.out = paste0("aligned_",basename(seqs_file)))
  
}

# create profile sequence 
create_profile = function(seq) {
  if(is.null(dim(seq))) {   # if only one sequence is given
    seq = matrix(seq,1,length(seq))
  }
  cols = ncol(seq)
  rows = nrow(seq)
  profile = matrix(0,5,cols)
  for(i in 1:cols) {
    for(j in 1:rows) {
      if(seq[j,i] == "A") {
        profile[1,i] = profile[1,i] + 1/rows
      } else if(seq[j,i] == "C"){
        profile[2,i] = profile[2,i] + 1/rows
      } else if(seq[j,i] == "G"){
        profile[3,i] = profile[3,i] + 1/rows
      } else if(seq[j,i] == "T"){
        profile[4,i] = profile[4,i] + 1/rows
      } else {
        profile[5,i] = profile[5,i] + 1/rows
      }
    }
  }
  return(profile)
}

# Needleman-Wunsch algorithm with profile matrices
NW = function(seq1,seq2){
  if(is.null(dim(seq1))) {   # if only one sequence is given
    seq1 = matrix(seq1,1,length(seq1))
  }
  
  if(is.null(dim(seq2))) {   # if only one sequence is given
    seq2 = matrix(seq2,1,length(seq2))
  }
  
  # create profile matrices
  profile1 = create_profile(seq1)
  profile2 = create_profile(seq2)

	# Length of sequences
	m = ncol(profile1)
	n = ncol(profile2)

	# Initialize score and traceback matrices
	V = matrix(-Inf,m+1,n+1)
	traceback = matrix("-",m+1,n+1)

	# Define scores
	gap = -1
	#gap.ext = -1
	match = 1
	mismatch = -1

	# Fill first values on D that are independent
	V[1,1] = 0

	# Fill first column of V
	for(i in 2:(m+1)) {
  	V[i,1] = gap*(i-1)
  	traceback[i,1] = 2
	}

	# Fill first row of V
	for(j in 2:(n+1)) {
  	V[1,j] = gap*(j-1)
  	traceback[1,j] = 1
	}

	# Fill dynamic programming table
	for(i in 2:(m+1)) {
  	for(j in 2:(n+1)) {
    	gap1 = V[i,j-1] + gap
    	gap2 = V[i-1,j] + gap
    	mis_match = V[i-1,j-1] + match_mismatch(i,j,profile1,profile2,match,mismatch)
    	
    	V[i,j]= max(gap1,gap2,mis_match)
    	traceback[i,j] = which.max(c(gap1,gap2,mis_match))
  	}
	}

  i = m+1
	j = n+1
	back_seq = c()
	
	# Retrieve backtracking sequence (gives us number of columns)
	while(i > 1 || j > 1) {
		if(traceback[i,j] == 1) {
			# Go left
		  back_seq = c(1,back_seq)
			j = j - 1
		} else if(traceback[i,j] == 2) {
			# Go up
		  back_seq = c(2,back_seq)
			i = i - 1
		} else {
			# Go diagonal
		  back_seq = c(3,back_seq)
			i = i - 1
			j = j - 1
		}
	}

	# Retrieve alignment using backtracing sequence
	aln1 = matrix(NA,nrow(seq1),length(back_seq))
	aln2 = matrix(NA,nrow(seq2),length(back_seq))
	i = j = 1 # keep track of positions without gaps on seq1, seq2 (i,j) respectively
	
	for(k in 1:length(back_seq)) {
		if(back_seq[k] == 1) {
			# Go left
			aln1[,k] = "-"
			aln2[,k] = seq2[,j]
			j = j + 1
		} else if(back_seq[k] == 2) {
			# Go up
			aln1[,k] = seq1[,i]
			aln2[,k] = "-"
			i = i + 1
		} else {
			# Go diagonal
			aln1[,k] = seq1[,i]
			aln2[,k] = seq2[,j]
			i = i + 1
			j = j + 1
		}
	}
	
	return(rbind(aln1,aln2))
}

# calculate match/mismatch score for dynamic programming diagonal movement
match_mismatch = function(i,j,prof1,prof2,match,mismatch) {
  m = 0
  for(k in 1:4) {
    for(l in 1:4) {
      m = m + prof1[k,i-1]*prof2[l,j-1]*ifelse(k==l,match,mismatch)
    }
  }
  return(m)
}

if(!interactive()) {
	ARGS = commandArgs(trailingOnly = TRUE)
	progressive_aln(ARGS[1])
}