clone_reducer
=============
<b>Author</b>: Michael R. McKain<br>
</br>
Version 1.0, October 2014
<br></br>
<b>Contact</b>: https://github.com/mrmckain
<h3>Description</h3>

Clone_reducer identifies clades in a gene tree that comprise a single accession/species with a bootstrap value greater than or equal to what the user chooses. The clade is reduced to a consensus sequence and a new alignment file is written. The consensus is based on majority rules.

Clone reducing script used in Estep et al. (2014). DOI will be provided when available.

<h4>Output</h4>

<b>Condensed Clones File</b>:<br></br>
	Final alignment file with sequences reduced.

<b>Logfile</b>:<br></br>
	Log of clones that were condensed given accession identifier.

<b>Removed Clones</b>:<br></br>
	Sequences of clones that were condensed.

<h4>Usage</h4><br></br>

perl clone_reducer.pl -alignment file -tree file [options] 

Options:

	-alignment 			Alignment file in FASTA format used to produce gene tree<br></br>
	-tree 				Newick tree file created from alignment
	-percent_length 		Percent in decimal format of unaligned sequence to alignment length [Default: 0.5]
	-bootstrap 			Minimum bootstrap value for clade to be considered for consensus sequence [Default: 50]
	-help 				Brief help message
	-man 				Full documentation
