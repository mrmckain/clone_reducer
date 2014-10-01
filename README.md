clone_reducer
=============
<h1><b>Author</b></h1>
Michael R. McKain, 2014.

Clone reducing script used in Estep et al. (2014). DOI will be provided when available.

Script looks at gene tree to identify clades with a bootstrap value of at least X that comprise the same accession/species. These clades are reduced (in the alignment) to a single consensus sequence. 

Output files include:

Logfile: Shows which sequences for what accessions were used to create a consensus sequence.

Reduced Alignment:

Removed Sequences: Full length sequences for those either used to make a consensus sequence or removed due to length.
