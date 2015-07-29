package edu.hawaii.mahdi.slimquick;


import com.beust.jcommander.Parameter;

public class JCommanderParms {
	 
	  @Parameter(names = "-in", description = "Input file name")
	  String inFileName;

	  @Parameter(names = "-nb_seqs", description = "Number of sequences in File (REMOVE THIS)")
	  Integer nbSequences;

	  @Parameter(names = "-kmer_size", description = "Kmer size, by default k=3 is used")
	  Integer kmerSize = 3;
	  
	  @Parameter(names = "-nb_trials", description = "Number of subsets trials")
	  Integer nbTrials = 128;

	  @Parameter(names = "-subset_size", description = "Number of kmers sampled in each trial")
	  Integer subsetSize = 27;
	  
	  @Parameter(names = "-partial_cleaning", description = "Number of subset before we start cleaning for singletons") 
	  Integer partialClean = 10;
	  
	  @Parameter(names = "--validate_clusters", description = "Compute cluster validation")
	  boolean validate = false;

	  @Parameter(names = "--correct_clustering", description = "True clusters of the input file. ")
	  String correctClusteringFile;
	  
	  @Parameter(names = "--nb_kmer_matches", description = "Min number of kmer matches between two sequence to be placed in same partition")
	  Integer MIN_NB_KMER_MATCHES = 220;

	  @Parameter(names = "--nb_sig_matches", description = "Min number of signature matches between two sequence to be placed in same partition")
	  Integer min_nb_com_sigs=4;

	  @Parameter(names = "--sub_sample_size", description = "Max nb seqs to consider in a large signature ")
	  Integer sub_sample_size=100;

	  @Parameter(names = "--max_sig_size", description = "Max number of seqeunces that map to a signature")
	  Integer max_sig_size=500;

	  @Parameter(names = { "-log", "-verbose" }, description = "Level of verbosity")
	  Integer verbose = 1;
	 
	  @Parameter(names = "-debug", description = "Debug mode")
	  boolean debug = false;
	  
	  @Parameter(names = "-help", help = true)
	  boolean help;

}
