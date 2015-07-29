package edu.hawaii.mahdi.slimquick;

import java.security.SignatureSpi;
import java.util.ArrayList;

import org.mapdb.*;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.io.BufferedReader;
import java.io.Console;
import java.io.FileReader;
import java.io.IOException;

import com.beust.jcommander.JCommander;
import com.google.common.collect.*;
import com.google.common.collect.Multiset.Entry;
import com.google.common.primitives.Ints;

import java.util.concurrent.ConcurrentNavigableMap;
import java.util.logging.Logger;


// TODO: When traversing the signatures, maintain a a hash table of 0/1 indicating whether seq X is similar to seq Y.
// TODO: make sure no sequences are left with -1. so far, removeSingletons can remove all the signatures of a sequence if it is alone in there.


public class SlimQuick {
	private final static Logger log = Logger.getLogger(SlimQuick.class.getName());
	
	

	
	/**
	 * Test whether sequences in a signature share less than MIN_NB_COM_SIGS of overlapping signatures 
	 * Major Limitation: Two very similar sequence (sharing > MIN_NB_COM_SIGS) can save a signature
	 * containing potentially thousand of other dissimilar sequence
	 * Perhaps remove sequences that share tons of similarity??
	 * @param sig	The signature to test
	 * @param seqs	The list containing all the Sequence objects
	 * @param signatures The map containing all the Signature objects
	 * @return true if the sequence contains sequences all with edges of ones.
	 */
	public static boolean findSignatureWithAllOneEdge(Signature sig, Sequence[] seqs, 
					Map<Double, Signature> signatures, JCommanderParms jct){
		//ex.of limitation sig = Sig1, MIN_NB_COM_SIGS=3
		// Seq1 => [Sig1, Sig2, Sig3, Sig4, ], Seq2 => [Sig1, Sig2, Sig3, Sig5], 
		// Seq3 => [Sig1, Sig2, Sig3,  ] Seq4 => [Sig4, Sig, Sig3 ]
		
		Multiset<Integer> neighborSigs = HashMultiset.create();
		
		
		// Collect all the signature associated with the sequences in sig.
		for(int i =0; i < sig.nbSequences; i++){
			Sequence seq = seqs[sig.sequences[i]];
			
			for (int seqSigId: seq.sigIds){
				if (seqSigId != -1){
					// add to the multiset increments the key by one
					neighborSigs.add(seqSigId);
				}
			}
		}
		// Do not count the current signature in the number of overlapping sigs.
		neighborSigs.remove(sig.signature);

		// TODO: How can we speed this up?
		int countNbNonOne=0;
		for (Entry<Integer> count: neighborSigs.entrySet()){
			if (count.getCount()>1){
				countNbNonOne++;
				if (countNbNonOne >=jct.min_nb_com_sigs){
					return false;
				}
			}
		}
		return true;
	}	
	/**
	 * Drop large cluster 
	 * @param seqs The list containing all the Sequence objects
	 * @param signatures The map containing all the Signature objects
	 * @param bands list of the bands for each trial
	 */
	public static void droppingLargeClusters(Sequence[] seqs, Map<Double, 
							Signature> signatures, int[][] bands, JCommanderParms jct){
		//TODO: move this to Signature class?
		List<Signature> largeSigs = new ArrayList<Signature>();
		for (Signature sig: signatures.values()){
			if (sig.nbSequences > jct.max_sig_size){

				largeSigs.add(sig);
			}
		}
		for (Signature sig: largeSigs){
			removeSignature(sig, seqs, signatures);
		}
	}
	/**
	 * compares seq1 to all the sequences whose ids are in tempList and computes the number 
	 * of matches sequences. A matches sequence is hard-coded as sequence having more than 
	 * MIN_NB_KMER_MATCHES, and less than 31 nonMatches and less than 31 sumNonMatches
	 * @param seq1	int id of the sequences of interest
	 * @param temp	List list of sequences to which we are comparing
	 * @param seqs	The list containing all the Sequence objects
	 * @return	number of sequences with matches in tempList
	 */	
	//TODO:keep this in memory to speed things up?
//	public static int nbKmerMatches2(int seq1, int[] tempList,  List<Sequence> seqs){
//		int nbSeqMatches=0;
//		Sequence mySeq1 = seqs.get(seq1);		
//		
//		for (int seq2: tempList){
//			if (seq1 == seq2)
//				continue;
//			
//			Sequence mySeq2 = seqs.get(seq2);
//
//			
//			int nbKmerMatches = 0;
//			int nonMatches =0;
//			int sumNonMatches=0;
//			List<Integer> mismatches= new ArrayList<Integer>(); // positions of mismatches
//
//		
//			for (int i=0; i< mySeq1.kmersVec.length; i++){
//				
//				if ((mySeq1.kmersVec[i] == mySeq2.kmersVec[i])){
//					// TODO make sure that the signature is not mostly 0's for those sequences in that cluster
//					nbKmerMatches++;
//				}else{
//					nonMatches++;
//					sumNonMatches += Math.abs(mySeq1.kmersVec[i] - mySeq2.kmersVec[i]);
//					mismatches.add(i);
//				}
//			}
//			
//			
//			
//			//TODO modify this hard-coded requirement to  a dynamic one
//			if ((nbKmerMatches >= MIN_NB_KMER_MATCHES) && (nonMatches < 31) && (sumNonMatches < 31) ){
//				nbSeqMatches++;
////				System.out.println("MATCH FOUND "+nbKmerMatches+"- Seq1:"+mySeq1.seqId+" Seq2:"+mySeq2.seqId
////						+" have "+nonMatches+" nonMatches "+sumNonMatches+" sumNonMatches"+ " and "+sumNonMatches_6k+" sumNonMatches_6k");
//			}else{
////				System.out.println("NO MATCHES "+nbKmerMatches+"- Seq1:"+mySeq1.seqId+" Seq2:"+mySeq2.seqId
////						+" have "+nonMatches+" nonMatches "+sumNonMatches+" sumNonMatches"+ " and "+sumNonMatches_6k+" sumNonMatches_6k");
//				
//			}
//		}
//		return nbSeqMatches;
//	}

	
	
	/**
	 * Determines whether two sequences are similar, i.e., both sequences have at least minSimilarSigs
	 * non-null signatures
	 * @param seq1Id
	 * @param seq2Id
	 * @param minSimilarSigs
	 * @param seqs
	 * @return
	 */
	
	//TODO: STORE THIS IN MEMEORY 
	public static boolean areSimilar(int seq1Id, int seq2Id, int minSimilarSigs, Sequence[] seqs){		
		Sequence seq1 = seqs[seq1Id];
		Sequence seq2 = seqs[seq2Id];
		int nbSimilarSigs = 0;
		for (int i=0; i<seq1.sigIds.length; i++){
			if ((seq1.sigIds[i]!=-1) && (seq2.sigIds[i] !=-1) && (seq1.sigIds[i]==seq2.sigIds[i])){
				nbSimilarSigs++;
			}
		}
		//TODO Move this in the For loop
		if (nbSimilarSigs >= minSimilarSigs){
			//System.out.println("Seq1 "+seq1.seqId+" seq2 "+seq2.seqId+" have "+nbSimilarSigs+" shared signatures");
			return true;
		}else{
			//System.out.println("Seq1 "+seq1.seqId+" seq2 "+seq2.seqId+" are not similar and have only "+nbSimilarSigs+" shared signatures");
			return false;
		}
	}
	
	
	
	/**
	 * For each sig, we make sure that each sequence has at least 1 other sequence to which it is similar
	 * in another signature, if not, the sequence is dropped from the signature
	 * sig is dropped if it has not other sequences
	 * @param bands
	 * @param seqs	The list containing all the Sequence objects
	 * @param signatures
	 * @param correct_clusters
	 */

	public static void dropAbherrentSigs(Sequence[] seqs, Map<Double, Signature> signatures, 
										int[] correct_clusters, JCommanderParms jct){
		
		System.out.println("Total number of signatures to process for abherrations:"+signatures.size());

		int totalProcessed =0;
		
		List<Integer> outliersList = new ArrayList<Integer>();
		List<Double> sigsToRemove=new ArrayList<Double>();

		for (Signature sig: signatures.values()){
			outliersList.clear();
			
			
			// if signature has more than MIN_NUM_COMPS
			// then randomly select a number of reference sequences to compare all the sequences against
			
			List<Integer> reps =  new ArrayList<Integer>(Arrays.asList(sig.sequences)).subList(0, sig.nbSequences);
			
			Multiset<Integer>  counts = HashMultiset.create();	// number of similar references for each sequence

			
			if (sig.nbSequences > jct.sub_sample_size){
				Collections.shuffle(reps);
				reps = reps.subList(0, jct.sub_sample_size);				
				for (int i =0; i < sig.nbSequences; i++){
					int seq1id = sig.sequences[i];
					for(int rep: reps){						
						// if seq1id and rep have at least min_nb_com_sigs, then increase their counts by 1
						if ((seq1id != rep) && (areSimilar(seq1id, rep, jct.min_nb_com_sigs, seqs) == true)){
							counts.add(seq1id);
							counts.add(rep);
						}
					}
				}
			}else{
				Integer[] seqList = sig.sequences;
				for (int i=0; i<sig.nbSequences-1; i++){
					for (int j=i+1; j<sig.nbSequences; j++){
						// if seqs i  and j have at least min_nb_com_sigs, then increase their counts by 1
						if (areSimilar(seqList[i], seqList[j], jct.min_nb_com_sigs, seqs) == true){
							counts.add(seqList[i]);
							counts.add(seqList[j]);
						}
					}						
				}
			}
			
			// If a sequence does not have at least n=1 (hard coded) other sequences
			// whith which it share at least min_nb_com_sigs, then remove it from the signature
			for (int i=0; i<sig.nbSequences; i++){
				int seqId = sig.sequences[i]; 
				//TODO make the percentage a user def. value Instead of hard code it
				// TODO should update the counts, after the reads is removed but cannot do it here to not slow things down
				if(	counts.count(seqId) < 1){
					//System.out.println("Sequence:"+seqId+" Does not have enough matches: REMOVED from sig: "+sig.signature);
					outliersList.add(seqId);
				}
			}	
			
			
			
//			for (Entry a : counts.entrySet()){
//				if (a.getCount() == 1){
//						//System.out.println("Total number of outliers is: "+ outliersList.size());
//						for (int i = 0; i< sig.nbSequences; i++ ){
//							Sequence s = seqs[sig.sequences[i]];
//							System.out.println("\t"+s.seqId+"\t"+Arrays.toString(s.sigIds));
//						}
//				}
//			}
			
			
			if(outliersList.size()>0){
				removeSignatureFromSeqs(sig, Ints.toArray(outliersList), seqs, signatures);
			}
			
			//if we removed all sequences form the signature then get rid of it
			if(sig.nbSequences<=1){
				sigsToRemove.add(sig.signature);
			}
			
			totalProcessed++;
			if (totalProcessed % 10000 ==0){
				System.out.println("\nNumber of signatures processed  for abherrations "+totalProcessed);
			}else if (totalProcessed % 100 ==0) {
				System.out.print(". ");
			}
		}
		System.out.println("Removing "+sigsToRemove.size()+" abherrant signatures");
		for (Double sig: sigsToRemove){
			signatures.remove(sig);
		}
	}
	

	
	/**
	 * reads a cd-hit formatted output file and return an array representing the clustering of the sequences
	 * each position in the array represent the cluster to which the reads at that position belongs
	 * @param fileName the correct clustering-formatted output file
	 * @param nbSequences in the file
	 * @return
	 * @throws IOException
	 */
    public static int[] getCorrectClustering(String fileName, int nbSequences) throws IOException{
    	int[] clusters = new int[nbSequences];
    
    	BufferedReader br = new BufferedReader(new FileReader(fileName));
        try {
            String line = "";
            int clust = -1;
            line = br.readLine();
            while (line != null) {
                if ( (line.length() > 7)  && line.substring(0, 7).equals("Cluster")){
                	clust = Integer.parseInt(line.split(" ")[1]);
                	
                }else{
                	clusters[Integer.parseInt(line)] = clust;
                }
                line = br.readLine();
            }
        } finally {
            br.close();
        }
        return clusters;
    }


    /**
     * remove all the signatures having only one sequence assigned to it.
     * @param seqs	The list containing all the Sequence objects
     * @param signatures	The map containing all the Signature objects
     */
	public static void removeSingletons(Sequence[] seqs, Map<Double, Signature> signatures){
		/* Sets all the signatures of size one of a sequence to null.
		 * */

		List<Signature> suspicousSigs = new ArrayList<Signature>();
		for(Signature sig: signatures.values()){
			//System.out.println("\t\t"+sig.nbSequences);
			if ((sig.nbSequences == 1)){
				suspicousSigs.add(sig);
			}
		}		
		for (Signature sig: suspicousSigs){
			removeSignature(sig, seqs, signatures);
		}
	}
	
	

	
	
	/**
	 * AAAAAAA
	 * @param seqs
	 * @param signatures
	 * @return
	 */
	
	public static int[] findPartitions(Sequence[] seqs, int nbSequences,  Map<Double, Signature> signatures){
		int[] clusters = new int[nbSequences];
		LinkedList<Integer> queue = new LinkedList<Integer>();
		Set<Integer> neighborsSet = new HashSet<Integer>();
		

		
		
		Map<Integer, Double> idToSig = new HashMap<Integer, Double>(signatures.size());
		for (Map.Entry<Double, Signature> sig: signatures.entrySet()){
			idToSig.put(sig.getValue().signatureId, sig.getKey());
		}
		
		// Set all values by default to -1
		for (int i=0; i< clusters.length; i++){
			clusters[i] = -1;
		}
		int nextClustId = -1;
		for (Sequence seq: seqs){	
			if (clusters[seq.seqId]==-1){
				// Assign it to its own cluster
				// Cannot belong to a previous cluster, or else it would have been labeled through another sequence
				nextClustId++;
				clusters[seq.seqId]= nextClustId;
			}else{
				continue;
			}
			// Add the sequence to a queue so as to inspect its children
			queue.push(seq.seqId);
			//System.out.println("Queued : "+ seq.seqId);

			while(queue.isEmpty()==false){
				int newSeqId = queue.removeLast();
				// For each signature that this sequence is involved in
				// get all the other sequences that are in the signatures
				for (int sigId: seqs[newSeqId].sigIds){
					if (sigId == -1){
						continue;
					}
					Signature sig = signatures.get(idToSig.get(sigId));

				
//				}catch( NullPointerException e){
//						System.out.println("Dying");
//						System.out.println("Working with seqId" + sigId);
//						System.out.println("ifToSig is: " + idToSig);
//						System.out.println("sigId is "+sigId + " signature is "+ idToSig.get(sigId));
//						
//						System.exit(0);
//						
//						
//					}
					
					
//					// once a signature is visited, or if null (had only one read in it)
					// we do not need to process it further
					if (sig == null || sig.visited==true)
							continue;

					// search for the position where the sequence 
					int myPosition = sig.getPosition(newSeqId);
					//int myPosition = 0;
					int seqSig;
					// We start at the position of the sequence since all sequences
					// before it would have been assigned already
					// since seqIds are incremental
					// we are sure any sequence before it has already a clusterId 
					for (int i = myPosition; i<sig.nbSequences; i++){
						seqSig = sig.sequences[i];
						if(clusters[seqSig] == -1){
							// if sequence was not assigned to a cluster
							// then add it to the queue to process other sequences in its cluster
							neighborsSet.add(seqSig);
						}
					}
					sig.visited=true;
				}
				// add the neighbors (sequences sharing this signature) to the queue 
				// so that we can add their neighbors
				// TODO: add this in the main loop?
				for(int seqId: neighborsSet){
					clusters[seqId] = clusters[seq.seqId];
					queue.addFirst(seqId);
					//System.out.println("Queued : "+seqId+" through: "+ seq.seqId);
				}
				neighborsSet.clear();
				//if ((queue.size() % 10000)==0)
				//System.out.println(" queue is: "+queue.size()+" in size");
			}
			
			if ((seq.seqId % 10000)==0)
				System.out.println(" X ");
		}
		System.out.println("");
		
		// System.out.println("seq: "+i+" is in cluster: "+clusters[i]);
		System.out.println("Number of clusters generated is: "+(nextClustId+1));
		return clusters;
	}

	/**
	 * Removes the signature from the signatures list and sets its pointers in the sequence to null
	 * @param sig signature to remove
	 * @param seqs	The list containing all the Sequence objects
	 * @param signatures	The map containing all the Signature objects
	 */
	public static void removeSignature(Signature sig, Sequence[] seqs, Map<Double, Signature> signatures){
		//System.out.println("remove sequence for sig"+ sig.signatureId+" which has "+sig.nbSequences+" sequences");
		
		
		for (int i =0; i < sig.nbSequences; i++){
			//REMOVE THIS CHECK AFTER VALIDATING THE seqId = seqs.get(seqId).seqId
			// MAHDI: Untested MOD AFTER CHANGE of seqs from ArrayList to [] 
			seqs[sig.sequences[i]].sigIds[sig.bandId]=-1;
		}
		
		signatures.remove(sig.signature);
	}
	/**
	 * Removes the signature from the signatures list and sets its pointers in the sequence to null
	 * @param sig signature to remove
	 * @param seqs	The list containing all the Sequence objects
	 * @param signatures	The map containing all the Signature objects
	 */
	public static void removeSignature(Double sig, Sequence[] seqs, Map<Double, Signature> signatures){
		//System.out.println("remove sequence for sig"+ sig.signatureId+" which has "+sig.nbSequences+" sequences");
		Signature signature = signatures.get(sig);
		for (int i =0; i < signature.nbSequences; i++){
			//REMOVE THIS CHECK AFTER VALIDATING THE seqId = seqs.get(seqId).seqId
			// MAHDI: Untested MOD AFTER CHANGE of seqs from ArrayList to []
			if (signature.sequences[i] != -1){
				seqs[signature.sequences[i]].sigIds[signature.bandId]=-1;
			}
		}
		signatures.remove(signature.signature);
	}
	
	
	
	/**
	 * remove a signature from a defined list of signatures. Same as removeSignature but restricted to seqIds
	 * @param sig	signature to remove
	 * @param seqIds list of seq ids to remove the seqeunces from
	 * @param seqs	The list containing all the Sequence objects
	 * @param signatures	The map containing all the Signature objects		
	 */
	public static void removeSignatureFromSeqs(Signature sig, int[] seqIds, Sequence[] seqs, Map<Double, Signature> signatures){
		
		//System.out.println("Getting ready to remove" + Arrays.toString(seqIds));
		for (int seqId: seqIds){
			//REMOVE THIS CHECK AFTER VALIDATING THE seqId = seqs.get(seqId).seqId
			if (seqId != seqs[seqId].seqId){
				System.out.println("ERROR IN removeSiganture: seqId != seqs.get(seqId).seqId ");
				System.exit(0);
			}
			int position = sig.getPosition(seqId);
			sig.sequences[position]= -1;
			seqs[seqId].sigIds[sig.bandId]=-1;
		}
		sig.nbSequences = sig.nbSequences - seqIds.length;
	}

	/**
	 * A signature is significant if it adds something new to the partial clustering
	 * i.e., assigns a sequence to a non- -1 cluster
	 * assigns two sequences previously in different clusters to the same cluster
	 * @param partialClustering
	 * @param sig
	 * @return
	 */
	public static boolean isSignificant(int[] partialClustering, Signature sig, int nextClusId){
		boolean significant = false;

//		System.out.println("before");
//		for (int i=0; i< sig.nbSequences; i++){
//			System.out.print(partialClustering[sig.sequences[i]]+" ");
//		}
//		System.out.print("\n");
//		
		
		int currentClust = partialClustering[sig.sequences[0]];
		// TODO: speed up by looking up only half the positions, i, and comparins i with i+1
		for (int i=0; i< sig.nbSequences; i++){
			int seqId = sig.sequences[i];
			if ((partialClustering[seqId] != currentClust) || (partialClustering[seqId] == -1) ){
				
				significant = true;
			}
			partialClustering[seqId] = nextClusId;
		} 
		
//		System.out.println("after");
//		for (int i=0; i< sig.nbSequences; i++){
//			System.out.print(partialClustering[sig.sequences[i]]+" ");
//		}
//		System.out.print("\n");
//
//		
//		System.out.println("significant is "+ significant);
		return significant;

	}
	
	
	public static void main(String[] args) throws IOException {
		JCommanderParms jct = new JCommanderParms();
		JCommander jCommander = new JCommander(jct, args);
		jCommander.setProgramName("SlimQuick");
	        if (jct.help) {
	        	jCommander.usage();
	            return;
	        }	

		
		
		int[] correct_clusters =null;
		if (jct.validate == true)
			correct_clusters =	getCorrectClustering(jct.correctClusteringFile, jct.nbSequences);
		
		
		long startTime = System.nanoTime();
		
		log.info("-- Reading Input File: " + jct.inFileName);
		Fasta myFasta = new Fasta(jct.inFileName, jct.nbSequences, jct.nbTrials);
		myFasta.parse(jct.nbSequences);
		System.out.println("-- Finished parsing sequences: " + myFasta.getNumSequences());
		
		
		Sequence[] seqs = myFasta.getSequences();		
	

		//Map<Double, Signature> signatures = new HashMap<Double, Signature>(jct.nbSequences * jct.nbTrials); // allocate max number of possible signatures
		
		DB db = DBMaker.memoryDB().cacheHashTableEnable().cacheSize(10000000).make();
		
		Map<Double, Signature> signatures = db.treeMap("test");
		int [][] bands = myFasta.getRandomSamples(jct.nbTrials, jct.subsetSize);

		Map<Double, Signature> tempSignatures = new HashMap<Double, Signature>(jct.nbSequences);
		Set<Double> newSignatures = new HashSet<Double>();
		
		
		log.info("Number of bands is: "+bands.length);
		
		
		long s = System.nanoTime();
		int subsetNumber = 0; // this is the band number

		
			
		// No more than than max_sig_size in a sig, if not, drop it at the end of the iteration
		Set<Double> largeSigs = new HashSet<Double>();
		int nextTempClsId=jct.nbSequences+1;
		// partial partitions to keep track of which signatures add value and which don't.
		int[] partialPartitions= new int[jct.nbSequences];
		for (int i=0; i< jct.nbSequences; i++){
			partialPartitions[i]=-1;
		}
		
		
		
		for (int[] band: bands){
			log.info("SubsetNumber: "+subsetNumber+"\t"+Arrays.toString(band));
			int nbSequencesProcessed = 0;
			
			for (Sequence seq: seqs){
				// if signature does not exist, create it and add it to seq
				// need to check because of sig collisions
				Double sig = myFasta.returnSignature(seq.kmersVec, band, subsetNumber);
				Signature signature = (Signature) tempSignatures.get(sig);
				
				
				if (signature == null){
					if (!largeSigs.contains(sig)){
						signature = new Signature(sig);
						signature.bandId=subsetNumber;
						tempSignatures.put(sig, signature);
						newSignatures.add(sig);
					}else{
						continue;
					}
				}
				if (signature.nbSequences < jct.max_sig_size){
						signature.addSequence(seq.seqId);
						seq.sigIds[subsetNumber] = signature.signatureId;
				}else{
					// remove the signature but also keep track of it so that we don't add to 
					// it again
					if (! largeSigs.contains(signature.signature)){
						removeSignature(signature, seqs, tempSignatures);
						largeSigs.add(signature.signature);
						
					}
				}

				
				
				// Drop large signatures
				
				nbSequencesProcessed++;
				if (nbSequencesProcessed % 100000 == 0)
					System.out.println("Processed : "+nbSequencesProcessed+" nb Sigs: " + (signatures.size() + tempSignatures.size()));
			}
			
			
			//clear the large cluster since there is no overlap in signatures across iterations	
			largeSigs.clear();
			
			
			
			System.out.println("Sorting out the utility of new signatures");
			for (Double tempSig: newSignatures){
					// does tempSig add significant info?
					// if so, add keep it, or else, drop it
				if (tempSignatures.containsKey(tempSig)){
					if (! isSignificant(partialPartitions, tempSignatures.get(tempSig), nextTempClsId)){
						removeSignature(tempSig, seqs, tempSignatures);
					}else{
						for(int i=0; i<tempSignatures.get(tempSig).nbSequences; i++){
							partialPartitions[i] =nextTempClsId; 
						}
						
					}
					nextTempClsId++;
				}
			}
			newSignatures.clear();
			System.out.println("Done sorting out the utility of new signatures");
				
			System.out.println("Removing");
			removeSingletons(seqs, tempSignatures);

			signatures.putAll(tempSignatures);
			System.out.println(" total number of signatures in this iteration is "+tempSignatures.size());
			tempSignatures.clear();
			
			
//			if (subsetNumber % jct.partialClean == 0){
//				System.out.println("removing Singletons");
//				removeSingletons(seqs, signatures);	
//			}
			
			subsetNumber++;	

		}
		
		
		// No need to do this since we remove do this in the loop independently for each signature
		// remove the singletons from the remaining bands
		// log.info("Dropping Singleton Signatures: "+(System.nanoTime() - startTime)/1000000000);
		// removeSingletons(seqs, signatures);
		
		// log.info("Dropping large signatures: "+(System.nanoTime() - startTime)/1000000000);
		// droppingLargeClusters(seqs, signatures, bands, jct); //cluster is a signature!!!!!
		// log.info("Identifying outliers in all signatures"+(System.nanoTime() - startTime)/1000000000);


//		log.info("Validating signatures of size 2 "+((System.nanoTime() - s)/1000000000));
//		List<Signature> sigsToRemove = new ArrayList<Signature>();
//		int nbTwos = 0;		
//		for(Signature sig: signatures.values()){
//			if ((sig.nbSequences==2)){
//				if(findSignatureWithAllOneEdge(sig, seqs, signatures, jct)==true){
//					nbTwos++;
//					sigsToRemove.add(sig);
//				}
//			}
//		}
//		for (Signature sig: sigsToRemove){
//			removeSignature(sig, seqs, signatures);
//		}
//		log.info("Removing "+nbTwos+" signatures of size 2");
		
		
		
		
		System.out.println(" total number of signatures before dropAnherrations is: "+signatures.size());
		dropAbherrentSigs(seqs, signatures, correct_clusters, jct);
		System.out.println(" total number of signatures After dropAnherrations is: "+signatures.size());

		System.exit(0);
		
		
	

		
		
		
		
		




		int[] new_partitions = findPartitions(seqs, jct.nbSequences, signatures);
		
		if (jct.validate == true){			
			System.out.println("Computing Rand index\n");
			System.out.println("Done");
		}
		
		log.info("Clustering is: \n "+ Arrays.toString(new_partitions));
		
		long endTime = System.nanoTime();
		long duration = (endTime - startTime)/1000000000; 
		
		log.info("Completed Successfully "+ duration+ "seconds");

	}
}
