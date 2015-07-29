package edu.hawaii.mahdi.slimquick;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;

public class Signature{
	
	public static int instanceCounter =0;
	
	public int nbSequences;
	
	public boolean visited= false;
	public double signature;
	
	public int signatureId=0;
	public Integer[] sequences;
	
	public int bandId; // the band number that was used to generate this signature
	
//	Comparator<Integer> c = new Comparator<Integer>() {
//	      public int compare(Integer s1, Integer s2) {
//	    	return (s1).compareTo(s2);
//	      }
//	      
//	};
	
	
	public Signature(double signature){
		this.signature = signature;
		//signatureId is its sequential number
		signatureId = instanceCounter; 
		instanceCounter++;
		sequences = new Integer[500];
		
		//System.out.println("this sig is "+signatureId+"max number of sigs is "+instanceCounter);
		nbSequences=0;

	}
	
//	public void addSequenceSign(int seq){
//		sequences[nbSequences] = seq;
//		nbSequences++;
//	}
	
	/**
	 * Finds the poisiton of the sequence id seq in the list of sequences 
	 * belonging to the current signature
	 * @param  url  an absolute URL giving the base location of the image
	 * @param  name the location of the image, relative to the url argument
	 * @return      the image at the specified URL
	 * @see         Image
	 */	
	public int getPosition(int seq){
		return java.util.Arrays.binarySearch(Arrays.copyOfRange(this.sequences, 0, this.nbSequences), seq);
	}
	
	public void addSequence(int seqId){
		this.sequences[this.nbSequences] = seqId;
		this.nbSequences++;
		
	}
	

}
