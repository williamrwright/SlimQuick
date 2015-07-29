package edu.hawaii.mahdi.slimquick;



public class Sequence {

	// Removed to save space .. this is not really needed beyong to create
	// the k-mers vector
	// String seqString;

	
	int seqId;
	byte[] kmersVec;
	
	//Signature[] signatures;
	int[] sigIds;
	

	
	int nbSubsets;
	
	public Sequence(int nbSubsets){
		//signatures = new Signature[nbSubsets];
		sigIds = new int[nbSubsets];
		for (int i=0; i<nbSubsets; i++){
			sigIds[i] = -1;	
		}
		
	}
	

	public void setId(int seqId) {
		this.seqId= seqId;
	}
	
//	public void setSeqString(String seqString) {
//		this.seqString = seqString;
//	}
	
	public void setKmersVec(byte[] kmersVec){
		this.kmersVec = kmersVec;
	}




	

	
}
