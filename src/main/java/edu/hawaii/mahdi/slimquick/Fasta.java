package edu.hawaii.mahdi.slimquick;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;

public class Fasta {
	private BufferedReader buffReader;

	
    private String nextline;
    private Sequence[] seqs;
    private Utils utils;
    private int nbSubsets;
    private int nbSequences;
    
    public Fasta(String filename, int nbSequences, int nbSubsets) {
   	
        try {
            this.buffReader = new BufferedReader(new InputStreamReader((InputStream) new FileInputStream(filename), Charset.defaultCharset()));
        }	
        catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        this.nextline = "";
        System.out.println("creating array for "+nbSequences+" sequences");
        this.seqs = new Sequence[nbSequences];
        this.utils = new Utils((byte) 4);
        this.nbSubsets = nbSubsets;
        this.nbSequences = nbSequences;
    }
   
    
       
    public void parse(int records) throws IOException {
    	
    	Sequence seq = null; 

    	String line;
        StringBuilder seqStr = null;
    
        if (!this.nextline.isEmpty()) {
            this.nextline = "";
            seq = new Sequence(this.nbSubsets);
            seq.setId(Integer.parseInt(this.nextline.substring(1)));
            seqStr = new StringBuilder();
        }
        int processed = 0;
        while ((line = this.buffReader.readLine()) != null) {
            if ((line = line.trim()).isEmpty()) continue;
            if (line.startsWith(">")) {
                if (seqStr != null) {
                    // seq.setSeqString(seqStr.toString());
                    // IMPORTANT ---- FOR testing purposes: reset this
                    seq.setKmersVec(this.utils.getKmers(seqStr.toString()));
                    seqs[processed] = seq;
                   
                    ++processed;
                    
                    if (processed % 10000 == 0)
                    	System.out.println(processed);
                }
                if (processed == records) {
                    this.nextline = line;
                    return;
                }
                seqStr = new StringBuilder();
                seq = new Sequence(nbSubsets);
                seq.setId(Integer.parseInt(line.substring(1)));
                continue;
            }
            seqStr.append(line);
        }
        if (seqStr != null) {
            this.nextline = "";
            //seq.setSeqString(seqStr.toString());
            // IMPORTANT ---- FOR testing purposes: reset this
            seq.setKmersVec(this.utils.getKmers(seqStr.toString()));
            seqs[processed] = seq;
        }
    }
    
    
    public double returnSignature(byte[] kmers, int [] band, int subsetNumber){
    	return utils.returnSignature(kmers, band, subsetNumber);
    	
    }
    
    public int [][] getRandomSamples(int nbSubsets, int subsetSize){
    	return utils.getRandomSamples(nbSubsets, subsetSize);
    }

    
    
    
    public int getNumSequences(){
    	return this.nbSequences;
    }
    	
    public Sequence[] getSequences(){
    	return this.seqs;
    }
	
}
