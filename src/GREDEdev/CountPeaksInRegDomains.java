package GREDEdev;

import java.util.ArrayList;
import GREDEdev.DataTypes.*;
import GREDEdev.IntervalST;
import GREDEdev.Interval1D;

public class CountPeaksInRegDomains {
	
	private ArrayList<GenomeRegion> targetDomains;
	private ArrayList<Peak> peaks;
	private ArrayList<IntervalST> RegDomainST;
	
	// create object
	public CountPeaksInRegDomains (String RegDomainName, String PeakName) {
		targetDomains = ReadTargetDomains (RegDomainName);
		peaks = ReadPeaks (PeakName);
	}
	
	// read reg domains files
	public ArrayList<GenomeRegion> ReadTargetDomains (String arg){
		ArrayList<GenomeRegion> regDomains = new ArrayList<GenomeRegion>();
		In in = new In(arg);
        while(!in.isEmpty()){
            String unsplitLine = in.readString();
            String[] items = unsplitLine.split("\t");
            String chr = items[0];
            int start = Integer.parseInt(items[1]);
            int end = Integer.parseInt(items[2]);
            String direction = items[3];
            String gene = items[4];
            GenomeRegion newRegion = new GenomeRegion(chr, start, end, direction, gene);
            regDomains.add(newRegion);
        }
		return regDomains;
	}
	
	// read peak files
	public ArrayList<Peak> ReadPeaks (String arg){
		ArrayList<Peak> peaks = new ArrayList<Peak>();
		In in = new In(arg);
        while(!in.isEmpty()){
            String unsplitLine = in.readString();
            String[] items = unsplitLine.split("\t");
            String id = items[0];
            String chr = items[1];
            int locus = Integer.parseInt(items[2]);
            float score = Float.parseFloat(items[3]);
            Peak newPeak = new Peak(id, chr, locus, score);
            peaks.add(newPeak);
        }
		return peaks;
	}
	
	
	// iteratively judge if the input peaks are within the input reg domains
	
	// 1. Adding targetDomains into RegDomainST indexed by chr numbers
	public void PreprocessTargetDomains () {
		RegDomainST = new ArrayList<IntervalST>();
		for (int i = 0; i < 21; i++) {
			RegDomainST.add(new IntervalST<String>());
		}
		for (int i = 0; i < targetDomains.size(); i++) {
			GenomeRegion currentDomain = targetDomains.get(i);
			int chrIndex = 0;
			if (currentDomain.getChr() == "chrX") {
				chrIndex = 20;
			}
			else if (currentDomain.getChr() == "chrY") {
				chrIndex = 21;
			}
			else {
				chrIndex = Integer.parseInt(currentDomain.getChr().replace("chr", ""));
			}			
			Interval1D newRegion = new Interval1D(currentDomain.getStart(), currentDomain.getEnd());
			RegDomainST.get(chrIndex).put(newRegion, currentDomain.getGene());
		}
	}
	
	
	// judge whether a peak is within any of the reg domains

}
