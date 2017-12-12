package GREDEdev;

import java.util.ArrayList;
import GREDEdev.DataTypes.*;

public class CountPeaksInRegDomains {
	
	private ArrayList<GenomeRegion> targetDomains;
	private ArrayList<Peak> peaks;
	
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
            GenomeRegion newRegion = new GenomeRegion();
            newRegion.chr = items[0];
            newRegion.start = Integer.parseInt(items[1]);
            newRegion.end = Integer.parseInt(items[2]);
            newRegion.direction = items[3];
            newRegion.gene = items[4];
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
            Peak newPeak = new Peak();
            newPeak.id = items[0];
            newPeak.chr = items[1];
            newPeak.locus = Integer.parseInt(items[2]);
            peaks.add(newPeak);
        }
		return peaks;
	}
	
	
	// iteratively judge if the input peaks are within the input reg domains
	
	
	// judge whether a peak is within any of the reg domains

}
