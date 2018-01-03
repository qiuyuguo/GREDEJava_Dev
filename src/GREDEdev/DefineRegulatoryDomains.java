package GREDEdev;

import java.util.ArrayList;
import GREDEdev.DataTypes.*;

public class DefineRegulatoryDomains {
    
    public static int[] gt;
    public ArrayList<GenomeRegion> RegDomains;
    public ArrayList<GenomeLocus> GenomeLoci;
    
    // create object: 
    // arg1 = genome table file
    // arg2 = TSS loci file
    public DefineRegulatoryDomains(String arg1, String arg2){
        GenomeTable(arg1);
        this.GenomeLoci = new ArrayList<GenomeLocus>();
        this.RegDomains = new ArrayList<GenomeRegion>();
        GenomeLoci = readGenomeLoci(arg2);
        buildRegDomains(GenomeLoci);
    }
    
    // construct genome table, i.e., a vector of individual chr length
    // input file: mm9.gt
    public void GenomeTable(String args){
        // ** should enable human genome (23 chr), too */
        gt = new int[21];
        In in = new In(args);
        int i = 0;
        while(!in.isEmpty()){
            String line = in.readLine();
            String[] items = line.split("\t");
            gt[i] = Integer.parseInt(items[1]);
            i++;
        }
    }
    
    // Read in TSS list and store it as array of GenomeLocus
    // input file: mm9.great3.0.txt
    public ArrayList<GenomeLocus> readGenomeLoci(String arg){
        ArrayList<GenomeLocus> GenomeLoci = new ArrayList<GenomeLocus>();
        In in = new In(arg);
        while(!in.isEmpty()){            
            String line = in.readLine();
            String[] items = line.split("\t");
            GenomeLocus item = new GenomeLocus(items[1], Integer.parseInt(items[2]), items[3], items[4]);
            GenomeLoci.add(item);
        }
        return GenomeLoci;
    }
    
    // calculate regulatory domain for items in array of GenomeLocus
    // this is based on the assumption that the gene list (gl) is pre-sorted
    public void buildRegDomains(ArrayList<GenomeLocus> gl){
        if (gl == null || gl.size() < 3){
            throw new java.lang.IllegalArgumentException();
        }
        for(int i = 1; i < gl.size()-1; i++){
            GenomeRegion item = generateRegDomain(gl.get(i), gl.get(i-1), gl.get(i+1));
            RegDomains.add(item);
        }
    }
    
    // generate regulatory domains:
    // basal range = +/- 1 kbp, extented range = 500 kb on both directions unless hitting basal range of the next gene
    // this process requires examination of 3 consecutive records from ref
    // ** should be able to take input so that to implement upstream, downsteam or both **/
    // ** should be able to take user-defined basal range and extented range */
    private GenomeRegion generateRegDomain(GenomeLocus gr, GenomeLocus grPrev, GenomeLocus grNext){
        int bR = 1000;
        int eR = 500000;
        String chr = gr.getChr();
        int start = gr.getLocus() - Math.max(bR, Math.min(Math.abs(grPrev.getLocus()+bR-gr.getLocus()), eR));
        int end = gr.getLocus() + Math.max(bR, Math.min(Math.abs(grNext.getLocus()-bR-gr.getLocus()), eR));
        String direction = gr.getDir();
        String gene = gr.getGene();
        GenomeRegion retPrimary = new GenomeRegion(chr, start, end, direction, gene);
        GenomeRegion ret = checkRange(retPrimary);
        return ret;
    }
    
    // check whether the given range is within the range of the chromosome 
    // adjust range based on the check and return the new range
    private GenomeRegion checkRange(GenomeRegion gr){
        
    		int chrN = findChr(gr.getChr());
        
    		String chr = "";
    		int start = 0;
    		int end = 0;
    		String direction = "";
    		String gene = "";
    		
    		if (gr.getStart() >= 0 && gr.getEnd() <= gt[chrN-1]){
            chr = revertChr(chrN);
            start = gr.getStart();
            end = gr.getEnd();
            direction = gr.getDir();
            gene = gr.getGene();
        }else if (gr.getStart() < 0 && gr.getEnd() <= gt[chrN-1]){
            chr = revertChr(chrN);
            start = 0;
            end = gr.getEnd();
            direction = gr.getDir();
            gene = gr.getGene();
        }else if (gr.getStart() >= 0 && gr.getEnd() > gt[chrN-1]){
            chr = revertChr(chrN);
            start = gr.getStart();
            end = gt[chrN-1];
            direction = gr.getDir();
            gene = gr.getGene();
        }
        
    		GenomeRegion ret = new GenomeRegion(chr, start, end, direction, gene);
        return ret;
    }
    
    // convert chr name to integer for array storage
    private int findChr(String chrStr){
        int chrInt;
        if (chrStr.equals("chrX")){
            chrInt = 20;
        }else if (chrStr.equals("chrY")){
            chrInt = 21;
        }else{
            //StdOut.println(chrStr);
            chrInt = Integer.parseInt(chrStr.replace("chr",""));
        }
        if (chrInt < 1 || chrInt > 21){
            throw new java.lang.IllegalArgumentException();
        }
        return chrInt;
    } 
    
    //convert chr integer to chr name for output
    private String revertChr(int chrInt){
        if (chrInt < 1 || chrInt > 21){
            throw new java.lang.IllegalArgumentException();
        }
        String chrStr;
        if (chrInt == 20){
            chrStr = "chrX";
        }else if (chrInt == 21){
            chrStr = "chrY";
        }else{
            chrStr = "chr"+Integer.toString(chrInt);
        }
        return chrStr;
    }   
    
    // print regulatory domains to file
    public void printRegDomains(){
        Out RDOut = new Out("mm9.RegDomains.txt");
        for (int i = 0; i < RegDomains.size(); i++){
            GenomeRegion item = RegDomains.get(i);
            String line = item.getChr() + "\t" + item.getStart() +"\t" + item.getEnd() + "\t" + item.getDir() + "\t" + item.getGene();
            RDOut.println(line);
        }
        RDOut.close();
    }
    
    public static void main(String[] args){
        DefineRegulatoryDomains DRDObj = new DefineRegulatoryDomains(args[0], args[1]);  
        DRDObj.printRegDomains();   
    }
    
    
}
