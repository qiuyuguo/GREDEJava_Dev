package GREDEdev;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;

import GREDEdev.DataTypes.*;

public class SelectRegDomains {
    
    private ArrayList<GenomeRegion> RefRegdomains;
    private ArrayList<String> GeneList;
    private RedBlackBST<String, GenomeRegion> RefRegDomainTable;
    private ArrayList<GenomeRegion> SelectedRegions;
    private ArrayList<GenomeRegion> RandomRegions;
    
    // Object construction
    // arg1: reg domain file (e.g., 'mm9.RegDomains.txt')
    // arg2: gene list for which the reg domains are drawn (e.g., 'tpm.NPC_0.vs.Cortex_0.neg.txt') 
    public SelectRegDomains (String arg1, String arg2){
        RefRegdomains = readRegDomains(arg1);
        GeneList = readGeneList(arg2);
        for (int i = 0; i < RefRegdomains.size(); i++){
            RefRegDomainTable.put(RefRegdomains.get(i).getGene(), RefRegdomains.get(i));
        }
    }
    
    // Read pre-calculated transcriptome-wide regulatory domains
    public ArrayList<GenomeRegion> readRegDomains(String arg){
    		RefRegdomains = new ArrayList<GenomeRegion>();
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
            RefRegdomains.add(newRegion);
        }
        return RefRegdomains;
    }
    
    // Read an input gene list
    public ArrayList<String> readGeneList(String arg){
        GeneList = new ArrayList<String>();
        In in = new In(arg);
        while(!in.isEmpty()){
            String item = in.readString();
            GeneList.add(item);
        }
        return GeneList;
    }
    
    // adjust the set of random regulatory domains so that there is no collision
    public void adjustRandomRegDomains(){
        RandomRegions = pickRandomRegDomains();
        Collections.sort(RandomRegions, new GenomeRegionComparator());
        for (int i = 1; i < RandomRegions.size(); i++){
            int start1 = RandomRegions.get(i-1).getStart();
            int end1 = RandomRegions.get(i-1).getEnd();
            int start2 = RandomRegions.get(i).getStart();
            int end2 = RandomRegions.get(i).getEnd();
            if(end2 - start1 < end1 - start1 + end2 - start2){
                int mid = (end1 + start2)/2;
                RandomRegions.get(i-1).setEnd(mid - 1);
                RandomRegions.get(i).setStart(mid);
            }
        }
    }
    
    // pick regulatory domains of a random set of genes matching the number of items in the input gene list
    public ArrayList<GenomeRegion> pickRandomRegDomains(){
        int RegCount = GeneList.size();
        int[] RandomIndexArray = new int[RegCount];
        for (int i = 0; i < RefRegdomains.size(); i++){
            RandomIndexArray[i] = i;
        }
        StdRandom.shuffle(RandomIndexArray);
        int[] RandomArrayIndexMatched = Arrays.copyOf(RandomIndexArray, RegCount);
        ArrayList<GenomeRegion> RandomArrayMatched = new ArrayList<GenomeRegion>();
        for (int i = 0; i < RandomArrayIndexMatched.length; i++){
            RandomArrayMatched.add(RefRegdomains.get(i));
        }
        return RandomArrayMatched;
    }
    
    // adjust the picked regulatory domains so that there is no collision
    public void adjustRegDomains(){
        SelectedRegions = pickRegDomains();
        Collections.sort(SelectedRegions, new GenomeRegionComparator());
        for (int i = 1; i < SelectedRegions.size(); i++){
            int start1 = SelectedRegions.get(i-1).getStart();
            int end1 = SelectedRegions.get(i-1).getEnd();
            int start2 = SelectedRegions.get(i).getStart();
            int end2 = SelectedRegions.get(i).getEnd();
            if(end2 - start1 < end1 - start1 + end2 - start2){
                int mid = (end1 + start2)/2;
                SelectedRegions.get(i-1).setEnd(mid - 1);
                SelectedRegions.get(i).setStart(mid);
            }
        }
    }
    
    // pick regulatory domains of given genes
    public ArrayList<GenomeRegion> pickRegDomains(){
        // match reg domains of given gene list
        // implemented by RedBlackBST
        SelectedRegions = new ArrayList<GenomeRegion>();
        for (int i = 0; i < GeneList.size(); i++){
            SelectedRegions.add(RefRegDomainTable.get(GeneList.get(i)));
        }
        return SelectedRegions;
    }
    
    
    // print both the selected reg domains and random reg domains with matched count
    public void printRegDomains(){
        Out selectedRDOut = new Out("target.RegDomains.txt");
        for (int i = 0; i < SelectedRegions.size(); i++){
            GenomeRegion item = SelectedRegions.get(i);
            String line = item.getChr() + "\t" + item.getStart() +"\t" + item.getEnd() + "\t" + item.getDir() + "\t" + item.getGene();
            selectedRDOut.println(line);
        }
        selectedRDOut.close();
        Out randomRDOut = new Out("random.RegDomains.txt");
        for (int i = 0; i < RandomRegions.size(); i++){
            GenomeRegion item = RandomRegions.get(i);
            String line = item.getChr() + "\t" + item.getStart() +"\t" + item.getEnd() + "\t" + item.getDir() + "\t" + item.getGene();
            randomRDOut.println(line);
        }
        randomRDOut.close();
    }
    
    // 
    public static void main(String[] args){
        SelectRegDomains SelectRegDomainsObj = new SelectRegDomains(args[0], args[1]);
        SelectRegDomainsObj.adjustRegDomains();
        SelectRegDomainsObj.adjustRandomRegDomains();
        SelectRegDomainsObj.printRegDomains();
    }
}
