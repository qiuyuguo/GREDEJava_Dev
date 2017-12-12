package GREDEdev;

import java.util.Comparator;

public class DataTypes {
    
    static class GenomeRegion{
        String chr;
        int start;
        int end;
        String direction;
        String gene;
    }
    
    static class GenomeLocus{
        String chr;
        int locus;
        String gene;
        String direction;
    }
    
    static class Peak{
    		String id;
        String chr;
        int locus;
        float score;
    }
    
    // comparator used to sort genome regions by 'start'
    static class GenomeRegionComparator implements Comparator<GenomeRegion>{
        public int compare(GenomeRegion g1, GenomeRegion g2){
            return g1.start -  g2.start;
        }
    }
    
}
