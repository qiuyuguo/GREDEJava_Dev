package GREDEdev;

import java.util.Comparator;

public class DataTypes {
    
    static class GenomeRegion{
        
    		private String chr;
        private int start;
        private int end;
        private String direction;
        private String gene;
        
        public GenomeRegion(String chr, int start, int end, String direction, String gene) {
        		this.chr = chr;
        		this.start = start;
        		this.end = end;
        		this.direction = direction;
        		this.gene = gene;
        }
        
        public String getChr() { return chr; }
        public int getStart() { return start; }
        public int getEnd() { return end; }
        public String getDir() { return direction; }
        public String getGene() { return gene; }
        
        public void setStart(int newStart) {
        		this.start = newStart;
        }
        public void setEnd(int newEnd) {
        		this.end = newEnd;
        }
        
    }
    
    static class GenomeLocus{
    		
    		private String chr;
    		private int locus;
    		private String gene;
    		private String direction;
    		
    		public GenomeLocus(String chr, int locus, String direction, String gene) {
        		this.chr = chr;
        		this.locus = locus;
        		this.direction = direction;
        		this.gene = gene;
        }
        
        public String getChr() { return chr; }
        public int getLocus() { return locus; }
        public String getDir() { return direction; }
        public String getGene() { return gene; }
        
    }
    
    static class Peak{
    	
    		private String id;
    		private String chr;
    		private int locus;
    		private float score;
    		
    		public Peak(String id, String chr, int locus, float score) {
    			this.id = id;
        		this.chr = chr;
        		this.locus = locus;
        		this.score = score;
        }
        
    		public String getId() { return id; }
        public String getChr() { return chr; }
        public int getLocus() { return locus; }
        public float getScore() { return score; }
        
    }
    
    // comparator used to sort genome regions by 'start'
    static class GenomeRegionComparator implements Comparator<GenomeRegion>{
        public int compare(GenomeRegion g1, GenomeRegion g2){
            return g1.getStart() -  g2.getStart();
        }
    }
    
}
