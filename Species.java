public class Species {
    private String name;          // A unique name associated with the species
    private String[] sequence;    // The biological sequence describing this species

    public Species(String name, String[] sequence) {
        this.name = name;
        this.sequence = sequence;
    }

    public String getName() {
        return this.name;
    }

    public String[] getSequence() {
        return this.sequence;
    }

    //  Returns the fraction of sequence elements
    //  that are different
    public static double distance(Species a, Species b) {
        String[] seq1 = a.getSequence();
        String[] seq2 = b.getSequence();
    
        if (seq1.length != seq2.length) {
            System.err.println("Error: Sequences must already be aligned");
            System.exit(5);
        } 
        
        int numDiffs = 0;
        for (int i = 0; i < seq1.length; i++) {
            if (!seq1[i].equals(seq2[i])) {
                numDiffs++;
            }
        }
        
        return ((double) numDiffs) / seq1.length;
    }
}
