/*
 *
 * Creates PhyloTree objects for each FASTA file in alignment list input file
 * Infers the hierarhical structure, computes all pairwise evolutionary 
 * distances, and reports other tree statistics
 * 
 * ----------------------------------------------------------------------------
 *
 * usage:
 *
 * java Driver fastaListFilename outputDir
 *
 * where the argument is
 * 
 *   fastaListFilename      a plaintext file with one line per FASTA alignment file
 *   outputDir              a directory where the trees and statistics will be written
 *
*/

import java.util.*;
import java.io.*;

public class Driver {
    private static final int PRINTING_DEPTH = 100;

    public static void main(String[] args) {
        if (args.length != 2) {
            System.err.println("Error: Wrong number of arguments. Need fasta list file and output directory.");
            System.exit(2);
        }
    
        String fastaListFilename = args[0];
        String outputDir = args[1];
        Scanner input = null;
        File inputFile = new File(fastaListFilename);

        try {
            input = new Scanner(inputFile);
        } 
	     catch (FileNotFoundException e) {
            System.err.println("Error: Unable to open file: " + fastaListFilename);
            System.exit(1);
        }

        int numFiles = 0;
        while (input.hasNext()) {
            String fastaFilename = input.next();
            numFiles++;
            System.out.println("\nLoading tree: " + numFiles);

            File fastaFile = new File(fastaFilename);

            PhyloTree tree = new PhyloTree(fastaFilename, PRINTING_DEPTH);

            File treeOutFile = new File(outputDir + "/" + fastaFile.getName() + ".tree");
            File distOutFile = new File(outputDir + "/" + fastaFile.getName() + ".distances");
            PrintStream treeOut = null;
            PrintStream distOut = null;

            try {
                treeOut = new PrintStream(treeOutFile);
                distOut = new PrintStream(distOutFile);
            } 
	         catch (FileNotFoundException e) {
                System.err.println("Error: Unable to open output file for writing" + e);
                System.exit(1);
            }

            System.out.print(tree);
            treeOut.print(tree.toTreeString());

            ArrayList<Species> speciesList = tree.getAllSpecies();
            java.text.DecimalFormat formatter = new java.text.DecimalFormat("0.00");

            if (speciesList != null) {
                for (int i = 0; i < speciesList.size(); i++) {
                     for (int j = 0; j < speciesList.size(); j++) {
                        String label1 = speciesList.get(i).getName();
                        String label2 = speciesList.get(j).getName();
                        distOut.println("EvDistance(" + label1 + "," + label2 + ") = " + formatter.format(tree.findEvolutionaryDistance(label1, label2))); 
                    }    
                }
            }

            System.out.println("# species is " + tree.countAllSpecies());
            System.out.println("Tree height is " + tree.getHeight());
            System.out.println("Weighted height is " + formatter.format(tree.getWeightedHeight()));
        }
    }
}
