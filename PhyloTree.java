import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Scanner;

/*
 * Defines a phylogenetic tree, which is a strictly binary tree 
 * that represents inferred hierarchical relationships between species
 * 
 * Iain Maryanow
 * Holden Matheson
 */

public class PhyloTree {
    private PhyloTreeNode overallRoot;    // The actual root of the overall tree
    private int printingDepth;            // How many spaces to indent the deepest 
                                          //    node when printing
    
    // Stores all the nodes by label, prevents having to traverse the tree
    // looking for each node
    private HashMap<String,PhyloTreeNode> hm = new HashMap<String,PhyloTreeNode>();


    // - A linked tree structure representing the inferred hierarchical
    //   species relationship has been created, and overallRoot points to
    //   the root of this tree
    public PhyloTree(String speciesFile, int printingDepth) {
    	this.printingDepth = printingDepth;
    	buildTree(loadSpeciesFile(speciesFile));
    }

    public PhyloTreeNode getOverallRoot() {
      return overallRoot;
    }

    public String toString() {
      return this.toString(overallRoot, 0.0, this.getWeightedHeight());
    }

    private String toString(PhyloTreeNode node, double weightedDepth, double maxDepth) {
      StringBuffer tree = new StringBuffer();
    	double k = (double)this.printingDepth * (weightedDepth / maxDepth);
    	
    	// Reverse the in-order traversal of the tree, building a string representation of the tree
    	if (node.getRightChild() != null) {
    	   tree.append(this.toString(node.getRightChild(), 
    						weightedNodeDepth(node.getRightChild()), 
    					   maxDepth));
    	}
    	
    	for(int i = 0; i < k; i++) {
    		tree.append('.');
    	}
    	tree.append(node.toString() + "\n");
    	
    	if (node.getLeftChild() != null) {
    		tree.append(this.toString(node.getLeftChild(), 
    						weightedNodeDepth(node.getLeftChild()),
    						maxDepth));
    	}
    		
      return tree.toString();
   }

   public String toTreeString() {
      return this.toTreeString(overallRoot);
   }

   private String toTreeString(PhyloTreeNode node) {
      StringBuffer tree = new StringBuffer();
    	
    	//Travels down the tree until it hits a leaf/species. 
      //Creates a tree string format from all species in tree.
    	if (node.isLeaf()){
    		tree.append(node.toString() + ":" + new java.text.DecimalFormat("0.00000").format(node.getParent().getDistanceToChild()));
    	}
    	else {
    		tree.append("(" + this.toTreeString(node.getRightChild()) + "," +
    					this.toTreeString(node.getLeftChild()) + ")");
    		if (node != this.overallRoot) {
    			tree.append(":" + new java.text.DecimalFormat("0.00000").format(node.getParent().getDistanceToChild()));
    		}
    	}
    	
      return tree.toString();
    }

    public int getHeight() {
        return this.nodeHeight(overallRoot);
    }

    public double getWeightedHeight() {
        return this.weightedNodeHeight(overallRoot);
    }

    public int countAllSpecies() {
        return overallRoot.getNumLeafs();
    }

    public java.util.ArrayList<Species> getAllSpecies() {
    	ArrayList<Species> aSpecies = new ArrayList<Species>();
    	getAllDescendantSpecies(overallRoot, aSpecies);
      
      return aSpecies;
    }

    public PhyloTreeNode findTreeNodeByLabel(String label) {
    	return hm.get(label);
    }

     public PhyloTreeNode findLeastCommonAncestor(String label1, String label2) {
        return findLeastCommonAncestor(hm.get(label1), hm.get(label2));
    }
    
    //  - If both nodes can be found: returns the sum of the weights 
    //    along the paths from their least common ancestor to each of
    //    the two nodes
   public double findEvolutionaryDistance(String label1, String label2) {   	 
      PhyloTreeNode p1, p2, p3, p4;
      p1 = hm.get(label1);
      p2 = hm.get(label2);
      double distance = 0.0;

      if (p1 == null || p2 == null) {
         return Double.POSITIVE_INFINITY;
      }
      else if (p1.equals(p2)) {
         return distance;
      }
      else {
         p3 = findLeastCommonAncestor(label1, label2);
        	
        	// If neither node is the least common ancestor(LCA)
        	// traverse from each to the LCA
        	if (!(p1.equals(p3) || p2.equals(p3))){
	        	p4 = p1.getParent();
	        	distance += p4.getDistanceToChild();
            
	        	while (!p3.equals(p4)){
	        		p4 = p4.getParent();
	        		distance += p4.getDistanceToChild();
	        	}

	        	p4 = p2.getParent();
	        	distance += p4.getDistanceToChild();
            
	        	while (!p3.equals(p4)){
	        		p4 = p4.getParent();
	        		distance += p4.getDistanceToChild();
	        	}
        	}
        	
        	// if p1 is the LCA, traverse from node 2
        	else if(p1.equals(p3)){
        		p4 = p2.getParent();
	        	distance += p4.getDistanceToChild();
            
	        	while (!p3.equals(p4)){
	        		p4 = p4.getParent();
	        		distance += p4.getDistanceToChild();
	        	}
        	}
        	
        	// if p2 is the LCA, traverse from node 1
        	else if(p2.equals(p3)){
        		p4 = p1.getParent();
	        	distance += p4.getDistanceToChild();
            
	        	while (!p3.equals(p4)){
	        		p4 = p4.getParent();
	        	   distance += p4.getDistanceToChild();
	         }
         }

         return distance;
      }
   }

   //  - Creates a linked tree structure representing the inferred hierarchical
   //    species relationship.
   private void buildTree(Species[] species) {
      ArrayList<PhyloTreeNode> forest= new ArrayList<PhyloTreeNode>();
    	MultiKeyMap<Double> distances = new MultiKeyMap<Double>();
    	PhyloTreeNode pNode, sNode1, sNode2;
    	double distance = 0.0;
    	
    	String[] s = new String[3];
    	
    	// For every species in the array, add a new node with them to the forest
    	// and store the node by label in a hash map
    	for(Species i : species){
    		pNode = new PhyloTreeNode(null, i);
    		forest.add(pNode);
    		hm.put(i.getName(), pNode);
    	}
      
    	// For every node in the forest, calculate the distance to 
      // every other node and place it in the MultiKeyMap distances
    	for(PhyloTreeNode i : forest){
    		for(PhyloTreeNode j : forest){
    			if (!i.equals(j) && (distances.get(i.getLabel(), j.getLabel()) == null)){
    				distances.put(i.getLabel(),j.getLabel(), Species.distance(i.getSpecies(), j.getSpecies()));
    			}
    		}
    	}
    	
    	// Until there is only the final node in the tree, keep creating 
      // new parent nodes from nodes that have the smallest distance to each other. 
    	while (forest.size() != 1){
    		s = distances.minDist().split("\\|");
         
    		sNode1 = findTreeNodeByLabel(s[0]);
    		sNode2 = findTreeNodeByLabel(s[1]);
         
    		if(s[0].compareTo(s[1]) < 0) {
    			pNode = new PhyloTreeNode(s[0] + "+" + s[1], null, sNode1, sNode2, 
                                     (Double.parseDouble(s[2])/2.0));
         }
    		else {
    			pNode = new PhyloTreeNode(s[1] + "+" + s[0], null, sNode2, sNode1, 
                                     (Double.parseDouble(s[2])/2.0));
    		}
         
         hm.put(pNode.getLabel(), pNode);
    		forest.remove(sNode1);
    		forest.remove(sNode2);
    		sNode1.setParent(pNode);
    		sNode2.setParent(pNode);
         
    		for(PhyloTreeNode i : forest){
    			if (i != null){
	    			distance = (double)(sNode1.getNumLeafs()) / 								
	    						  (double)(sNode1.getNumLeafs() + sNode2.getNumLeafs()) *  		
	    						  distances.get(i.getLabel(), s[0]) + 							
	 					   	  (double)(sNode2.getNumLeafs()) / 
	 					   	  (double)(sNode1.getNumLeafs() + sNode2.getNumLeafs()) * 
	 					   	  distances.get(i.getLabel(), s[1]);
                          
	    			distances.put(pNode.getLabel(), i.getLabel(), distance); 					
					distances.remove(i.getLabel(), s[0]);
					distances.remove(i.getLabel(), s[1]);
    			}
    		}
         
    		distances.remove(s[0], s[1]);
    		forest.add(pNode);
    	}
      
    	overallRoot = forest.get(0);
    }

   public static int nodeDepth(PhyloTreeNode node) {
      int count = 0;
      
    	if (node != null) {
    	   PhyloTreeNode parent = node.getParent();
         
        	while (parent != null){
        	   count++;
        	   parent = parent.getParent();
        	}
        	return count;	
      }
      
      return -1;
   }
    
   public static int nodeHeight(PhyloTreeNode node) {
      if (node != null){
         int l,r;
         
        	if (node.isLeaf()) {
        		return 0;
         }
        	else {
        		l = nodeHeight(node.getLeftChild());
        		r = nodeHeight(node.getRightChild());
        	   return 1 + (l >= r ? l : r); 
         }
      }
        
      return -1;
   }

   public static double weightedNodeHeight(PhyloTreeNode node) {
      if (node != null) {
         double l,r;
         
        	if (node.isLeaf()) {
        		return 0.0;
         }
        	else{
        		l = weightedNodeHeight(node.getLeftChild()) + node.getDistanceToChild();
        		r = weightedNodeHeight(node.getRightChild()) + node.getDistanceToChild();
        		return (l >= r ? l : r); 
        	}
      }
      
      return Double.NEGATIVE_INFINITY;
   }

   public static Species[] loadSpeciesFile(String filename) {
      boolean newSpecies = true;
    	String line = "";
    	String[] aLine;
    	String name = "";
    	StringBuffer sequence = new StringBuffer(); 
    	ArrayList<Species> aSpecies = new ArrayList<Species>();
    	
    	try {
	      Scanner fs = new Scanner(new File(filename));
         
	      if (fs.hasNextLine()) {
	   	   line = fs.nextLine();
	      }
         
		   // For all lines in the file, when a '>' is found, 
         //a new species will begin, the name being the last element in the current line.
	      while (newSpecies && fs.hasNextLine()) {
	   	   sequence.setLength(0);
	   	   newSpecies = false;
   		   aLine = line.split("\\|");
   		   name = aLine[aLine.length-1];
   				
   				
   		   // The species sequence will be all lines concatenated before another '>' is found.
   		   while(!newSpecies && fs.hasNextLine()){
		      line = fs.nextLine();
	            if(line.charAt(0) == '>'){
				      newSpecies = true;
				      String[] seq = new String[sequence.length()];
                  
				      for(int i = 0; i < sequence.length(); i++){
					      seq[i]= Character.toString(sequence.charAt(i));
				      }
			         aSpecies.add(new Species(name, seq));
					}
					else{
						sequence.append(line);
					}
   		   }
	      }
         
	   	String[] seq = new String[sequence.length()];
			for(int i = 0; i < sequence.length(); i++){
				seq[i]= Character.toString(sequence.charAt(i));
			}
			aSpecies.add(new Species(name, seq));
	   		
	   	fs.close();
	   }
	   catch (FileNotFoundException e){
	      System.err.println("Error: Unable to open file " + filename);
	      System.exit(1);
	   }
	   catch (NoSuchElementException e){
	      System.err.println("Error: Unable to parse file " + filename);
	      System.exit(2);
	   }
      
    	return (Species[])aSpecies.toArray(new Species[aSpecies.size()]);
   }

/*    // getAllDescendantSpecies
    // Pre-conditions:
    //    - node points to a node in a phylogenetic tree structure
    //    - descendants is a non-null reference variable to an empty arraylist object
    // Post-conditions:
    //    - descendants is populated with all species in the subtree rooted at node
    //      in in-/pre-/post-order (they are equivalent here)*/
    private static void getAllDescendantSpecies(PhyloTreeNode node,java.util.ArrayList<Species> descendants) {
       // Recursively finds all species from node to bottom of tree.
       if(node != null){
        	if (node.isLeaf())
        		descendants.add(node.getSpecies());
        	else{
        		getAllDescendantSpecies(node.getLeftChild(), descendants);
        		getAllDescendantSpecies(node.getRightChild(), descendants); 
        	}
        }
    }
/*    // findTreeNodeByLabel
    // Pre-conditions:
    //    - node points to a node in a phylogenetic tree structure
    //    - label is the label of a tree node that you intend to locate
    // Post-conditions:
    //    - If no node with the label exists in the subtree, return null
    //    - Else: return the PhyloTreeNode with the specified label 
    // Notes:
    //    - Assumes labels are unique in the tree*/
    private static PhyloTreeNode findTreeNodeByLabel(PhyloTreeNode node, String label) {
    	
    	// Method is not used, we store the nodes in a HashMap by label to increase lookup speed.
    	// Traverses the tree, looking for the node with the target label.
        if (node.isLeaf()){
        	if(node.getLabel().equals(label))
        		return node;
        	else
        		return null;
        }
        else{
        	PhyloTreeNode l, r;
        	l = findTreeNodeByLabel(node.getLeftChild(), label);
        	r = findTreeNodeByLabel(node.getRightChild(), label);
        	return l == null ? r : l;
        }
    }

/*    // findLeastCommonAncestor
    // Pre-conditions:
    //    - node1 and node2 point to nodes in the phylogenetic tree
    // Post-conditions:
    //    - If node1 or node2 are null, return null
    //    - Else: returns the PhyloTreeNode of their common ancestor 
    //      with the largest depth*/
     /*private static PhyloTreeNode findLeastCommonAncestor(PhyloTreeNode node1, PhyloTreeNode node2) {
    	 
    	// Recursively finds all ancestors which node1 and node2 share and then returns the deepest (least common)
        if(node1 == null || node2 == null)
    	 return null;
        else{
        	PhyloTreeNode l,r;
        	if(node1.equals(node2))
        		return node1;
        	l = findLeastCommonAncestor(node1.getParent(), node2);
        	r = findLeastCommonAncestor(node1, node2.getParent());
        	if (l == null)
        		return r;
        	else if (r == null)
        		return l;
        	else
        		return nodeDepth(l) > nodeDepth(r) ? l : r;
        }
    }*/
    
    private PhyloTreeNode findLeastCommonAncestor(PhyloTreeNode node1, PhyloTreeNode node2) {
      return findLeastCommonAncestor(overallRoot, node1, node2);
    }
     
   private static PhyloTreeNode findLeastCommonAncestor(PhyloTreeNode root, PhyloTreeNode node1, PhyloTreeNode node2) {
      if (root == null) {
         return null;
      }
      
      if (root.equals(node1) || root.equals(node2)) {
         return root;
      }
      
      PhyloTreeNode left = findLeastCommonAncestor(root.getLeftChild(), node1, node2);
      PhyloTreeNode right = findLeastCommonAncestor(root.getRightChild(), node1, node2);
      
      if (left != null && right != null) {
         return root;
      }
      
      return (left != null) ? left : right;
   } 
    
     /*    // weightedNodeDepth
     // Pre-conditions:
     //    - node is null or the root of tree (possibly subtree)
     // Post-conditions:
     //    - If null: returns -1
     //    - Else: returns the depth of the node within the overall tree*/
     public static double weightedNodeDepth(PhyloTreeNode node) {
    	 
    	// Walks up from node to overallRoot summing the weights along the way. 
        double i = 0;
     	if (node != null){
     		PhyloTreeNode parent = node.getParent();
         	while (parent != null){
         	i += parent.getDistanceToChild();
         	parent = parent.getParent();
         	}
         	return i;	
         }
         else
         	return -1;
     }
}
