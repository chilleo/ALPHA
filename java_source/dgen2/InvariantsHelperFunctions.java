package edu.rice.cs.bioinfo.programs.phylonet.algos.network.test;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

public class InvariantsHelperFunctions {


    /**
     * @param topo         Topology which must include outgroupName as a leaf
     * @param outgroupName
     * @return A list of trees, which do not include the outgroup
     */
    public static List<Tree> getAllTreesMinusOutgroup(Network<Double> topo, String outgroupName) {
        // If this method throws an index out of bounds exception, ensure that topo contains outgroupName

        String[] theTaxa = new String[topo.getLeafCount() - 1]; //{"A", "B", "C", "D"};
        int i = 0;
        for (NetNode leaf : topo.getLeaves()) {
            if (!leaf.getName().equals(outgroupName)) {
                theTaxa[i++] = leaf.getName();
            }
        }
        Arrays.sort(theTaxa);

        List<Tree> oldTrees = Trees.generateAllBinaryTrees(theTaxa);
        List<Tree> trees = new ArrayList<>();

        for (Tree t : oldTrees) {
            // Recreate all trees so that node IDs will start at 0
            // Otherwise something goes wrong in GeneTreeProbability
            trees.add(Trees.readTree(t.toNewick()));
        }

        return trees;
    }

    /**
     * @param topo         Topology which must include outgroupName as a leaf
     * @param outgroupName
     * @return A list of trees, which include the outgroup
     */
    public static List<Tree> getAllTreesWithOutgroup(Network<Double> topo, String outgroupName) {
        List<Tree> withoutOutgroup = getAllTreesMinusOutgroup(topo, outgroupName);
        List<Tree> withOutgroup = new ArrayList<>();
        for (Tree t : withoutOutgroup) {
            withOutgroup.add(addOutgroup(t, outgroupName));
        }
        return withOutgroup;
    }

    // Check for fuzzy equality of (double) floating-point values
    public static boolean closeEnough(double x, double y, InvariantsLitmusTest.equalityCheckType whichWayToCheckEquality) {

        if (whichWayToCheckEquality == InvariantsLitmusTest.equalityCheckType.DGEN) {
            //10^-14 literal difference
            double epsilon = 0.00000000000001;
            return (Math.abs(x - y) <= epsilon);
        } else if (whichWayToCheckEquality == InvariantsLitmusTest.equalityCheckType.FuzzyLoose) {
            //10^-4
            double epsilon = 0.0001;
            double differenceMagnitude = Math.abs(x - y);
            double averageMagnitude = (Math.abs(x) + Math.abs(y)) / 2.0;
            return differenceMagnitude < epsilon * averageMagnitude;

        } else if (whichWayToCheckEquality == InvariantsLitmusTest.equalityCheckType.FuzzyStrict) {
            //10^-6
            double epsilon = 0.000001;
            double differenceMagnitude = Math.abs(x - y);
            double averageMagnitude = (Math.abs(x) + Math.abs(y)) / 2.0;
            return differenceMagnitude < epsilon * averageMagnitude;

        } else {
            //default to dgen calc
            //10^-14 literal difference
            double epsilon = 0.00000000000001;
            return (Math.abs(x - y) <= epsilon);
        }

        //will probably change this back to just checking that they are within a small threshold at some point (abs(x-y) < tiny num)
        //done now with option for different ones
    }

    // Find indices of (trees with) equal probabilities
    public static List<List<Integer>> findInvariantGroups(List<Double> probabilities, InvariantsLitmusTest.equalityCheckType whichEqualityType) {
        List<List<Integer>> groups = new ArrayList<>();

        for (int i = 0; i < probabilities.size(); i++) {
            boolean foundAGroup = false;
            for (List<Integer> group : groups) {
                if (closeEnough(probabilities.get(i), probabilities.get(group.get(0)), whichEqualityType)) {
                    group.add(i);
                    foundAGroup = true;
                    break;
                }
            }
            if (!foundAGroup) {
                List<Integer> newGroup = new ArrayList<>();
                newGroup.add(i);
                groups.add(newGroup);
            }
        }

        return groups;
    }

    public static List<List<Integer>> updateInvariantGroups(List<List<Integer>> oldGroups, List<Double> probabilities, InvariantsLitmusTest.equalityCheckType whichEqualityCheck, String debugType) {
        List<List<Integer>> newGroups = new ArrayList<>();
        for (List<Integer> group : oldGroups) {
            while (group.size() > 0) {
                List<Integer> newGroup = new ArrayList<>();
                newGroup.add(group.get(0));
                group.remove(0);
                for (Integer i : group) {
                    if (closeEnough(probabilities.get(i), probabilities.get(newGroup.get(0)), whichEqualityCheck)) {
                        newGroup.add(i);
                    }else{
                        //a refinement has happened, lets print to screen to see whats going on right now for how many random runs we need (why are the refinements happening)
                        //Sayln(debugType+" "+i+" has been removed from the group containing "+newGroup.get(0));
                        //Sayln(debugType+" "+i+" probability: "+probabilities.get(i));
                        //Sayln(debugType+" "+newGroup.get(0)+" probability: "+probabilities.get(newGroup.get(0)));
                    }
                }
                for (Integer i : newGroup) {
                    if (group.contains(i)) {
                        group.remove(group.indexOf(i));
                    }
                }
                newGroups.add(newGroup);
            }
        }
        return newGroups;
    }


    public static List<Network> getAllParametrizedNetworks(Network topology, int numRandomSettings) {
        List<Network> parametrizedNetworks = new ArrayList<>();

        for (int i = 0; i < numRandomSettings; i++) {
            Network<Double> newNet = randomizeBranchLengths(topology);
            parametrizedNetworks.add(newNet);
        }
        return parametrizedNetworks;
    }

    public static Network randomizeBranchLengths(Network topology) {
        double[] branchLengths = {0.5, 1, 2, 4};
        Network<Double> newNet = topology.clone();

        for (NetNode<Double> node : newNet.getTreeNodes()) {
            for (NetNode parent : node.getParents()) {
                double newLength;
                newLength = branchLengths[(int) (Math.random() * branchLengths.length)];
//                newLength = 4.0;
//                newLength = Math.random() * 4; //10;
                node.setParentDistance(parent, newLength);

            }
        }
        for (NetNode<Double> hybridNode : newNet.getNetworkNodes()) {
            for (NetNode parent : hybridNode.getParents()) {
                double newLength;
                newLength = 0.0;
//                newLength = branchLengths[(int) (Math.random() * branchLengths.length)];
//                newLength = NetNode.NO_DISTANCE; // THIS IS NEGATIVE INFINITY, DO NOT USE
                hybridNode.setParentDistance(parent, newLength);
            }
        }
        return newNet;
    }

    public static List<Network> getAllParametrizedNetworksRandom(Network topology, int numRandomSettings) {
        List<Network> parametrizedNetworks = new ArrayList<>();

        for (int i = 0; i < numRandomSettings; i++) {
            Network<Double> newNet = randomizeBranchLengthsUniformBounded(topology);
            parametrizedNetworks.add(newNet);
        }
        return parametrizedNetworks;
    }

    public static Network randomizeBranchLengthsUniformBounded(Network topology) {
        Network<Double> newNet = topology.clone();

        //current min and max values are hard coded (ten and 10^-6  right now)
        double minValue = 0.000001;
        double maxValue = 10.0;

        for (NetNode<Double> node : newNet.getTreeNodes()) {
            for (NetNode parent : node.getParents()) {
                double newLength = minValue + (maxValue - minValue) * (new Random()).nextDouble();
                node.setParentDistance(parent, newLength);
            }
        }
        for (NetNode<Double> hybridNode : newNet.getNetworkNodes()) {
            for (NetNode parent : hybridNode.getParents()) {
                double newLength;
                newLength = 0.0;
//                newLength = branchLengths[(int) (Math.random() * branchLengths.length)];
//                newLength = NetNode.NO_DISTANCE; // THIS IS NEGATIVE INFINITY, DO NOT USE
                hybridNode.setParentDistance(parent, newLength);
            }
        }
        return newNet;
    }


    public static Tree addOutgroup(Tree tree, String outgroupName) {
        TMutableNode oldRoot = (TMutableNode) tree.getRoot();
        STITree newTree = (STITree) Trees.readTree(";");
        STITree outgroup = (STITree) Trees.readTree(outgroupName);

        newTree.getRoot().adoptChild(oldRoot);
        newTree.getRoot().adoptChild(outgroup.getRoot());

        //TODO adopt child not working correctly in updating tree properties like set of nodes etc, adding this inefficient line to deal with it for now
        newTree = (STITree) Trees.readTree(newTree.toNewickWD());
        return newTree;
    }

    public static ArrayList<String> GetSitePatternsForTreeOrderedAlphabeticallyByTaxaName(STITree sayMyPatterns, String outgroupName) {

        //debugging print statement
        //System.out.println("DOING TREE "+sayMyPatterns);

        ArrayList<HashMap<String,String>> returnMe = new ArrayList<HashMap<String,String>>();

        //make sure to root it here
        sayMyPatterns.rerootTreeAtNode(sayMyPatterns.getNode(outgroupName));

        Map<STITreeCluster,TNode> theClusters = sayMyPatterns.getClusters(sayMyPatterns.getLeaves());

        for(TNode patternThisCluster : theClusters.values()){

            //if its a leaf node we dont need it and if its everything we dont need it, otherwise we make a pattern

            if(patternThisCluster.getLeafCount() == 1 || patternThisCluster.getLeafCount() == sayMyPatterns.getLeafCount()-1){
                //do nothing (testing string below)
                //System.out.println("NOT COUNTING "+ patternThisCluster.toString()+" W #LEAVES EQUAL " +patternThisCluster.getLeafCount());
            }else{
                //System.out.println("WE ARE COUNTING "+ patternThisCluster.toString()+" W #LEAVES EQUAL " +patternThisCluster.getLeafCount());

                HashMap<String,String> taxaToABValue = new HashMap<String,String>();
                for(TNode bLeaf : patternThisCluster.getLeaves()){
                    taxaToABValue.put(bLeaf.getName(),"B");
                }
                for(String anyLeaf : sayMyPatterns.getLeaves()){
                    if(!taxaToABValue.containsKey(anyLeaf)){
                        taxaToABValue.put(anyLeaf,"A");
                    }
                }
                returnMe.add(taxaToABValue);

            }

        }

        //sayMyPatterns.getCl

        //if i wanna return list of strings of patterns like i was originally thinking, extract from list of hashmaps
        ArrayList<String> returnMeStrings = new ArrayList<String>();
        ArrayList<String> sortedTaxa = new ArrayList<String>();
        for(String anyLeaf2 : sayMyPatterns.getLeaves()){
            sortedTaxa.add(anyLeaf2);
        }
        Collections.sort(sortedTaxa);
        for(HashMap onePattern : returnMe){
            StringBuilder thePatternAlphSorted = new StringBuilder();
            for(String oneSiteTaxa : sortedTaxa){
                thePatternAlphSorted.append(onePattern.get(oneSiteTaxa));
            }
            returnMeStrings.add(thePatternAlphSorted.toString());
        }



        //print out patterns for testing
        //debug printing the full site patterns list
        /*
        System.out.println("For Tree "+sayMyPatterns+" the site patterns are: ");
        for(HashMap onePattern : returnMe){
            System.out.print("PTRN: ");
            for(Object oneSiteTaxa : onePattern.entrySet()){
                System.out.print("{"+((Map.Entry)oneSiteTaxa).getKey()+",");
                System.out.print(((Map.Entry)oneSiteTaxa).getValue()+"}");
            }
            System.out.print("\n");

            System.out.print("PTRN(ALPH SORT, CONCAT STRING): ");
            StringBuilder thePatternAlphSortedToPrint = new StringBuilder();
            for(String oneSiteTaxa : sortedTaxa){
                System.out.print("{"+oneSiteTaxa+",");
                System.out.print(onePattern.get(oneSiteTaxa)+"}");
                thePatternAlphSortedToPrint.append(onePattern.get(oneSiteTaxa));
            }
            System.out.print("\n" +thePatternAlphSortedToPrint+"\n");

        }
        System.out.print("\n");
        */

        //can return in either format from here now, hashmaps of taxa to A/B site patterns for just list of strings of the patterns (where taxa are alphabetically sorted)
        //the hashmaps are in returnMe
        //the list of strings are in returnMeStrings

        return returnMeStrings;
    }



    public static Comparator<List<Integer>> sortListByFirstElement = new Comparator<List<Integer>>() {
        public int compare(List<Integer> l1, List<Integer> l2) {
            return l1.get(0) - l2.get(0);
        }
    };

    public static void Say(String sayThis){
        System.out.print(sayThis);
    }

    public static void Sayln(String sayThis){
        System.out.println(sayThis);
    }

}