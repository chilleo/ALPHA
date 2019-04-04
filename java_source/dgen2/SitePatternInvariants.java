package edu.rice.cs.bioinfo.programs.phylonet.algos.network.test;

import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.*;

import static edu.rice.cs.bioinfo.programs.phylonet.algos.network.test.InvariantsHelperFunctions.*;
import static edu.rice.cs.bioinfo.programs.phylonet.algos.network.test.InvariantsHelperFunctions.sortListByFirstElement;

public class SitePatternInvariants{



    //the site patterns (alphabetical)
    ArrayList<String> sortedSitePatterns;

    //the taxa sorted alphabetically (from left to right this is how the site patterns are arranged inside their strings
    ArrayList<String> sortedTaxa;

    // Invariant groups which are consistent across parameter settings
    //the invariants, the actual invariant groups
    List<List<Integer>> sitePatternInvariants;



    /**
     *
     * @param theTreeOrNetwork
     * @param theOutgroupName
     * @param whichEqualityCheck
     */
    public SitePatternInvariants(Network<Double> theTreeOrNetwork, String theOutgroupName, InvariantsLitmusTest.equalityCheckType whichEqualityCheck, int numberOfRandomSettings){

        List<Tree> trees = getAllTreesWithOutgroup(theTreeOrNetwork, theOutgroupName);

        //init a couple things
        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        boolean firstSetting = true;
        //the gt probabilities
        List<Double> probabilities = new ArrayList<>();

        //make a bunch of tree/networks with random BL
        List<Network> parameterSettings = getAllParametrizedNetworksRandom(theTreeOrNetwork, numberOfRandomSettings); // 100 random settings

        //initialize
        //invariants of site patterns and site pattern objects needed, more initializitions
        //the invariants, the actual invariant groups
        sitePatternInvariants = new ArrayList<>();



        //loop over trees/networks with rand BLs to find invariants
        for (Network netWithParameters: parameterSettings) {

            //initialize a couple things here that need to reset every iteration of the loop
            //the probabilities (sum of prob gt|st/sn)
            HashMap<String,Double> sitePatternProbabilities = new HashMap<>();
            //ordered probabilities (keeping everything ordered to maintain the connection between prob and site pattern
            List<Double> sitePatternProbsAlphabetical = new ArrayList<>();


            //first, get probs of gts|st/sn
            double[] probsArray = new double[trees.size()];
            try {
                gtp.calculateGTDistribution(netWithParameters, trees, null, probsArray);
            }catch(Exception e){
                int investigate = 0;
            }

            probabilities = new ArrayList<>();
            for (double d: probsArray) {
                probabilities.add(d);
            }

            //now, loop through trees, get all site patterns for a tree, and add probs to hashmap total probabilities
            for(int iIterateTrees = 0; iIterateTrees < trees.size(); iIterateTrees++){
                Tree iTree = trees.get(iIterateTrees);
                double iTreeProbability = probsArray[iIterateTrees];

                for(String iSitePattern: GetSitePatternsForTreeOrderedAlphabeticallyByTaxaName((STITree) iTree,theOutgroupName)){
                    //add value to hashmap value
                    if(!sitePatternProbabilities.containsKey(iSitePattern)) {
                        sitePatternProbabilities.put(iSitePattern,iTreeProbability);
                    }else{
                        sitePatternProbabilities.put(iSitePattern,(sitePatternProbabilities.get(iSitePattern)+iTreeProbability));
                    }
                }
            }
            //turn hashmap into ordered list for findInvariantGroups
            //a.k.a. sitePatternProbsAlphabetical = ordered alphabetically sitePatternProbabilities
            sortedSitePatterns = new ArrayList(sitePatternProbabilities.keySet());
            Collections.sort(sortedSitePatterns);
            for(String sp : sortedSitePatterns){
                sitePatternProbsAlphabetical.add(sitePatternProbabilities.get(sp));
            }



            if (firstSetting) {
                firstSetting = false;
                sitePatternInvariants = findInvariantGroups(sitePatternProbsAlphabetical,whichEqualityCheck);
            } else {
                sitePatternInvariants = updateInvariantGroups(sitePatternInvariants,sitePatternProbsAlphabetical,whichEqualityCheck,"SitePattern");
            }
//            System.out.println(probabilities);
        }

        sitePatternInvariants.sort(sortListByFirstElement);

        //once done, sitePatternInvariants has been created and filled, and sortedSitePatterns has the patterns.
        //lastly, the letters are ordered for the alphabetical list of taxa names from left to right

        //make sorted string array of taxa
        sortedTaxa = new ArrayList<String>();
        for(String anyLeaf2 : trees.get(0).getLeaves()){
            sortedTaxa.add(anyLeaf2);
        }
        Collections.sort(sortedTaxa);

    }


    public void PrintInvariants(){

        StringBuilder resultsString = new StringBuilder();

        Sayln("Invariants for the following taxa :");
        for(int i = 0; i < sortedTaxa.size(); i++){
            Say(sortedTaxa.get(i)+",");
        }
        Sayln("\n\nThe site patterns:");

        int i = 0;
        for (String thePattern: sortedSitePatterns) {
            resultsString.append( " " + i++ + ": " + thePattern);
        }

        resultsString.append("\ninvariants from random sampling and standard likelihood:");

        resultsString.append("\n" + sitePatternInvariants);

        Sayln(resultsString.toString());



    }

    public void PrintSitePatternsOnly(){
        int i = 0;
        for (String thePattern: sortedSitePatterns) {
            Say( " " + i++ + ": " + thePattern);
        }
    }

    public void PrintOrderedTaxaOnly(){
        Say(String.valueOf(sortedTaxa));
    }


    public List<List<Integer>> GetRefinedInvariants(List<Integer> oneTreeInvariant) {

        //jus to make my life easier ima put the one tree invariants into a set and just can call contains function
        Set<Integer> oneTreeInvariantsSet = new HashSet<Integer>(oneTreeInvariant);

        //make a copy of the invariants with extras removed
        List<List<Integer>> refinedInvariants = new ArrayList<List<Integer>>();
        for(List<Integer> copyMyContents : this.sitePatternInvariants){

            ArrayList<Integer> copyIntoHere = new ArrayList<Integer>();

            for(Integer copyMe : copyMyContents){
                //only add the thing in to copy it if it exists inside of oneTreeInvariants
                if(oneTreeInvariantsSet.contains(copyMe)) {
                    copyIntoHere.add(copyMe);
                }
            }

            if(!copyIntoHere.isEmpty())
                refinedInvariants.add(copyIntoHere);
        }



        return refinedInvariants;

    }

    public String toStringThisInvariantShowPatterns(List<Integer> oneInvariant) {

        StringBuilder returnMe = new StringBuilder();
        returnMe.append("[");
        int lastElementTracker = 0;
        for(Integer thisPattern : oneInvariant){
            String thePattern = this.sortedSitePatterns.get(thisPattern);
            returnMe.append("'"+thePattern+"'");

            //dont add final comma
            if(lastElementTracker != oneInvariant.size()-1)
                returnMe.append(',');

            lastElementTracker++;
        }

        returnMe.append("]");
        return returnMe.toString();
    }
}


