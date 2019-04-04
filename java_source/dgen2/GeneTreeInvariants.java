package edu.rice.cs.bioinfo.programs.phylonet.algos.network.test;

import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static edu.rice.cs.bioinfo.programs.phylonet.algos.network.test.InvariantsHelperFunctions.*;

public class GeneTreeInvariants{

    Network<Double> theTreeOrNetwork;
    String theOutgroupName;

    List<Tree> trees;

    // Invariant groups which are consistent across parameter settings
    List<List<Integer>> consistentInvariantGroups = new ArrayList<>();

    //objects for debugging, should be able to get rid of any and comment them out at basically any time
    List<Double> debugProbabilities;

    /**
     * Constructor. Creates all the field values and finds all the invariants
     * for the given inputs and such.
     * @param givenTreeOrNetwork
     * @param givenOutgroupName
     * @param whichEqualityCheck
     * @param numberOfRandomSettings
     */
    public GeneTreeInvariants(Network<Double> givenTreeOrNetwork, String givenOutgroupName,InvariantsLitmusTest.equalityCheckType whichEqualityCheck, int numberOfRandomSettings){

        theTreeOrNetwork = givenTreeOrNetwork;
        theOutgroupName = givenOutgroupName;

        //make all gts
        trees =  getAllTreesWithOutgroup(theTreeOrNetwork, theOutgroupName);

        //initialize some stuff
        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        boolean firstSetting = true;
        List<Double> probabilities = new ArrayList<>();



        //make a bunch of trees/networks with random BL settings
        List<Network> parameterSettings = getAllParametrizedNetworksRandom(theTreeOrNetwork, numberOfRandomSettings); // 100 random settings

        //for each of these randomly parameterized networks, do the invariants check
        for (Network netWithParameters: parameterSettings) {

            //calc prob of gt give st/sn
            double[] probsArray = new double[trees.size()];
            gtp.calculateGTDistribution(netWithParameters, trees, null, probsArray);
            probabilities = new ArrayList<>();
            for (double d: probsArray) {
                probabilities.add(d);
            }

            //if first make invariant groupings, otherwise update invariant groupings
            if (firstSetting) {
                firstSetting = false;
                consistentInvariantGroups = findInvariantGroups(probabilities,whichEqualityCheck);

            } else {
                consistentInvariantGroups = updateInvariantGroups(consistentInvariantGroups, probabilities,whichEqualityCheck,"Tree");
            }
//            System.out.println(probabilities);
            debugProbabilities = probabilities;
        }

        consistentInvariantGroups.sort(sortListByFirstElement);
        //debug printing the 'correct' invariant
        //System.out.println(consistentInvariantGroups);


    }




    public void PrintInvariants(){

        StringBuilder resultsString = new StringBuilder();
        resultsString.append("\n////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n");

        resultsString.append("\ntrees:");
        int i = 0;
        for (Tree t: trees) {
            resultsString.append(" " + i++ + ": " + t);
        }
        resultsString.append("\ninvariants:");
        resultsString.append("\n" + consistentInvariantGroups);

        System.out.println(resultsString);

    }

    public void PrintGeneTreesOnly(){
        int i = 0;
        for (Tree t: trees) {
            Say(" " + i++ + ": " + t);
        }
    }

    public void GenerateExampleProbabilities(){
        //make one random network
        List<Network> parameterSettings = getAllParametrizedNetworksRandom(theTreeOrNetwork, 1); // 100 random settings
        Network netWithParameters = parameterSettings.get(0);

        //init
        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        List<Double> probabilities = new ArrayList<>();

        //calc prob of gt give st/sn
        double[] probsArray = new double[trees.size()];
        gtp.calculateGTDistribution(netWithParameters, trees, null, probsArray);
        probabilities = new ArrayList<>();
        for (double d: probsArray) {
            probabilities.add(d);
        }

        StringBuilder resultsString = new StringBuilder();
        //print probabilities
        resultsString.append("\nexample classic probs:\n");
        int i = 0;
        for (Double d: probabilities) {
            resultsString.append( i++ + ": " + d + ", ");
        }
        Say(resultsString.toString());

    }


}
