package edu.rice.cs.bioinfo.programs.phylonet.algos.network.test;

import static edu.rice.cs.bioinfo.programs.phylonet.algos.network.test.InvariantsLitmusTest.Call_Create_Network_Helper;
import static edu.rice.cs.bioinfo.programs.phylonet.algos.network.test.InvariantsLitmusTest.GenerateDgenStatistic;

public class DGEN2 {

    public static void main(String[] args) {


        /*
        String treeString = "(((P1,P2),(P3,P4)),O);";
        double introgressionProb = 0.9;
        String introgressionList = "[('P1', 'P3')]";

        String outGroupName = "O";
        String saveStatHere = "/Users/leo/rice/res/data/dgen/tmp/fig3/stat_5taxp1p3_fig3.txt";
        int numberOfRandomRuns = 100;

        //gonna use the python helper for now (which is what ill tell users to use if they want help with simple networks
        String networkString = Call_Create_Network_Helper(treeString,introgressionList,introgressionProb);
*/

        String treeString = args[0]; //etc do for all
        String networkString = args[1];
        String outGroupName = args[2];
        String saveStatHere = args[3];
        int numberOfRandomRuns = Integer.parseInt(args[4]);

        //call the java code to actually create the statistic
        GenerateDgenStatistic(treeString,networkString,outGroupName,saveStatHere,numberOfRandomRuns);

        int debugTesting = 0;

    }

}
