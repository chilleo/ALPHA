package edu.rice.cs.bioinfo.programs.phylonet.algos.network.test;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityIntegrated;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITreeCluster;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.omg.PortableInterceptor.SYSTEM_EXCEPTION;
import sun.nio.ch.Net;

import java.io.*;
import java.lang.reflect.Array;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.IntStream;

import static edu.rice.cs.bioinfo.programs.phylonet.LeoScript_FullPaperCwR.GetOutputFromProgramNoHangFromDir;
import static edu.rice.cs.bioinfo.programs.phylonet.algos.network.test.InvariantsHelperFunctions.Sayln;

/**
 * Created by hunter on 9/5/18.
 */
public class InvariantsLitmusTest {

    private static Comparator<List<Integer>> sortListByFirstElement = new Comparator<List<Integer>>() {
        public int compare(List<Integer> l1, List<Integer> l2) {
            return l1.get(0) - l2.get(0);
        }
    };

    public static void main(String[] args) {



        DoWholePaper();

        /*
        LeoTesting();

        Network<Double> treeForAbba = Networks.readNetwork("(O, (G, (H, C)));");
        Network treeForAbbaWithoutOG = Networks.readNetwork("(C, (B, A));");
//        Network netForAbba = Networks.readNetwork();
        //abba baba tree case
        runTest(treeForAbba, "O");

        List<Tree> testTrees = getAllTreesMinusOutgroup(treeForAbba, "O");
        testTrees = getAllTreesWithOutgroup(treeForAbba, "O");

        //abba baba network case <- the values do not look right, needs to be looked into why
        //Network<Double> netForAbba = Networks.readNetwork("(O:1.0, ((H:1.0,(C:1.0)#H1:1::0.3):1.0,G:1.0,#H1:1::0.7));");
        double inheritanceProbability = 0.1;
        Network<Double> netForAbba = Networks.readNetwork("(O, ((G, #H1:0::" + inheritanceProbability + "), (H, (C)#H1:0::" + (1 - inheritanceProbability) + ")));");
        runTest(netForAbba,"O");

        //dfoil fig tree case
        Network<Double> treeForDFOIL = Networks.readNetwork("(((P1, P2), (P3, P4)), O);");
        //dfoil fig tree case
        runTest(treeForDFOIL, "O");


        //dfoil fig network case
        //add reticultion between P1 and P3 or P1 and P4 (can use whatever way hunter has been using on his own)
        inheritanceProbability = 0.1;
        Network<Double> netForDFOIL = Networks.readNetwork("(((P1, (P2)#H1:0::" + inheritanceProbability + "), ((P3, #H1:0::" + (1 - inheritanceProbability) + "), P4)), O);");
        //dfoil fig net case
        runTest(netForDFOIL, "O");
        */


        //leo note - not really going past here for now

//        Network<Double> netForFullTest = Networks.readNetwork("(A, (B, (C, (D, E))));");
//        for (Tree t: getAllTreesMinusOutgroup(netForFullTest,"A")) {
//            Network netForTest = Networks.readNetwork(t.toNewick());
//            System.out.println();
//            runTest(netForTest,"A");
//        }
//        Network<Double> balancedNet = Networks.readNetwork("((A,C)hnode2,(B,D)hnode1);");
//        System.out.println(getAllTrees(balancedNet));
//        runTest(balancedNet);
    }

    private static void DoWholePaper() {

        //CichlidMakeBlocks(); //testing 10k and 100k first chunk and realign and seeing much better outgroup alignment

        //MakeAbbaBaba(); //the sanity check and also makes the d statistic in case i need it

        //RunExperiment1Fig3(); // DONE FOR NOW (fig3) - and rerun with block -em
        //RunExperiment2Fig4(); //DONE FOR NOW - and rerun with block -em
        //RunExperimentFig7(); //DONE - reran with -em start and explicit end
        //RunExperimentCichlids(); //DONE //MAKING THIS. for stat gen and single pass run (verbose) 3:24-3:26,

        int nextrundebasjkdlfjkdf = 0;




    }

    private static void RunExperimentCichlids() {

        //*******************************
        //******** start from first line and edit onwards once i get the actual file i need
        //*******************************

        //had to change it to a phylip? (also there are many letters besides actg (Y etc))

        //delete everything in the folder storing the results
        //rm -rf /path/to/directory/*
        try {


            //String cmdPath = "/Users/leo/rice/res/alpha/alpha/CommandLineFiles";
            GetOutputFromProgramNoHangFromDir("rm -rf /Users/leo/rice/res/data/dgen/tmp/figC/*",null);

            GetOutputFromProgramNoHangFromDir("rm /Users/leo/rice/res/data/dgen/tmp/figC/fig_figC.txt",null);
            //GetOutputFromProgramNoHangFromDir("echo \"D or p Value, Statistic, Significant\" > /Users/leo/rice/res/data/dgen/tmp/fig7/fig_fig7.txt",null);

        } catch (IOException e) {
            e.printStackTrace();
        }

        String treeString = "(((Cfus,(Cdec,Ceja)),(Cmam,Cgui)),Sgui);";
        double introgressionProb = 0.9;
        String introgressionList = "[('Cmam', 'Ceja')]";

        String outGroupName = "Sgui";
        String saveStatHere = "/Users/leo/rice/res/data/dgen/tmp/figC/stat_cichlid.txt";
        int numberOfRandomRuns = 100;

        //gonna use the python helper for now (which is what ill tell users to use if they want help with simple networks
        String networkString = Call_Create_Network_Helper(treeString,introgressionList,introgressionProb);

        //call the java code to actually create the statistic
        //MOVED TO TRY BLOCK TO BE RIGHT NEXT TO OTHER MAIN THINGS HAPPENING
        //GenerateDgenStatistic(treeString,networkString,outGroupName,saveStatHere,numberOfRandomRuns);

        //and then finally call python code to run the statistic on some actual data

        //migRate modulation, the first part
        //String simPath;

        //simPath = "/Users/leo/rice/res/data/dgen/simulations/4taxa/withGeneFlow/multipleReticulations/emStartAndEnd_3to2_4to3_6tax/";

        //create list of files to pass into dgen run
        //StringBuilder filesList = new StringBuilder();
        //filesList.append("/Users/leo/rice/res/data/cichlid/alignment/cichlid6tax.phylip-sequential.txt");

        String plotName = "/Users/leo/rice/res/data/dgen/tmp/figC/plot_figC";
        String resultsName = "/Users/leo/rice/res/data/dgen/tmp/figC/results_figC.txt";

        String dgenCommand; // = "["+filesList+"]";

        //newer test, test running newest version at cmd line
        // WORKS!
        dgenCommand = "python -c \"from RunDGEN import *; run_saved_dgen('"+saveStatHere+"', ['/Users/leo/rice/res/data/cichlid/alignment/cichlid6tax.phylip-sequential.txt'], verbose=True, plot='"+plotName+"', meta='"+"Dgen"+"')\" > "+resultsName;

        //command to run window based analysis
        int windowSize = 10000;
        int windowOffset = 500;
        String plotNameWindows = "/Users/leo/rice/res/data/dgen/tmp/figC/plot_figC_windows";
        String resultsNameWindows = "/Users/leo/rice/res/data/dgen/tmp/figC/results_figC_windows.txt";

        String dgenCommandWindows = "python -c \"from RunDGEN import *; run_saved_dgen('"+saveStatHere+"', ['/Users/leo/rice/res/data/cichlid/alignment/cichlid6tax.phylip-sequential.txt'], window_size="+windowSize+", window_offset="+windowOffset+", verbose=True, plot='"+plotNameWindows+"', meta='"+"Dgen"+"')\" > "+resultsNameWindows;

        //System.out.println(dgenCommand);

        String dResults;
        String dResultsWindows;

        //doing now also for 100kbp first chunk unaligned and aligned w mafft
        String dResults100kbpUnaligned;
        String dResults100kbpAligned;
        String plotNameWindowsUnal = "/Users/leo/rice/res/data/dgen/tmp/figC/plot_figC_windowsUnal";
        String resultsNameWindowsUnal = "/Users/leo/rice/res/data/dgen/tmp/figC/results_figC_windowsUnal.txt";
        String plotNameWindowsAlig = "/Users/leo/rice/res/data/dgen/tmp/figC/plot_figC_windowsAlig";
        String resultsNameWindowsAlig = "/Users/leo/rice/res/data/dgen/tmp/figC/results_figC_windowsAlig.txt";

        String dgenCommandWindows100kbpUnaligned = "python -c \"from RunDGEN import *; run_saved_dgen('"+saveStatHere+"', ['/Users/leo/rice/res/data/cichlid/alignment/cichlidMediumUnaligned.phylip-sequential.txt'], window_size="+windowSize+", window_offset="+windowOffset+", verbose=True, plot='"+plotNameWindowsUnal+"', meta='"+"Dgen"+"')\" > "+resultsNameWindowsUnal;
        String dgenCommandWindows100kbpAligned = "python -c \"from RunDGEN import *; run_saved_dgen('"+saveStatHere+"', ['/Users/leo/rice/res/data/cichlid/alignment/cichlidMediumAligned.phylip-sequential.txt'], window_size="+windowSize+", window_offset="+windowOffset+", verbose=True, plot='"+plotNameWindowsAlig+"', meta='"+"Dgen"+"')\" > "+resultsNameWindowsAlig;

        //doing now also also for 200kbp first chunk unaligned and aligned w mafft
        String dResultsUnaligned200k;
        String dResultsAligned200k;
        String plotNameWindowsUnal200k = "/Users/leo/rice/res/data/dgen/tmp/figC/plot_figC_windowsUnal200k";
        String resultsNameWindowsUnal200k = "/Users/leo/rice/res/data/dgen/tmp/figC/results_figC_windowsUnal200k.txt";
        String plotNameWindowsAlig200k = "/Users/leo/rice/res/data/dgen/tmp/figC/plot_figC_windowsAlig200k";
        String resultsNameWindowsAlig200k = "/Users/leo/rice/res/data/dgen/tmp/figC/results_figC_windowsAlig200k.txt";

        String dgenCommandWindowsUnaligned200k = "python -c \"from RunDGEN import *; run_saved_dgen('"+saveStatHere+"', ['/Users/leo/rice/res/data/cichlid/alignment/cichlidLarge.phylip-sequential.txt'], window_size="+windowSize+", window_offset="+windowOffset+", verbose=True, plot='"+plotNameWindowsUnal200k+"', meta='"+"Dgen"+"')\" > "+resultsNameWindowsUnal200k;
        String dgenCommandWindowsAligned200k = "python -c \"from RunDGEN import *; run_saved_dgen('"+saveStatHere+"', ['/Users/leo/rice/res/data/cichlid/alignment/cichlidLargeAligned.phylip-sequential.txt'], window_size="+windowSize+", window_offset="+windowOffset+", verbose=True, plot='"+plotNameWindowsAlig200k+"', meta='"+"Dgen"+"')\" > "+resultsNameWindowsAlig200k;


        //System.out.println(dgenCommand);

        try {
            String cmdPath = "/Users/leo/rice/res/alpha/alpha/CommandLineFiles";

            //started gen, results, and 2x 100kbp results at 2:24, ended at

            //generate the stat
            GenerateDgenStatistic(treeString,networkString,outGroupName,saveStatHere,numberOfRandomRuns);

            //run stat on whole alignment
            dResults = GetOutputFromProgramNoHangFromDir(dgenCommand,cmdPath);

            //run stat on whole alignment sent by them
            //dResultsWindows = GetOutputFromProgramNoHangFromDir(dgenCommandWindows,cmdPath);

            //run on unal and alig first 100kbp
            dResults100kbpUnaligned = GetOutputFromProgramNoHangFromDir(dgenCommandWindows100kbpUnaligned,cmdPath);
            dResults100kbpAligned = GetOutputFromProgramNoHangFromDir(dgenCommandWindows100kbpAligned,cmdPath);

            //run on
            dResultsUnaligned200k = GetOutputFromProgramNoHangFromDir(dgenCommandWindowsUnaligned200k,cmdPath);
            dResultsAligned200k = GetOutputFromProgramNoHangFromDir(dgenCommandWindowsAligned200k,cmdPath);

            int stop = 0;




        } catch (IOException e) {
            e.printStackTrace();
        }



    }



    private static void RunExperimentFig7() {

            //delete everything in the folder storing the results
            //rm -rf /path/to/directory/*
            try {
                //String cmdPath = "/Users/leo/rice/res/alpha/alpha/CommandLineFiles";
                GetOutputFromProgramNoHangFromDir("rm -rf /Users/leo/rice/res/data/dgen/tmp/fig7/*",null);

                GetOutputFromProgramNoHangFromDir("rm /Users/leo/rice/res/data/dgen/tmp/fig7/fig_fig7.txt",null);
                GetOutputFromProgramNoHangFromDir("echo \"D or p Value, Statistic, Significant\" > /Users/leo/rice/res/data/dgen/tmp/fig7/fig_fig7.txt",null);

            } catch (IOException e) {
                e.printStackTrace();
            }

            String treeString = "((((P1,P2),(A1,A2)),P3),O);";
            double introgressionProb = 0.9;
            String introgressionList = "[('P3', 'P2'),('O', 'P3')]";

            String outGroupName = "O";
            String saveStatHere = "/Users/leo/rice/res/data/dgen/tmp/fig7/stat_6taxTwoRetic_fig7.txt";
            int numberOfRandomRuns = 100;

            //gonna use the python helper for now (which is what ill tell users to use if they want help with simple networks
            String networkString = Call_Create_Network_Helper(treeString,introgressionList,introgressionProb);

            //call the java code to actually create the statistic
            GenerateDgenStatistic(treeString,networkString,outGroupName,saveStatHere,numberOfRandomRuns);

            //and then finally call python code to run the statistic on some actual data

            //migRate modulation, the first part
                String simPath;

                //simPath = "/Users/leo/rice/res/data/dgen/simulations/4taxa/withGeneFlow/multipleReticulations/3to2_4to3_6tax/";
                simPath = "/Users/leo/rice/res/data/dgen/simulations/4taxa/withGeneFlow/multipleReticulations/emStartAndEnd_3to2_4to3_6tax/";

        //create list of files to pass into dgen run
                StringBuilder filesList = new StringBuilder();
                for(int simNum = 1; simNum<11; simNum++){
                    if(filesList.length()<1)
                        filesList.append("'"+simPath+"sim"+simNum+"/seqfileNamed'");
                    else
                        filesList.append(",'"+simPath+"sim"+simNum+"/seqfileNamed'");
                }

                String plotName = "/Users/leo/rice/res/data/dgen/tmp/fig7/plot_dgen_6taxTwoRetic";
                String resultsName = "/Users/leo/rice/res/data/dgen/tmp/fig7/results_dgen_6taxTwoRetic.txt";

                String dgenCommand = "["+filesList+"]";

                //newer test, test running newest version at cmd line
                // WORKS!
                dgenCommand = "python -c \"from RunDGEN import *; run_saved_dgen('"+saveStatHere+"', ["+filesList.toString()+"], plot='"+plotName+"', meta='"+"Dgen"+"')\" > "+resultsName;

                String dResults;

                //System.out.println(dgenCommand);

                try {
                    String cmdPath = "/Users/leo/rice/res/alpha/alpha/CommandLineFiles";
                    dResults = GetOutputFromProgramNoHangFromDir(dgenCommand,cmdPath);

                    //make fig file from two result types
                    GetOutputFromProgramNoHangFromDir("cat "+plotName+"_0.txt >> /Users/leo/rice/res/data/dgen/tmp/fig7/fig_fig7.txt",null);
                    //previously run d stat results
                    GetOutputFromProgramNoHangFromDir("cat /Users/leo/rice/res/ALPHA/wabiFinal-ALPHA-master/CommandLineFiles/plot_6TaxaS5_0.txt >> /Users/leo/rice/res/data/dgen/tmp/fig7/fig_fig7.txt",null);
                    //swap wording to Dstat in the old results file
                    GetOutputFromProgramNoHangFromDir("sed -i -- 's/6TaxaS5/Dstat/g' fig_fig7.txt","/Users/leo/rice/res/data/dgen/tmp/fig7");



                } catch (IOException e) {
                    e.printStackTrace();
                }



    }



    private static void MakeAbbaBaba() {

        String treeString = "(((P1,P2),P3),ZZ);";
        double introgressionProb = 0.9;
        String introgressionList = "[('P2', 'P3')]";

        String outGroupName = "ZZ";
        String saveStatHere = "/Users/leo/rice/res/data/dgen/tmp/AbbaBaba/stat_AbbaBabaZZ.txt";
        int numberOfRandomRuns = 100;

        //gonna use the python helper for now (which is what ill tell users to use if they want help with simple networks
        String networkString = Call_Create_Network_Helper(treeString,introgressionList,introgressionProb);

        //call the java code to actually create the statistic
        GenerateDgenStatistic(treeString,networkString,outGroupName,saveStatHere,numberOfRandomRuns);

        //lets make another with O as outgroup just to have around
         treeString = "(((P1,P2),P3),O);";
         introgressionProb = 0.9;
         introgressionList = "[('P2', 'P3')]";

         outGroupName = "O";
         saveStatHere = "/Users/leo/rice/res/data/dgen/tmp/AbbaBaba/stat_AbbaBabaO.txt";
        //int numberOfRandomRuns = 100;

        //gonna use the python helper for now (which is what ill tell users to use if they want help with simple networks
         networkString = Call_Create_Network_Helper(treeString,introgressionList,introgressionProb);

        //call the java code to actually create the statistic
        GenerateDgenStatistic(treeString,networkString,outGroupName,saveStatHere,numberOfRandomRuns);


    }

    private static void RunExperiment1Fig3() {

        try {
            //String cmdPath = "/Users/leo/rice/res/alpha/alpha/CommandLineFiles";
            GetOutputFromProgramNoHangFromDir("rm -rf /Users/leo/rice/res/data/dgen/tmp/fig3/*",null);

        //create initial file for the final storage of info for R plotting
            GetOutputFromProgramNoHangFromDir("rm /Users/leo/rice/res/data/dgen/tmp/fig3/fig_fig3.txt",null);
            GetOutputFromProgramNoHangFromDir("echo \"p Value, Migration Rate, Significant\" > /Users/leo/rice/res/data/dgen/tmp/fig3/fig_fig3.txt",null);

        } catch (IOException e) {
            e.printStackTrace();
        }

        String treeString = "(((P1,P2),(P3,P4)),O);";
        double introgressionProb = 0.9;
        String introgressionList = "[('P1', 'P3')]";

        String outGroupName = "O";
        String saveStatHere = "/Users/leo/rice/res/data/dgen/tmp/fig3/stat_5taxp1p3_fig3.txt";
        int numberOfRandomRuns = 100;

        //gonna use the python helper for now (which is what ill tell users to use if they want help with simple networks
        String networkString = Call_Create_Network_Helper(treeString,introgressionList,introgressionProb);

        //currently, strange error. lines i used to test and find out there was an extra paranthesis that apparently causes problems
        //networkString = "(((((P1)#H1:0::0.9),P2),((#H1:0::0.1,P3),P4)),ZZ);";
        //String networkStringLessP = "((((P1)#H1:0::0.9,P2),((#H1:0::0.1,P3),P4)),ZZ);";
        //worked before...nearly same
        //double inhP = 0.9;
        //String networkString2 = "(((P1, (P2)#H1:0::" + inhP + "), ((P3, #H1:0::" + (1 - inhP) + "), P4)), ZZ);";
        //String networkString3 = "(((P2, (P1)#H1:0::" + inhP + "), ((P3, #H1:0::" + (1 - inhP) + "), P4)), ZZ);";
        //dendroscope compare
        //(((((P1)#H1:0),P2),((#H1:0,P3),P4)),ZZ);
        //(((P2, (P1)#H1:0), ((P3, #H1:0), P4)), ZZ);

        //call the java code to actually create the statistic
        GenerateDgenStatistic(treeString,networkString,outGroupName,saveStatHere,numberOfRandomRuns);

        //and then finally call python code to run the statistic on some actual data

        for(double migRate : new double[]{0.0, 0.01, 0.05,0.1,0.25,0.5}){
            String simPath;
            if(migRate == 0.0)
                simPath = "/Users/leo/rice/res/data/dgen/simulations/5taxa/noGeneFlow/dFoilStandard50kbp/";
            else
                //simPath = "/Users/leo/rice/res/data/dgen/simulations/5taxa/withGeneFlow/dFoilStandard50kbpM"+migRate+"/";
                //this fixes the possibility of ongoing migration not being totally correct (that the box ends when things speciate vs possibly going a lil bit higher for ancestral migration
                simPath = "/Users/leo/rice/res/data/dgen/simulations/5taxa/withGeneFlow/emStartAndEnd_dFoilStandard50kbpM"+migRate+"/";

            //create list of files to pass into dgen run
            StringBuilder filesList = new StringBuilder();
            for(int simNum = 1; simNum<11; simNum++){
                if(filesList.length()<1)
                    filesList.append("'"+simPath+"sim"+simNum+"/seqfile'");
                else
                    filesList.append(",'"+simPath+"sim"+simNum+"/seqfile'");
            }

            //String statName = "/Users/leo/rice/res/data/dgen/tmp/stat_5taxp1p3_fig3_O.txt";
            String plotName = "/Users/leo/rice/res/data/dgen/tmp/fig3/plot_dgen_dFoilStandard50kbpM"+migRate;
            String resultsName = "/Users/leo/rice/res/data/dgen/tmp/fig3/results_dgen_dFoilStandard50kbpM"+migRate+".txt";

            String dgenCommand = "["+filesList+"]";

            //newer test, test running newest version at cmd line
            //this one works
            //dgenCommand = "python -c \"from RunDGEN import *; run_saved_dgen('/Users/leo/rice/res/data/dgen/tmp/testStatFoil2and3iThinkO.txt', ['/Users/leo/rice/res/data/dgen/simulations/5taxa/noGeneFlow/dFoilStandard50kbp/sim1/seqfile','/Users/leo/rice/res/data/dgen/simulations/5taxa/noGeneFlow/dFoilStandard50kbp/sim2/seqfile'], plot='/Users/leo/rice/res/data/dgen/tmp/plot_dgen_dFoilStandard50kbp', meta='0.0')\" > "+resultsName;
            // WORKS!
            dgenCommand = "python -c \"from RunDGEN import *; run_saved_dgen('"+saveStatHere+"', ["+filesList.toString()+"], plot='"+plotName+"', meta='"+migRate+"')\" > "+resultsName;


            //test running the old version at the command line (works)
            //String statName = "/Users/leo/rice/res/ALPHA/ALPHA-master104/CommandLineFiles/stat_5taxP1toP3.txt";
            //String plotName = "/Users/leo/rice/res/data/dgen/tmp/plot_dFoilStandard50kbp";
            //String resultsName = "/Users/leo/rice/res/data/dgen/tmp/results_dFoilStandard50kbp.txt";
            //dgenCommand = "python -c \"from CalculateGeneralizedDStatistic import *; calculate_generalized(["+filesList.toString()+"], statistic='"+statName+"',plot='"+plotName+"',meta='"+migRate+"')\" > "+resultsName;
            String dResults;

            try {
                String cmdPath = "/Users/leo/rice/res/alpha/alpha/CommandLineFiles";
                dResults = GetOutputFromProgramNoHangFromDir(dgenCommand,cmdPath);

                //add results to fig file for R plotting
                //System.out.println("cat "+plotName+" >> /Users/leo/rice/res/data/dgen/tmp/fig3/fig_fig3.txt");
                GetOutputFromProgramNoHangFromDir("cat "+plotName+"_0.txt >> /Users/leo/rice/res/data/dgen/tmp/fig3/fig_fig3.txt",null);

            } catch (IOException e) {
                e.printStackTrace();
            }


        }

    }

    private static void RunExperiment2Fig4() {

        //delete everything in the folder storing the results
        //rm -rf /path/to/directory/*
        try {
            //String cmdPath = "/Users/leo/rice/res/alpha/alpha/CommandLineFiles";
            GetOutputFromProgramNoHangFromDir("rm -rf /Users/leo/rice/res/data/dgen/tmp/fig4/*",null);

            //create initial files for the final storage of info for R plotting
            GetOutputFromProgramNoHangFromDir("rm /Users/leo/rice/res/data/dgen/tmp/fig4/fig_fig4rate.txt",null);
            GetOutputFromProgramNoHangFromDir("echo \"p Value, Migration Rate, Significant\" > /Users/leo/rice/res/data/dgen/tmp/fig4/fig_fig4rate.txt",null);
            GetOutputFromProgramNoHangFromDir("rm /Users/leo/rice/res/data/dgen/tmp/fig4/fig_fig4time.txt",null);
            GetOutputFromProgramNoHangFromDir("echo \"p Value, Migration Time, Significant\" > /Users/leo/rice/res/data/dgen/tmp/fig4/fig_fig4time.txt",null);

        } catch (IOException e) {
            e.printStackTrace();
        }

        String treeString = "(((P1,P2),(P3,(P4,P5))),O);";
        double introgressionProb = 0.9;
        String introgressionList = "[('P1', 'P3')]";

        String outGroupName = "O";
        String saveStatHere = "/Users/leo/rice/res/data/dgen/tmp/fig4/stat_6taxp1p3_fig4.txt";
        int numberOfRandomRuns = 100;

        //gonna use the python helper for now (which is what ill tell users to use if they want help with simple networks
        String networkString = Call_Create_Network_Helper(treeString,introgressionList,introgressionProb);

        //call the java code to actually create the statistic
        GenerateDgenStatistic(treeString,networkString,outGroupName,saveStatHere,numberOfRandomRuns);

        //and then finally call python code to run the statistic on some actual data

        //migRate modulation, the first part

        for(double migRate : new double[]{0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001}) {
            String simPath;

            DecimalFormat df = new DecimalFormat("#");
            df.setMaximumFractionDigits(12);
            //System.out.println(df.format(migRate));
            String migRateString = "0" + df.format(migRate);

            //added emStartAndEnd to folder name and remade those ones to ensure the -ems end at the right time instead of leaving it up to ms
            simPath = "/Users/leo/rice/res/data/dgen/simulations/6taxa/withGeneFlow/StrengthTrend/emStartAndEnd_M"+migRateString+"/";

            //create list of files to pass into dgen run
            StringBuilder filesList = new StringBuilder();
            for(int simNum = 1; simNum<21; simNum++){
                if(filesList.length()<1)
                    filesList.append("'"+simPath+"sim"+simNum+"/seqfileNamed'");
                else
                    filesList.append(",'"+simPath+"sim"+simNum+"/seqfileNamed'");
            }

            String plotName = "/Users/leo/rice/res/data/dgen/tmp/fig4/plot_dgen_6taxMigTrend"+migRateString;
            String resultsName = "/Users/leo/rice/res/data/dgen/tmp/fig4/results_dgen_6taxMigTrend"+migRateString+".txt";

            String dgenCommand = "["+filesList+"]";

            //newer test, test running newest version at cmd line
             // WORKS!
            dgenCommand = "python -c \"from RunDGEN import *; run_saved_dgen('"+saveStatHere+"', ["+filesList.toString()+"], plot='"+plotName+"', meta='"+migRateString+"')\" > "+resultsName;

            String dResults;

            try {
                String cmdPath = "/Users/leo/rice/res/alpha/alpha/CommandLineFiles";
                dResults = GetOutputFromProgramNoHangFromDir(dgenCommand,cmdPath);

                //add results to fig file for R plotting
                GetOutputFromProgramNoHangFromDir("cat "+plotName+"_0.txt >> /Users/leo/rice/res/data/dgen/tmp/fig4/fig_fig4rate.txt",null);


            } catch (IOException e) {
                e.printStackTrace();
            }


        }

        //migTime modulation, the second part

        for(double migTime : new double[]{0.1, 0.3, 0.5, 0.7, 0.9}) {
            String simPath;

            //DecimalFormat df = new DecimalFormat("#");
            //df.setMaximumFractionDigits(12);
            //String migRateString = "0" + df.format(migRate);

            //added emStartAndEnd to folder name and remade those ones to ensure the -ems end at the right time instead of leaving it up to ms
            simPath = "/Users/leo/rice/res/data/dgen/simulations/6taxa/withGeneFlow/TimeTrend/emStartAndEnd_M0.0000002t"+migTime+"t/";



            //create list of files to pass into dgen run
            StringBuilder filesList = new StringBuilder();
            for(int simNum = 1; simNum<21; simNum++){
                if(filesList.length()<1)
                    filesList.append("'"+simPath+"sim"+simNum+"/seqfileNamed'");
                else
                    filesList.append(",'"+simPath+"sim"+simNum+"/seqfileNamed'");
            }

            String plotName = "/Users/leo/rice/res/data/dgen/tmp/fig4/plot_dgen_6taxTimeTrend"+migTime;
            String resultsName = "/Users/leo/rice/res/data/dgen/tmp/fig4/results_dgen_6taxTimeTrend"+migTime+".txt";

            String dgenCommand = "["+filesList+"]";

            //newer test, test running newest version at cmd line
            // WORKS!
            dgenCommand = "python -c \"from RunDGEN import *; run_saved_dgen('"+saveStatHere+"', ["+filesList.toString()+"], plot='"+plotName+"', meta='"+migTime+"')\" > "+resultsName;

            String dResults;

            try {
                String cmdPath = "/Users/leo/rice/res/alpha/alpha/CommandLineFiles";
                dResults = GetOutputFromProgramNoHangFromDir(dgenCommand,cmdPath);

                //add results to fig file for R plotting
                GetOutputFromProgramNoHangFromDir("cat "+plotName+"_0.txt >> /Users/leo/rice/res/data/dgen/tmp/fig4/fig_fig4time.txt",null);


            } catch (IOException e) {
                e.printStackTrace();
            }


        }

    }


    public static void CichlidMakeBlocks(){

        int smallBlockSize = 10000;
        int mediumBlockSize = 100000;
        int largeBlockSize = 200000;

        try{

            //the original phylip file
            BufferedReader cichlidOriginal = new BufferedReader(new FileReader("/Users/leo/rice/res/data/cichlid/alignment/cichlid6tax.phylip-sequential.txt"));

            //the new blocked ones
            BufferedWriter writeSmall = new BufferedWriter(new FileWriter("/Users/leo/rice/res/data/cichlid/alignment/cichlid6taxSmall.phylip-sequential.txt"));
            BufferedWriter writeMedium = new BufferedWriter(new FileWriter("/Users/leo/rice/res/data/cichlid/alignment/cichlid6taxMedium.phylip-sequential.txt"));
            BufferedWriter writeLarge = new BufferedWriter(new FileWriter("/Users/leo/rice/res/data/cichlid/alignment/cichlid6taxLarge.phylip-sequential.txt"));


            String firstLine = cichlidOriginal.readLine();
            writeSmall.write(firstLine.replace("14528710",String.valueOf(smallBlockSize))+"\n");
            writeMedium.write(firstLine.replace("14528710",String.valueOf(mediumBlockSize))+"\n");
            writeLarge.write(firstLine.replace("14528710",String.valueOf(largeBlockSize))+"\n");

            for(int iCich =1; iCich < 7; iCich++) {
                System.out.println(iCich);
                String thisLine = cichlidOriginal.readLine();

                //1 bp test
                writeSmall.write(thisLine.substring(0,10+smallBlockSize)+"\n");
                //2 bp test
                writeMedium.write(thisLine.substring(0,10+mediumBlockSize)+"\n");
                writeLarge.write(thisLine.substring(0,10+largeBlockSize)+"\n");
            }


            cichlidOriginal.close();
            writeSmall.close();
            writeMedium.close();
            writeLarge.close();


        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    public static String Call_Create_Network_Helper(String treeString, String introgressionList, double introgressionProb) {

        //python helper function creates network string
        String cmdString = "python -c \"from RunDGEN import *; print Create_Network_Helper('"+treeString+"',"+introgressionList+","+introgressionProb+")\"";
        String cmdPath = "/Users/leo/rice/res/alpha/alpha/CommandLineFiles";
        String cmdResults = "";
        try {
            cmdResults = GetOutputFromProgramNoHangFromDir(cmdString, cmdPath);
        } catch (IOException e) {
            e.printStackTrace();
        }

        String networkString = cmdResults.replace("\n","");

        return networkString;

    }

    /**
     *
     * @param topology is the tree or network topology (with inheritance probability), which contains outgroup
     * @param outgroupName
     */
    public static void runTest(Network topology, String outgroupName) {

        List<Tree> trees = getAllTreesWithOutgroup(topology, outgroupName);

        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
//        GeneTreeProbability gtp = new GeneTreeProbability();
        GeneTreeProbabilityIntegrated gtpi = new GeneTreeProbabilityIntegrated();
        double blParam = 1.0;
        double blMax = Double.POSITIVE_INFINITY;
        gtpi.setBranchLengthExponentialPrior(blParam, blMax);

        List<Network> parameterSettings = getAllParametrizedNetworks(topology, 1); // 100 random settings
        boolean firstSetting = true;
        // Invariant groups which are consistent across parameter settings
        List<List<Integer>> consistentInvariantGroups = new ArrayList<>();

        List<Double> probabilities = new ArrayList<>();
        for (Network netWithParameters: parameterSettings) {
            // use this line when gtp is a GeneTreeProbability
//            probabilities = gtp.calculateGTDistribution(netWithParameters, trees, null, false);
            // use these [6] lines when gtp is a GeneTreeProbabilityYF
            double[] probsArray = new double[trees.size()];
            gtp.calculateGTDistribution(netWithParameters, trees, null, probsArray);
            probabilities = new ArrayList<>();
            for (double d: probsArray) {
                probabilities.add(d);
            }

            if (firstSetting) {
                firstSetting = false;
                consistentInvariantGroups = findInvariantGroups(probabilities,equalityCheckType.DGEN);
            } else {
                consistentInvariantGroups = updateInvariantGroups(consistentInvariantGroups, probabilities,equalityCheckType.DGEN);
            }
//            System.out.println(probabilities);
        }

        List<Double> integratedProbs = gtpi.calculateGTDistribution(topology, trees, null, false);

        List<List<Integer>> integratedInvariantGroups = findInvariantGroups(integratedProbs,equalityCheckType.DGEN);

        // Strip internal node names for convenience
        for (Object onode: topology.getTreeNodes()) {
            NetNode node = (NetNode)onode;
            if (!node.isLeaf()) node.setName("");
        }
        for (Object onode: topology.getNetworkNodes()) {
            NetNode node = (NetNode)onode;
            if (!node.isLeaf()) node.setName("");
        }

        String resultsString = "\n////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n";
        resultsString += "\ntesting: " + topology.toString();
        resultsString += "\ntrees:";
        int i = 0;
        for (Tree t: trees) {
            resultsString += " " + i++ + ": " + t;
        }
        resultsString += "\ninvariants from random sampling and standard likelihood:";
        consistentInvariantGroups.sort(sortListByFirstElement);
        resultsString += "\n" + consistentInvariantGroups;
        resultsString += "\ninvariants from integrated likelihood:";
        resultsString += "\n" + integratedInvariantGroups;

        resultsString += "\n\nexample parameter setting:\n";
        resultsString += parameterSettings.get(parameterSettings.size() - 1);

        resultsString += "\nexample classic probs:\n";
        i = 0;
        for (Double d: probabilities) {
            resultsString += i++ + ": " + d + ", ";
        }

        resultsString += "\nintegrated probs:\n";
        i = 0;
        for (Double d: integratedProbs) {
            resultsString += i++ + ": " + d + ", ";
        }

        System.out.println(resultsString);
    }

    /**
     *
     * @param topo Topology which must include outgroupName as a leaf
     * @param outgroupName
     * @return A list of trees, which do not include the outgroup
     */
    private static List<Tree> getAllTreesMinusOutgroup(Network<Double> topo, String outgroupName) {
        // If this method throws an index out of bounds exception, ensure that topo contains outgroupName

        String[] theTaxa = new String[topo.getLeafCount()-1]; //{"A", "B", "C", "D"};
        int i = 0;
        for (NetNode leaf: topo.getLeaves()) {
            if(!leaf.getName().equals(outgroupName)){
                theTaxa[i++] = leaf.getName();
            }
        }
        Arrays.sort(theTaxa);

        List<Tree> oldTrees = Trees.generateAllBinaryTrees(theTaxa);
        List<Tree> trees = new ArrayList<>();

        for (Tree t: oldTrees) {
            // Recreate all trees so that node IDs will start at 0
            // Otherwise something goes wrong in GeneTreeProbability
            trees.add(Trees.readTree(t.toNewick()));
        }

        return trees;
    }

    /**
     *
     * @param topo Topology which must include outgroupName as a leaf
     * @param outgroupName
     * @return A list of trees, which include the outgroup
     */
    private static List<Tree> getAllTreesWithOutgroup(Network<Double> topo, String outgroupName) {
        List<Tree> withoutOutgroup = getAllTreesMinusOutgroup(topo, outgroupName);
        List<Tree> withOutgroup = new ArrayList<>();
        for (Tree t: withoutOutgroup) {
            withOutgroup.add(addOutgroup(t, outgroupName));
        }
        return withOutgroup;
    }

    // Check for fuzzy equality of (double) floating-point values
    private static boolean closeEnough(double x, double y,equalityCheckType whichWayToCheckEquality) {

        if(whichWayToCheckEquality == equalityCheckType.DGEN){
            //10^-14 literal difference
            double epsilon = 0.00000000000001;
            return (Math.abs(x-y)<= epsilon);
        }else if(whichWayToCheckEquality == equalityCheckType.FuzzyLoose){
            //10^-4
            double epsilon = 0.0001;
            double differenceMagnitude = Math.abs(x - y);
            double averageMagnitude = (Math.abs(x) + Math.abs(y)) / 2.0;
            return differenceMagnitude < epsilon * averageMagnitude;

        }else if(whichWayToCheckEquality == equalityCheckType.FuzzyStrict){
            //10^-6
            double epsilon = 0.000001;
            double differenceMagnitude = Math.abs(x - y);
            double averageMagnitude = (Math.abs(x) + Math.abs(y)) / 2.0;
            return differenceMagnitude < epsilon * averageMagnitude;

        }else{
            //default to dgen calc
            //10^-14 literal difference
            double epsilon = 0.00000000000001;
            return (Math.abs(x-y)<= epsilon);
        }

        //will probably change this back to just checking that they are within a small threshold at some point (abs(x-y) < tiny num)
        //done now with option for different ones
    }

    // Find indices of (trees with) equal probabilities
    private static List<List<Integer>> findInvariantGroups(List<Double> probabilities, equalityCheckType whichEqualityType) {
        List<List<Integer>> groups = new ArrayList<>();

        for (int i = 0; i < probabilities.size(); i++) {
            boolean foundAGroup = false;
            for (List<Integer> group: groups) {
                if (closeEnough(probabilities.get(i), probabilities.get(group.get(0)),whichEqualityType)) {
                    group.add(i);
                    foundAGroup = true;
                    break;
                }
            }
            if (! foundAGroup) {
                List<Integer> newGroup = new ArrayList<>();
                newGroup.add(i);
                groups.add(newGroup);
            }
        }

        return groups;
    }

    private static List<List<Integer>> updateInvariantGroups(List<List<Integer>> oldGroups, List<Double> probabilities, equalityCheckType whichEqualityCheck) {
        List<List<Integer>> newGroups = new ArrayList<>();
        for (List<Integer> group: oldGroups) {
            while (group.size() > 0) {
                List<Integer> newGroup = new ArrayList<>();
                newGroup.add(group.get(0));
                group.remove(0);
                for (Integer i: group) {
                    if (closeEnough(probabilities.get(i), probabilities.get(newGroup.get(0)),whichEqualityCheck)) {
                        newGroup.add(i);
                    }
                }
                for (Integer i: newGroup) {
                    if (group.contains(i)) {
                        group.remove(group.indexOf(i));
                    }
                }
                newGroups.add(newGroup);
            }
        }
        return newGroups;
    }



    private static List<Network> getAllParametrizedNetworks(Network topology, int numRandomSettings) {
        List<Network> parametrizedNetworks = new ArrayList<>();

        for (int i = 0; i < numRandomSettings; i++) {
            Network<Double> newNet = randomizeBranchLengths(topology);
            parametrizedNetworks.add(newNet);
        }
        return parametrizedNetworks;
    }

    private static Network randomizeBranchLengths(Network topology) {
        double[] branchLengths = {0.5, 1, 2, 4};
        Network<Double> newNet = topology.clone();

        for (NetNode<Double> node: newNet.getTreeNodes()) {
            for (NetNode parent: node.getParents()) {
                double newLength;
                newLength = branchLengths[(int) (Math.random() * branchLengths.length)];
//                newLength = 4.0;
//                newLength = Math.random() * 4; //10;
                node.setParentDistance(parent, newLength);

            }
        }
        for (NetNode<Double> hybridNode: newNet.getNetworkNodes()) {
            for (NetNode parent: hybridNode.getParents()) {
                double newLength;
                newLength = 0.0;
//                newLength = branchLengths[(int) (Math.random() * branchLengths.length)];
//                newLength = NetNode.NO_DISTANCE; // THIS IS NEGATIVE INFINITY, DO NOT USE
                hybridNode.setParentDistance(parent, newLength);
            }
        }
        return newNet;
    }

    private static List<Network> getAllParametrizedNetworksRandom(Network topology, int numRandomSettings) {
        List<Network> parametrizedNetworks = new ArrayList<>();

        for (int i = 0; i < numRandomSettings; i++) {
            Network<Double> newNet = randomizeBranchLengthsUniformBounded(topology);
            parametrizedNetworks.add(newNet);
        }
        return parametrizedNetworks;
    }

    private static Network randomizeBranchLengthsUniformBounded(Network topology) {
        Network<Double> newNet = topology.clone();

        //current min and max values are hard coded (ten and 10^-6  right now)
        double minValue = 0.000001;
        double maxValue = 10.0;

        for (NetNode<Double> node: newNet.getTreeNodes()) {
            for (NetNode parent: node.getParents()) {
                double newLength = minValue + (maxValue-minValue)*(new Random()).nextDouble();
                node.setParentDistance(parent, newLength);
            }
        }
        for (NetNode<Double> hybridNode: newNet.getNetworkNodes()) {
            for (NetNode parent: hybridNode.getParents()) {
                double newLength;
                newLength = 0.0;
//                newLength = branchLengths[(int) (Math.random() * branchLengths.length)];
//                newLength = NetNode.NO_DISTANCE; // THIS IS NEGATIVE INFINITY, DO NOT USE
                hybridNode.setParentDistance(parent, newLength);
            }
        }
        return newNet;
    }


    private static Tree addOutgroup(Tree tree, String outgroupName) {
        TMutableNode oldRoot = (TMutableNode) tree.getRoot();
        STITree newTree = (STITree) Trees.readTree(";");
        STITree outgroup = (STITree) Trees.readTree(outgroupName);

        newTree.getRoot().adoptChild(oldRoot);
        newTree.getRoot().adoptChild(outgroup.getRoot());

        //TODO adopt child not working correctly in updating tree properties like set of nodes etc, adding this inefficient line to deal with it for now
        newTree = (STITree) Trees.readTree(newTree.toNewickWD());
        return newTree;
    }



    private static void FindGeneTreeInvariantsBracket(Tree treeHypothesis, Network networkHypothesis, String outgroupName){

    }

    private static void FindGeneTreeInvariantsIntegrated(Tree treeHypothesis, Network networkHypothesis, String outgroupName){

    }

    private static void TestDStatisticCase(){

    }


    private static ArrayList<String> GetSitePatternsForTreeOrderedAlphabeticallyByTaxaName(STITree sayMyPatterns, String outgroupName) {

        System.out.println("DOING TREE "+sayMyPatterns);

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

        //if i wanna return list of strings of patterns like i was originally thinking
        ArrayList<String> returnMeStrings = new ArrayList<String>();

        //debug printing the full site patterns list
        System.out.println("For Tree "+sayMyPatterns+" the site patterns are: ");

        //make sorted string array of taxa
        ArrayList<String> sortedTaxa = new ArrayList<String>();
        for(String anyLeaf2 : sayMyPatterns.getLeaves()){
            sortedTaxa.add(anyLeaf2);
        }
        Collections.sort(sortedTaxa);

        //print out patterns for testing
        for(HashMap onePattern : returnMe){
            System.out.print("PTRN: ");
            for(Object oneSiteTaxa : onePattern.entrySet()){
                System.out.print("{"+((Map.Entry)oneSiteTaxa).getKey()+",");
                System.out.print(((Map.Entry)oneSiteTaxa).getValue()+"}");
            }
            System.out.print("\n");

            System.out.print("PTRN(ALPH SORT, CONCAT STRING): ");
            StringBuilder thePatternAlphSorted = new StringBuilder();
            for(String oneSiteTaxa : sortedTaxa){
                System.out.print("{"+oneSiteTaxa+",");
                System.out.print(onePattern.get(oneSiteTaxa)+"}");
                thePatternAlphSorted.append(onePattern.get(oneSiteTaxa));
            }
            System.out.print("\n" +thePatternAlphSorted+"\n");
            returnMeStrings.add(thePatternAlphSorted.toString());

        }
        System.out.print("\n");

        //can return in either format from here now, hashmaps of taxa to A/B site patterns for just list of strings of the patterns (where taxa are alphabetically sorted)
        //the hashmaps are in returnMe
        //the list of strings are in returnMeStrings

        return returnMeStrings;
    }


    private static void TestSitePatternGen() {

        String outgroupName = "O";

        Network<Double> treeForAbba = Networks.readNetwork("(O, (G, (H, C)));");
        List<Tree> trees = getAllTreesWithOutgroup(treeForAbba,"O");
        for(Tree sayMyPatterns : trees){
            ArrayList<String> sitePatterns = GetSitePatternsForTreeOrderedAlphabeticallyByTaxaName((STITree) sayMyPatterns, outgroupName);
        }

        //dfoil fig tree case
        Network<Double> treeForDFOIL = Networks.readNetwork("(((P1,P2),(P3,P4)),O);");
        //dfoil fig tree case
        trees = getAllTreesWithOutgroup(treeForDFOIL,"O");
        for(Tree sayMyPatterns : trees){
            ArrayList<String> sitePatterns = GetSitePatternsForTreeOrderedAlphabeticallyByTaxaName((STITree) sayMyPatterns, outgroupName); GetSitePatternsForTreeOrderedAlphabeticallyByTaxaName((STITree) sayMyPatterns, outgroupName);
        }


    }


    private static void LeoTesting(){

        TestFinalStatisticCreation();

        //TestFigureDFOIL(); //DONE i think, this was me finding out the dfoil paper error

        //TestSitePatternGen(); //tests and shows generation of site patterns for trees


        //tests about how many random pulls are necessary to see the GT invariants
        TestRandomRequirementsClean();

        TestRandomRequirements();

        Network<Double> treeForAbba = Networks.readNetwork("(O, (G, (H, C)));");
        Network treeForAbbaWithoutOG = Networks.readNetwork("(C, (B, A));");
//        Network netForAbba = Networks.readNetwork();
        //abba baba tree case
        runTest(treeForAbba, "O");

        //abba baba network case <- the values do not look right, needs to be looked into why
        Network<Double> netForAbba = Networks.readNetwork("(O:1.0, ((H:1.0,(C:1.0)#H1:0::0.3):1.0,G:1.0,#H1:0::0.7));");
        runTest(netForAbba,"O");

        //dfoil fig tree case
        Network<Double> treeForDFOIL = Networks.readNetwork("(((P1,P2),(P3,P4)),O);");
        //dfoil fig tree case
        runTest(treeForDFOIL, "O");


        //dfoil fig network case
        //add reticultion between P1 and P3 or P1 and P4 (can use whatever way hunter has been using on his own)
        Network<Double> netForDFOIL = Networks.readNetwork("(((P1,P2),(P3,P4)),O);");
        //dfoil fig tree case
        runTest(treeForDFOIL, "O");



        //leo note - not really going past here for now

    }

    private static void TestFinalStatisticCreation() {

        //older test case ill try again in a bit
        //Network<Double> whatsWrongTree = Networks.readNetwork("(ZZ,(((P1,P2),P4),P3));");
        //runTestRandomCompareClean(whatsWrongTree, "ZZ", equalityCheckType.DGEN);
        //Network<Double> test123 = Networks.readNetwork("(((((P1)#H1:0::0.9),P2),((#H1:0::0.1,P3),P4)),O);"); //was testing that the network strings generated by python look right. they do
        //int debugStopPoint = 0;








        /*
        int countRtest = netForAbba.getReticulationCount();
        Iterable<NetNode<Double>> testMe = netForAbba.getNetworkNodes();

        for(NetNode<Double> aNetNode : testMe){
            double someNewProbability = 0.7;
            aNetNode.setParentProbability(aNetNode.getParents().iterator().next(),someNewProbability);
            aNetNode.setParentProbability(aNetNode.getParents().iterator().next(),1-someNewProbability);
            int didThisWork = 0;

        }
        */

        //can set parent probabilities inside testMe

        //simplest case first

        //checking into the simple cases for site patterns invariants
        String treeString = "(ZZ, (G, (H, C)));";
        //abba baba network
        double inheritanceProbability = 0.1;
        String networkString = "(ZZ, ((G, #H1:0::" + inheritanceProbability + "), (H, (C)#H1:0::" + (1 - inheritanceProbability) + ")));";
        String outGroupName = "ZZ";
        String saveStatHere = "/Users/leo/rice/res/data/dgen/tmp/testStat.txt";
        int numberOfRandomRuns = 100;
        GenerateDgenStatistic(treeString,networkString,outGroupName,saveStatHere,numberOfRandomRuns);



        //dfoil fig tree case
        treeString = "(((P1, P2), (P3, P4)), ZZ);";
        //dfoil fig network case
        //add reticultion between P1 and P3 or P1 and P4 (can use whatever way hunter has been using on his own)
        inheritanceProbability = 0.1;
        networkString = "(((P1, (P2)#H1:0::" + inheritanceProbability + "), ((P3, #H1:0::" + (1 - inheritanceProbability) + "), P4)), ZZ);";
        saveStatHere = "/Users/leo/rice/res/data/dgen/tmp/testStatFoil2and3iThink.txt";
        numberOfRandomRuns = 100;
        GenerateDgenStatistic(treeString,networkString,outGroupName,saveStatHere,numberOfRandomRuns);

        //add another taxa to make it 6
        treeString = "((((P1, P2), (P3, P4)), P5), ZZ);";
        //dfoil fig network case
        //add reticultion between P1 and P3 or P1 and P4 (can use whatever way hunter has been using on his own)
        inheritanceProbability = 0.1;
        networkString = "((((P1, (P2)#H1:0::" + inheritanceProbability + "), ((P3, #H1:0::" + (1 - inheritanceProbability) + "), P4)), P5), ZZ);";
        saveStatHere = "/Users/leo/rice/res/data/dgen/tmp/testStat6tax.txt";
        numberOfRandomRuns = 100;
        GenerateDgenStatistic(treeString,networkString,outGroupName,saveStatHere,numberOfRandomRuns);







        //going to higher taxa random trees:
        for(int numTaxMinusOut = 4; numTaxMinusOut < 7; numTaxMinusOut++){
            //ArrayList<String> randomTaxaList = new ArrayList<>();
            String[] randomTaxaList = new String[numTaxMinusOut];
            for(int taxNumLabel = 1; taxNumLabel<=numTaxMinusOut; taxNumLabel++)
                randomTaxaList[taxNumLabel-1] = "P"+taxNumLabel;

            //make a random tree
            Tree theRandomTree = Trees.generateRandomTree(randomTaxaList);
            theRandomTree = addOutgroup(theRandomTree,"ZZ");
            Network<NetNodeInfo> tryingRandomTreesAndNetworks = Networks.readNetwork(theRandomTree.toNewick()); //have to swap to network - swap to cast if theres a cast

            System.out.println("\n\nInvestigating the"+(numTaxMinusOut+1)+" tax random tree "+tryingRandomTreesAndNetworks.toString());

            //cant store everything in console

            //for 8 taxa ended up getting:
            //SameAsDGENEquality|GT|0.0_0.0_0.0_0.0_0.01_0.03_0.05_0.11_
            //SameAsDGENEquality|SP|0.07_0.09_0.19_0.47_0.71_0.99__1.0__1.0_
            //takes a pretty long time still for 8 taxa

            runTestRandomCompareClean(tryingRandomTreesAndNetworks, "ZZ", equalityCheckType.DGEN);

            runTestRandomCompareClean(tryingRandomTreesAndNetworks, "ZZ", equalityCheckType.FuzzyLoose);
            runTestRandomCompareClean(tryingRandomTreesAndNetworks, "ZZ", equalityCheckType.FuzzyStrict);

            int debugStop = 0;

            //jiafan recommends some kind of order walk for the tree to set heights (the Random Tree)
            //make random network by adding random retic. have to add times to tree since jiafan used em in his random retic adder
            //for(NetNode<NetNodeInfo> giveMeValue : tryingRandomTreesAndNetworks.getTreeNodes()){
            //    giveMeValue.getData().storeHeight(((BniNetNode<NetNodeInfo>)giveMeValue).);
            //}
            /*
            JiafanNetworkUtils.addRandomReticulation(tryingRandomTreesAndNetworks,0.9);
            System.out.println("Investigating the "+(numTaxMinusOut+1)+" tax random network "+tryingRandomTreesAndNetworks.toString());

            runTestRandomCompareClean(tryingRandomTreesAndNetworks, "ZZ", equalityCheckType.FuzzyLoose);
            runTestRandomCompareClean(tryingRandomTreesAndNetworks, "ZZ", equalityCheckType.FuzzyStrict);
            runTestRandomCompareClean(tryingRandomTreesAndNetworks, "ZZ", equalityCheckType.DGEN);
            */
        }


        //going to higher taxa random one retic networks:


        int stopHereCurrently = 0;

    }

    public static void GenerateDgenStatistic(String treeString, String networkString, String outGroupName, String saveStatHere, int numberOfRandomRuns) {

        //generate invariants for the tree case and network case
        //tree first
        Network<Double> dTree = Networks.readNetwork(treeString);
        SitePatternInvariants treeInvariants = new SitePatternInvariants(dTree,outGroupName,equalityCheckType.DGEN,numberOfRandomRuns);

        //generate invariants for the network case
        Network<Double> dNetwork = Networks.readNetwork(networkString);
        SitePatternInvariants networkInvariants = new SitePatternInvariants(dNetwork,outGroupName,equalityCheckType.DGEN,numberOfRandomRuns);


        //print out a file for each amalgamed results, one big one to maybe make plotting easier, and the same but for averages for each
        try{
            BufferedWriter statWriter = new BufferedWriter(new FileWriter(saveStatHere));


            //write out ordered taxa e.g. Taxa: ['P4', 'P3', 'P2', 'P1', 'O']
            statWriter.write("Taxa: [");
            int dropLastComma = 0;
            for(String aTaxaName : treeInvariants.sortedTaxa){
                statWriter.write('\'' + aTaxaName + '\'');
                if(dropLastComma != treeInvariants.sortedTaxa.size() - 1){
                    statWriter.write(",");
                }
            }
            statWriter.write("]\n");

            //write out species tree and network
            statWriter.write("Species Tree: "+treeString+"\n");
            statWriter.write("Species Network: "+networkString+"\n");
            statWriter.write("Outgroup: "+outGroupName+"\n");

            //do all the necessary work for new dgen stat and write all the info out to a file
            //for every tree invariant, if it has been broken up, we are interested in it
            for(List<Integer> oneTreeInvariant : treeInvariants.sitePatternInvariants){

                //first, check to see if all the SPs in this invariant are also all in the same SP invariant in the network case (they could be in a new invariant with more potentially? so have to check it this way i think. can just grab the SPs outta their new invariants and then check equality)
                List<List<Integer>> correspondingNetworkInvariants = networkInvariants.GetRefinedInvariants(oneTreeInvariant);


                //not a great way of doing but just gonna check string comparison to check for equality(they should all be sorted to it should be fine)
                //setting to true to debug for now
                //new way for equality check, will actuall check everything (if net inv. is size 1 also because otherwise it for sure is not the same
                //nvm, lazy way of checking string of first element should be fine
                boolean sameInvariant = true;
                if(!oneTreeInvariant.toString().equals(correspondingNetworkInvariants.get(0).toString()))
                    sameInvariant = false;
                //if(correspondingNetworkInvariants.size() == 1){
                //    for(a : correspondingNetworkInvariants.get(0)){
                //
                //    }
                //}


                if(!sameInvariant){
                    //System.out.println(oneTreeInvariant.toString());
                    //System.out.println(correspondingNetworkInvariants.get(0).toString());
                //if(!oneTreeInvariant.toString().equals(correspondingNetworkInvariants.toString())){

                    //invariants on individual lines
                    //statWriter.write("Tree Invariants: "+treeInvariants.toStringThisInvariantShowPatterns(oneTreeInvariant)+"\n");
                    //for(List<Integer> oneNetInvariant : correspondingNetworkInvariants) {
                    //    statWriter.write("CorrespondingNetInvariant Group: "+networkInvariants.toStringThisInvariantShowPatterns(oneNetInvariant)+"\n");
                    //}

                    //tree inv. followed by net invariants on one line
                    statWriter.write("Tree Invariant Followed By Network Invariants: ["+treeInvariants.toStringThisInvariantShowPatterns(oneTreeInvariant));
                    for(List<Integer> oneNetInvariant : correspondingNetworkInvariants) {
                        statWriter.write(","+networkInvariants.toStringThisInvariantShowPatterns(oneNetInvariant));
                    }
                    statWriter.write("]\n");

                }


            }







            //text for final all file
            ArrayList<StringBuilder> writeAllText = new ArrayList<StringBuilder>();
            StringBuilder averagesText = new StringBuilder();


            //close writers
            statWriter.close();
        }catch(Exception e){
            e.printStackTrace();
        }


    }


    private static void TestRandomRequirements(){

        //simplest case first

        //checking into the simple cases for site patterns invariants

        Network<Double> treeForAbba = Networks.readNetwork("(O, (G, (H, C)));");
        //abba baba tree case
        runTestRandomCompare(treeForAbba, "O", equalityCheckType.DGEN);

        //do it for a few cases and then some random cases increasing in # of taxa

        //abba baba network
        double inheritanceProbability = 0.1;
        Network<Double> netForAbba = Networks.readNetwork("(O, ((G, #H1:0::" + inheritanceProbability + "), (H, (C)#H1:0::" + (1 - inheritanceProbability) + ")));");
        runTestRandomCompare(netForAbba,"O", equalityCheckType.DGEN);

        //dfoil fig tree case
        Network<Double> treeForDFOIL = Networks.readNetwork("(((P1, P2), (P3, P4)), O);");
        runTestRandomCompare(treeForDFOIL, "O", equalityCheckType.DGEN);


        //dfoil fig network case
        //add reticultion between P1 and P3 or P1 and P4 (can use whatever way hunter has been using on his own)
         inheritanceProbability = 0.1;
        Network<Double> netForDFOIL = Networks.readNetwork("(((P1, (P2)#H1:0::" + inheritanceProbability + "), ((P3, #H1:0::" + (1 - inheritanceProbability) + "), P4)), O);");
        //runTestRandomCompare(netForDFOIL, "O", equalityCheckType.FuzzyLoose);
        //runTestRandomCompare(netForDFOIL, "O", equalityCheckType.FuzzyStrict);
        runTestRandomCompare(netForDFOIL, "O", equalityCheckType.DGEN);

        int debstoppoint = 0;






        //going to higher taxa random trees:
        for(int numTaxMinusOut = 4; numTaxMinusOut < 8; numTaxMinusOut++){
            //ArrayList<String> randomTaxaList = new ArrayList<>();
            String[] randomTaxaList = new String[numTaxMinusOut];
            for(int taxNumLabel = 1; taxNumLabel<=numTaxMinusOut; taxNumLabel++)
                randomTaxaList[taxNumLabel-1] = "P"+taxNumLabel;

            //make a random tree
            Tree theRandomTree = Trees.generateRandomTree(randomTaxaList);
            theRandomTree = addOutgroup(theRandomTree,"O");
            Network<NetNodeInfo> tryingRandomTreesAndNetworks = Networks.readNetwork(theRandomTree.toNewick()); //have to swap to network - swap to cast if theres a cast

            System.out.println("\n\nInvestigating the"+(numTaxMinusOut+1)+" tax random tree "+tryingRandomTreesAndNetworks.toString());

            runTestRandomCompare(tryingRandomTreesAndNetworks, "O", equalityCheckType.FuzzyLoose);
            runTestRandomCompare(tryingRandomTreesAndNetworks, "O", equalityCheckType.FuzzyStrict);
            runTestRandomCompare(tryingRandomTreesAndNetworks, "O", equalityCheckType.DGEN);

            //jiafan recommends some kind of order walk for the tree to set heights (the Random Tree)
            //make random network by adding random retic. have to add times to tree since jiafan used em in his random retic adder
            //for(NetNode<NetNodeInfo> giveMeValue : tryingRandomTreesAndNetworks.getTreeNodes()){
            //    giveMeValue.getData().storeHeight(((BniNetNode<NetNodeInfo>)giveMeValue).);
            //}
            /*
            JiafanNetworkUtils.addRandomReticulation(tryingRandomTreesAndNetworks,0.9);
            System.out.println("Investigating the "+(numTaxMinusOut+1)+" tax random network "+tryingRandomTreesAndNetworks.toString());

            runTestRandomCompare(tryingRandomTreesAndNetworks, "O", equalityCheckType.FuzzyLoose);
            runTestRandomCompare(tryingRandomTreesAndNetworks, "O", equalityCheckType.FuzzyStrict);
            runTestRandomCompare(tryingRandomTreesAndNetworks, "O", equalityCheckType.DGEN);
            */
        }


        //going to higher taxa random one retic networks:


        int stopHereCurrently = 0;


    }

    private static void TestRandomRequirementsClean(){

        //slightly more complex case first now for testing. the tree (ZZ,(((P1,P2),P4),P3)); seems to be a good one that
        //is not always getting done right on the first try
        //these first lines give a good glimpse into what is going on w gt invariants not all being found in a nice
        //isolated setting. just go to the update invariants function in InvariantsHelperFunctions.java and show the debug messages
        //when u update an invariant group and look at the results of everything
        Network<Double> whatsWrongTree = Networks.readNetwork("(ZZ,(((P1,P2),P4),P3));");
        runTestRandomCompareClean(whatsWrongTree, "ZZ", equalityCheckType.DGEN);
        int debugStopPoint = 0;







        //simplest case first

        //checking into the simple cases for site patterns invariants

        Network<Double> treeForAbba = Networks.readNetwork("(ZZ, (G, (H, C)));");
        //abba baba tree case
        runTestRandomCompareClean(treeForAbba, "ZZ", equalityCheckType.DGEN);

        //do it for a few cases and then some random cases increasing in # of taxa

        //abba baba network
        double inheritanceProbability = 0.1;
        Network<Double> netForAbba = Networks.readNetwork("(ZZ, ((G, #H1:0::" + inheritanceProbability + "), (H, (C)#H1:0::" + (1 - inheritanceProbability) + ")));");
        runTestRandomCompareClean(netForAbba,"ZZ", equalityCheckType.DGEN);

        //dfoil fig tree case
        Network<Double> treeForDFOIL = Networks.readNetwork("(((P1, P2), (P3, P4)), ZZ);");
        System.out.println("\n\nInvestigating "+treeForDFOIL.toString());

        runTestRandomCompareClean(treeForDFOIL, "ZZ", equalityCheckType.DGEN);


        //dfoil fig network case
        //add reticultion between P1 and P3 or P1 and P4 (can use whatever way hunter has been using on his own)
        inheritanceProbability = 0.1;
        Network<Double> netForDFOIL = Networks.readNetwork("(((P1, (P2)#H1:0::" + inheritanceProbability + "), ((P3, #H1:0::" + (1 - inheritanceProbability) + "), P4)), ZZ);");
        System.out.println("\n\nInvestigating "+netForDFOIL.toString());
        //runTestRandomCompare(netForDFOIL, "O", equalityCheckType.FuzzyLoose);
        //runTestRandomCompare(netForDFOIL, "O", equalityCheckType.FuzzyStrict);
        runTestRandomCompareClean(netForDFOIL, "ZZ", equalityCheckType.DGEN);

        int debstoppoint = 0;






        //going to higher taxa random trees:
        for(int numTaxMinusOut = 4; numTaxMinusOut < 7; numTaxMinusOut++){
            //ArrayList<String> randomTaxaList = new ArrayList<>();
            String[] randomTaxaList = new String[numTaxMinusOut];
            for(int taxNumLabel = 1; taxNumLabel<=numTaxMinusOut; taxNumLabel++)
                randomTaxaList[taxNumLabel-1] = "P"+taxNumLabel;

            //make a random tree
            Tree theRandomTree = Trees.generateRandomTree(randomTaxaList);
            theRandomTree = addOutgroup(theRandomTree,"ZZ");
            Network<NetNodeInfo> tryingRandomTreesAndNetworks = Networks.readNetwork(theRandomTree.toNewick()); //have to swap to network - swap to cast if theres a cast

            System.out.println("\n\nInvestigating the"+(numTaxMinusOut+1)+" tax random tree "+tryingRandomTreesAndNetworks.toString());

            //cant store everything in console

            //for 8 taxa ended up getting:
            //SameAsDGENEquality|GT|0.0_0.0_0.0_0.0_0.01_0.03_0.05_0.11_
            //SameAsDGENEquality|SP|0.07_0.09_0.19_0.47_0.71_0.99__1.0__1.0_
            //takes a pretty long time still for 8 taxa

            runTestRandomCompareClean(tryingRandomTreesAndNetworks, "ZZ", equalityCheckType.DGEN);

            runTestRandomCompareClean(tryingRandomTreesAndNetworks, "ZZ", equalityCheckType.FuzzyLoose);
            runTestRandomCompareClean(tryingRandomTreesAndNetworks, "ZZ", equalityCheckType.FuzzyStrict);

            int debugStop = 0;

            //jiafan recommends some kind of order walk for the tree to set heights (the Random Tree)
            //make random network by adding random retic. have to add times to tree since jiafan used em in his random retic adder
            //for(NetNode<NetNodeInfo> giveMeValue : tryingRandomTreesAndNetworks.getTreeNodes()){
            //    giveMeValue.getData().storeHeight(((BniNetNode<NetNodeInfo>)giveMeValue).);
            //}
            /*
            JiafanNetworkUtils.addRandomReticulation(tryingRandomTreesAndNetworks,0.9);
            System.out.println("Investigating the "+(numTaxMinusOut+1)+" tax random network "+tryingRandomTreesAndNetworks.toString());

            runTestRandomCompareClean(tryingRandomTreesAndNetworks, "ZZ", equalityCheckType.FuzzyLoose);
            runTestRandomCompareClean(tryingRandomTreesAndNetworks, "ZZ", equalityCheckType.FuzzyStrict);
            runTestRandomCompareClean(tryingRandomTreesAndNetworks, "ZZ", equalityCheckType.DGEN);
            */
        }


        //going to higher taxa random one retic networks:


        int stopHereCurrently = 0;


    }

    private static void TestFigureDFOIL(){

        //at this point, this has basically been "solved"
        //there are two trees in the big tree grouping (the two symmetric-y trees)
        //whos probability is 2x the prob of all the others because once everything coalesces above
        //the root, these two trees have two orderings both of which give their tree whereas
        //all the others have one ordering so they are in their own group which dfoil messed up on
        //we may or may not point all this stuff out in writing and whereas perhaps there are more mistakes
        //to be found when going to the network case, for now will not look further into it

        //dfoil fig tree case
        Network<Double> treeForDFOIL = Networks.readNetwork("(((P1, P2), (P3, P4)), ZZ);");
        System.out.println("\n\nInvestigating "+treeForDFOIL.toString());

        GeneTreeInvariants gtiTree = new GeneTreeInvariants(treeForDFOIL,"ZZ",equalityCheckType.DGEN,100);
        gtiTree.PrintInvariants();
        gtiTree.GenerateExampleProbabilities();

        //dfoil fig network case
        //add reticultion between P1 and P3 or P1 and P4 (can use whatever way hunter has been using on his own)
        double inheritanceProbability = 0.1;
        //Network<Double> netForDFOIL = Networks.readNetwork("(((P1, (P2)#H1:0::" + inheritanceProbability + "), ((P3, #H1:0::" + (1 - inheritanceProbability) + "), P4)), ZZ);");
        //(((P1, (P2)#H1:0::0.1), ((P3, #H1:0::0.9), P4)), ZZ);

        //previous one was between p2 and p3 (incorrect for whats in the paper fig), this is for p1 vs p3
        Network<Double> netForDFOIL = Networks.readNetwork("((((P1)#H1:0::" + inheritanceProbability + ", P2), ((P3, #H1:0::" + (1 - inheritanceProbability) + "), P4)), ZZ);");
        //dendroscope string: ((((P1)#H1:0, P2), ((P3, #H1:0), P4)), ZZ);



        System.out.println("\n\nInvestigating "+netForDFOIL.toString());

        GeneTreeInvariants gtiNet = new GeneTreeInvariants(netForDFOIL,"ZZ",equalityCheckType.DGEN,100);
        gtiNet.PrintInvariants();
        gtiNet.GenerateExampleProbabilities();

        int debstoppoint = 0;

    }



    public static void runTestRandomCompare(Network topology, String outgroupName, equalityCheckType whichEqualityCheck) {

        List<Tree> trees = getAllTreesWithOutgroup(topology, outgroupName);

        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
        //GeneTreeProbability gtp = new GeneTreeProbability();


        List<Network> parameterSettings = getAllParametrizedNetworksRandom(topology, 100); // 100 random settings
        boolean firstSetting = true;
        // Invariant groups which are consistent across parameter settings
        List<List<Integer>> consistentInvariantGroups = new ArrayList<>();

        //invariants of site patterns and site pattern objects needed
        List<List<Integer>> sitePatternInvariants = new ArrayList<>();
        HashMap<String,Double> sitePatternProbabilities = new HashMap<>();
        //double[] sitePatternProbsArrayAlphabetical = new double[trees.size()];
        List<Double> sitePatternProbsAlphabetical = new ArrayList<>();


        List<Double> probabilities = new ArrayList<>();

        //open up parameter settings and make sure it looks right - also prolly will need to change to exponential decaying random instead of uniform
        for (Network netWithParameters: parameterSettings) {
            // use this line when gtp is a GeneTreeProbability
//            probabilities = gtp.calculateGTDistribution(netWithParameters, trees, null, false);
            // use these [6] lines when gtp is a GeneTreeProbabilityYF
            double[] probsArray = new double[trees.size()];
            gtp.calculateGTDistribution(netWithParameters, trees, null, probsArray);
            probabilities = new ArrayList<>();
            for (double d: probsArray) {
                probabilities.add(d);
            }

            //loop through trees, get all site patterns for a tree, and add probs to hashmap total probabilities
            for(int iIterateTrees = 0; iIterateTrees < trees.size(); iIterateTrees++){

                Tree iTree = trees.get(iIterateTrees);
                double iTreeProbability = probsArray[iIterateTrees];

                for(String iSitePattern: GetSitePatternsForTreeOrderedAlphabeticallyByTaxaName((STITree) iTree,"O")){
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
            ArrayList<String> sortedSitePatterns = new ArrayList(sitePatternProbabilities.keySet());
            Collections.sort(sortedSitePatterns);
            for(String sp : sortedSitePatterns){
                sitePatternProbsAlphabetical.add(sitePatternProbabilities.get(sp));
            }



            if (firstSetting) {
                firstSetting = false;
                consistentInvariantGroups = findInvariantGroups(probabilities,whichEqualityCheck);
                sitePatternInvariants = findInvariantGroups(sitePatternProbsAlphabetical,whichEqualityCheck);
            } else {
                consistentInvariantGroups = updateInvariantGroups(consistentInvariantGroups, probabilities,whichEqualityCheck);
                sitePatternInvariants = updateInvariantGroups(sitePatternInvariants,sitePatternProbsAlphabetical,whichEqualityCheck);
            }
//            System.out.println(probabilities);
        }

        consistentInvariantGroups.sort(sortListByFirstElement);
        //debug printing the 'correct' invariant
        //System.out.println(consistentInvariantGroups);

        //DO A SECOND RUN OF 1x or 2x or whatever random settings instead of a bunch and see how many times the invariants match the 'true' ones\

        int totalNumberOfRuns = 100;

        int ranSetMin = 1;
        int ranSetMax = 10;
        System.out.print("\n");
        if(whichEqualityCheck == equalityCheckType.FuzzyLoose)
            System.out.print("Loose_FuzzEquality|");
        if(whichEqualityCheck == equalityCheckType.FuzzyStrict)
            System.out.print("StrictFuzzEquality|");
        if(whichEqualityCheck == equalityCheckType.DGEN)
            System.out.print("SameAsDGENEquality|");

        System.out.print("NumRandomSettings: "+ranSetMin+" to "+ranSetMax+" Fraction in agreement: ");

        for(int numRandomSettings = ranSetMin; numRandomSettings<ranSetMax; numRandomSettings++) {
            //System.out.println("Doing invariants with " +numRandomSettings+" random BL setting");

            int numCorrect = 0;
            for(int numRuns = 0; numRuns<totalNumberOfRuns; numRuns++){


                parameterSettings = getAllParametrizedNetworksRandom(topology, numRandomSettings); // 100 random settings
                firstSetting = true;
                List<List<Integer>> consistentInvariantGroupsAttempt = new ArrayList<>();
                probabilities = new ArrayList<>();

                //open up parameter settings and make sure it looks right
                for (Network netWithParameters : parameterSettings) {
                    double[] probsArray = new double[trees.size()];
                    gtp.calculateGTDistribution(netWithParameters, trees, null, probsArray);
                    probabilities = new ArrayList<>();
                    for (double d : probsArray) {
                        probabilities.add(d);
                    }

                    if (firstSetting) {
                        firstSetting = false;
                        consistentInvariantGroupsAttempt = findInvariantGroups(probabilities,whichEqualityCheck);
                    } else {
                        consistentInvariantGroupsAttempt = updateInvariantGroups(consistentInvariantGroupsAttempt, probabilities,whichEqualityCheck);
                    }
//            System.out.println(probabilities);
                }

                //invariants found lets check and see if they found the 'right' ones or not
                consistentInvariantGroupsAttempt.sort(sortListByFirstElement);
                if(consistentInvariantGroups.equals(consistentInvariantGroupsAttempt)){
                    numCorrect++;
                }

                //debug the found invariant
                //System.out.println(consistentInvariantGroupsAttempt);

                //debug print trees and probs.
                /*
                System.out.println("\n\nexample parameter setting:\n");
                System.out.println(parameterSettings.get(parameterSettings.size() - 1));

                System.out.println("\nexample classic probs:\n");
                int iDebug = 0;
                for (Double d: probabilities) {
                    System.out.println(iDebug + ": " + d + ", "+trees.get(iDebug));
                    iDebug++;
                }*/

                int debBP = 0;
                debBP++;

            }

            //take care for int division
            double fractionCorrectInvariants = (double)numCorrect/(double)totalNumberOfRuns;
            //System.out.println("NumRandomSettings: "+numRandomSettings+" Fraction of times the invariants found were 'correct': "+fractionCorrectInvariants);
            if (fractionCorrectInvariants == 1.0) {
                System.out.print("_");
            }
            System.out.print(fractionCorrectInvariants+"_");

        }




        //big debug code chunk

        // Strip internal node names for convenience
        for (Object onode: topology.getTreeNodes()) {
            NetNode node = (NetNode)onode;
            if (!node.isLeaf()) node.setName("");
        }
        for (Object onode: topology.getNetworkNodes()) {
            NetNode node = (NetNode)onode;
            if (!node.isLeaf()) node.setName("");
        }

        String resultsString = "\n////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n";
        resultsString += "\ntesting: " + topology.toString();
        resultsString += "\ntrees:";
        int i = 0;
        for (Tree t: trees) {
            resultsString += " " + i++ + ": " + t;
        }
        resultsString += "\ninvariants from random sampling and standard likelihood:";
        consistentInvariantGroups.sort(sortListByFirstElement);
        resultsString += "\n" + consistentInvariantGroups;
        //resultsString += "\ninvariants from integrated likelihood:";
        //resultsString += "\n" + integratedInvariantGroups;

        resultsString += "\n\nexample parameter setting:\n";
        resultsString += parameterSettings.get(parameterSettings.size() - 1);

        resultsString += "\nexample classic probs:\n";
        i = 0;
        for (Double d: probabilities) {
            resultsString += i++ + ": " + d + ", ";
        }

        //resultsString += "\nintegrated probs:\n";
        //i = 0;
        //for (Double d: integratedProbs) {
        //    resultsString += i++ + ": " + d + ", ";
        //}

        resultsString += "\nsite patterns:";
        i = 0;
        ArrayList<String> sortedSitePatterns = new ArrayList(sitePatternProbabilities.keySet());
        Collections.sort(sortedSitePatterns);
        for (String t: sortedSitePatterns) {
            resultsString += " " + i++ + ": " + t;
        }
        resultsString += "\ninvariants from random sampling and standard likelihood:";
        sitePatternInvariants.sort(sortListByFirstElement);
        resultsString += "\n" + consistentInvariantGroups;
        //resultsString += "\ninvariants from integrated likelihood:";
        //resultsString += "\n" + integratedInvariantGroups;

        //resultsString += "\n\nexample parameter setting:\n";
        //resultsString += parameterSettings.get(parameterSettings.size() - 1);

        resultsString += "\nexample classic probs:\n";
        i = 0;
        for (Double d: sitePatternProbsAlphabetical) {
            resultsString += i++ + ": " + d + ", ";
        }

        System.out.println(resultsString);
        //end big debug code chunk


    }


    public static void runTestRandomCompareClean(Network topology, String outgroupName, equalityCheckType whichEqualityCheck) {


        //experiment settings
        int ranSetMin = 1;
        int ranSetMax = 129;
        int totalNumberOfRuns = 100;
        StringBuilder resultStringGT = new StringBuilder();


        //get tree and site pattern 'baseline' aka 100 random settings invariants
        int numberRunsForBaseline = 500;
        GeneTreeInvariants baselineInvariantsGT = new GeneTreeInvariants(topology,outgroupName,whichEqualityCheck,numberRunsForBaseline);
        SitePatternInvariants baselineInvariantsSP = new SitePatternInvariants(topology,outgroupName,whichEqualityCheck,numberRunsForBaseline);



        //debug printing
        /*
        Sayln("Baseline GT and SP Invariants:");
        baselineInvariantsGT.PrintGeneTreesOnly();
        Sayln("");
        Sayln(String.valueOf(baselineInvariantsGT.consistentInvariantGroups));
        baselineInvariantsSP.PrintOrderedTaxaOnly();
        Sayln("");
        baselineInvariantsSP.PrintSitePatternsOnly();
        Sayln("");
        Sayln(String.valueOf(baselineInvariantsSP.sitePatternInvariants));
        */
        System.out.print("\nNumRandomSettings: "+ranSetMin+" to "+ranSetMax+" (each column doubles starting at 1) Fraction in agreement: \n\n");


        if(whichEqualityCheck == equalityCheckType.FuzzyLoose)
            resultStringGT.append("Loose_FuzzEquality|");
        if(whichEqualityCheck == equalityCheckType.FuzzyStrict)
            resultStringGT.append("StrictFuzzEquality|");
        if(whichEqualityCheck == equalityCheckType.DGEN)
            resultStringGT.append("SameAsDGENEquality|");

        StringBuilder resultStringSP = new StringBuilder();
        resultStringSP.append(resultStringGT.toString() + "SP|");
        resultStringGT.append("GT|");


        for(int numRandomSettings = ranSetMin; numRandomSettings<ranSetMax; numRandomSettings = numRandomSettings*2) {
            //System.out.println("Doing invariants with " +numRandomSettings+" random BL setting");

            int numCorrectGT = 0;
            int numCorrectSP = 0;

            //debug printing
            //Sayln("Experiment Gt and Sp Invariants");

            //Sayln("************** DOING RUNS FOR "+numRandomSettings+" RANDOM SETTINGS ***************");

            for(int numRuns = 0; numRuns<totalNumberOfRuns; numRuns++) {


                //GeneTreeInvariants experimentInvariantsGT = baselineInvariantsGT; //temporary for debugging
                GeneTreeInvariants experimentInvariantsGT = new GeneTreeInvariants(topology,outgroupName,whichEqualityCheck,numRandomSettings);
                SitePatternInvariants experimentInvariantsSP = new SitePatternInvariants(topology,outgroupName,whichEqualityCheck,numRandomSettings);


                //invariants found lets check and see if they found the 'right' ones or not
                if (baselineInvariantsGT.consistentInvariantGroups.equals(experimentInvariantsGT.consistentInvariantGroups)) {
                    numCorrectGT++;
                }else{
                    //lets take a look and see what the cause is for it not matching
                    //Sayln("Experiment Gt and Sp Invariants");
                    //System.out.println(experimentInvariantsGT.consistentInvariantGroups + "\n");

                    int checkMeOut = 0;
                }
                if (baselineInvariantsSP.sitePatternInvariants.equals(experimentInvariantsSP.sitePatternInvariants)) {
                    numCorrectSP++;
                }


                //debug printing
                //System.out.println(experimentInvariantsGT.consistentInvariantGroups + "\n" + experimentInvariantsSP.sitePatternInvariants);
                //System.out.println(experimentInvariantsSP.sitePatternInvariants);


            }

            //take care for int division
            double fractionCorrectInvariantsGT = (double)numCorrectGT/(double)totalNumberOfRuns;
            if (fractionCorrectInvariantsGT == 1.0) {
                resultStringGT.append("_");
            }
            resultStringGT.append(fractionCorrectInvariantsGT+"_");

            //and do for SP too
            double fractionCorrectInvariantsSP = (double)numCorrectSP/(double)totalNumberOfRuns;
            if (fractionCorrectInvariantsSP == 1.0) {
                resultStringSP.append("_");
            }
            resultStringSP.append(fractionCorrectInvariantsSP+"_");




        }


        System.out.println(resultStringGT.toString());
        System.out.println(resultStringSP.toString());



    }



    enum equalityCheckType
    {
        DGEN, FuzzyLoose, FuzzyStrict;
    }



/*

public class InvariantsLitmusTest {

    public static void main(String[] args) {

        TestSitePatternGen();

        Network<Double> treeForAbba = Networks.readNetwork("(O, (G, (H, C)));");
        Network treeForAbbaWithoutOG = Networks.readNetwork("(C, (B, A));");
//        Network netForAbba = Networks.readNetwork();
        //abba baba tree case
        runTest(treeForAbba, "O");

        //abba baba network case <- the values do not look right, needs to be looked into why
        Network<Double> netForAbba = Networks.readNetwork("(O:1.0, ((H:1.0,(C:1.0)#H1:0::0.3):1.0,G:1.0,#H1:0::0.7));");
        runTest(netForAbba,"O");

        //dfoil fig tree case
        Network<Double> treeForDFOIL = Networks.readNetwork("(((P1,P2),(P3,P4)),O);");
        //dfoil fig tree case
        runTest(treeForDFOIL, "O");


        //dfoil fig network case
        //add reticultion between P1 and P3 or P1 and P4 (can use whatever way hunter has been using on his own)
        Network<Double> netForDFOIL = Networks.readNetwork("(((P1,P2),(P3,P4)),O);");
        //dfoil fig tree case
        runTest(treeForDFOIL, "O");



        //leo note - not really going past here for now

        Network<Double> netForFullTest = Networks.readNetwork("(A, (B, (C, (D, E))));");
        for (Tree t: getAllTrees(netForFullTest,"A")) {
            Network netForTest = Networks.readNetwork(t.toNewick());
            runTest(netForTest,"A");
        }
        Network<Double> balancedNet = Networks.readNetwork("((A,C)hnode2,(B,D)hnode1);");
//        System.out.println(getAllTrees(balancedNet));
//        runTest(balancedNet);
    }



    public static void runTest(Network topology, String outgroupName) {

        List<Tree> trees = getAllTreesWithOutgroup(topology, outgroupName);


        for(Tree sayMyPatterns : trees){
            ArrayList<String> sitePatterns = GetSitePatternsForTreeOrderedAlphabeticallyByTaxaName((STITree) sayMyPatterns, outgroupName); GetSitePatternsForTreeOrderedAlphabeticallyByTaxaName((STITree) sayMyPatterns, outgroupName);
        }


        topology = addOutgroup(topology, outgroupName);

        GeneTreeProbabilityYF gtp = new GeneTreeProbabilityYF();
//        GeneTreeProbability gtp = new GeneTreeProbability();

        GeneTreeProbabilityIntegrated gtpi = new GeneTreeProbabilityIntegrated();
        double blParam = 1.0;
        double blMax = Double.POSITIVE_INFINITY;
        gtpi.setBranchLengthExponentialPrior(blParam, blMax);

        List<Network> parameterSettings = getAllParametrizedNetworks(topology);
        boolean firstSetting = true;
//        List<Tuple<Integer, Integer>> consistentPairInvariants = new ArrayList<>();
        List<List<Integer>> consistentInvariantGroups = new ArrayList<>();

        List<Double> probabilities = new ArrayList<>();
        for (Network netWithParameters: parameterSettings) {
//            probabilities = gtp.calculateGTDistribution(netWithParameters, trees, null, false);
            double[] probsArray = new double[trees.size()];
            gtp.calculateGTDistribution(netWithParameters, trees, null, probsArray);
            probabilities = new ArrayList<>();
            for (double d: probsArray) {
                probabilities.add(d);
            }

            if (firstSetting) {
                firstSetting = false;
                consistentInvariantGroups = findInvariantGroups(probabilities);
            } else {
                consistentInvariantGroups = updateInvariantGroups(consistentInvariantGroups, probabilities);
            }
//            System.out.println(probabilities);
        }

        List<Double> integratedProbs = gtpi.calculateGTDistribution(topology, trees, null, false);

        List<List<Integer>> integratedInvariantGroups = findInvariantGroups(integratedProbs);

//        Networks.removeAllParameters(topology);
        String netString = topology.toString();
        String resultsString = "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n";
        resultsString += "\ntesting: " + netString;
        resultsString += "\ntrees:";
        int i = 0;
        for (Tree t: trees) {
            resultsString += " " + i++ + ": " + t;
        }
        resultsString += "\ninvariants from random sampling and standard likelihood:";
        resultsString += "\n" + consistentInvariantGroups;
        resultsString += "\ninvariants from integrated likelihood:";
        resultsString += "\n" + integratedInvariantGroups;
//        resultsString += "\n\nclassic probs:\n";
//        int i = 0;
//        for (Double d: probabilities) {
//            resultsString += i++ + ": " + d + ", ";
//        }
        resultsString += "\nintegrated probs:\n";
        i = 0;
        for (Double d: integratedProbs) {
            resultsString += i++ + ": " + d + ", ";
        }
        System.out.println(resultsString);
    }

    private static List<Tree> getAllTrees(Network<Double> topo,String outgroupName) {

        String[] theTaxa = new String[topo.getLeafCount()-1]; //{"A", "B", "C", "D"};
        int i = 0;
        for (NetNode leaf: topo.getLeaves()) {
            if(!leaf.getName().equals(outgroupName)){
                theTaxa[i++] = leaf.getName();
            }
        }

        List<Tree> oldTrees = Trees.generateAllBinaryTrees(theTaxa);
        List<Tree> trees = new ArrayList<>();

        for (Tree t: oldTrees) {
            // Recreate all trees so that node IDs will start at 0
            // Otherwise something goes wrong in the probability calculation
            trees.add(Trees.readTree(t.toNewick()));
        }

        return trees;
    }

    private static List<Tree> getAllTreesWithOutgroup(Network<Double> topo, String outgroupName) {
        List<Tree> withoutOutgroup = getAllTrees(topo,outgroupName);
        List<Tree> withOutgroup = new ArrayList<>();
        for (Tree t: withoutOutgroup) {
            withOutgroup.add(addOutgroup(t, outgroupName));
        }
        return withOutgroup;
    }

    // Check for fuzzy equality of (double) floating-point values
    private static boolean closeEnough(double x, double y) {
        double epsilon = 0.0001;
        double differenceMagnitude = Math.abs(x - y);
        double averageMagnitude = (Math.abs(x) + Math.abs(y)) / 2.0;
        return differenceMagnitude < epsilon * averageMagnitude;
    }

    // Find indices of (trees with) equal probabilities
    private static List<List<Integer>> findInvariantGroups(List<Double> probabilities) {
        List<List<Integer>> groups = new ArrayList<>();

        for (int i = 0; i < probabilities.size(); i++) {
            boolean foundAGroup = false;
            for (List<Integer> group: groups) {
                if (closeEnough(probabilities.get(i), probabilities.get(group.get(0)))) {
                    group.add(i);
                    foundAGroup = true;
                    break;
                }
            }
            if (! foundAGroup) {
                List<Integer> newGroup = new ArrayList<>();
                newGroup.add(i);
                groups.add(newGroup);
            }
        }

        return groups;
    }

    // todo: preserve ordering
    private static List<List<Integer>> updateInvariantGroups(List<List<Integer>> oldGroups, List<Double> probabilities) {
        List<List<Integer>> newGroups = new ArrayList<>();
        for (List<Integer> group: oldGroups) {
            while (group.size() > 0) {
                List<Integer> newGroup = new ArrayList<>();
                newGroup.add(group.get(0));
                group.remove(0);
                for (Integer i: group) {
                    if (closeEnough(probabilities.get(i), probabilities.get(newGroup.get(0)))) {
                        newGroup.add(i);
                    }
                }
                for (Integer i: newGroup) {
                    if (group.contains(i)) {
                        group.remove(group.indexOf(i));
                    }
                }
                newGroups.add(newGroup);
            }
        }
        return newGroups;
    }

    private static List<Network> getAllParametrizedNetworks(Network topology) {
        List<Network> parametrizedNetworks = new ArrayList<>();

        int numRandomSettings = 100;
        for (int i = 0; i < numRandomSettings; i++) {
            Network<Double> newNet = randomizeBranchLengths(topology);
            parametrizedNetworks.add(newNet);
        }
        return parametrizedNetworks;
    }

    private static Network randomizeBranchLengths(Network topology) {
        double[] branchLengths = {0.5, 1, 2, 4};
        Network<Double> newNet = topology.clone();

        for (NetNode<Double> node: newNet.getTreeNodes()) {
            double newLength = branchLengths[(int) (Math.random() * branchLengths.length)];
//            newLength = 4.0;
//            double newLength = Math.random() * 4; //10;
            for (NetNode parent: node.getParents()) {
                node.setParentDistance(parent, newLength);
            }
        }
        for (NetNode<Double> hybridNode: newNet.getNetworkNodes()) {

            for (NetNode parent: hybridNode.getParents()) {
                double newLength = branchLengths[(int) (Math.random() * branchLengths.length)];
//                newLength = 0.0;
//                newLength = NetNode.NO_DISTANCE; // THIS IS NEGATIVE INFINITY, DO NOT USE

                hybridNode.setParentDistance(parent, newLength);
            }
        }
        return newNet;
    }

    // Refactor these two for more DRY?
    // todo: implement
    private static Network addOutgroup(Network net, String outgroupName) {
        Network withOutgroup = net.clone();
        NetNode newRoot;
//        withOutgroup.resetRoot(newRoot);
        return net;
    }

    // Refactor these two for more DRY?
    // todo: implement
    private static Tree addOutgroup(Tree tree, String outgroupName) {
        TNode oldRoot = ((STITree)tree).getRoot();
        //Tree newRoot = Trees.readTree("newRoot;");
        Tree testMe = Trees.readTree(";");
        Tree outgroup = Trees.readTree(outgroupName);

        //((STITree)newRoot).getRoot().adoptChild((TMutableNode) oldRoot);
        //((STITree)newRoot).getRoot().adoptChild((TMutableNode) outgroup.getRoot());
        //((STITree)testMe).getRoot().adoptChild();
        ((STITree)testMe).getRoot().adoptChild((TMutableNode) oldRoot);
        ((STITree)testMe).getRoot().adoptChild((TMutableNode) outgroup.getRoot());



//        try {
//            outgroup = new STITree(outgroupName);
//        } catch(Exception e) {
//            e.printStackTrace();
//        }

//        ((STITree)newRoot).getRoot().adoptChild(outgroup.getRoot());

//                ((STITree)tree).getRoot().adoptChild(); //adopt old tree as child



        testMe = Trees.readTree(testMe.toNewickWD());
        return testMe;
    }



    private static void FindGeneTreeInvariantsBracket(Tree treeHypothesis, Network networkHypothesis, String outgroupName){

    }

    private static void FindGeneTreeInvariantsIntegrated(Tree treeHypothesis, Network networkHypothesis, String outgroupName){

    }

    private static void TestDStatisticCase(){

    }

    private static ArrayList<String> GetSitePatternsForTreeOrderedAlphabeticallyByTaxaName(STITree sayMyPatterns, String outgroupName) {

        System.out.println("DOING TREE "+sayMyPatterns);

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

        //if i wanna return list of strings of patterns like i was originally thinking
        ArrayList<String> returnMeStrings = new ArrayList<String>();

        //debug printing the full site patterns list
        System.out.println("For Tree "+sayMyPatterns+" the site patterns are: ");

        //make sorted string array of taxa
        ArrayList<String> sortedTaxa = new ArrayList<String>();
        for(String anyLeaf2 : sayMyPatterns.getLeaves()){
            sortedTaxa.add(anyLeaf2);
        }
        Collections.sort(sortedTaxa);

        //print out patterns for testing
        for(HashMap onePattern : returnMe){
            System.out.print("PTRN: ");
            for(Object oneSiteTaxa : onePattern.entrySet()){
                System.out.print("{"+((Map.Entry)oneSiteTaxa).getKey()+",");
                System.out.print(((Map.Entry)oneSiteTaxa).getValue()+"}");
            }
            System.out.print("\n");

            System.out.print("PTRN(ALPH SORT, CONCAT STRING): ");
            StringBuilder thePatternAlphSorted = new StringBuilder();
            for(String oneSiteTaxa : sortedTaxa){
                System.out.print("{"+oneSiteTaxa+",");
                System.out.print(onePattern.get(oneSiteTaxa)+"}");
                thePatternAlphSorted.append(onePattern.get(oneSiteTaxa));
            }
            System.out.print("\n" +thePatternAlphSorted+"\n");
            returnMeStrings.add(thePatternAlphSorted.toString());

        }
        System.out.print("\n");

        //can return in either format from here now, hashmaps of taxa to A/B site patterns for just list of strings of the patterns (where taxa are alphabetically sorted)
        //the hashmaps are in returnMe
        //the list of strings are in returnMeStrings

        return null;
    }


    private static void TestSitePatternGen() {

        String outgroupName = "O";

        Network<Double> treeForAbba = Networks.readNetwork("(O, (G, (H, C)));");
        List<Tree> trees = getAllTreesWithOutgroup(treeForAbba,"O");
        for(Tree sayMyPatterns : trees){
            ArrayList<String> sitePatterns = GetSitePatternsForTreeOrderedAlphabeticallyByTaxaName((STITree) sayMyPatterns, outgroupName);
        }

        //dfoil fig tree case
        Network<Double> treeForDFOIL = Networks.readNetwork("(((P1,P2),(P3,P4)),O);");
        //dfoil fig tree case
        trees = getAllTreesWithOutgroup(treeForDFOIL,"O");
        for(Tree sayMyPatterns : trees){
            ArrayList<String> sitePatterns = GetSitePatternsForTreeOrderedAlphabeticallyByTaxaName((STITree) sayMyPatterns, outgroupName); GetSitePatternsForTreeOrderedAlphabeticallyByTaxaName((STITree) sayMyPatterns, outgroupName);
        }


    }

}

*/




}











