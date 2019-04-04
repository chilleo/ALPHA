package edu.rice.cs.bioinfo.programs.phylonet;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.summary.RFDistanceBL;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.apache.commons.io.FileUtils;
import org.jblas.util.Random;
import org.jgrapht.*;
import org.jgrapht.alg.shortestpath.DijkstraShortestPath;
import org.jgrapht.graph.*;

import java.io.*;
import java.lang.reflect.Array;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Created by leo on 6/1/16.
 */
public class LeoScript_FullPaperCwR {
    private static final boolean DEBUG = true;
    private static String globalTimeWriter;

    public static void main(String[] args) {

        DoWholePaper();

        //TestRFCalc();

    }




    private static void DoWholePaper() {

        //DoExperiment2("/Users/leo/rice/tmp", "Experiment2_Test", 1);


        //DoExperiment1("/Users/leo/rice/tmp", "Experiment1_Test", 1);

        //DONE
        //DoExperimentGeneral("/Users/leo/rice/tmp", "dfoilstd", 3, "ms 5 1 -t 14000 -T -r 20000 10000 -I 5 1 1 1 1 1 -ej 1.0 2 1 -ej 1.0 4 3 -ej 1.2 3 1 -ej 1.5 5 1", "seq-gen -z " + Random.nextInt(Integer.MAX_VALUE) + " -mHKY -l 10000 -s .028 -p 10000", false, false, true);

        //start - 15 taxa, 3 per population - done
        //DoExperimentGeneral("/Users/leo/rice/tmp", "dfoilstd5x3", 3, "ms 15 1 -t 14000 -T -r 20000 10000 -I 5 3 3 3 3 3 -ej 1.0 2 1 -ej 1.0 4 3 -ej 1.2 3 1 -ej 1.5 5 1", "seq-gen -z " + Random.nextInt(Integer.MAX_VALUE) + " -mHKY -l 10000 -s .028 -p 10000", false, false, false);

        //start - 100 taxa, 20 per population - done x1
        //DoExperimentGeneral("/Users/leo/rice/tmp", "dfoilstd5x20", 1, "ms 100 1 -t 14000 -T -r 20000 10000 -I 5 20 20 20 20 20 -ej 1.0 2 1 -ej 1.0 4 3 -ej 1.2 3 1 -ej 1.5 5 1", "seq-gen -z " + Random.nextInt(Integer.MAX_VALUE) + " -mHKY -l 10000 -s .028 -p 10000", false, false, false);

        //first experiment with Branch Lengths - not really done yet
        //DoExperimentGeneral("/Users/leo/rice/tmp", "dfoilstdBL", 1, "ms 5 1 -t 14000 -T -r 20000 10000 -I 5 1 1 1 1 1 -ej 1.0 2 1 -ej 1.0 4 3 -ej 1.2 3 1 -ej 1.5 5 1", "seq-gen -z " + Random.nextInt(Integer.MAX_VALUE) + " -mHKY -l 10000 -s .028 -p 10000", false, false, false);

        //first experiment based on bio stuff, primates, with example auto command generation-DONE
        //DoExperimentGeneral("/Users/leo/rice/tmp", "primatesWind10x", 10, MakeSampleStringPrimates(10000, 1.0*Math.pow(10,-8), 1), MakeSeqGenString(10000, 2.5*Math.pow(10,-8), 40000), true, true, true);

        //DOING WINDOW CENTERING FOR PREDICTIONS!! (TEST) first experiment based on bio stuff, primates, with example auto command generation-DONE
        //DoExperimentGeneral("/Users/leo/rice/tmp", "primatesWindCen10x", 10, MakeSampleStringPrimates(10000, 1.0*Math.pow(10,-8), 1), MakeSeqGenString(10000, 2.5*Math.pow(10,-8), 40000), true, true, true);

        //comparing centered and not centered, different window size and window offset numbers
        //GETTING MEMORY ERRORS - CAN SWAP SHORTEST PATH CALC TO SOMETHING MUCH MORE EFFICIENT?
        //WHAT TAKES ALL THE MEM THO? IF ITS THE TREES I CAN PROLLY REMOVE EM FROM MEMORY AS I ADD STUFF TO GRAPH AS STRINGS ETC
        //DoExperimentGeneral("/Users/leo/rice/tmp", "primatesNew5x", 10, MakeSampleStringPrimates(10000, 1.0*Math.pow(10,-8), 1), MakeSeqGenString(10000, 2.5*Math.pow(10,-8), 40000), true, true, true);

        //checking on optimal RF dist based on window BS trees alone assuming perfect ones were picked out
        //1xB is a full finished run. 1xC ima do testing on why i cant get a single 0RF dist tree to happen
        //1xC done, with fixed rf dists (unrooted). took ~35 minutes.
        //doing twenty runs now which should finish in the early AM (limited window size/ offset combos to make it go a lot faster)
        //JK - canceled that to add the smoothwinds calcs back in, will time one run in 1xD and then launch a multi run
        //1xD took ~42 minutes, so will do 15 for now
        //debug run
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "primatesOptimal1xD", 1, MakeSampleStringPrimates(10000, 1.0*Math.pow(10,-8), 1, 40000.0), MakeSeqGenString(10000, 2.5*Math.pow(10,-8), 40000.0), false, false, true);

        //finished a 15xD run with primates
        //doing 15xE (had to fix population size, it was doing branch lengths that were half as long as they were supposed to be (basically was doing branch lengths that were slightly off, fixed now)
        //15xE done. and 15xF newest accurate results
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "primatesOptimal15xF", 15, MakeSampleStringPrimates(10000, 1.0*Math.pow(10,-8), 1, 40000.0), MakeSeqGenString(10000, 2.5*Math.pow(10,-8), 40000.0), true, true, true);

        //human sub population tree
        //15xE done. and 15xF newest accurate results
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "humansOptimal15xF", 15, MakeSampleStringHumans(10000, 1.1*Math.pow(10,-8), 1, 20000.0), MakeSeqGenString(10000, 1.25*Math.pow(10,-8), 20000.0), true, true, true);

        //mosquito tree. mut. rate based on fly and pop size kinda random from random pape for now of 10k. human recomb rate for now
        //15xF done, did with 2x longer branch lengths tho, doing 15xFv2 with fixed bl and also going down lower in offset
        //mosq tree figs2d trial, still comparing lessmem and normal
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "mosqOptimal15xFdubChMemFigS2D", 15, MakeSampleStringMosq(10000, 1.1*Math.pow(10,-8), 1, 10000.0), MakeSeqGenString(10000, 6.2*Math.pow(10,-8), 10000.0), true, true, true);

        //higher tax
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "mosqOptimal15xFdubChMemFigS2D12tax", 15, MakeSampleStringMosq(10000, 1.1*Math.pow(10,-8), 2, 10000.0), MakeSeqGenString(10000, 6.2*Math.pow(10,-8), 10000.0), true, true, true);
        //did both of these
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "mosqOptimal8xFdubChMemFigS2D24tax", 8, MakeSampleStringMosq(10000, 1.1*Math.pow(10,-8), 4, 10000.0), MakeSeqGenString(10000, 6.2*Math.pow(10,-8), 10000.0), true, true, true);
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "mosqOptimal8xFdubChMemFigS2D48tax", 8, MakeSampleStringMosq(10000, 1.1*Math.pow(10,-8), 8, 10000.0), MakeSeqGenString(10000, 6.2*Math.pow(10,-8), 10000.0), true, true, true);

        //feb 2019 onwards

        //just doin a fresh run of stuff before now on feb 5, 2019 (these basically just copy the lines above
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "new_mosqOptimal15xFdubChMemFigS2D12tax", 15, MakeSampleStringMosq(10000, 1.1*Math.pow(10,-8), 2, 10000.0), MakeSeqGenString(10000, 6.2*Math.pow(10,-8), 10000.0), true, true, true);
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "new_mosqOptimal8xFdubChMemFigS2D24tax", 8, MakeSampleStringMosq(10000, 1.1*Math.pow(10,-8), 4, 10000.0), MakeSeqGenString(10000, 6.2*Math.pow(10,-8), 10000.0), true, true, true);
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "new_mosqOptimal8xFdubChMemFigS2D48tax", 8, MakeSampleStringMosq(10000, 1.1*Math.pow(10,-8), 8, 10000.0), MakeSeqGenString(10000, 6.2*Math.pow(10,-8), 10000.0), true, true, true);

        //adding birds in
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "new_birdOptimal10xFdubChMem4tax", 10, MakeSampleStringBirds(10000, Math.pow(10.0,-8.0), 1, 1000000), MakeSeqGenString(10000, 1.7*Math.pow(10,-9), 1000000.0), true, true, true);
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "new_birdOptimal10xFdubChMem8tax", 10, MakeSampleStringBirds(10000, Math.pow(10.0,-8.0), 1, 1000000), MakeSeqGenString(10000, 1.7*Math.pow(10,-9), 1000000.0), true, true, true);
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "new_birdOptimal10xFdubChMem12tax", 10, MakeSampleStringBirds(10000, Math.pow(10.0,-8.0), 1, 1000000), MakeSeqGenString(10000, 1.7*Math.pow(10,-9), 1000000.0), true, true, true);

        //rerunning recent run of primates and human sub pop
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "new_primatesOptimal15xF_11tax", 15, MakeSampleStringPrimates(10000, 1.0*Math.pow(10,-8), 1, 40000.0), MakeSeqGenString(10000, 2.5*Math.pow(10,-8), 40000.0), true, true, true);

        //human sub population tree
       // DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "new_humansOptimal15xF_4tax", 15, MakeSampleStringHumans(10
        // 000, 1.1*Math.pow(10,-8), 1, 20000.0), MakeSeqGenString(10000, 1.25*Math.pow(10,-8), 20000.0), true, true, true);

        //adding butterflies in
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "new_butterflyOptimal10xFdubChMem4tax", 10, MakeSampleStringButterfly(10000, Math.pow(10.0,-8.0), 1, 2000000), MakeSeqGenString(10000, 2.9*Math.pow(10,-9), 2000000.0), true, true, true);
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "new_butterflyOptimal10xFdubChMem8tax", 10, MakeSampleStringButterfly(10000, Math.pow(10.0,-8.0), 2, 2000000), MakeSeqGenString(10000, 2.9*Math.pow(10,-9), 2000000.0), true, true, true);
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "new_butterflyOptimal10xFdubChMem12tax", 10, MakeSampleStringButterfly(10000, Math.pow(10.0,-8.0), 3, 2000000), MakeSeqGenString(10000, 2.9*Math.pow(10,-9), 2000000.0), true, true, true);

        //artificially inflate distance for primates
        //RENAME THIS ONE!!!!
        //DoOptimalityExperimentGeneral("/Users/leo/rice/tmp", "new_primatesOptimal15xF_11tax_10xMuts", 10, MakeSampleStringPrimates(10000, 1.0*Math.pow(10,-8), 1, 40000.0), MakeSeqGenString(10000, 10*2.5*Math.pow(10,-8), 40000.0), true, true, true);

        //script for doing everything in pape before plotting in R
        DoAllFinalPaperResults();
        //run R script after running above to generate plots


        AnalyzeHCGData("/Users/leo/rice/tmp","realHCGData",true,true,true);

        AnalyzeMosqData("/Users/leo/rice/tmp","realMosqData",true,true,true);



    }

    private static void DoAllFinalPaperResults() {


        //setting everything up for the 6 tax mosquito example. and then doing 12 and 18. how long it takes..?
        // done? DONE
        //DoOptimalityExperimentGeneralFinal("/Users/leo/rice/tmp", "final_mosq_6tax_20x", 20, MakeSampleStringMosq(10000, 1.1*Math.pow(10,-8), 1, 10000.0), MakeSeqGenString(10000, 6.2*Math.pow(10,-8), 10000.0), true, true, true);
        //DoOptimalityExperimentGeneralFinal("/Users/leo/rice/tmp", "final_mosq_12tax_20x", 20, MakeSampleStringMosq(10000, 1.1*Math.pow(10,-8), 2, 10000.0), MakeSeqGenString(10000, 6.2*Math.pow(10,-8), 10000.0), true, true, true);
        //DoOptimalityExperimentGeneralFinal("/Users/leo/rice/tmp", "final_mosq_18tax_20x", 20, MakeSampleStringMosq(10000, 1.1*Math.pow(10,-8), 3, 10000.0), MakeSeqGenString(10000, 6.2*Math.pow(10,-8), 10000.0), true, true, true);

        //butterfly
        //DoOptimalityExperimentGeneralFinal("/Users/leo/rice/tmp", "final_butterfly_4tax_20x", 20, MakeSampleStringButterfly(10000, Math.pow(10.0,-8.0), 1, 2000000), MakeSeqGenString(10000, 2.9*Math.pow(10,-9), 2000000.0), true, true, true);
        //DoOptimalityExperimentGeneralFinal("/Users/leo/rice/tmp", "final_butterfly_8tax_20x", 20, MakeSampleStringButterfly(10000, Math.pow(10.0,-8.0), 2, 2000000), MakeSeqGenString(10000, 2.9*Math.pow(10,-9), 2000000.0), true, true, true);
        //DoOptimalityExperimentGeneralFinal("/Users/leo/rice/tmp", "final_butterfly_12tax_20x", 20, MakeSampleStringButterfly(10000, Math.pow(10.0,-8.0), 3, 2000000), MakeSeqGenString(10000, 2.9*Math.pow(10,-9), 2000000.0), true, true, true);

        //bird
        //DoOptimalityExperimentGeneralFinal("/Users/leo/rice/tmp", "final_bird_4tax_20x", 20, MakeSampleStringBirds(10000, Math.pow(10.0,-8.0), 1, 1000000), MakeSeqGenString(10000, 1.7*Math.pow(10,-9), 1000000.0), true, true, true);
        //DoOptimalityExperimentGeneralFinal("/Users/leo/rice/tmp", "final_bird_8tax_20x", 20, MakeSampleStringBirds(10000, Math.pow(10.0,-8.0), 2, 1000000), MakeSeqGenString(10000, 1.7*Math.pow(10,-9), 1000000.0), true, true, true);
        //DoOptimalityExperimentGeneralFinal("/Users/leo/rice/tmp", "final_bird_12tax_20x", 20, MakeSampleStringBirds(10000, Math.pow(10.0,-8.0), 3, 1000000), MakeSeqGenString(10000, 1.7*Math.pow(10,-9), 1000000.0), true, true, true);

        //human subpop
        //DoOptimalityExperimentGeneralFinal("/Users/leo/rice/tmp", "final_humans_4tax_20x", 20, MakeSampleStringHumans(10000, 1.1*Math.pow(10,-8), 1, 20000.0), MakeSeqGenString(10000, 1.25*Math.pow(10,-8), 20000.0), true, true, true);

        //primates 11tax
        //DoOptimalityExperimentGeneralFinal("/Users/leo/rice/tmp", "final_primates_11tax_20x", 20, MakeSampleStringPrimates(10000, 1.0*Math.pow(10,-8), 1, 40000.0), MakeSeqGenString(10000, 2.5*Math.pow(10,-8), 40000.0), true, true, true);


        int debugStopHereCheck = 0;


    }

    private static void DoOptimalityExperimentGeneralFinal(String whereToDoExperiment, String experimentName, int numberOfDuplicateRuns, String msCmd, String seqCmd, boolean runNewSims, boolean runSlowMethods, boolean runFastMethods){

        //params
        int genomeLength = 10000;


        //commands
        //msCmd = "ms 15 1 -t 15 -r 10 4000 -T";
        //seqCmd = "seq-gen -z " + Random.nextInt(Integer.MAX_VALUE) + " -mHKY -l 4000 -s .00544 -p 4000";


        String experimentFolder = whereToDoExperiment+"/"+experimentName;
        if(runNewSims){DeleteFolder(experimentFolder);} //just to be safe, if im doing new sims i for sure gotta rerun everything. someimtes things intelligently dont rerun
        MakeFolder(experimentFolder); //create experiment folder


        for(int i=1; i<=numberOfDuplicateRuns; i++) {

            String baseFolder = experimentFolder+"/Run"+i; //base folder where the exp. happens
            MakeFolder(baseFolder); //create this run's folder
            MakeFolder(baseFolder + "/results"); //create folder for results
            MakeFolder(baseFolder + "/timingsFolder"); //create folder for timings results
            globalTimeWriter = baseFolder+"/timings.txt"; //older way of storing timings, for now i will just leave it

            //variables for wallclock timings of rent, argw, and getting raxml trees + doing SW to plot
            long rentSeconds;
            long argwSeconds;
            long raxmlSeconds;
            long swSeconds;

            //a good reference for one of the last times i was doing this is
            //in the file /Users/leo/rice/res/code/bmc17/start/doAll.sh

            if(runNewSims) {
                //run ms and seq gen
                DPrint("starting ms");
                RunMS(baseFolder, msCmd);
                DPrint("done w ms. starting seq gen");
                RunSeqGen(baseFolder, seqCmd);
                DPrint("done w seq gen");
            }

            if (runFastMethods) {
                //run RENT+, prolly first on infinite sites and then on jc, calculate score(s)

                FormatToRENT(baseFolder);
                Date rentStart = new Date();
                RunRENT(baseFolder);
                rentSeconds = ((new Date()).getTime() - rentStart.getTime())/1000;
                RecordTiming(baseFolder + "/timingsFolder/rent",rentSeconds);
                ProcessRENT(baseFolder, genomeLength);
                compareTwoMSFiles(baseFolder+"/ms/trueTrees",baseFolder+"/rent/finalResultsRENT",baseFolder+"/results/rent",genomeLength);



                //trying less settings for final figs. 3 window sizes and 3 offsets

                for(int iSize = 1000; iSize<=3000; iSize=iSize+1000){
                    for(int jOff = 800; jOff>=200; jOff=jOff/2) {


                        //for now only go down to offset of 100 to make it go faster (not any longer)
                        //in the future would like to check out 50, 25, and maybe even 12 (cast int of 12.5, add cast to division)


                        if(jOff > iSize){continue;}

                        System.out.println("wind " + iSize + " off "+ jOff);
                        TPrint("wind " + iSize + " off "+ jOff);

                        //can do this one time
                        Windower(baseFolder, iSize, jOff);
                        Date raxmlStart = new Date();
                        RunWindowedRAXML(baseFolder, iSize, jOff);
                        raxmlSeconds = ((new Date()).getTime() - raxmlStart.getTime())/1000;
                        RecordTiming(baseFolder + "/timingsFolder/raxml" + iSize + "wind"+ jOff + "off",raxmlSeconds);

                        //this part is 2x one centered one without

                        boolean doCentering = true;
                        String centeringString = "-CNTRD";

                        //optimality check thing (not doing any more for now - prolly a bit too misleading honestly cuz of how everything works and with it allowing switching every time ms switches)
                        //compareMSFileVsOptimal(baseFolder + "/ms/trueTrees", baseFolder, baseFolder + "/results/OPTraxml" +centeringString+ iSize + "bp" + jOff + "os", genomeLength,iSize, jOff, doCentering);

                        ProcessRAXML(baseFolder, iSize, jOff,genomeLength,doCentering);
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/raxml"+centeringString + iSize + "bp" + jOff + "os" + "/treefileRAXML", baseFolder + "/results/raxml"+centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //run our method, calc score(s)
                        //memory errors for now

                        //mem saver so i can actually run low offsets
                        //seems to be working correctly now - but not running it this second, now i am lets see if this works
                        //if i check notes in timing i should now make sure the rf dists for shortest paths now match below since i
                        //moved this code up one line
                        Date swStart = new Date();
                        SmoothWindsLessMem(baseFolder, iSize, jOff, genomeLength, false, doCentering);
                        swSeconds = ((new Date()).getTime() - swStart.getTime())/1000;
                        RecordTiming(baseFolder + "/timingsFolder/sw" + iSize + "wind"+ jOff + "off",swSeconds);
                        RecordTiming(baseFolder + "/timingsFolder/raxml+sw" + iSize + "wind"+ jOff + "off",raxmlSeconds+swSeconds); //combine both raxml and sw
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //above worked but below (smoothwinds function, with "Exception in thread "main" java.lang.OutOfMemoryError: GC overhead limit exceeded") failed when setting was...500window size and 25 offset
                        System.out.println("less mem done ws "+iSize+" off "+jOff);


                        //run our method (WITH BL), calc score(s)
                        //SmoothWinds(baseFolder, iSize, jOff, genomeLength, true, doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsBL" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsBL" + centeringString + iSize + "bp" + jOff + "os", genomeLength);


                        //FIRST HALF WAS WITH CENTERING, THIS HALF IS WITHOUT CENTERING
                        //centering seems pretty much always better, not doing this any more

                        /*
                        doCentering = false;
                        centeringString = "";

                        //this does the optimality comparison (not doing any more for now - prolly a bit too misleading honestly cuz of how everything works and with it allowing switching every time ms switches)
                        compareMSFileVsOptimal(baseFolder + "/ms/trueTrees", baseFolder, baseFolder + "/results/OPTraxml"+centeringString + iSize + "bp" + jOff + "os", genomeLength,iSize, jOff, doCentering);


                        ProcessRAXML(baseFolder, iSize, jOff,genomeLength,doCentering);
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/raxml"+centeringString + iSize + "bp" + jOff + "os" + "/treefileRAXML", baseFolder + "/results/raxml" +centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //run our method, calc score(s)
                        //memory errors for now

                        SmoothWindsLessMem(baseFolder, iSize, jOff, genomeLength, false, doCentering);
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //before, before i made less mem one
                        //SmoothWinds(baseFolder, iSize, jOff,genomeLength, false,doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwinds"+centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwinds"+centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //run our method (WITH BL), calc score(s)
                        //SmoothWinds(baseFolder, iSize, jOff,genomeLength,true,doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsBL"+centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsBL"+centeringString + iSize + "bp" + jOff + "os", genomeLength);
                        */


                    }
                }
            }






            if(runSlowMethods) {
                //run ARGweaver, prolly first on infinite sites and then on jc, calculate score(s)
                //seems to work but is hanging? also might have to clear out directory when i run this

                //running argw takes forever, testing right now so will not keep re running it
                Date argwStart = new Date();
                RunARGW(baseFolder);
                argwSeconds = ((new Date()).getTime() - argwStart.getTime())/1000;
                RecordTiming(baseFolder + "/timingsFolder/argw",argwSeconds);

                ProcessARGW(baseFolder);
                compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/argw/finalLocalTrees", baseFolder + "/results/argw", genomeLength);

            }

            //record timings


            //calc scores


        }

        //amalgamate results from all runs
        AmalgamateResults(experimentFolder, numberOfDuplicateRuns);
        AmalgamateTimingsResults(experimentFolder, numberOfDuplicateRuns);

    }

    private static void RecordTiming(String saveHere, long timeInSeconds) {
        try{
            BufferedWriter writeResults = new BufferedWriter(new FileWriter(saveHere));
            StringBuilder resultsText = new StringBuilder();

            resultsText.append("Seconds:" + timeInSeconds + "\n");

            writeResults.write(resultsText.toString());
            writeResults.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    private static String MakeSampleStringPrimates(int seqLength, double recombRate, int samplesPerPop, double popSize) {

        //for each biologically inspired ms sim, the overall structure is hardcoded based on which bio example it is.
        //here for primates we have an 11 taxon tree from https://www.nature.com/articles/nature12228, with 20y gen time and 40k constant fixed popSize
        double genTime = 20.0;
        //double popSize = 40000.0; //using 40k

        //im actually not using the mutations output from ms for anything, so i dont need the -t theta

        //for recombination paramarameter, as taken from the docs for ms:
        //"To include crossing-over (recombination) in the model,  one uses the-rswitch and specifies the cross-over rate parameter,ρwhich is 4N0r, whereristhe probability of cross-over per generation between the ends of the locus beingsimulated.  In addition, the user specifies the numbers of sites, nsites, betweenwhich recombination can occur.  One should think of nsites as the number ofbase pairs in the locus.  For example, to simulate samples of size 15 for a locuswhich is 2,501 base pairs long with cross-over probability between adjacent basepairs of 10−8per generation, and assumingN0= 106one enters:ms 15 1000 -t 10.004 -r 100.0 2501In this case,ρ= 100, since the cross-over probability between the two ends ofthe locus is 10−8(2501−1) = 2.5×10−5, and thusρ= 4×106×2.5×10−5= 100."
        double paramR = 4 * popSize * seqLength * recombRate; //recomb rate, as mentioned in ms doc, is going to be in "with cross-over probability between adjacent basepairs of 10−8per generation"

        //-ejs are in units of "In each case the first parameter following the switch is the timewhen  the  demographic  change  occurred,  measured  from  the  present  in  unitsof  4N0generations."

        //here for primates i have manually done all the calcs, which end up with the following -ejs, these all used fixed pop size of 40k and gen time of 20:
        /*
        **** SPLITS for ms (in chron. order) ****

        (note- time for -ej is (years/genTime)/(4*40000) or x/160000)

        //had to correct these ones below, were off by a factor
        //fix in the final result because i calculate it fully from pop size (80k/genTime)/(4*popSize) etc

        8 Cross River and 7 Western Lowland Gorilla split 80,000ya

        -ej 0.125 8 7

        87(now 7) and 9 Eastern Lowland Gorilla split 150,000ya

        -ej 0.234375 9 7

        3 Eastern Chimp and 2 Central Chimp split 175,000ya

        -ej 0.2734375 3 2

        5 western chimp and 4 nigeria-cameroon chimp split 289,000ya

        -ej 0.4515625 5 4

        5432 chimps split 419,000ya

        -ej 0.6546875 4 2

        10 bornean orangutan and 11 sumatran split 485,000ya

        -ej 0.7578125 11 10

        5432 chimps and 6 bonobo split 871,000ya

        -ej 1.3609375 6 2

        1 human and 65432 chimp/bonobo split 3,765,000ya

        -ej 5.8828125 2 1

        654321 humans etc and 789 gorillas split 5,619,000ya

        -ej 8.7796875 7 1

        987654321 humans etc and 1011 orangutans split 11,138,000ya

        -ej 17.403125 10 1
        */

        //currently not using inf. sites for anything, but just in case for later ill add the -t in even tho i really dont need it right now
        //using mut. rate of 2.5 * 10 -8
        //String infSitesLine = "-t " + 4*popSize*(2.5*Math.pow(10,-8));
        //String infSitesLine =

        //giving the following command:

        return "ms "+11*samplesPerPop+" 1 -T -r "+paramR+" "+seqLength+" -I 11 "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" -ej "+
                (80000.0/genTime)/(4.0*popSize)+" 8 7 -ej "+
                (150000.0/genTime)/(4.0*popSize)+" 9 7 -ej "+
                (175000.0/genTime)/(4.0*popSize)+" 3 2 -ej "+
                (289000.0/genTime)/(4.0*popSize)+" 5 4 -ej "+
                (419000.0/genTime)/(4.0*popSize)+" 4 2 -ej "+
                (485000.0/genTime)/(4.0*popSize)+" 11 10 -ej "+
                (871000.0/genTime)/(4.0*popSize)+" 6 2 -ej "+
                (3765000.0/genTime)/(4.0*popSize)+" 2 1 -ej "+
                (5619000.0/genTime)/(4.0*popSize)+" 7 1 -ej "+
                (11138000.0/genTime)/(4.0*popSize)+" 10 1";

    }

    private static String MakeSampleStringHumans(int seqLength, double recombRate, int samplesPerPop, double popSize) {

        //for each biologically inspired ms sim, the overall structure is hardcoded based on which bio example it is.
        //here for mosq we have a 4 taxon tree from http://science.sciencemag.org/content/early/2017/09/27/science.aao6266/tab-figures-data and the sup. material
        double genTime = 30.0; //specifically mentioned in SM 8.2, along w 1.25 e -8 mutation rate
        //double popSize = 20000.0; //from S22

        //im actually not using the mutations output from ms for anything, so i dont need the -t theta

        //for recombination paramarameter, as taken from the docs for ms:
        //"To include crossing-over (recombination) in the model,  one uses the-rswitch and specifies the cross-over rate parameter,ρwhich is 4N0r, whereristhe probability of cross-over per generation between the ends of the locus beingsimulated.  In addition, the user specifies the numbers of sites, nsites, betweenwhich recombination can occur.  One should think of nsites as the number ofbase pairs in the locus.  For example, to simulate samples of size 15 for a locuswhich is 2,501 base pairs long with cross-over probability between adjacent basepairs of 10−8per generation, and assumingN0= 106one enters:ms 15 1000 -t 10.004 -r 100.0 2501In this case,ρ= 100, since the cross-over probability between the two ends ofthe locus is 10−8(2501−1) = 2.5×10−5, and thusρ= 4×106×2.5×10−5= 100."
        //here we use recomb rate of 1.25*10-8 * 0.88 according to the science paper (the say the ratio of recomb. rate to mutation rate is 0.88)
        //this equals to 1.1 * 10^-8
        double paramR = 4 * popSize * seqLength * recombRate; //recomb rate, as mentioned in ms doc, is going to be in "with cross-over probability between adjacent basepairs of 10−8per generation"


        //-ejs are in units of "In each case the first parameter following the switch is the timewhen  the  demographic  change  occurred,  measured  from  the  present  in  unitsof  4N0generations."

        //here for primates i have manually done all the calcs, which end up with the following -ejs, these all used fixed pop size of 40k and gen time of 20:
        /*
        **** SPLITS for ms (in chron. order) ****

        (note- time for -ej is (years/genTime)/(4*40000) or x/160000)


        FOR HUMAN SUBPOPS:

         non african (4) splits from east african (3):
         75k (75000.0) years ago

         non african/ east african (3) splits from west african (2):
         150k (150000.0) years ago

         non african/east african / west african (2) splits from central african (1):
         250k (250000.0) years ago


        */

        //currently not using inf. sites for anything, but just in case for later ill add the -t in even tho i really dont need it right now
        //using mut. rate of 2.5 * 10 -8
        //String infSitesLine = "-t " + 4*popSize*(2.5*Math.pow(10,-8));

        //giving the following command:

        return "ms "+4*samplesPerPop+" 1 -T -r "+paramR+" "+seqLength+" -I 4 "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" -ej "+
                (75000.0/genTime)/(4.0*popSize)+" 4 3 -ej "+
                (150000.0/genTime)/(4.0*popSize)+" 3 2 -ej "+
                (250000.0/genTime)/(4.0*popSize)+" 2 1";
    }


    private static String MakeSampleStringMosq(int seqLength, double recombRate, int samplesPerPop, double popSize) {

        //for each biologically inspired ms sim, the overall structure is hardcoded based on which bio example it is.
        //here for mosq we have a 6 taxon tree from fig S2 A from dingqiaos paper in mol. ecol. : https://www.readcube.com/articles/supplement?doi=10.1111%2Fmec.13544&index=0&ssl=1&st=374a48a2062267a185f5ce4b5b21507f&preview=1
        //double genTime = 30.0; //specifically mentioned in SM 8.2, along w 1.25 e -8 mutation rate
        //double popSize = 10000.0; //found a random paper with mosquito ne estimates of 10k and some of like 1 and 2k, using 10k for now

        //im actually not using the mutations output from ms for anything, so i dont need the -t theta

        //for recombination paramarameter, as taken from the docs for ms:
        //"To include crossing-over (recombination) in the model,  one uses the-rswitch and specifies the cross-over rate parameter,ρwhich is 4N0r, whereristhe probability of cross-over per generation between the ends of the locus beingsimulated.  In addition, the user specifies the numbers of sites, nsites, betweenwhich recombination can occur.  One should think of nsites as the number ofbase pairs in the locus.  For example, to simulate samples of size 15 for a locuswhich is 2,501 base pairs long with cross-over probability between adjacent basepairs of 10−8per generation, and assumingN0= 106one enters:ms 15 1000 -t 10.004 -r 100.0 2501In this case,ρ= 100, since the cross-over probability between the two ends ofthe locus is 10−8(2501−1) = 2.5×10−5, and thusρ= 4×106×2.5×10−5= 100."
        //here we use recomb rate of 1.25*10-8 * 0.88 according to the science paper (the say the ratio of recomb. rate to mutation rate is 0.88)
        //this equals to 1.1 * 10^-8
        double paramR = 4 * popSize * seqLength * recombRate; //recomb rate, as mentioned in ms doc, is going to be in "with cross-over probability between adjacent basepairs of 10−8per generation"


        //-ejs are in units of "In each case the first parameter following the switch is the timewhen  the  demographic  change  occurred,  measured  from  the  present  in  unitsof  4N0generations."

        //here for primates i have manually done all the calcs, which end up with the following -ejs, these all used fixed pop size of 40k and gen time of 20:
        /*
        **** SPLITS for ms (in chron. order) ****

        (note- time for -ej is (years/genTime)/(4*40000) or x/160000)


        FOR MOSQ POPS:

        Im taking these values from dingqiaos paper where it was already in units of ?

        (i think the units were 2N for diploid so im going to go ahead and just keep them to be 4n haploid for ms)
           from the paper: "Thebranch lengths of the phylogenetic network are given incoalescent units, where one unit equals 2Negenera-tions, where Neis the effective population size. "
           //the part from the paper wasnt actually on the fig but it should be fair to say the fig is in 2n units

        //everything NEEDS TO BE SCALED FROM 2N generations (diploid) to 4N generations (diploid) by dividing by two!

         col splits from gam:
         at time 0.6 (this was the only value not in the paper fig, luay went ahead and just told me to use it and to make the tree ultrametric to get the other times)
         0.6/2 = 0.3

        those split from ara at time:
        1.56
        1.56/2=0.78

        those split from qua at time:
        3.53
        3.53/2=1.765

        those split from mel at time:
        4.13
        4.13/2=2.065

        those split from mer at time:
        4.36
        4.36/2=2.18

        FOR FIG S2 (D) tree we have:

        col splits from gam:
         5.93 - 2.12 - 1.16
        (5.93 - 2.12 - 1.16)/2 = 1.325

        those split from ara at time:
        5.93-2.12
        (5.93-2.12)/2= 1.905

        those split from qua at time:
        5.93
        5.93/2= 2.965

        those split from mel at time:
        5.93+0.79
        5.93+0.79/2= 3.36

        those split from mer at time:
        5.93 + 0.79 + 2.24
        (5.93 + 0.79 + 2.24)/2= 4.48

        */

        //currently not using inf. sites for anything, but just in case for later ill add the -t in even tho i really dont need it right now
        //using mut. rate of 2.5 * 10 -8
        //String infSitesLine = "-t " + 4*popSize*(2.5*Math.pow(10,-8));

        //giving the following command:

        //not correct, but can test with for debugging and scoring with longer branches etc
        //this was the Figure S2 (A) tree with a 0.6 between bottom leaves, before correcting for the scaling factor of 2
        /*
        return "ms "+6*samplesPerPop+" 1 -T -r "+paramR+" "+seqLength+" -I 6 "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" -ej "+
                "0.6"+" 6 5 -ej "+
                "1.56"+" 5 4 -ej "+
                "3.53"+" 4 3 -ej "+
                "4.13"+" 3 2 -ej "+
                "4.36"+" 2 1";
        */

        //The fig S2 (A) tree, correctly scaled, (assuming ultrametric and 0.6/2 for bottom leaves)
        /*
        return "ms "+6*samplesPerPop+" 1 -T -r "+paramR+" "+seqLength+" -I 6 "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" -ej "+
                "0.3"+" 6 5 -ej "+
                "0.78"+" 5 4 -ej "+
                "1.765"+" 4 3 -ej "+
                "2.065"+" 3 2 -ej "+
                "2.18"+" 2 1";
        */

        //based on fig s2d as highlighted in some of the comments above
        return "ms "+6*samplesPerPop+" 1 -T -r "+paramR+" "+seqLength+" -I 6 "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" -ej "+
                "1.325"+" 6 5 -ej "+
                "1.905"+" 5 4 -ej "+
                "2.965"+" 4 3 -ej "+
                "3.36"+" 3 2 -ej "+
                "4.48"+" 2 1";





    }

    private static String MakeSampleStringBirds(int seqLength, double recombRate, int samplesPerPop, double popSize) {

        //for each biologically inspired ms sim, the overall structure is hardcoded based on which bio example it is.
        //here for mosq we have a 6 taxon tree from fig S2 A from dingqiaos paper in mol. ecol. : https://www.readcube.com/articles/supplement?doi=10.1111%2Fmec.13544&index=0&ssl=1&st=374a48a2062267a185f5ce4b5b21507f&preview=1
        //double genTime = 30.0; //specifically mentioned in SM 8.2, along w 1.25 e -8 mutation rate
        //double popSize = 10000.0; //found a random paper with mosquito ne estimates of 10k and some of like 1 and 2k, using 10k for now

        //double birdScaling = .0068; //for bird using 1.7*10^-9 as per https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2873279/. in it they used units of sites per year but also assumed one year gene time so it should be good for sites/gen. so 4*1000000* 1.7e-9

        //double recombinationRate = Math.pow(10.0,-8.0);

        //int popSize = 1000000; //based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2873279/, several Ne's but i used the one stated to be sort of an ancestral one (10^6)
        double genTime = 1.0; // based again on https://www.sciencedirect.com/science/article/pii/S1055790384710402 gallus gallus is AFB 12 months, most others 12 months though some others like 36 is semi common, prolly could modulate this around 1,2, and 3 for 12 24 and 36 months. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2873279/ also uses one year as gen time

        int firstDivTimeColumbAndPasserea = 67000000; //from jarvis
        int secondDivTimeGallo = 87000000; //from jarvis
        int thirdDivTimePalaeo = 100000000; //from jarvis

        double firstSplitTimeInYears = firstDivTimeColumbAndPasserea;
        double secondSplitTimeInYears = secondDivTimeGallo;
        double thirdSplitTimeInYears = thirdDivTimePalaeo;

        //im actually not using the mutations output from ms for anything, so i dont need the -t theta

        //for recombination paramarameter, as taken from the docs for ms:
        //"To include crossing-over (recombination) in the model,  one uses the-rswitch and specifies the cross-over rate parameter,ρwhich is 4N0r, whereristhe probability of cross-over per generation between the ends of the locus beingsimulated.  In addition, the user specifies the numbers of sites, nsites, betweenwhich recombination can occur.  One should think of nsites as the number ofbase pairs in the locus.  For example, to simulate samples of size 15 for a locuswhich is 2,501 base pairs long with cross-over probability between adjacent basepairs of 10−8per generation, and assumingN0= 106one enters:ms 15 1000 -t 10.004 -r 100.0 2501In this case,ρ= 100, since the cross-over probability between the two ends ofthe locus is 10−8(2501−1) = 2.5×10−5, and thusρ= 4×106×2.5×10−5= 100."
        //here we use recomb rate of 1.25*10-8 * 0.88 according to the science paper (the say the ratio of recomb. rate to mutation rate is 0.88)
        //this equals to 1.1 * 10^-8
        double paramR = 4 * popSize * seqLength * recombRate; //recomb rate, as mentioned in ms doc, is going to be in "with cross-over probability between adjacent basepairs of 10−8per generation"


        //-ejs are in units of "In each case the first parameter following the switch is the timewhen  the  demographic  change  occurred,  measured  from  the  present  in  unitsof  4N0generations."


        //giving the following command:

        return "ms "+4*samplesPerPop+" 1 -T -r "+paramR+" "+seqLength+" -I 4 "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" -ej "+
                (firstSplitTimeInYears/genTime)/(4.0*popSize)+" 4 3 -ej "+
                (secondSplitTimeInYears/genTime)/(4.0*popSize)+" 3 2 -ej "+
                (thirdSplitTimeInYears/genTime)/(4.0*popSize)+" 2 1";





    }

    private static String MakeSampleStringButterfly(int seqLength, double recombRate, int samplesPerPop, double popSize) {


        //mut rate inside here = double butterflyScaling = 0.0232; //for buttefly itll be 4*2000000* 2.9e-9 from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4271535/



        //double recombinationRate = Math.pow(10.0,-8.0);
        //String msCmd = MakeButerflyString();
        //String msCmd = MakeSampleStringMaker(50000,recombinationRate,1,40000,20,3765000,5619000,11138000);

        //int popSize = 2000000; //based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4271535/
        double genTime = 0.25; // based again on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4271535/ which said four gens per year

        double paramR = 4 * popSize * seqLength * recombRate; //recomb rate, as mentioned in ms doc, is going to be in "with cross-over probability between adjacent basepairs of 10−8per generation"


        int firstDivTimeMelpFromMelp = 100000; //from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0036464
        int secondDivTimeMelpFromTim = 1300000; //from S1 in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4769579/
        int thirdDivTimeSilv = 2100000; //from S1 in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4769579/

        double firstSplitTimeInYears = firstDivTimeMelpFromMelp;
        double secondSplitTimeInYears = secondDivTimeMelpFromTim;
        double thirdSplitTimeInYears = thirdDivTimeSilv;

        //2.5 mya estimate for tim and melp in here that im not using now, from https://journals.plos.org/plosbiology/article/figure?id=10.1371/journal.pbio.1002353.g005

        return "ms "+4*samplesPerPop+" 1 -T -r "+paramR+" "+seqLength+" -I 4 "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" "+samplesPerPop+" -ej "+
                (firstSplitTimeInYears/genTime)/(4.0*popSize)+" 4 3 -ej "+
                (secondSplitTimeInYears/genTime)/(4.0*popSize)+" 3 2 -ej "+
                (thirdSplitTimeInYears/genTime)/(4.0*popSize)+" 2 1";

    }




    private static String MakeSeqGenString(int seqLength, double mutationRate, double popSize) {

        //convert mutation rate to -s param
        //ms outputs BL as "The branch lengths are in units of 4N0 generations."
        //so multiply by 4N to pop the actual # of generations back out and then multiply generations by mut rate to get to what seq gen needs -
        //seq gen says "This option allows the user to set a value with which to scale the branch lengths in order to make them equal the expected number of substitutions per site for each branch."
        //so mut rate will be in units of expected number of substitutions per site / generation, or just substitutions per site per generation
        double paramS = 4 * popSize * mutationRate;

        //for 2.5 * 10-8, we have both:
        /*
        mutation rate based on (2.5 x 10-8) from (http://johnhawks.net/weblog/reviews/genomics/variation/human-mutation-rate-review-2010.html)
        and again we have http://www.genetics.org/content/genetics/156/1/297.full.pdf including the 2.5e-8 and says it could be up to 3.4e-8
        */

        //glue it all together. also worthy of note, can generate the random int here and it will flow through to both the phylip and fasta command runs in the experiment running code (it resuses the random seed, which has to happen)
        return "seq-gen -z " + Random.nextInt(Integer.MAX_VALUE) + " -mHKY -l "+seqLength+" -s "+paramS+" -p "+seqLength;
    }

    private static void DoExperimentGeneral(String whereToDoExperiment, String experimentName, int numberOfDuplicateRuns, String msCmd, String seqCmd, boolean runNewSims, boolean runSlowMethods, boolean runFastMethods){

        //params
        int genomeLength = 10000;


        //commands
        //msCmd = "ms 15 1 -t 15 -r 10 4000 -T";
        //seqCmd = "seq-gen -z " + Random.nextInt(Integer.MAX_VALUE) + " -mHKY -l 4000 -s .00544 -p 4000";


        String experimentFolder = whereToDoExperiment+"/"+experimentName;
        if(runNewSims){DeleteFolder(experimentFolder);} //just to be safe, if im doing new sims i for sure gotta rerun everything. someimtes things intelligently dont rerun
        MakeFolder(experimentFolder); //create experiment folder


        for(int i=1; i<=numberOfDuplicateRuns; i++) {

            String baseFolder = experimentFolder+"/Run"+i; //base folder where the exp. happens
            MakeFolder(baseFolder); //create this run's folder
            MakeFolder(baseFolder + "/results"); //create folder for results
            globalTimeWriter = baseFolder+"/timings.txt";

            //a good reference for one of the last times i was doing this is
            //in the file /Users/leo/rice/res/code/bmc17/start/doAll.sh

            if(runNewSims) {
                //run ms and seq gen
                DPrint("starting ms");
                RunMS(baseFolder, msCmd);
                DPrint("done w ms. starting seq gen");
                RunSeqGen(baseFolder, seqCmd);
                DPrint("done w seq gen");
            }

            if (runFastMethods) {
                //run RENT+, prolly first on infinite sites and then on jc, calculate score(s)

                FormatToRENT(baseFolder);
                RunRENT(baseFolder);
                ProcessRENT(baseFolder, genomeLength);
                compareTwoMSFiles(baseFolder+"/ms/trueTrees",baseFolder+"/rent/finalResultsRENT",baseFolder+"/results/rent",genomeLength);


                //debug block try 500bp/500bp raxml
                /*
                int windowSizeD=1600;
                int windowOffsetD=400;
                Windower(baseFolder, windowSizeD, windowOffsetD);
                RunWindowedRAXML(baseFolder, windowSizeD, windowOffsetD);
                ProcessRAXML(baseFolder, windowSizeD,genomeLength, windowOffsetD);
                compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/raxml" + windowSizeD + "bp" + windowOffsetD + "os" + "/treefileRAXML", baseFolder + "/results/raxml" + windowSizeD + "bp" + windowOffsetD + "os", genomeLength);
                */


                //try aaaaall the wind/offsets
                for(int iSize = 500; iSize<=3000; iSize=iSize+500){
                    for(int jOff = 50; jOff<=800; jOff=jOff*2) {

                        //window size 500 and offset of 25 gave me java heap space error for a 10k sequence
                        //when messing with the smooth winds graph creation :/ starting w 50 now..

                        //debug a fast run first
                        //iSize = 2000;
                        //jOff = 1000;


                        if(jOff > iSize){continue;}

                        System.out.println("wind " + iSize + " off "+ jOff);
                        TPrint("wind " + iSize + " off "+ jOff);

                        //can do this one time
                        Windower(baseFolder, iSize, jOff);
                        RunWindowedRAXML(baseFolder, iSize, jOff);

                        //this part is 2x one centered one without

                        boolean doCentering = true;
                        String centeringString = "-CNTRD";

                        ProcessRAXML(baseFolder, iSize, jOff,genomeLength,doCentering);
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/raxml"+centeringString + iSize + "bp" + jOff + "os" + "/treefileRAXML", baseFolder + "/results/raxml"+centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //run our method, calc score(s)
                        SmoothWinds(baseFolder, iSize, jOff, genomeLength, false, doCentering);
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwinds" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwinds" + centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //run our method (WITH BL), calc score(s)
                        SmoothWinds(baseFolder, iSize, jOff, genomeLength, true, doCentering);
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsBL" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsBL" + centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //FIRST HALF WAS WITH CENTERING, THIS HALF IS WITHOUT CENTERING

                        doCentering = false;
                        centeringString = "";

                        ProcessRAXML(baseFolder, iSize, jOff,genomeLength,doCentering);
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/raxml"+centeringString + iSize + "bp" + jOff + "os" + "/treefileRAXML", baseFolder + "/results/raxml"+centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //run our method, calc score(s)
                        SmoothWinds(baseFolder, iSize, jOff,genomeLength, false,doCentering);
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwinds"+centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwinds"+centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //run our method (WITH BL), calc score(s)
                        SmoothWinds(baseFolder, iSize, jOff,genomeLength,true,doCentering);
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsBL"+centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsBL"+centeringString + iSize + "bp" + jOff + "os", genomeLength);



                    }
                }




                //not doing these two now
                /*
                int windowSize = 1000;
                int windowOffset = 100;
                //run window passing, calc score(s)
                Windower(baseFolder, windowSize, windowOffset);

                //run raxml on each window
                RunWindowedRAXML(baseFolder, windowSize, windowOffset);
                //turn windows into an ms style gene tree file
                ProcessRAXML(baseFolder, windowSize, windowOffset,genomeLength);
                compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/raxml" + windowSize + "bp" + windowOffset + "os" + "/treefileRAXML", baseFolder + "/results/raxml" + windowSize + "bp" + windowOffset + "os", genomeLength);


                //lets try 500bp/500bp raxml
                windowSize=500;
                windowOffset=500;
                Windower(baseFolder, windowSize, windowOffset);
                RunWindowedRAXML(baseFolder, windowSize, windowOffset);
                ProcessRAXML(baseFolder, windowSize, windowOffset,genomeLength);
                compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/raxml" + windowSize + "bp" + windowOffset + "os" + "/treefileRAXML", baseFolder + "/results/raxml" + windowSize + "bp" + windowOffset + "os", genomeLength);


                //run our method, calc score(s)
                SmoothWinds(baseFolder, windowSize, windowOffset,genomeLength, false);
                compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwinds" + windowSize + "bp" + windowOffset + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwinds" + windowSize + "bp" + windowOffset + "os", genomeLength);

                //run our method (WITH BL), calc score(s)
                SmoothWinds(baseFolder, windowSize, windowOffset,genomeLength,true);
                compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsBL" + windowSize + "bp" + windowOffset + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsBL" + windowSize + "bp" + windowOffset + "os", genomeLength);
                */


            }


            if(runSlowMethods) {
                //run ARGweaver, prolly first on infinite sites and then on jc, calculate score(s)
                //seems to work but is hanging? also might have to clear out directory when i run this

                //running argw takes forever, testing right now so will not keep re running it
                RunARGW(baseFolder);
                ProcessARGW(baseFolder);
                compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/argw/finalLocalTrees", baseFolder + "/results/argw", genomeLength);

            }

            //calc scores


            //test
            //compareTwoMSFiles(baseFolder+"/ms/trueTrees",baseFolder+"/badTest.txt",baseFolder+"/results/badTest",genomeLength);


            int stopDebug = 0;

        }

        //amalgamate results from all runs
        AmalgamateResults(experimentFolder, numberOfDuplicateRuns);

    }

    private static void DoOptimalityExperimentGeneral(String whereToDoExperiment, String experimentName, int numberOfDuplicateRuns, String msCmd, String seqCmd, boolean runNewSims, boolean runSlowMethods, boolean runFastMethods){

        //params
        int genomeLength = 10000;


        //commands
        //msCmd = "ms 15 1 -t 15 -r 10 4000 -T";
        //seqCmd = "seq-gen -z " + Random.nextInt(Integer.MAX_VALUE) + " -mHKY -l 4000 -s .00544 -p 4000";


        String experimentFolder = whereToDoExperiment+"/"+experimentName;
        if(runNewSims){DeleteFolder(experimentFolder);} //just to be safe, if im doing new sims i for sure gotta rerun everything. someimtes things intelligently dont rerun
        MakeFolder(experimentFolder); //create experiment folder


        for(int i=1; i<=numberOfDuplicateRuns; i++) {

            String baseFolder = experimentFolder+"/Run"+i; //base folder where the exp. happens
            MakeFolder(baseFolder); //create this run's folder
            MakeFolder(baseFolder + "/results"); //create folder for results
            globalTimeWriter = baseFolder+"/timings.txt";

            //a good reference for one of the last times i was doing this is
            //in the file /Users/leo/rice/res/code/bmc17/start/doAll.sh

            if(runNewSims) {
                //run ms and seq gen
                DPrint("starting ms");
                RunMS(baseFolder, msCmd);
                DPrint("done w ms. starting seq gen");
                RunSeqGen(baseFolder, seqCmd);
                DPrint("done w seq gen");
            }

            if (runFastMethods) {
                //run RENT+, prolly first on infinite sites and then on jc, calculate score(s)

                FormatToRENT(baseFolder);
                RunRENT(baseFolder);
                ProcessRENT(baseFolder, genomeLength);
                compareTwoMSFiles(baseFolder+"/ms/trueTrees",baseFolder+"/rent/finalResultsRENT",baseFolder+"/results/rent",genomeLength);


                //debug block try 500bp/500bp raxml
                /*
                int windowSizeD=1600;
                int windowOffsetD=400;
                Windower(baseFolder, windowSizeD, windowOffsetD);
                RunWindowedRAXML(baseFolder, windowSizeD, windowOffsetD);
                ProcessRAXML(baseFolder, windowSizeD,genomeLength, windowOffsetD);
                compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/raxml" + windowSizeD + "bp" + windowOffsetD + "os" + "/treefileRAXML", baseFolder + "/results/raxml" + windowSizeD + "bp" + windowOffsetD + "os", genomeLength);
                */


                //try aaaaall the wind/offsets
                //slower set, more exhaustive
                //for(int iSize = 500; iSize<=3000; iSize=iSize+500){
                //    for(int jOff = 800; jOff>=100; jOff=jOff/2) {

                for(int iSize = 1000; iSize<=2000; iSize=iSize+500){
                    for(int jOff = 800; jOff>=100; jOff=jOff/2) {


                                //for now only go down to offset of 100 to make it go faster
                                //in the future would like to check out 50, 25, and maybe even 12 (cast int of 12.5, add cast to division)


                                //window size 500 and offset of 25 gave me java heap space error for a 10k sequence
                                //when messing with the smooth winds graph creation :/ starting w 50 now..

                                //debug a fast run first
                                //iSize = 2000;
                                //jOff = 1000;


                                if(jOff > iSize){continue;}

                                System.out.println("wind " + iSize + " off "+ jOff);
                                TPrint("wind " + iSize + " off "+ jOff);

                                //can do this one time
                                Windower(baseFolder, iSize, jOff);
                                RunWindowedRAXML(baseFolder, iSize, jOff);

                                //this part is 2x one centered one without

                                boolean doCentering = true;
                                String centeringString = "-CNTRD";

                                //optimality check thing (not doing any more for now - prolly a bit too misleading honestly cuz of how everything works and with it allowing switching every time ms switches)
                                //compareMSFileVsOptimal(baseFolder + "/ms/trueTrees", baseFolder, baseFolder + "/results/OPTraxml" +centeringString+ iSize + "bp" + jOff + "os", genomeLength,iSize, jOff, doCentering);

                                ProcessRAXML(baseFolder, iSize, jOff,genomeLength,doCentering);
                                compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/raxml"+centeringString + iSize + "bp" + jOff + "os" + "/treefileRAXML", baseFolder + "/results/raxml"+centeringString + iSize + "bp" + jOff + "os", genomeLength);

                                //run our method, calc score(s)
                                //memory errors for now

                        //mem saver so i can actually run low offsets
                        //seems to be working correctly now - but not running it this second, now i am lets see if this works
                        //if i check notes in timing i should now make sure the rf dists for shortest paths now match below since i
                        //moved this code up one line
                        SmoothWindsLessMem(baseFolder, iSize, jOff, genomeLength, false, doCentering);
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //above worked but below (smoothwinds function, with "Exception in thread "main" java.lang.OutOfMemoryError: GC overhead limit exceeded") failed when setting was...500window size and 25 offset
                        System.out.println("less mem done ws "+iSize+" off "+jOff);

                        //too much mem for small offsets etc
                        //SmoothWinds(baseFolder, iSize, jOff, genomeLength, false, doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwinds" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwinds" + centeringString + iSize + "bp" + jOff + "os", genomeLength);
                        //System.out.println("done ws "+iSize+" off "+jOff);

                        //mem saver so i can actually run low offsets
                        //seems to be working correctly now - but not running it this second
                        /*
                        SmoothWindsLessMem(baseFolder, iSize, jOff, genomeLength, false, doCentering);
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os", genomeLength);
                        */

                        //run our method (WITH BL), calc score(s)
                        //SmoothWinds(baseFolder, iSize, jOff, genomeLength, true, doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsBL" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsBL" + centeringString + iSize + "bp" + jOff + "os", genomeLength);


                                //FIRST HALF WAS WITH CENTERING, THIS HALF IS WITHOUT CENTERING

                                doCentering = false;
                                centeringString = "";

                                //this does the optimality comparison (not doing any more for now - prolly a bit too misleading honestly cuz of how everything works and with it allowing switching every time ms switches)
                                compareMSFileVsOptimal(baseFolder + "/ms/trueTrees", baseFolder, baseFolder + "/results/OPTraxml"+centeringString + iSize + "bp" + jOff + "os", genomeLength,iSize, jOff, doCentering);


                                ProcessRAXML(baseFolder, iSize, jOff,genomeLength,doCentering);
                                compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/raxml"+centeringString + iSize + "bp" + jOff + "os" + "/treefileRAXML", baseFolder + "/results/raxml" +centeringString + iSize + "bp" + jOff + "os", genomeLength);

                                //run our method, calc score(s)
                                //memory errors for now

                        SmoothWindsLessMem(baseFolder, iSize, jOff, genomeLength, false, doCentering);
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //before, before i made less mem one
                        //SmoothWinds(baseFolder, iSize, jOff,genomeLength, false,doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwinds"+centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwinds"+centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //run our method (WITH BL), calc score(s)
                        //SmoothWinds(baseFolder, iSize, jOff,genomeLength,true,doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsBL"+centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsBL"+centeringString + iSize + "bp" + jOff + "os", genomeLength);



                            }
                        }
                    }






            if(runSlowMethods) {
                //run ARGweaver, prolly first on infinite sites and then on jc, calculate score(s)
                //seems to work but is hanging? also might have to clear out directory when i run this

                //running argw takes forever, testing right now so will not keep re running it
                RunARGW(baseFolder);
                ProcessARGW(baseFolder);
                compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/argw/finalLocalTrees", baseFolder + "/results/argw", genomeLength);

            }

            //calc scores


            //test
            //compareTwoMSFiles(baseFolder+"/ms/trueTrees",baseFolder+"/badTest.txt",baseFolder+"/results/badTest",genomeLength);


            int stopDebug = 0;

        }

        //amalgamate results from all runs
        AmalgamateResults(experimentFolder, numberOfDuplicateRuns);

    }



    private static void AnalyzeHCGData(String whereToDoExperiment, String experimentName, boolean runRENT, boolean runSW, boolean runARGW){


        String experimentFolder = whereToDoExperiment+"/"+experimentName;
        //if(runNewSims){DeleteFolder(experimentFolder);} //just to be safe, if im doing new sims i for sure gotta rerun everything. someimtes things intelligently dont rerun
        MakeFolder(experimentFolder); //create experiment folder

        //hardcoded for now
        String whereTheDataIs = "/Users/leo/rice/res/data/hobolth/target106";
        //first gonna do this stuff for a single file
        String analyzeThisFile = "/Users/leo/rice/res/data/hobolth/target106/output.2.fas";
        String analyzeThisFileName = analyzeThisFile.split("/")[analyzeThisFile.split("/").length-1];

        //int thisFileGenomeLength = Integer.parseInt(analyzeThisFile.split("\\.")[2]);
        //for hcg i am just punching this in
        int thisFileGenomeLength = 32609;

        String baseFolder = experimentFolder; //base folder where the exp. happens
        //MakeFolder(baseFolder); //create this run's folder
        MakeFolder(baseFolder + "/results"); //create folder for results
        globalTimeWriter = baseFolder+"/timings.txt";

        //a good reference for one of the last times i was doing this is
        //in the file /Users/leo/rice/res/code/bmc17/start/doAll.sh



        if (runRENT) {
            //run RENT+, prolly first on infinite sites and then on jc, calculate score(s)

            FormatRunProcessRENTSingleFile(baseFolder,thisFileGenomeLength,analyzeThisFile);

            //FormatToRENT(baseFolder);
            //RunRENT(baseFolder);
            //ProcessRENT(baseFolder, thisFileLength);
            //compareTwoMSFiles(baseFolder+"/ms/trueTrees",baseFolder+"/rent/finalResultsRENT",baseFolder+"/results/rent",genomeLength);
        }



        if(runARGW) {
            //run ARGweaver, prolly first on infinite sites and then on jc, calculate score(s)
            //seems to work but is hanging? also might have to clear out directory when i run this

            //running argw takes forever, testing right now so will not keep re running it
            RunARGWSingleFile(baseFolder,analyzeThisFile,true);
            ProcessARGWSingleFile(baseFolder,analyzeThisFile);
            //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/argw/finalLocalTrees", baseFolder + "/results/argw", genomeLength);

        }


        if(runSW){
            //try aaaaall the wind/offsets
            //slower set, more exhaustive
            //for(int iSize = 500; iSize<=3000; iSize=iSize+500){
            //    for(int jOff = 800; jOff>=100; jOff=jOff/2) {

            //region SW Stuff

                        //for now only go down to offset of 100 to make it go faster
                        //in the future would like to check out 50, 25, and maybe even 12 (cast int of 12.5, add cast to division)


                        //window size 500 and offset of 25 gave me java heap space error for a 10k sequence
                        //when messing with the smooth winds graph creation :/ starting w 50 now..

                        //debug a fast run first
                        //iSize = 2000;
                        //jOff = 1000;


                        int iSize = 1000;
                        int jOff = 50;
                        int genomeLength = thisFileGenomeLength;

                        //can do this one time
                        Windower(baseFolder, iSize, jOff);
                        RunWindowedRAXML(baseFolder, iSize, jOff);

                        //this part is 2x one centered one without

                        boolean doCentering = true;
                        String centeringString = "-CNTRD";

                        //optimality check thing (not doing any more for now - prolly a bit too misleading honestly cuz of how everything works and with it allowing switching every time ms switches)
                        //compareMSFileVsOptimal(baseFolder + "/ms/trueTrees", baseFolder, baseFolder + "/results/OPTraxml" +centeringString+ iSize + "bp" + jOff + "os", genomeLength,iSize, jOff, doCentering);

                        ProcessRAXML(baseFolder, iSize, jOff,genomeLength,doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/raxml"+centeringString + iSize + "bp" + jOff + "os" + "/treefileRAXML", baseFolder + "/results/raxml"+centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //run our method, calc score(s)
                        //memory errors for now

                        //mem saver so i can actually run low offsets
                        //seems to be working correctly now - but not running it this second, now i am lets see if this works
                        //if i check notes in timing i should now make sure the rf dists for shortest paths now match below since i
                        //moved this code up one line
                        SmoothWindsLessMem(baseFolder, iSize, jOff, genomeLength, false, doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //above worked but below (smoothwinds function, with "Exception in thread "main" java.lang.OutOfMemoryError: GC overhead limit exceeded") failed when setting was...500window size and 25 offset
                        System.out.println("less mem done ws "+iSize+" off "+jOff);

                        //too much mem for small offsets etc
                        //SmoothWinds(baseFolder, iSize, jOff, genomeLength, false, doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwinds" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwinds" + centeringString + iSize + "bp" + jOff + "os", genomeLength);
                        //System.out.println("done ws "+iSize+" off "+jOff);

                        //mem saver so i can actually run low offsets
                        //seems to be working correctly now - but not running it this second

                        //SmoothWindsLessMem(baseFolder, iSize, jOff, genomeLength, false, doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os", genomeLength);


                        //run our method (WITH BL), calc score(s)
                        //SmoothWinds(baseFolder, iSize, jOff, genomeLength, true, doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsBL" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsBL" + centeringString + iSize + "bp" + jOff + "os", genomeLength);


                        //FIRST HALF WAS WITH CENTERING, THIS HALF IS WITHOUT CENTERING

                        //doCentering = false;
                        //centeringString = "";
                        //ProcessRAXML(baseFolder, iSize, jOff,genomeLength,doCentering);
                        //SmoothWindsLessMem(baseFolder, iSize, jOff, genomeLength, false, doCentering);


            //endregion SW Stuff

        }

        //calc scores


        //test
        //compareTwoMSFiles(baseFolder+"/ms/trueTrees",baseFolder+"/badTest.txt",baseFolder+"/results/badTest",genomeLength);
        compareTwoMSFiles(baseFolder+"/rent/finished/"+analyzeThisFileName+".finalTrees",baseFolder+"/argw/"+analyzeThisFileName+"/finalLocalTrees",baseFolder+"/COMPARE"+analyzeThisFileName,thisFileGenomeLength);


        int stopDebug = 0;



        //amalgamate results from all runs
        //AmalgamateResults(experimentFolder, numberOfDuplicateRuns);

    }

    private static void AnalyzeMosqData(String whereToDoExperiment, String experimentName, boolean runRENT, boolean runSW, boolean runARGW){


        String experimentFolder = whereToDoExperiment+"/"+experimentName;
        //if(runNewSims){DeleteFolder(experimentFolder);} //just to be safe, if im doing new sims i for sure gotta rerun everything. someimtes things intelligently dont rerun
        MakeFolder(experimentFolder); //create experiment folder

        //hardcoded for now
        String whereTheDataIs = "/Users/leo/rice/res/data/mosquitoDing/gene/2L";
        //first gonna do this stuff for a single file
        String analyzeThisFile = "/Users/leo/rice/res/data/mosquitoDing/gene/2L/2L.1316545.10000.fa";
        String analyzeThisFileName = analyzeThisFile.split("/")[analyzeThisFile.split("/").length-1];
        int thisFileGenomeLength = Integer.parseInt(analyzeThisFile.split("\\.")[2]);


            String baseFolder = experimentFolder; //base folder where the exp. happens
            //MakeFolder(baseFolder); //create this run's folder
            MakeFolder(baseFolder + "/results"); //create folder for results
            globalTimeWriter = baseFolder+"/timings.txt";

            //a good reference for one of the last times i was doing this is
            //in the file /Users/leo/rice/res/code/bmc17/start/doAll.sh



            if (runRENT) {
                //run RENT+, prolly first on infinite sites and then on jc, calculate score(s)

                FormatRunProcessRENTSingleFile(baseFolder,thisFileGenomeLength,analyzeThisFile);

                //FormatToRENT(baseFolder);
                //RunRENT(baseFolder);
                //ProcessRENT(baseFolder, thisFileLength);
                //compareTwoMSFiles(baseFolder+"/ms/trueTrees",baseFolder+"/rent/finalResultsRENT",baseFolder+"/results/rent",genomeLength);
            }



            if(runARGW) {
                //run ARGweaver, prolly first on infinite sites and then on jc, calculate score(s)
                //seems to work but is hanging? also might have to clear out directory when i run this

                //running argw takes forever, testing right now so will not keep re running it
                RunARGWSingleFile(baseFolder,analyzeThisFile,true);
                ProcessARGWSingleFile(baseFolder,analyzeThisFile);
                //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/argw/finalLocalTrees", baseFolder + "/results/argw", genomeLength);

            }


            if(runSW){
                //try aaaaall the wind/offsets
                //slower set, more exhaustive
                //for(int iSize = 500; iSize<=3000; iSize=iSize+500){
                //    for(int jOff = 800; jOff>=100; jOff=jOff/2) {

                //region SW Stuff
                /*
                for(int iSize = 1000; iSize<=2000; iSize=iSize+500){
                    for(int jOff = 800; jOff>=100; jOff=jOff/2) {


                        //for now only go down to offset of 100 to make it go faster
                        //in the future would like to check out 50, 25, and maybe even 12 (cast int of 12.5, add cast to division)


                        //window size 500 and offset of 25 gave me java heap space error for a 10k sequence
                        //when messing with the smooth winds graph creation :/ starting w 50 now..

                        //debug a fast run first
                        //iSize = 2000;
                        //jOff = 1000;


                        if(jOff > iSize){continue;}

                        System.out.println("wind " + iSize + " off "+ jOff);
                        TPrint("wind " + iSize + " off "+ jOff);

                        //can do this one time
                        Windower(baseFolder, iSize, jOff);
                        RunWindowedRAXML(baseFolder, iSize, jOff);

                        //this part is 2x one centered one without

                        boolean doCentering = true;
                        String centeringString = "-CNTRD";

                        //optimality check thing (not doing any more for now - prolly a bit too misleading honestly cuz of how everything works and with it allowing switching every time ms switches)
                        //compareMSFileVsOptimal(baseFolder + "/ms/trueTrees", baseFolder, baseFolder + "/results/OPTraxml" +centeringString+ iSize + "bp" + jOff + "os", genomeLength,iSize, jOff, doCentering);

                        ProcessRAXML(baseFolder, iSize, jOff,genomeLength,doCentering);
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/raxml"+centeringString + iSize + "bp" + jOff + "os" + "/treefileRAXML", baseFolder + "/results/raxml"+centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //run our method, calc score(s)
                        //memory errors for now

                        //mem saver so i can actually run low offsets
                        //seems to be working correctly now - but not running it this second, now i am lets see if this works
                        //if i check notes in timing i should now make sure the rf dists for shortest paths now match below since i
                        //moved this code up one line
                        SmoothWindsLessMem(baseFolder, iSize, jOff, genomeLength, false, doCentering);
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //above worked but below (smoothwinds function, with "Exception in thread "main" java.lang.OutOfMemoryError: GC overhead limit exceeded") failed when setting was...500window size and 25 offset
                        System.out.println("less mem done ws "+iSize+" off "+jOff);

                        //too much mem for small offsets etc
                        //SmoothWinds(baseFolder, iSize, jOff, genomeLength, false, doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwinds" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwinds" + centeringString + iSize + "bp" + jOff + "os", genomeLength);
                        //System.out.println("done ws "+iSize+" off "+jOff);

                        //mem saver so i can actually run low offsets
                        //seems to be working correctly now - but not running it this second

                        //SmoothWindsLessMem(baseFolder, iSize, jOff, genomeLength, false, doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os", genomeLength);


                        //run our method (WITH BL), calc score(s)
                        //SmoothWinds(baseFolder, iSize, jOff, genomeLength, true, doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsBL" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsBL" + centeringString + iSize + "bp" + jOff + "os", genomeLength);


                        //FIRST HALF WAS WITH CENTERING, THIS HALF IS WITHOUT CENTERING

                        doCentering = false;
                        centeringString = "";

                        //this does the optimality comparison (not doing any more for now - prolly a bit too misleading honestly cuz of how everything works and with it allowing switching every time ms switches)
                        compareMSFileVsOptimal(baseFolder + "/ms/trueTrees", baseFolder, baseFolder + "/results/OPTraxml"+centeringString + iSize + "bp" + jOff + "os", genomeLength,iSize, jOff, doCentering);


                        ProcessRAXML(baseFolder, iSize, jOff,genomeLength,doCentering);
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/raxml"+centeringString + iSize + "bp" + jOff + "os" + "/treefileRAXML", baseFolder + "/results/raxml" +centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //run our method, calc score(s)
                        //memory errors for now

                        SmoothWindsLessMem(baseFolder, iSize, jOff, genomeLength, false, doCentering);
                        compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsLessMem" + centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //before, before i made less mem one
                        //SmoothWinds(baseFolder, iSize, jOff,genomeLength, false,doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwinds"+centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwinds"+centeringString + iSize + "bp" + jOff + "os", genomeLength);

                        //run our method (WITH BL), calc score(s)
                        //SmoothWinds(baseFolder, iSize, jOff,genomeLength,true,doCentering);
                        //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwindsBL"+centeringString + iSize + "bp" + jOff + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwindsBL"+centeringString + iSize + "bp" + jOff + "os", genomeLength);



                    }
                }
                */
                //endregion SW Stuff

            }

            //calc scores


            //test
            //compareTwoMSFiles(baseFolder+"/ms/trueTrees",baseFolder+"/badTest.txt",baseFolder+"/results/badTest",genomeLength);
        compareTwoMSFiles(baseFolder+"/rent/finished/"+analyzeThisFileName+".finalTrees",baseFolder+"/argw/"+analyzeThisFileName+"/finalLocalTrees",baseFolder+"/COMPARE"+analyzeThisFileName,thisFileGenomeLength);


        int stopDebug = 0;



        //amalgamate results from all runs
        //AmalgamateResults(experimentFolder, numberOfDuplicateRuns);

    }

    private static void FormatRunProcessRENTSingleFile(String baseFolder, int genomeLength, String analyzeThisFile) {

        MakeFolder(baseFolder+"/rent");

        FormatToRENTSingleFile(baseFolder,analyzeThisFile);
        RunRENTSingleFile(baseFolder,analyzeThisFile);
        ProcessRENTSingleFile(baseFolder,genomeLength,analyzeThisFile);

    }


    private static void DPrint(String s) {
        if(DEBUG){
            //System.out.println(new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.SSS").format(new Date())+": "+s);
            System.out.println(new SimpleDateFormat("HH:mm:ss").format(new Date())+"| "+s);
        }else{
            //do nothing
        }
    }

    private static void TPrint(String s) {

        File file = new File(globalTimeWriter);

        try {
            if (!file.exists()) {
                //write line to file

                BufferedWriter output = new BufferedWriter(new FileWriter(globalTimeWriter));
                output.append(new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date()) + "| " + s+"\n");
                output.close();

                //System.out.println(new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.SSS").format(new Date())+": "+s);
                //System.out.println(new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date()) + "| " + s);

            } else {
                //append line to file
                BufferedWriter output = new BufferedWriter(new FileWriter(globalTimeWriter, true));
                output.append(new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date()) + "| " + s+"\n");
                output.close();
                //System.out.println(new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date()) + "| " + s);


            }
        }catch(Exception e){
            e.printStackTrace();
        }


    }



    private static void AmalgamateResults(String experimentFolder, int runs) {

        //there will always be a Run1, grab the names of the final results files that were generated from the Run1 results folder (there will be one of these in each run
        // (from https://stackoverflow.com/questions/5694385/getting-the-filenames-of-all-files-in-a-folder)
        List<String> resultsNames = new ArrayList<String>();
        HashMap<String,ArrayList<Double>> resultsAmalg = new HashMap<String,ArrayList<Double>>(); //this will gather all the results
        HashMap<String,Double> resultsTotals = new HashMap<String,Double>(); //this will gather all the results


        //File[] files = folder.listFiles((dir, name) -> !name.equals(".DS_Store")); - nice line from https://stackoverflow.com/questions/30486404/java-list-files-in-a-folder-avoiding-ds-store
        //to avoid the dstupid ds store file on macs
        File[] resultsFiles = new File(experimentFolder+"/Run1/results").listFiles((dir, name) -> !name.equals(".DS_Store")); //If this pathname does not denote a directory, then listFiles() returns null.
        for (File file : resultsFiles) {
            if (file.isFile() && !file.getName().contains("RFvsNumSites")) {
                resultsNames.add(file.getName());
                resultsAmalg.put(file.getName(), new ArrayList());
                resultsTotals.put(file.getName(),0.0);
            }

        }

        //now amalgamate for each run and each result file in that run
        for(int i=1; i<=runs; i++) {
            for(String thisResult : resultsNames){

                //grab the RF value
                try {
                    BufferedReader brRF = new BufferedReader(new FileReader(experimentFolder+"/Run"+i+"/results/"+thisResult));
                    double thisRF = Double.parseDouble(brRF.readLine().split(":")[1]);
                    resultsAmalg.get(thisResult).add(thisRF);
                    resultsTotals.put(thisResult,resultsTotals.get(thisResult) + thisRF);
                } catch (Exception e) {
                    e.printStackTrace();
                }


            }
        }

        //print out a file for each amalgamed results, one big one to maybe make plotting easier, and the same but for averages for each
        try {
            MakeFolder(experimentFolder + "/FinalResults");

            BufferedWriter writeAll = new BufferedWriter(new FileWriter(experimentFolder + "/FinalResults/all"));
            BufferedWriter writeAllAverage = new BufferedWriter(new FileWriter(experimentFolder + "/FinalResults/ave_all"));


            //text for final all file
            ArrayList<StringBuilder> writeAllText = new ArrayList<StringBuilder>();
            StringBuilder averagesText = new StringBuilder();
            IntStream.rangeClosed(0, runs).forEach(i -> writeAllText.add(new StringBuilder())); // populate empties

            //singles and build bigger files
            for (String thisResult : resultsNames) {

                //the average rf dist for this method across all runs
                String thisResultAverageRF = String.valueOf(resultsTotals.get(thisResult) / (double) runs);

                //add to writeAllText first line, and values to the rest of the lines, and also add a line for averages
                writeAllText.get(0).append(thisResult + ",");
                IntStream.rangeClosed(1, runs).forEach(i -> writeAllText.get(i).append(resultsAmalg.get(thisResult).get(i - 1) + ",")); // populate
                averagesText.append(thisResultAverageRF + ",");

                //writers
                BufferedWriter writeThisResult = new BufferedWriter(new FileWriter(experimentFolder + "/FinalResults/" + thisResult));
                BufferedWriter writeThisResultAverage = new BufferedWriter(new FileWriter(experimentFolder + "/FinalResults/ave_" + thisResult));

                //write ave
                writeThisResultAverage.write(thisResultAverageRF);

                //write list of all
                for (double singleValueRF : resultsAmalg.get(thisResult)) {
                    writeThisResult.write(singleValueRF + "\n");
                }

                //close writers
                writeThisResult.close();
                writeThisResultAverage.close();

            }

            //write bigger ones. first the full one then the averages file which is easy
            for (StringBuilder writeAllTextLine : writeAllText) {
                writeAll.write(writeAllTextLine.toString() + "\n");
            }

            //i made it so i can just piggy back off first line of the last file and then just add average line
            writeAllAverage.write(writeAllText.get(0).toString() + "\n");
            writeAllAverage.write(averagesText.toString());

            //close writers
            writeAll.close();
            writeAllAverage.close();

            //******
            //writing the stuff formatted for R
            // had to add this dumb stuff at the very end so its poorly coded
            // r required this style of formatting to make the right kinds of figs
            //******

            BufferedWriter rFormattingAll = new BufferedWriter(new FileWriter(experimentFolder + "/FinalResults/r_all"));
            BufferedWriter rFormattingAllAverage = new BufferedWriter(new FileWriter(experimentFolder + "/FinalResults/r_ave_all"));

            //titles
            rFormattingAll.write("Method,Accuracy\n");
            rFormattingAllAverage.write("Method,Accuracy\n");

            //write lines
            String[] methods = writeAllText.get(0).toString().split(","); //has a comma at end
            String[] averages = averagesText.toString().split(","); //has a comma at end
            for (int methodIndex = 0; methodIndex < methods.length-1; methodIndex++){
                //for writing all info
                for (int runIndex = 1; runIndex <= runs; runIndex++) {
                    String[] thisRunInfo = writeAllText.get(runIndex).toString().split(","); //has a comma at end
                    rFormattingAll.write(methods[methodIndex]+","+thisRunInfo[methodIndex]+"\n");
                }
                //just writes average
                rFormattingAllAverage.write(methods[methodIndex]+","+averages[methodIndex]+"\n");
        }

            rFormattingAll.close();
            rFormattingAllAverage.close();


        }catch(Exception e){
            e.printStackTrace();
        }

        //int doStuffWithResults = 0;

    }


    private static void AmalgamateTimingsResults(String experimentFolder, int runs) {

        //there will always be a Run1, grab the names of the final results files that were generated from the Run1 results folder (there will be one of these in each run
        // (from https://stackoverflow.com/questions/5694385/getting-the-filenames-of-all-files-in-a-folder)
        List<String> resultsNames = new ArrayList<String>();
        HashMap<String,ArrayList<Double>> resultsAmalg = new HashMap<String,ArrayList<Double>>(); //this will gather all the results
        HashMap<String,Double> resultsTotals = new HashMap<String,Double>(); //this will gather all the results


        //File[] files = folder.listFiles((dir, name) -> !name.equals(".DS_Store")); - nice line from https://stackoverflow.com/questions/30486404/java-list-files-in-a-folder-avoiding-ds-store
        //to avoid the dstupid ds store file on macs
        File[] resultsFiles = new File(experimentFolder+"/Run1/timingsFolder").listFiles((dir, name) -> !name.equals(".DS_Store")); //If this pathname does not denote a directory, then listFiles() returns null.
        for (File file : resultsFiles) {
            if (file.isFile()) {
                resultsNames.add(file.getName());
                resultsAmalg.put(file.getName(), new ArrayList());
                resultsTotals.put(file.getName(),0.0);
            }

        }

        //now amalgamate for each run and each result file in that run
        for(int i=1; i<=runs; i++) {
            for(String thisResult : resultsNames){

                //grab the RF value
                try {
                    BufferedReader brRF = new BufferedReader(new FileReader(experimentFolder+"/Run"+i+"/timingsFolder/"+thisResult));
                    double thisSeconds = Double.parseDouble(brRF.readLine().split(":")[1]);
                    resultsAmalg.get(thisResult).add(thisSeconds);
                    resultsTotals.put(thisResult,resultsTotals.get(thisResult) + thisSeconds);
                } catch (Exception e) {
                    e.printStackTrace();
                }


            }
        }

        //print out a file for each amalgamed results, one big one to maybe make plotting easier, and the same but for averages for each
        try{
            MakeFolder(experimentFolder+"/FinalTimings");

            BufferedWriter writeAll = new BufferedWriter(new FileWriter(experimentFolder+"/FinalTimings/all"));
            BufferedWriter writeAllAverage = new BufferedWriter(new FileWriter(experimentFolder+"/FinalTimings/ave_all"));

            //text for final all file
            ArrayList<StringBuilder> writeAllText = new ArrayList<StringBuilder>();
            StringBuilder averagesText = new StringBuilder();
            IntStream.rangeClosed(0, runs).forEach(i -> writeAllText.add(new StringBuilder())); // populate empties

            //singles and build bigger files
            for(String thisResult : resultsNames){

                //the average rf dist for this method across all runs
                String thisResultAverageSeconds = String.valueOf(resultsTotals.get(thisResult) / (double) runs);

                //add to writeAllText first line, and values to the rest of the lines, and also add a line for averages
                writeAllText.get(0).append(thisResult+",");
                IntStream.rangeClosed(1, runs).forEach(i -> writeAllText.get(i).append(resultsAmalg.get(thisResult).get(i-1) + ",")); // populate
                averagesText.append(thisResultAverageSeconds+",");

                //writers
                BufferedWriter writeThisResult = new BufferedWriter(new FileWriter(experimentFolder+"/FinalResults/"+thisResult));
                BufferedWriter writeThisResultAverage = new BufferedWriter(new FileWriter(experimentFolder+"/FinalResults/ave_"+thisResult));

                //write ave
                writeThisResultAverage.write(thisResultAverageSeconds);

                //write list of all
                for(double singleValueRF : resultsAmalg.get(thisResult)){
                    writeThisResult.write(singleValueRF+"\n");
                }

                //close writers
                writeThisResult.close();
                writeThisResultAverage.close();
            }

            //write bigger ones. first the full one then the averages file which is easy
            for(StringBuilder writeAllTextLine : writeAllText) {
                writeAll.write(writeAllTextLine.toString()+"\n");
            }

            //i made it so i can just piggy back off first line of the last file and then just add average line
            writeAllAverage.write(writeAllText.get(0).toString()+"\n");
            writeAllAverage.write(averagesText.toString());

            //close writers
            writeAll.close();
            writeAllAverage.close();

            //******
            //writing the stuff formatted for R
            // had to add this dumb stuff at the very end so its poorly coded
            // r required this style of formatting to make the right kinds of figs
            //******

            BufferedWriter rFormattingAll = new BufferedWriter(new FileWriter(experimentFolder + "/FinalTimings/r_all"));
            BufferedWriter rFormattingAllAverage = new BufferedWriter(new FileWriter(experimentFolder + "/FinalTimings/r_ave_all"));

            //titles
            rFormattingAll.write("Method,Time\n");
            rFormattingAllAverage.write("Method,Time\n");

            //write lines
            String[] methods = writeAllText.get(0).toString().split(","); //has a comma at end
            String[] averages = averagesText.toString().split(","); //has a comma at end
            for (int methodIndex = 0; methodIndex < methods.length-1; methodIndex++){
                //for writing all info
                for (int runIndex = 1; runIndex <= runs; runIndex++) {
                    String[] thisRunInfo = writeAllText.get(runIndex).toString().split(","); //has a comma at end
                    rFormattingAll.write(methods[methodIndex]+","+thisRunInfo[methodIndex]+"\n");
                }
                //just writes average
                rFormattingAllAverage.write(methods[methodIndex]+","+averages[methodIndex]+"\n");
            }

            rFormattingAll.close();
            rFormattingAllAverage.close();




        }catch(Exception e){
            e.printStackTrace();
        }

        //int doStuffWithResults = 0;

    }








    //region OLD CODE
    /*
    private static void DoExperiment1(String whereToDoExperiment, String experimentName, int numberOfDuplicateRuns){

        //params
        int genomeLength = 4000;
        int windowSize = 1000;
        int windowOffset = 100;


        System.out.println("exp1!");
        String experimentFolder = whereToDoExperiment+"/"+experimentName;
        MakeFolder(experimentFolder); //create experiment folder




        for(int i=1; i<=numberOfDuplicateRuns; i++) {

            String baseFolder = experimentFolder+"/Run"+i; //base folder where the exp. happens
            MakeFolder(baseFolder); //create this run's folder
            MakeFolder(baseFolder + "/results"); //create folder for results

            //a good reference for one of the last times i was doing this is
            //in the file /Users/leo/rice/res/code/bmc17/start/doAll.sh

            //run ms and seq gen
            RunMS(baseFolder, "ms 15 1 -t 15 -r 10 4000 -T");

            RunSeqGen(baseFolder, "seq-gen -z " + Random.nextInt(Integer.MAX_VALUE) + " -mHKY -l 4000 -s .00544 -p 4000");

            //run RENT+, prolly first on infinite sites and then on jc, calculate score(s)
            FormatToRENT(baseFolder);
            RunRENT(baseFolder);
            ProcessRENT(baseFolder, genomeLength);

            //run ARGweaver, prolly first on infinite sites and then on jc, calculate score(s)
            //seems to work but is hanging? also might have to clear out directory when i run this

            //running argw takes forever, testing right now so will not keep re running it
            //RunARGW(baseFolder);
            ProcessARGW(baseFolder);

            //run window passing, calc score(s)
            Windower(baseFolder, windowSize, windowOffset);

            //run raxml on each window
            RunWindowedRAXML(baseFolder, windowSize, windowOffset);
            //turn windows into an ms style gene tree file
            ProcessRAXML(baseFolder, windowSize, windowOffset);
            compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/raxml" + windowSize + "bp" + windowOffset + "os" + "/treefileRAXML", baseFolder+"/results/raxml"+windowSize+"bp"+windowOffset+"os",genomeLength);


            //lets try 500bp/500bp raxml
            windowSize=500;
            windowOffset=500;
            Windower(baseFolder, windowSize, windowOffset);
            RunWindowedRAXML(baseFolder, windowSize, windowOffset);
            ProcessRAXML(baseFolder, windowSize, windowOffset);
            compareTwoMSFiles(baseFolder+"/ms/trueTrees",baseFolder+"/raxml"+windowSize+"bp"+windowOffset+"os"+"/treefileRAXML",baseFolder+"/results/raxml"+windowSize+"bp"+windowOffset+"os",genomeLength);


            //run our method, calc score(s)
            SmoothWinds(baseFolder, windowSize, windowOffset,false);
            compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/smoothwinds" + windowSize + "bp" + windowOffset + "os" + "/treefileSmoothWinds", baseFolder + "/results/smoothwinds" + windowSize + "bp" + windowOffset + "os", genomeLength);


            //calc scores

            compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/argw/finalLocalTrees", baseFolder + "/results/argw", genomeLength);
            compareTwoMSFiles(baseFolder+"/ms/trueTrees",baseFolder+"/rent/finalResultsRENT",baseFolder+"/results/rent",genomeLength);


            int stopDebug = 0;

        }


        //after all duplicate runs, amalgamate results:
        AmalgamateResults(experimentFolder,numberOfDuplicateRuns);




    }

    private static void DoExperiment2(String whereToDoExperiment, String experimentName, int numberOfDuplicateRuns){

        //params
        int genomeLength = 40;
        int windowSize = 1000;
        int windowOffset = 100;


        System.out.println("exp2!");
        String experimentFolder = whereToDoExperiment+"/"+experimentName;
        MakeFolder(experimentFolder); //create experiment folder




        for(int i=1; i<=numberOfDuplicateRuns; i++) {

            String baseFolder = experimentFolder+"/Run"+i; //base folder where the exp. happens
            MakeFolder(baseFolder); //create this run's folder

            //a good reference for one of the last times i was doing this is
            //in the file /Users/leo/rice/res/code/bmc17/start/doAll.sh

            //run ms and seq gen
            RunMS(baseFolder, "ms 15 1 -t 15 -r 10 40 -T");

            //RunSeqGen(baseFolder,"seq-gen -z "+Random.nextInt(Integer.MAX_VALUE)+" -mHKY -l 40 -s .00544 -p 40");
            //high mutation rate
            RunSeqGen(baseFolder, "seq-gen -z " + Random.nextInt(Integer.MAX_VALUE) + " -mHKY -l 40 -s 1.00544 -p 40");

            //run RENT+, prolly first on infinite sites and then on jc, calculate score(s)
            FormatToRENT(baseFolder);
            RunRENT(baseFolder);
            ProcessRENT(baseFolder,genomeLength);

            //run ARGweaver, prolly first on infinite sites and then on jc, calculate score(s)
            //seems to work but is hanging? also might have to clear out directory when i run this

            //running argw takes forever, testing right now so will not keep re running it
            //RunARGW(baseFolder);
            //ProcessARGW(baseFolder);

            //run window passing, calc score(s)
            //Windower(baseFolder,windowSize,windowOffset);

            //run raxml on each window
            //RunWindowedRAXML(baseFolder,windowSize,windowOffset);

            //turn windows into an ms style gene tree file
            //ProcessRAXML(baseFolder, windowSize, windowOffset);

            //run our method, calc score(s)


            //calc scores
            MakeFolder(baseFolder + "/results");
            //compareTwoMSFiles(baseFolder + "/ms/trueTrees", baseFolder + "/argw/finalLocalTrees",baseFolder+"/results/argw",genomeLength);
            //compareTwoMSFiles(baseFolder+"/ms/trueTrees",baseFolder+"/raxml/treefileRAXML",baseFolder+"/results/raxml",genomeLength);
            compareTwoMSFiles(baseFolder+"/ms/trueTrees",baseFolder+"/rent/finalResultsRENT",baseFolder+"/results/rent",genomeLength);

            int stopDebug = 0;

        }

    }
    */
    //endregion OLD CODE





    private static void RunMS(String whereToRun, String commandString){

        //add \ms folder to where we are running all this stuff
        whereToRun = whereToRun+"/ms";

        MakeFolder(whereToRun); //create ms output folder

        //String commandString = "ms 15 1 -t 15 -r 10 4000 -T";
        System.out.println(commandString);

        //String commandResults = RunCommandLineCommand(commandString);
        //System.out.println(commandResults);

        //write output to file
        try {
            TPrint("start|"+commandString); //timing info in case i need it later
            String commandResults = GetOutputFromProgramNoHangFromDir(commandString, null);
            TPrint("end|"+commandString); //timing info in case i need it later
            //System.out.println(commandResults);


            BufferedWriter out = new BufferedWriter(new FileWriter(whereToRun+"/fullOutput"));
            BufferedWriter cmd = new BufferedWriter(new FileWriter(whereToRun+"/commandUsed"));
            BufferedWriter trees = new BufferedWriter(new FileWriter(whereToRun+"/trueTrees"));
            BufferedWriter isMuts = new BufferedWriter(new FileWriter(whereToRun+"/infSitesMutations"));

            //save full output
            out.write(commandResults);

            //save command
            cmd.write(commandResults.split("\n")[0]);

            //save treefile
            trees.write(commandResults.split("//\n")[1].split("segsites")[0]);

            //save infinite sites SNPs?
            //gotta wrap in an if statement now cuz im not necessarily generating the inf sites mutations any more
            if(commandResults.contains("positions")) {
                isMuts.write(commandResults.split("positions: ")[1]);
            }

            //close the writers
            out.close();
            cmd.close();
            trees.close();
            isMuts.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

    }



    private static void RunRENT(String whereToRun){

        //add \ms folder to where we are running all this stuff
        String whereTheSimsAre = whereToRun+"/ms";
        String whereToRunRENT = whereToRun+"/rent";
        String whereTheSeqsAre = whereToRun+"/seq";

        MakeFolder(whereToRunRENT); //create ms output folder

        //sadly, was forced to handle multiple rent input files, so everything is in a big outter loop
        File dir = new File(whereTheSeqsAre);
        String[] extensions = new String[] { "rnt" };
        List<File> files = (List<File>) FileUtils.listFiles(dir, extensions, true);
        for (File file : files) {

            //command string, one below runs their method and one runs the method with their in house correctness checker
            //old testing commands
            //String commandString = "java -jar /Users/leo/rice/res/path/RentPlus.jar " + whereTheSimsAre + "/infSitesMutations";
            //String commandStringWithCheck = "java -jar /Users/leo/rice/res/path/RentPlus.jar " + whereTheSimsAre + "/infSitesMutations" + " " + whereTheSimsAre + "/fullOutput";



            //old way, can potentially hang cuz of weirdness
            //String commandResults = RunCommandLineCommand(commandString);
            //System.out.println(commandResults);
            //String commandResultsWithCheck = RunCommandLineCommand(commandStringWithCheck);

            //write output to file
            try {
                //new real commands
                String commandString = "java -jar /Users/leo/rice/res/path/RentPlus.jar " + file.getCanonicalPath();
                //String commandStringWithCheck = "java -jar /Users/leo/rice/res/path/RentPlus.jar " + whereTheSimsAre + "/infSitesMutations" + " " + whereTheSimsAre + "/fullOutput";

                TPrint("start|"+commandString); //timing info in case i need it later
                String commandResults = GetOutputFromProgramNoHangFromDir(commandString, null);
                TPrint("end|"+commandString); //timing info in case i need it later
                System.out.println(commandResults);

                //String commandResultsWithCheck = GetOutputFromProgramNoHangFromDir(commandStringWithCheck, null);


                BufferedWriter out = new BufferedWriter(new FileWriter(file.getCanonicalPath().replace(".rnt",".fullOutput")));
                BufferedWriter trees = new BufferedWriter(new FileWriter(file.getCanonicalPath().replace(".rnt",".rentTrees")));
                //BufferedWriter out = new BufferedWriter(new FileWriter(whereToRunRENT + "/fullOutput"));
                //BufferedWriter trees = new BufferedWriter(new FileWriter(whereToRunRENT + "/rentTrees"));
                //BufferedWriter outWithCheck = new BufferedWriter(new FileWriter(whereToRunRENT + "/fullOutputWithCheck"));


                //save full output
                out.write(commandResults);

                //save just trees
                String[] treeLines = commandResults.split("At");
                for (int i = 1; i < treeLines.length; i++) {
                    trees.write("At" + treeLines[i].split("\n")[0] + "\n");
                }

                //outWithCheck.write(commandResultsWithCheck);


                //close the writers
                out.close();
                trees.close();
                //outWithCheck.close();

            } catch (Exception e) {
                e.printStackTrace();
            }
        }

    }

    private static void RunRENTSingleFile(String whereToRun, String analyzeThisFilePath){

        String analyzeThisFileName = analyzeThisFilePath.split("/")[analyzeThisFilePath.split("/").length-1];


        //add \ms folder to where we are running all this stuff
        //String whereTheSimsAre = whereToRun+"/ms";
        //String whereToRunRENT = whereToRun+"/rent";
        String whereTheSeqsAre = whereToRun+"/rent/sites/"+analyzeThisFileName;

        //MakeFolder(whereToRun+"/rent"); //create ms output folder, should already be done and can comment out but w

        //sadly, was forced to handle multiple rent input files, so everything is in a big outter loop
        File dir = new File(whereTheSeqsAre);
        String[] extensions = new String[] { "rnt" };
        List<File> files = (List<File>) FileUtils.listFiles(dir, extensions, true);
        for (File file : files) {

            //write output to file
            try {
                //new real commands
                String commandString = "java -jar /Users/leo/rice/res/path/RentPlus.jar " + file.getCanonicalPath();
                //String commandStringWithCheck = "java -jar /Users/leo/rice/res/path/RentPlus.jar " + whereTheSimsAre + "/infSitesMutations" + " " + whereTheSimsAre + "/fullOutput";

                TPrint("start|"+commandString); //timing info in case i need it later
                String commandResults = GetOutputFromProgramNoHangFromDir(commandString, null);
                TPrint("end|"+commandString); //timing info in case i need it later
                System.out.println(commandResults);

                //String commandResultsWithCheck = GetOutputFromProgramNoHangFromDir(commandStringWithCheck, null);


                BufferedWriter out = new BufferedWriter(new FileWriter(file.getCanonicalPath().replace(".rnt",".fullOutput")));
                BufferedWriter trees = new BufferedWriter(new FileWriter(file.getCanonicalPath().replace(".rnt",".rentTrees")));
                //BufferedWriter out = new BufferedWriter(new FileWriter(whereToRunRENT + "/fullOutput"));
                //BufferedWriter trees = new BufferedWriter(new FileWriter(whereToRunRENT + "/rentTrees"));
                //BufferedWriter outWithCheck = new BufferedWriter(new FileWriter(whereToRunRENT + "/fullOutputWithCheck"));


                //save full output
                out.write(commandResults);

                //save just trees
                String[] treeLines = commandResults.split("At");
                for (int i = 1; i < treeLines.length; i++) {
                    trees.write("At" + treeLines[i].split("\n")[0] + "\n");
                }

                //outWithCheck.write(commandResultsWithCheck);


                //close the writers
                out.close();
                trees.close();
                //outWithCheck.close();

            } catch (Exception e) {
                e.printStackTrace();
            }
        }

    }


    private static void ProcessRENT(String whereToRun,int genomeLength){

        //add \ms folder to where we are running all this stuff
        String whereToRunRENT = whereToRun+"/rent";
        String whereTheSeqsAre = whereToRun+"/seq";


        //SNAGGED FROM LEO SCRIPT TEMP (rentplusprep3) - WITH A MAJOR FIX FOR HOW MANY BP TO BE WRITTEN IN THE MS STYLE REPRESENATION

        try {

            File f = new File(whereTheSeqsAre);
            File[] matchingFiles = f.listFiles(new FilenameFilter() {
                public boolean accept(File dir, String name) {
                    return name.endsWith(".rnt.trees");
                }
            });


            int howManyBPSoFar = 0;
            StringBuilder theStuffToWrite = new StringBuilder();

            for (int concatI = 0; concatI < matchingFiles.length; concatI++) {

                //first line treated slightly differently (before and after that point are considered that tree, others are just that point and beyond
                boolean firstLine = true;

                BufferedReader br = new BufferedReader(new FileReader(whereTheSeqsAre + "/sitesRENT" + String.valueOf(concatI) + ".rnt.trees"));

                String line;

                int nextPosition = 0;
                int currentPosition = 0;
                String nextTree = "";
                String currentTree;

                while ((line = br.readLine()) != null) {

                    String[] theLineArray = line.split("\t");
                    int whereAmI = Integer.parseInt(theLineArray[0]);

                    //String whereAmI = theLineArray[0].replace("[","").replace("]","");
                    String theTreeString = theLineArray[1];

                    //ordering of names and everything is all messed up, have to add a bunch of dumb stuff to fix it like with argweaver
                    //region ---- fix names section
                    //step1 - get name correspondences
                    BufferedReader namesRENT = new BufferedReader(new FileReader(whereTheSeqsAre + "/namesRENT"));
                    HashMap<String,String> namesMap = new HashMap<String,String>();
                    String nameLine;
                    int startName = 1;
                    while ((nameLine = namesRENT.readLine()) != null) {
                        namesMap.put(String.valueOf(startName),nameLine);
                        startName++;
                    }
                    //step2 - edit nodes to attach correspondences (avoids duplicate node names)
                    Tree treeRename = Trees.readTree(theTreeString);

                    //add temp name
                    for(TNode renameMe : treeRename.getNodes()){

                        //only messing with leaves
                        if(renameMe.getChildCount()==0) {
                            //rename each node based on name mapping in smc file
                            String currentNodeName = renameMe.getName();
                            String newNodeName = "temp_"+namesMap.get(currentNodeName);
                            ((STINode) renameMe).setName(newNodeName);
                            //System.out.println("renamed "+currentNodeName+" to "+newNodeName);
                        }else{
                            //annoying thing necessary for argweaver, doesnt matter here
                            //((STINode) renameMe).setName("i_"+renameMe.getName());
                        }
                    }

                    //then go back and remove temp_> from names since i had to do that to avoid duplicate naming issues
                    for(TNode renameMe : treeRename.getNodes()){
                        //only messing with leaves
                        if(renameMe.getChildCount()==0) {
                            //rename each node based on name mapping in smc file
                            String currentNodeName = renameMe.getName();
                            String newNodeName = currentNodeName.substring(6);
                            ((STINode) renameMe).setName(newNodeName);
                        }
                    }

                    theTreeString = treeRename.toNewick();
                    //endregion ---- end name fixing

                    //
                    //int msWhereAmI = whereAmI - howManyBPSoFar;
                    //howManyBPSoFar += msWhereAmI;

                    //looking one ahead, with special case for first line
                    if(firstLine) {
                        //currentPosition = nextPosition;
                        currentTree = theTreeString;
                        nextPosition = whereAmI;
                        nextTree = theTreeString;
                        firstLine = false;
                    }else {
                        currentPosition = nextPosition;
                        currentTree = nextTree;
                        nextPosition = whereAmI;
                        nextTree = theTreeString;
                    }


                    //System.out.println("["+msWhereAmI+"]"+theTree+";");
                    //theStuffToWrite.append("[" + msWhereAmI + "]" + theTreeString + "\n");
                    theStuffToWrite.append("[" + (nextPosition-currentPosition) + "]" + currentTree + "\n");



                }

                //after a null was hit, we are at the end, put one last line:
                theStuffToWrite.append("[" + (genomeLength-nextPosition) + "]" + nextTree + "\n");


                br.close();

            }

            //write the stuff out
            PrintWriter writer = new PrintWriter(whereToRunRENT + "/finalResultsRENT", "UTF-8");
            writer.write(theStuffToWrite.toString());
            writer.close();

            //int deb = 0;

        }catch(Exception e){
            e.printStackTrace();
        }


    }


    private static void ProcessRENTSingleFile(String whereToRun,int genomeLength,String analyzeThisFilePath){

        String analyzeThisFileName = analyzeThisFilePath.split("/")[analyzeThisFilePath.split("/").length-1];

        //add \ms folder to where we are running all this stuff
        //String whereToRunRENT = whereToRun+"/rent";
        String whereTheSeqsAre = whereToRun+"/rent/sites/"+analyzeThisFileName;


        //SNAGGED FROM LEO SCRIPT TEMP (rentplusprep3) - WITH A MAJOR FIX FOR HOW MANY BP TO BE WRITTEN IN THE MS STYLE REPRESENATION

        try {

            File f = new File(whereTheSeqsAre);
            File[] matchingFiles = f.listFiles(new FilenameFilter() {
                public boolean accept(File dir, String name) {
                    return name.endsWith(".rnt.trees");
                }
            });


            int howManyBPSoFar = 0;
            StringBuilder theStuffToWrite = new StringBuilder();

            for (int concatI = 0; concatI < matchingFiles.length; concatI++) {

                //first line treated slightly differently (before and after that point are considered that tree, others are just that point and beyond
                boolean firstLine = true;

                BufferedReader br = new BufferedReader(new FileReader(whereTheSeqsAre + "/" + String.valueOf(concatI) + ".rnt.trees"));

                String line;

                int nextPosition = 0;
                int currentPosition = 0;
                String nextTree = "";
                String currentTree;

                while ((line = br.readLine()) != null) {

                    String[] theLineArray = line.split("\t");
                    int whereAmI = Integer.parseInt(theLineArray[0]);

                    //String whereAmI = theLineArray[0].replace("[","").replace("]","");
                    String theTreeString = theLineArray[1];

                    //ordering of names and everything is all messed up, have to add a bunch of dumb stuff to fix it like with argweaver
                    //region ---- fix names section
                    //step1 - get name correspondences
                    BufferedReader namesRENT = new BufferedReader(new FileReader(whereToRun + "/rent/names/" + analyzeThisFileName));
                    HashMap<String,String> namesMap = new HashMap<String,String>();
                    String nameLine;
                    int startName = 1;
                    while ((nameLine = namesRENT.readLine()) != null) {
                        namesMap.put(String.valueOf(startName),nameLine);
                        startName++;
                    }
                    //step2 - edit nodes to attach correspondences (avoids duplicate node names)
                    Tree treeRename = Trees.readTree(theTreeString);

                    //add temp name
                    for(TNode renameMe : treeRename.getNodes()){

                        //only messing with leaves
                        if(renameMe.getChildCount()==0) {
                            //rename each node based on name mapping in smc file
                            String currentNodeName = renameMe.getName();
                            String newNodeName = "temp_"+namesMap.get(currentNodeName);
                            ((STINode) renameMe).setName(newNodeName);
                            //System.out.println("renamed "+currentNodeName+" to "+newNodeName);
                        }else{
                            //annoying thing necessary for argweaver, doesnt matter here
                            //((STINode) renameMe).setName("i_"+renameMe.getName());
                        }
                    }

                    //then go back and remove temp_> from names since i had to do that to avoid duplicate naming issues
                    for(TNode renameMe : treeRename.getNodes()){
                        //only messing with leaves
                        if(renameMe.getChildCount()==0) {
                            //rename each node based on name mapping in smc file
                            String currentNodeName = renameMe.getName();
                            String newNodeName = currentNodeName.substring(6);
                            ((STINode) renameMe).setName(newNodeName);
                        }
                    }

                    theTreeString = treeRename.toNewick();
                    //endregion ---- end name fixing

                    //
                    //int msWhereAmI = whereAmI - howManyBPSoFar;
                    //howManyBPSoFar += msWhereAmI;

                    //looking one ahead, with special case for first line
                    if(firstLine) {
                        //currentPosition = nextPosition;
                        currentTree = theTreeString;
                        nextPosition = whereAmI;
                        nextTree = theTreeString;
                        firstLine = false;
                    }else {
                        currentPosition = nextPosition;
                        currentTree = nextTree;
                        nextPosition = whereAmI;
                        nextTree = theTreeString;
                    }


                    //System.out.println("["+msWhereAmI+"]"+theTree+";");
                    //theStuffToWrite.append("[" + msWhereAmI + "]" + theTreeString + "\n");
                    theStuffToWrite.append("[" + (nextPosition-currentPosition) + "]" + currentTree + "\n");



                }

                //after a null was hit, we are at the end, put one last line:
                theStuffToWrite.append("[" + (genomeLength-nextPosition) + "]" + nextTree + "\n");


                br.close();

            }

            //write the stuff out
            MakeFolder(whereToRun+"/rent/finished");
            PrintWriter writer = new PrintWriter(whereToRun + "/rent/finished/" + analyzeThisFileName+".finalTrees", "UTF-8");
            writer.write(theStuffToWrite.toString());
            writer.close();

            //int deb = 0;

        }catch(Exception e){
            e.printStackTrace();
        }


    }



    private static void RunSeqGen(String whereToRun, String commandString){

        //add \ms folder to where we are running all this stuff
        String whereTheSimsAre = whereToRun+"/ms";
        String whereToRunSeqGen = whereToRun+"/seq";

        MakeFolder(whereToRunSeqGen); //create ms output folder

        //using seq gen itself to do fasta AND phyli[, requiring me to run it in fixed seed across two runs, so generate a seed:
        int seed = Random.nextInt(Integer.MAX_VALUE);

        //command string, double check this later, also apparently fasta output is possible now too but cant find the command
        commandString = commandString +" <"+whereTheSimsAre+"/trueTrees";
        //DOWNLOAD NEW VERS FOR FASTA OUTPUT
        //String commandStringFasta = "seq-gen -of -z "+seed+ " -mHKY -l 4000 -s .00544 -p 4000 <"+whereTheSimsAre+"/trueTrees";
        String commandStringFasta = commandString.replace("seq-gen","seq-gen -of");

        System.out.println(commandString);

        //FOR SOME REASON THE UPDATED VERSION OF RUN COMMAND WONT WORK HERE SINCE THERE IS NOT A PROPER WAIT HAPPENING
        //SO THAT THE FILES ARE READY (cannot use 'GetOutputFromProgramNoHangFromDir')
        //String commandResults = RunCommandLineCommand(commandString);
        //System.out.println(commandResults);
        //String commandResultsFasta = RunCommandLineCommand(commandStringFasta);

        //write output to file
        try {
            //swapped above lines to these so it doesnt randomly hang
            String commandResults = GetOutputFromProgramNoHangFromDir(commandString, null);
            System.out.println(commandResults);
            String commandResultsFasta = GetOutputFromProgramNoHangFromDir(commandStringFasta, null);

            //remove top part from the two outputs so im only saving the fasta and phylip parts
            StringBuilder phylipFile = new StringBuilder();
            StringBuilder fastaFile = new StringBuilder();
            //do for phylip
            boolean pastIntro = false;
            for(String line : commandResults.split("\n")){
                if(pastIntro){
                    phylipFile.append(line+"\n");
                }
                if(line.contains("Time taken")){
                    pastIntro = true;
                }
            }
            //do for fasta
            pastIntro = false;
            for(String line : commandResultsFasta.split("\n")){
                if(pastIntro){
                    fastaFile.append(line+"\n");
                }
                if(line.contains("Time taken")){
                    pastIntro = true;
                }
            }

            //write out everything

            BufferedWriter out = new BufferedWriter(new FileWriter(whereToRunSeqGen+"/phylip"));
            BufferedWriter out2 = new BufferedWriter(new FileWriter(whereToRunSeqGen+"/fasta"));
            BufferedWriter cmd = new BufferedWriter(new FileWriter(whereToRunSeqGen+"/commands"));


            //save full output
            out.write(phylipFile.toString());
            out2.write(fastaFile.toString());

            //save commands for reference
            cmd.write(commandString +"\n"+commandStringFasta);


            //close the writers
            out.close();
            out2.close();
            cmd.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

    }


    private static void FormatToRENT(String whereToRun){

        //add \ms folder to where we are running all this stuff
        String whereSeqGenIs = whereToRun+"/seq";



        //input fasta and write output to rent sites file
        try {

            //COPIED AND PASTED THIS STUFF FROM MY LEOSCRIPT TEMP

            //BufferedReader br = new BufferedReader(new FileReader("/Users/leo/rice/res/tmp/sim/seqfile4800"));

            //BufferedReader br = new BufferedReader(new FileReader("/Users/leo/rice/res/code/scripts/phylonet-hmm/results2/simWF1200/rentPlus/seqfileWF1200Formatted2"));
            BufferedReader br = new BufferedReader(new FileReader(whereSeqGenIs+"/fasta"));

            StringBuilder firstLineWithPositions = new StringBuilder();


            String line;
            String name;
            StringBuilder seq;
            int seqLength = 0;
            String nextName = null;
            HashMap<String,String> theSeqs = new HashMap<String,String>();
            while ((line = br.readLine()) != null) {

                if(nextName == null) {
                    name = new String(line);
                    seq = new StringBuilder();
                }else{
                    name = nextName;
                    seq = new StringBuilder();
                    seq.append(line);
                }

                //read seq

                while ((line = br.readLine()) != null && line.charAt(0) != '>') {
                    seq.append(line);
                }
                nextName = line;
                theSeqs.put(name,seq.toString());
                seqLength = seq.length(); //(should all be the same anyways)
                //double checking they are all the same
                //System.out.println(seq.length());


            }


            //data
            //seguqnce lengths is stored in seqLength
            int numberOfSequences = theSeqs.size();
            HashMap<String,StringBuilder> theNewSeqs = new HashMap<String,StringBuilder>();

            int stop=0;


            //write out the region where this window was taken from just for convenience
            //PrintWriter writer2 = new PrintWriter(fourthArg + "/w" + windowName + ".region", "UTF-8");


            for(int currentPos = 0; currentPos < seqLength; currentPos++) {

                //set to see how many unique letters and see if its a biallelic site
                HashSet<Character> howManyUnique = new HashSet<Character>();

                for(HashMap.Entry<String, String> entry : theSeqs.entrySet()) {
                    String key = entry.getKey();
                    String value = entry.getValue();

                    howManyUnique.add(value.charAt(currentPos));

                }

                if(howManyUnique.size() == 2){
                    //do stuff
                    //System.out.println("biallelic SNP at site" + currentPos);
                    firstLineWithPositions.append(currentPos+" ");
                    //get orangutan letter, that will be the 0, others will be 1
                    //String debug1 = theSeqs.get(">OrOut ");
                    //char theZeroLetter = theSeqs.get(">OrOut ").charAt(currentPos);
                    //not doing this any more, just arbitrarily assign, take first letter as 'theZeroLetter'
                    char theZeroLetter = howManyUnique.iterator().next();

                    //iterate through the seqs, if the letter is theZeroLetter, output a 0, else output a 1
                    //the outgroup should always just be 0s (not sure if this matters but doing it this way for now) - update, still dont htink it matters and doing it arbitrarily
                    for(HashMap.Entry<String, String> entry : theSeqs.entrySet()) {
                        String key = entry.getKey();
                        String value = entry.getValue();
                        char thisLetter = value.charAt(currentPos);

                        if(thisLetter==theZeroLetter){

                            if(theNewSeqs.get(key) == null){
                                //theNewSeqs.get(key). = new StringBuilder();
                                theNewSeqs.put(key,(new StringBuilder()));
                            }

                            theNewSeqs.get(key).append("0");
                        }else{
                            if(theNewSeqs.get(key) == null){
                                //theNewSeqs.get(key). = new StringBuilder();
                                theNewSeqs.put(key,(new StringBuilder()));
                            }
                            theNewSeqs.get(key).append("1");
                        }

                        //howManyUnique.add(value.charAt(currentPos));

                    }







                }else{
                    //dont do stuff

                }


            }



            int howManyFiles = 0;
            int howManyToPrintAtATime = 50000;
            String[] firstLineWithPositionsArray = firstLineWithPositions.toString().split(" ");
            int howManySNPs = firstLineWithPositionsArray.length;
            int howManyHaveBeenPrinted = 0;

            //print out names in order since i will need later to turn back to MS format tree file
            PrintWriter writer2 = new PrintWriter(whereSeqGenIs+"/namesRENT", "UTF-8");
            for (HashMap.Entry<String, StringBuilder> entry : theNewSeqs.entrySet()) {
                String key = entry.getKey();
                //StringBuilder value = entry.getValue();
                writer2.println(key);
            }
            writer2.close();


            while(howManyHaveBeenPrinted<howManySNPs) {

                int howManyToPrint = Math.min(howManyToPrintAtATime,(howManySNPs-howManyHaveBeenPrinted));

                PrintWriter writer = new PrintWriter(whereSeqGenIs+"/sitesRENT"+String.valueOf(howManyFiles)+".rnt", "UTF-8");
                //names in order of their seqs in the sites file so i can make sure the mapping is done properly later
                howManyFiles++;


                for (int firstLineI = howManyHaveBeenPrinted; firstLineI < (howManyHaveBeenPrinted+howManyToPrint); firstLineI++) {
                    writer.print(firstLineWithPositionsArray[firstLineI] + " ");
                }
                writer.print("\n");
                //writer.println(firstLineWithPositions);
                for (HashMap.Entry<String, StringBuilder> entry : theNewSeqs.entrySet()) {
                    String key = entry.getKey();
                    StringBuilder value = entry.getValue();
                    writer.println(value.toString().substring(howManyHaveBeenPrinted, (howManyHaveBeenPrinted+howManyToPrint)));
                }


                //if i later want to make sure i know which was which i can do it, just make an array or whatever and output in another file which line was which instead of the above stuff
                howManyHaveBeenPrinted += howManyToPrint;

                writer.close();
            }



        } catch (Exception e) {
            e.printStackTrace();
        }

    }


    private static void FormatToRENTSingleFile(String whereToRun, String analyzeThisFilePath){

        String analyzeThisFileName = analyzeThisFilePath.split("/")[analyzeThisFilePath.split("/").length-1];

        //input fasta and write output to rent sites file
        try {

            //COPIED AND PASTED THIS STUFF FROM MY LEOSCRIPT TEMP

            //BufferedReader br = new BufferedReader(new FileReader("/Users/leo/rice/res/tmp/sim/seqfile4800"));

            //BufferedReader br = new BufferedReader(new FileReader("/Users/leo/rice/res/code/scripts/phylonet-hmm/results2/simWF1200/rentPlus/seqfileWF1200Formatted2"));
            BufferedReader br = new BufferedReader(new FileReader(analyzeThisFilePath));

            StringBuilder firstLineWithPositions = new StringBuilder();


            String line;
            String name;
            StringBuilder seq;
            int seqLength = 0;
            String nextName = null;
            HashMap<String,String> theSeqs = new HashMap<String,String>();
            while ((line = br.readLine()) != null) {

                if(nextName == null) {
                    name = new String(line);
                    seq = new StringBuilder();
                }else{
                    name = nextName;
                    seq = new StringBuilder();
                    seq.append(line);
                }

                //read seq

                while ((line = br.readLine()) != null && line.charAt(0) != '>') {
                    seq.append(line);
                }
                nextName = line;
                theSeqs.put(name,seq.toString());
                seqLength = seq.length(); //(should all be the same anyways)
                //double checking they are all the same
                //System.out.println(seq.length());


            }


            //data
            //seguqnce lengths is stored in seqLength
            int numberOfSequences = theSeqs.size();
            HashMap<String,StringBuilder> theNewSeqs = new HashMap<String,StringBuilder>();

            int stop=0;


            //write out the region where this window was taken from just for convenience
            //PrintWriter writer2 = new PrintWriter(fourthArg + "/w" + windowName + ".region", "UTF-8");


            for(int currentPos = 0; currentPos < seqLength; currentPos++) {

                //set to see how many unique letters and see if its a biallelic site
                HashSet<Character> howManyUnique = new HashSet<Character>();

                for(HashMap.Entry<String, String> entry : theSeqs.entrySet()) {
                    String key = entry.getKey();
                    String value = entry.getValue();

                    howManyUnique.add(value.charAt(currentPos));

                }

                if(howManyUnique.size() == 2){
                    //do stuff
                    //System.out.println("biallelic SNP at site" + currentPos);
                    firstLineWithPositions.append(currentPos+" ");
                    //get orangutan letter, that will be the 0, others will be 1
                    //String debug1 = theSeqs.get(">OrOut ");
                    //char theZeroLetter = theSeqs.get(">OrOut ").charAt(currentPos);
                    //not doing this any more, just arbitrarily assign, take first letter as 'theZeroLetter'
                    char theZeroLetter = howManyUnique.iterator().next();

                    //iterate through the seqs, if the letter is theZeroLetter, output a 0, else output a 1
                    //the outgroup should always just be 0s (not sure if this matters but doing it this way for now) - update, still dont htink it matters and doing it arbitrarily
                    for(HashMap.Entry<String, String> entry : theSeqs.entrySet()) {
                        String key = entry.getKey();
                        String value = entry.getValue();
                        char thisLetter = value.charAt(currentPos);

                        if(thisLetter==theZeroLetter){

                            if(theNewSeqs.get(key) == null){
                                //theNewSeqs.get(key). = new StringBuilder();
                                theNewSeqs.put(key,(new StringBuilder()));
                            }

                            theNewSeqs.get(key).append("0");
                        }else{
                            if(theNewSeqs.get(key) == null){
                                //theNewSeqs.get(key). = new StringBuilder();
                                theNewSeqs.put(key,(new StringBuilder()));
                            }
                            theNewSeqs.get(key).append("1");
                        }

                        //howManyUnique.add(value.charAt(currentPos));

                    }







                }else{
                    //dont do stuff

                }


            }



            int howManyFiles = 0;
            int howManyToPrintAtATime = 50000;
            String[] firstLineWithPositionsArray = firstLineWithPositions.toString().split(" ");
            int howManySNPs = firstLineWithPositionsArray.length;
            int howManyHaveBeenPrinted = 0;

            //print out names in order since i will need later to turn back to MS format tree file
            MakeFolder(whereToRun+"/rent/names");
            PrintWriter writer2 = new PrintWriter(whereToRun+"/rent/names/"+analyzeThisFileName, "UTF-8");
            for (HashMap.Entry<String, StringBuilder> entry : theNewSeqs.entrySet()) {
                String key = entry.getKey();
                //StringBuilder value = entry.getValue();
                writer2.println(key);
            }
            writer2.close();


            while(howManyHaveBeenPrinted<howManySNPs) {

                int howManyToPrint = Math.min(howManyToPrintAtATime,(howManySNPs-howManyHaveBeenPrinted));

                MakeFolder(whereToRun+"/rent/sites");
                MakeFolder(whereToRun+"/rent/sites/"+analyzeThisFileName);
                PrintWriter writer = new PrintWriter(whereToRun+"/rent/sites/"+analyzeThisFileName+"/"+String.valueOf(howManyFiles)+".rnt", "UTF-8");
                //names in order of their seqs in the sites file so i can make sure the mapping is done properly later
                howManyFiles++;


                for (int firstLineI = howManyHaveBeenPrinted; firstLineI < (howManyHaveBeenPrinted+howManyToPrint); firstLineI++) {
                    writer.print(firstLineWithPositionsArray[firstLineI] + " ");
                }
                writer.print("\n");
                //writer.println(firstLineWithPositions);
                for (HashMap.Entry<String, StringBuilder> entry : theNewSeqs.entrySet()) {
                    String key = entry.getKey();
                    StringBuilder value = entry.getValue();
                    writer.println(value.toString().substring(howManyHaveBeenPrinted, (howManyHaveBeenPrinted+howManyToPrint)));
                }


                //if i later want to make sure i know which was which i can do it, just make an array or whatever and output in another file which line was which instead of the above stuff
                howManyHaveBeenPrinted += howManyToPrint;

                writer.close();
            }



        } catch (Exception e) {
            e.printStackTrace();
        }

    }





    //this and the useful delete function from the nice example at http://roufid.com/how-to-delete-folder-recursively-in-java/#native_java
    private static void DeleteFolder(String folderToDelete) {

        File file = new File(folderToDelete);

        if(!file.exists()){
            return;
        }

        try {
            //Deleting the directory recursively.
            delete(file);
            //System.out.println("Directory has been deleted recursively !");
        } catch (IOException e) {
            System.out.println("Problem occurs when deleting the directory : " + folderToDelete);
            e.printStackTrace();
        }

    }

    private static void delete(File file) throws IOException {

        for (File childFile : file.listFiles()) {

            if (childFile.isDirectory()) {
                delete(childFile);
            } else {
                if (!childFile.delete()) {
                    throw new IOException();
                }
            }
        }

        if (!file.delete()) {
            throw new IOException();
        }
    }



    private static void RunARGW(String whereToRun){

        //add \ms folder to where we are running all this stuff
        String whereTheSimsAre = whereToRun+"/seq";
        String whereToRunARGW = whereToRun+"/argw";

        DeleteFolder(whereToRunARGW); //for argw i gotta remove everything if its gonna run again
        MakeFolder(whereToRunARGW); //create ms output folder



        //command string, one below runs their method and one runs the method with their in house correctness checker
        //String[] commandStrings = {"cd "+whereToRunARGW,"/Users/leo/rice/res/tools/mdrasmus-argweaver-fee3d32/bin/arg-sample -f "+whereTheSimsAre+"/fasta"};
        //String[] commandStrings = {"cd "+whereToRunARGW,"touch test"};
        //String commandString = "touch test";
        String commandString = "/Users/leo/rice/res/tools/mdrasmus-argweaver-fee3d32/bin/arg-sample -f "+whereTheSimsAre+"/fasta";
        System.out.println(commandString);


        //write output to file
        try {

            //for argweaver we have to erase all the files to get it to run again
            FileUtils.cleanDirectory(new File(whereToRunARGW));

            //now run the actual command
            //String commandResults = RunCommandLineCommands(commandStrings);
            //String commandResults = RunCommandLineCommandFromDir(commandString, whereToRunARGW);
            TPrint("start|"+commandString); //timing info in case i need it later
            String commandResults = GetOutputFromProgramNoHangFromDir(commandString, whereToRunARGW);
            TPrint("end|"+commandString); //timing info in case i need it later
            System.out.println(commandResults);

            BufferedWriter out = new BufferedWriter(new FileWriter(whereToRunARGW+"/fullOutput"));
            BufferedWriter cmd = new BufferedWriter(new FileWriter(whereToRunARGW+"/command"));


            //save full output
            out.write(commandResults);

            //save commands ran (lambda makes it concise...minus forcing an extra try catch block >.<
            /*
            Arrays.stream(commandStrings)
                    .forEach(e -> {
                        try {
                            cmd.write(e + "\n");
                        } catch (IOException e1) {
                            e1.printStackTrace();
                        }
                    });
            */

            //close the writers
            out.close();
            cmd.close();


        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private static void RunARGWSingleFile(String whereToRun, String analyzeThisFilePath, boolean replaceDashesWithNs){

        //id replaceNs = true, swaps dashes for Ns (recommended by maintainer, dashes not supported)
        String optionalReplaceDashString = "";

        if(replaceDashesWithNs){
            try {
                optionalReplaceDashString = ".NoDash";
                BufferedWriter outN = new BufferedWriter(new FileWriter(analyzeThisFilePath+".NoDash"));
                BufferedReader brN = new BufferedReader(new FileReader(analyzeThisFilePath));

                String line;
                while ((line = brN.readLine()) != null) {
                    outN.write(line.replace("-","N")+"\n");
                }

                outN.close();
                brN.close();
            }catch(Exception e){
                e.printStackTrace();
            }
        }

        String analyzeThisFileName = analyzeThisFilePath.split("/")[analyzeThisFilePath.split("/").length-1];

        //add \ms folder to where we are running all this stuff
        //String whereTheSimsAre = whereToRun+"/seq";
        //String whereToRunARGW = whereToRun+"/argw";

        //DeleteFolder(whereToRunARGW); //for argw i gotta remove everything if its gonna run again
        MakeFolder(whereToRun+"/argw"); //create ms output folder

        DeleteFolder(whereToRun+"/argw/"+analyzeThisFileName); //for argw i gotta remove everything if its gonna run again
        MakeFolder(whereToRun+"/argw/"+analyzeThisFileName); //create ms output folder

        String whereToRunARGW = whereToRun+"/argw/"+analyzeThisFileName;


        //command string, one below runs their method and one runs the method with their in house correctness checker
        //String[] commandStrings = {"cd "+whereToRunARGW,"/Users/leo/rice/res/tools/mdrasmus-argweaver-fee3d32/bin/arg-sample -f "+whereTheSimsAre+"/fasta"};
        //String[] commandStrings = {"cd "+whereToRunARGW,"touch test"};
        //String commandString = "touch test";
        String commandString = "/Users/leo/rice/res/tools/mdrasmus-argweaver-fee3d32/bin/arg-sample -f "+analyzeThisFilePath+optionalReplaceDashString;
        System.out.println(commandString);


        //write output to file
        try {

            //for argweaver we have to erase all the files to get it to run again
            FileUtils.cleanDirectory(new File(whereToRunARGW));

            //now run the actual command
            //String commandResults = RunCommandLineCommands(commandStrings);
            //String commandResults = RunCommandLineCommandFromDir(commandString, whereToRunARGW);
            TPrint("start|"+commandString); //timing info in case i need it later
            String commandResults = GetOutputFromProgramNoHangFromDir(commandString, whereToRunARGW);
            TPrint("end|"+commandString); //timing info in case i need it later
            System.out.println(commandResults);

            BufferedWriter out = new BufferedWriter(new FileWriter(whereToRunARGW+"/fullOutput"));
            BufferedWriter cmd = new BufferedWriter(new FileWriter(whereToRunARGW+"/command"));


            //save full output
            out.write(commandResults);

            //save commands ran (lambda makes it concise...minus forcing an extra try catch block >.<
            /*
            Arrays.stream(commandStrings)
                    .forEach(e -> {
                        try {
                            cmd.write(e + "\n");
                        } catch (IOException e1) {
                            e1.printStackTrace();
                        }
                    });
            */

            //close the writers
            out.close();
            cmd.close();

            //after done running, if had to do dash swap, can delete .NoDash file
            if(replaceDashesWithNs){
                (new File(analyzeThisFilePath+optionalReplaceDashString)).delete();
            }

        } catch (Exception e) {
            e.printStackTrace();
        }

    }



    private static void ProcessARGW(String whereToRun){


        String whereToRunARGW = whereToRun+"/argw";

        String highestLikelihoodARG = getHighestLikelihoodArg(whereToRunARGW);

        System.out.println(highestLikelihoodARG);

        //unzip the arg with the highest likelihood
        String commandResults = RunCommandLineCommand("gzip -d " + highestLikelihoodARG);

        //turn the ARG into an ms style tree file
        argweaverToMS(highestLikelihoodARG,whereToRunARGW+"/finalLocalTrees");




        //command string, one below runs their method and one runs the method with their in house correctness checker
        //String[] commandStrings = {"cd "+whereToRunARGW,"/Users/leo/rice/res/tools/mdrasmus-argweaver-fee3d32/bin/arg-sample -f "+whereTheSimsAre+"/fasta"};
        //String[] commandStrings = {"cd "+whereToRunARGW,"touch test"};
        //String commandString = "touch test";
        //String commandString = "/Users/leo/rice/res/tools/mdrasmus-argweaver-fee3d32/bin/arg-sample -f "+whereTheSimsAre+"/fasta";
        //System.out.println(commandString);
        //now run the actual command
        //String commandResults = RunCommandLineCommands(commandStrings);
        //String commandResults = RunCommandLineCommandFromDir(commandString, whereToRunARGW);
        //String commandResults = GetOutputFromProgramNoHangFromDir(commandString, whereToRunARGW);
        //System.out.println(commandResults);

    }


    private static void ProcessARGWSingleFile(String whereToRun, String analyzeThisFilePath){

        String analyzeThisFileName = analyzeThisFilePath.split("/")[analyzeThisFilePath.split("/").length-1];


        String whereToRunARGW = whereToRun+"/argw/"+analyzeThisFileName;

        String highestLikelihoodARG = getHighestLikelihoodArg(whereToRunARGW);

        System.out.println(highestLikelihoodARG);

        //unzip the arg with the highest likelihood
        String commandResults = RunCommandLineCommand("gzip -d " + highestLikelihoodARG);

        //turn the ARG into an ms style tree file
        argweaverToMS(highestLikelihoodARG,whereToRunARGW+"/finalLocalTrees");

    }




    public static void Windower(String whereToRun, int windowSize, int windowOffset) {

        //copied and pasted basically from leo script temp, operates on fasta files

        String whereTheSimsAre = whereToRun+"/seq";
        String whereToMakeWindows = whereToRun+"/windows"+windowSize+"bp"+windowOffset+"os";

        MakeFolder(whereToMakeWindows); //create ms output folder

        try {
            //BufferedReader br = new BufferedReader(new FileReader("/Users/leo/rice/res/tmp/sim/seqfile4800"));
            BufferedReader br = new BufferedReader(new FileReader(whereTheSimsAre+"/fasta"));

            String line;
            String name;
            StringBuilder seq;
            int seqLength = 0;
            String nextName = null;
            HashMap<String,String> theSeqs = new HashMap<String,String>();
            while ((line = br.readLine()) != null) {

                if(nextName == null) {
                    name = new String(line);
                    seq = new StringBuilder();
                }else{
                    name = nextName;
                    seq = new StringBuilder();
                    seq.append(line);
                }

                //read seq

                while ((line = br.readLine()) != null && line.charAt(0) != '>') {
                    seq.append(line);
                }
                nextName = line;
                //seq = br.readLine();
                theSeqs.put(name,seq.toString());
                seqLength = seq.length(); //(should all be the same anyways)
                //double checking they are all the same
                //System.out.println(seq.length());

            }


            int windowName = 0;
            int basePairsDone = 0;
            for(int currentPos = 0; basePairsDone < seqLength; currentPos += windowOffset) {

                //handles last window if theres overrun
                if(basePairsDone < seqLength && currentPos+windowSize >= seqLength){
                    windowSize = seqLength-currentPos;
                }

                PrintWriter writer = new PrintWriter(whereToMakeWindows + "/w" + windowName + ".fasta", "UTF-8");
                //write out the region where this window was taken from just for convenience
                PrintWriter writer2 = new PrintWriter(whereToMakeWindows + "/w" + windowName + ".region", "UTF-8");

                for(HashMap.Entry<String, String> entry : theSeqs.entrySet()) {
                    String key = entry.getKey();
                    String value = entry.getValue();

                    writer.println(key);
                    if(currentPos+windowSize < seqLength) {
                        //print in lengths of 60...need to change eventually probably
                        writer.println(value.substring(currentPos, currentPos + windowSize));

                    }else{
                        //print in lengths of 60...need to change eventually probably
                        writer.println(value.substring(currentPos, seqLength));
                    }
                    writer.println();

                }


                //writer.println("blah");
                writer2.println(currentPos + " to " + (currentPos+windowSize));

                writer.close();
                writer2.close();

                windowName++;

                //tick up basePairsDone
                basePairsDone = currentPos+windowSize;
            }







        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e){
            e.printStackTrace();
        }



    }


    public static void WindowerWithSeqPath(String whereToRun, String sequencePath, int windowSize, int windowOffset) {

        //copied and pasted basically from leo script temp, operates on fasta files

        //String whereTheSimsAre = whereToRun+"/seq";
        String whereToMakeWindows = whereToRun+"/windows"+windowSize+"bp"+windowOffset+"os";

        MakeFolder(whereToMakeWindows); //create ms output folder

        try {
            //BufferedReader br = new BufferedReader(new FileReader("/Users/leo/rice/res/tmp/sim/seqfile4800"));
            BufferedReader br = new BufferedReader(new FileReader(sequencePath));

            String line;
            String name;
            StringBuilder seq;
            int seqLength = 0;
            String nextName = null;
            HashMap<String,String> theSeqs = new HashMap<String,String>();
            while ((line = br.readLine()) != null) {

                if(nextName == null) {
                    name = new String(line);
                    seq = new StringBuilder();
                }else{
                    name = nextName;
                    seq = new StringBuilder();
                    seq.append(line);
                }

                //read seq

                while ((line = br.readLine()) != null && line.charAt(0) != '>') {
                    seq.append(line);
                }
                nextName = line;
                //seq = br.readLine();
                theSeqs.put(name,seq.toString());
                seqLength = seq.length(); //(should all be the same anyways)
                //double checking they are all the same
                //System.out.println(seq.length());

            }


            int windowName = 0;
            int basePairsDone = 0;
            for(int currentPos = 0; basePairsDone < seqLength; currentPos += windowOffset) {

                //handles last window if theres overrun
                if(basePairsDone < seqLength && currentPos+windowSize >= seqLength){
                    windowSize = seqLength-currentPos;
                }

                PrintWriter writer = new PrintWriter(whereToMakeWindows + "/w" + windowName + ".fasta", "UTF-8");
                //write out the region where this window was taken from just for convenience
                PrintWriter writer2 = new PrintWriter(whereToMakeWindows + "/w" + windowName + ".region", "UTF-8");

                for(HashMap.Entry<String, String> entry : theSeqs.entrySet()) {
                    String key = entry.getKey();
                    String value = entry.getValue();

                    writer.println(key);
                    if(currentPos+windowSize < seqLength) {
                        //print in lengths of 60...need to change eventually probably
                        writer.println(value.substring(currentPos, currentPos + windowSize));

                    }else{
                        //print in lengths of 60...need to change eventually probably
                        writer.println(value.substring(currentPos, seqLength));
                    }
                    writer.println();

                }


                //writer.println("blah");
                writer2.println(currentPos + " to " + (currentPos+windowSize));

                writer.close();
                writer2.close();

                windowName++;

                //tick up basePairsDone
                basePairsDone = currentPos+windowSize;
            }







        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e){
            e.printStackTrace();
        }



    }



    public static void RunWindowedRAXML(String whereToRun, int windowSize, int windowOffset){

        TPrint("start|"+"raxml"+windowSize+"bp"+windowOffset+"os"); //timing info in case i need it later

        //add \ms folder to where we are running all this stuff
        String whereTheWindowsAre = whereToRun+"/windows"+windowSize+"bp"+windowOffset+"os";
        String whereToRunRAXML = whereToRun+"/raxml"+windowSize+"bp"+windowOffset+"os";


        MakeFolder(whereToRunRAXML); //create ms output folder


        //iterate through all fastas, example from http://www.avajava.com/tutorials/lessons/how-do-i-get-all-files-with-certain-extensions-in-a-directory-including-subdirectories.html
        File dir = new File(whereTheWindowsAre);
        String[] extensions = new String[] { "fasta" };
        List<File> files = (List<File>) FileUtils.listFiles(dir, extensions, true);
        for (File file : files) {

            //do everything here for a particular file
            try {

                //System.out.println("file: " + file.getCanonicalPath());

                //command string, one below runs their method and one runs the method with their in house correctness checker
                String commandString = "raxmlHPC -f a -x12345 -p 12345 -# 100 -k -m GTRGAMMA -s "+file.getName()+" -n "+file.getName();

                //String commandResults = RunCommandLineCommandFromDir(commandString, whereTheWindowsAre);
                String commandResults = GetOutputFromProgramNoHangFromDir(commandString, whereTheWindowsAre);
                //System.out.println(commandResults);
                //System.out.println(commandString);


                BufferedWriter out = new BufferedWriter(new FileWriter(whereToRunRAXML + "/fullOutput"));

                //save full output
                out.write(commandResults);

                //close the writers
                out.close();

            } catch (Exception e) {
                e.printStackTrace();
            }

        }


        TPrint("end|" + "raxml" + windowSize + "bp" + windowOffset + "os"); //timing info in case i need it later
    }

    public static void RunWindowedRAXMLNoTPrint(String whereToRun, int windowSize, int windowOffset){

        //TPrint("start|"+"raxml"+windowSize+"bp"+windowOffset+"os"); //timing info in case i need it later

        //add \ms folder to where we are running all this stuff
        String whereTheWindowsAre = whereToRun+"/windows"+windowSize+"bp"+windowOffset+"os";
        String whereToRunRAXML = whereToRun+"/raxml"+windowSize+"bp"+windowOffset+"os";


        MakeFolder(whereToRunRAXML); //create ms output folder


        //iterate through all fastas, example from http://www.avajava.com/tutorials/lessons/how-do-i-get-all-files-with-certain-extensions-in-a-directory-including-subdirectories.html
        File dir = new File(whereTheWindowsAre);
        String[] extensions = new String[] { "fasta" };
        List<File> files = (List<File>) FileUtils.listFiles(dir, extensions, true);
        for (File file : files) {

            //do everything here for a particular file
            try {

                //System.out.println("file: " + file.getCanonicalPath());

                //command string, one below runs their method and one runs the method with their in house correctness checker
                String commandString = "raxmlHPC -f a -x12345 -p 12345 -# 100 -k -m GTRGAMMA -s "+file.getName()+" -n "+file.getName();

                //String commandResults = RunCommandLineCommandFromDir(commandString, whereTheWindowsAre);
                String commandResults = GetOutputFromProgramNoHangFromDir(commandString, whereTheWindowsAre);
                //System.out.println(commandResults);
                //System.out.println(commandString);


                BufferedWriter out = new BufferedWriter(new FileWriter(whereToRunRAXML + "/fullOutput"));

                //save full output
                out.write(commandResults);

                //close the writers
                out.close();

            } catch (Exception e) {
                e.printStackTrace();
            }

        }


        //TPrint("end|" + "raxml" + windowSize + "bp" + windowOffset + "os"); //timing info in case i need it later
    }


    private static void ProcessRAXML(String baseFolder, int windowSize, int windowOffset, int genomeLength, boolean doCentering) {

        String centeringLabel = "";
        if(doCentering){
            centeringLabel="-CNTRD";
        }

        //add \ms folder to where we are running all this stuff
        String whereTheWindowsAre = baseFolder+"/windows"+windowSize+"bp"+windowOffset+"os";
        String whereToRunRAXML = baseFolder+"/raxml"+centeringLabel+windowSize+"bp"+windowOffset+"os";

        MakeFolder(whereToRunRAXML); //create ms output folder if for some reason it wasnt made



        //how many 'raxml best trees' are there?
        File dir = new File(whereTheWindowsAre);
        File [] files = dir.listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.startsWith("RAxML_bestTree");
            }
        });
        int howManyWindows = files.length;

        try {

            BufferedWriter out = new BufferedWriter(new FileWriter(whereToRunRAXML + "/treefileRAXML"));

            ArrayList<Integer> stretchLengths = GetStretchLengths(windowSize, windowOffset, genomeLength, doCentering);

            for (int i = 0; i < howManyWindows; i++) {


                //grab tree from file
                BufferedReader bestTreeFile = new BufferedReader(new FileReader(whereTheWindowsAre+"/RAxML_bestTree.w"+i+".fasta"));
                String bestTreeString = "";
                String line;
                while ((line = bestTreeFile.readLine()) != null) {
                    bestTreeString+=line;
                }

                //if Im using offset do the more complicated version where the last stretch is longer
                out.write("[" + stretchLengths.get(i)+"]"+bestTreeString+"\n"); //now just doing this in the stretch length getter function
                /*
                if(i == howManyWindows-1){
                    out.write("["+windowSize+"]"+bestTreeString+"\n"); //last tree is longer (in the small offset case)
                }
                else{
                    out.write("[" + windowOffset+"]"+bestTreeString+"\n");
                }*/

            }

            //close the writers
            out.close();

        } catch(Exception e){
            e.printStackTrace();
        }

    }

    private static ArrayList<Integer> GetStretchLengths(int windowSize, int windowOffset, int genomeLength, boolean doCentering) {

        //boolean doCentering = true;

        ArrayList<Integer> returnMe = new ArrayList<Integer>();

        int basePairsDone = 0;
        int currentPos = 0;
        while(basePairsDone < genomeLength) {

            int thisStretch = 0;

            //if doing centering, i just need to push the first window forward so the rest (minus the last) are outputting an offset# of bp located at their center
            //if not doing centering this first if will be ingnored and itll just do an offset for the first one (the last else), like normal
            //also, this only matters when offset<windowsize, but this equation also happens to evaluate to offset in that case, so it doesnt matter
            if(doCentering == true && currentPos==0){
                thisStretch = windowOffset + ((windowSize-windowOffset)/2); //first one is a bit more complicated
                //tick current position up a bit more than it would usually be ticked up by (always ticked up by just offset before)
                //currentPos+=((windowSize-windowOffset)/2);
            } else if(basePairsDone+windowOffset >= genomeLength){
                //handles last window if theres overrun, and works if there isnt. (just make sure to grab the rest of the bp for the last window)
                windowSize = genomeLength-currentPos;
                thisStretch = windowSize; // done, whatever is left is stretch
            }else{
                thisStretch = windowOffset; //not done, offset is stretch
            }

            returnMe.add(thisStretch);
            currentPos+=thisStretch;

            //tick up basePairsDone
            if(basePairsDone==0){
                basePairsDone = windowSize;}
            else{
                basePairsDone += windowOffset;
            }
        }

        return returnMe;
    }


    private static void SmoothWinds(String whereToRun, int windowSize, int windowOffset, int genomeLength, boolean useBL, boolean doCentering){


        //add \ms folder to where we are running all this stuff
        String whereTheWindowsAre = whereToRun+"/windows"+windowSize+"bp"+windowOffset+"os";
        String whereIsRAXML = whereToRun+"/raxml"+windowSize+"bp"+windowOffset+"os";
        String whereToRunSmoothWinds;

        String centeringLabel = "";
        if(doCentering){
            centeringLabel = "-CNTRD";
        }

        if(useBL){
            TPrint("start|"+"smoothwindsBL"+centeringLabel+windowSize+"bp"+windowOffset+"os"); //timing info in case i need it later
            whereToRunSmoothWinds = whereToRun + "/smoothwindsBL"+centeringLabel + windowSize + "bp" + windowOffset + "os";
        }else {
            TPrint("start|"+"smoothwinds"+centeringLabel+windowSize+"bp"+windowOffset+"os"); //timing info in case i need it later
            whereToRunSmoothWinds = whereToRun + "/smoothwinds"+centeringLabel + windowSize + "bp" + windowOffset + "os";
        }

        //this will be populated after reading in all the bootstrap files, to be made into a graph after
        ArrayList<HashSet<String>> orderedWindowBootstrapTrees = new ArrayList<HashSet<String>>();

        MakeFolder(whereToRunSmoothWinds); //create ms output folder

        //how many 'raxml bootstrap files' are there?
        File dir = new File(whereTheWindowsAre);
        File [] files = dir.listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.startsWith("RAxML_bootstrap");
            }
        });
        int howManyWindows = files.length;


        //loop through bootstrap files
        try {

            //BufferedWriter out = new BufferedWriter(new FileWriter(whereToRunSmoothWinds + "/treefileSmoothWinds"));

            for (int i = 0; i < howManyWindows; i++) {

                HashSet<String> theseTrees = new HashSet<String>();

                //grab trees from file
                BufferedReader bsTreeFile = new BufferedReader(new FileReader(whereTheWindowsAre+"/RAxML_bootstrap.w"+i+".fasta"));
                //String bestTreeString = "";
                String line;
                while ((line = bsTreeFile.readLine()) != null) {
                    theseTrees.add(line);
                    //bestTreeString+=line;
                }

                orderedWindowBootstrapTrees.add(theseTrees);

            }

            //close the writers
            //out.close();

        } catch(Exception e){
            e.printStackTrace();
        }


        //make graph out of my arraylist of hashsets based on tree similarity formulation

        // constructs a directed graph, each tree has an int appended in the node name starting with 0 on the left most trees
        org.jgrapht.Graph<String, DefaultWeightedEdge> graphBootTrees = new SimpleDirectedWeightedGraph<String, DefaultWeightedEdge>(DefaultWeightedEdge.class);

        //add a start node and the first layer of vertices
        graphBootTrees.addVertex("start");
        for(String firstSetTree : orderedWindowBootstrapTrees.get(0)){
            graphBootTrees.addVertex(firstSetTree+0);
            //add edge from start to this one
            DefaultWeightedEdge thisEdge = graphBootTrees.addEdge("start", firstSetTree+0);
            graphBootTrees.setEdgeWeight(thisEdge, 0.0);
        }

        //now add second and onwards layers of vertices with edges from the last layer to this layer, and end node attached to last layer
        graphBootTrees.addVertex("end");
        for(int currentTreeIndex = 1; currentTreeIndex<orderedWindowBootstrapTrees.size(); currentTreeIndex++){

            for(String currentSetTree : orderedWindowBootstrapTrees.get(currentTreeIndex)){
                //add a vertex for it
                graphBootTrees.addVertex(currentSetTree+currentTreeIndex);
                //add weighted edges before all the ones directly before this one to this one
                for(String previousSetTree : orderedWindowBootstrapTrees.get(currentTreeIndex-1)){
                    //add edge from previous to current tree node based on edge weight calculator function
                    DefaultWeightedEdge thisEdge = graphBootTrees.addEdge(previousSetTree+(currentTreeIndex-1), currentSetTree+currentTreeIndex);
                    if(useBL){
                        graphBootTrees.setEdgeWeight(thisEdge, CalculateEdgeWeightBLV1(previousSetTree, currentSetTree));
                    }else {
                        graphBootTrees.setEdgeWeight(thisEdge, CalculateEdgeWeightV1(previousSetTree, currentSetTree));
                    }
                }

                //if this is the last layer add a connection to the end node
                if(currentTreeIndex==(orderedWindowBootstrapTrees.size()-1)){
                    DefaultWeightedEdge thisEdge = graphBootTrees.addEdge(currentSetTree+currentTreeIndex,"end");
                    graphBootTrees.setEdgeWeight(thisEdge, 0.0);
                }

            }

        }


        //get shortest path

        //these two lines are the only out of a ton of internet examples that actually compiles and works 0.0
        DijkstraShortestPath<String, DefaultWeightedEdge> runDijkstra = new DijkstraShortestPath<String, DefaultWeightedEdge>(graphBootTrees);
        List<String> shortestPath = runDijkstra.getPath("start","end").getVertexList();
        double shortestPathWeight = runDijkstra.getPath("start","end").getWeight();

        //print the ordered list of vertices on shortest path
        //System.out.println("Shortest path from start to end:");
        //System.out.println(shortestPath);
        //debugging for new formulation to save memory
        System.out.println("Shortest path weight: "+shortestPathWeight);
        TPrint("w graph implement, total path RF " + shortestPathWeight);
        //List<DefaultWeightedEdge> sp;

        //these two lines are the only out of a ton of internet examples that actually compiles and works 0.0
        //DijkstraShortestPath<String, DefaultWeightedEdge> runDijkstra = new DijkstraShortestPath<String, DefaultWeightedEdge>(graph);
        //List<String> shortestPath = runDijkstra.getPath("vertex1","vertex5").getVertexList();

        //print the ordered list of vertices on shortest path
        //System.out.println(shortestPath);

        int look = 0;

        //print the trees out
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(whereToRunSmoothWinds + "/treefileSmoothWinds"));

            //start at 1 to skip 'start', and end one short to skip 'end'
            //howManyWindows = (shortestPath.size()-2); //already populated
            ArrayList<Integer> stretchLengths = GetStretchLengths(windowSize,windowOffset,genomeLength, doCentering);
            for (int i = 0; i < howManyWindows; i++) {

                //grab tree from file
                String bestTreeString = shortestPath.get(i+1).split(";")[0] + ";"; //strip end num off

                //if Im using offset do the more complicated version where the last stretch is longer
                out.write("[" + stretchLengths.get(i) + "]" + bestTreeString + "\n"); //now just doing this in stretch length getter function
                /*
                if (i == howManyWindows - 1) {
                    out.write("[" + windowSize + "]" + bestTreeString + "\n"); //last tree is longer (in the small offset case)
                } else {
                    out.write("[" + windowOffset + "]" + bestTreeString + "\n");
                }*/

            }

            //close the writers
            out.close();
        }catch(Exception e){
            e.printStackTrace();
        }


    int checkk=0;





        TPrint("end|"+"smoothwinds"+windowSize+"bp"+windowOffset+"os"); //timing info in case i need it later
    }

    public static void SmoothWindsLessMem(String whereToRun, int windowSize, int windowOffset, int genomeLength, boolean useBL, boolean doCentering){


        //add \ms folder to where we are running all this stuff
        String whereTheWindowsAre = whereToRun+"/windows"+windowSize+"bp"+windowOffset+"os";
        String whereIsRAXML = whereToRun+"/raxml"+windowSize+"bp"+windowOffset+"os";
        String whereToRunSmoothWinds;

        String centeringLabel = "";
        if(doCentering){
            centeringLabel = "-CNTRD";
        }

        if(useBL){
            TPrint("start|"+"smoothwindsLessMemBL"+centeringLabel+windowSize+"bp"+windowOffset+"os"); //timing info in case i need it later
            whereToRunSmoothWinds = whereToRun + "/smoothwindsLessMemBL"+centeringLabel + windowSize + "bp" + windowOffset + "os";
        }else {
            TPrint("start|"+"smoothwindsLessMem"+centeringLabel+windowSize+"bp"+windowOffset+"os"); //timing info in case i need it later
            whereToRunSmoothWinds = whereToRun + "/smoothwindsLessMem"+centeringLabel + windowSize + "bp" + windowOffset + "os";
        }

        //this will be populated after reading in all the bootstrap files, to be made into a graph after
        //ArrayList<HashSet<String>> orderedWindowBootstrapTrees = new ArrayList<HashSet<String>>();

        MakeFolder(whereToRunSmoothWinds); //create ms output folder

        //how many 'raxml bootstrap files' are there?
        File dir = new File(whereTheWindowsAre);
        File [] files = dir.listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.startsWith("RAxML_bootstrap");
            }
        });
        int howManyWindows = files.length;

        ArrayList<WindowPathGraphNode> previousNodes = new ArrayList<WindowPathGraphNode>();
        ArrayList<WindowPathGraphNode> theseNodes;

        //loop through bootstrap files, im now going to do everything inside this loop so im only doing one 'layer' at a time
        try {

            //BufferedWriter out = new BufferedWriter(new FileWriter(whereToRunSmoothWinds + "/treefileSmoothWinds"));



            for (int i = 0; i < howManyWindows; i++) {

                HashSet<String> theseTrees = new HashSet<String>();

                //grab trees from file
                BufferedReader bsTreeFile = new BufferedReader(new FileReader(whereTheWindowsAre+"/RAxML_bootstrap.w"+i+".fasta"));
                //String bestTreeString = "";
                String line;
                while ((line = bsTreeFile.readLine()) != null) {
                    theseTrees.add(line);
                    //bestTreeString+=line;
                }

                // THIS IF STATEMENT HAS ALL THE GRAPH PROCESSING STUFF, FOLLOWED BY GRABBING THE BEST PATH
                //process the very first trees slightly differently
                if(i == 0){
                    //init the nodes and add to previous nodes list
                    previousNodes = new ArrayList<WindowPathGraphNode>();
                    for(String singleInitialTree : theseTrees){
                        WindowPathGraphNode makeThisTreeNode = new WindowPathGraphNode(singleInitialTree);
                        makeThisTreeNode.InitFirstLayerNode(); //sets best score to 0 from -1.0
                        makeThisTreeNode.addSelfToOwnBestPath();
                        previousNodes.add(makeThisTreeNode);
                    }
                }else{
                    //create this set of window path graph nodes and compare each one to all the previous nodes
                    theseNodes = new ArrayList<WindowPathGraphNode>();
                    for(String singleCurrentTree : theseTrees){
                        WindowPathGraphNode thisTreeNode = new WindowPathGraphNode(singleCurrentTree);

                        //compare to all previous nodes
                        for(WindowPathGraphNode previousNode : previousNodes){
                            thisTreeNode.ProcessThisPathRF(previousNode);
                        }

                        //after, a best path to this node has been found, add self to own best path
                        thisTreeNode.addSelfToOwnBestPath();
                        theseNodes.add(thisTreeNode);
                    }

                    //this layer has been process, wipe out previous nodes (make previous nodes = these nodes)
                    previousNodes = theseNodes;
                }

                //orderedWindowBootstrapTrees.add(theseTrees);

            }

            //close the writers
            //out.close();


        } catch(Exception e){
            e.printStackTrace();
        }



        //once its done, the nodes should be populated with best paths to each. lets grab the best of the best for the single best path
        ArrayList<String> singleBestPath = previousNodes.get(0).bestPathToThisTree; //init arbitrarily to first entry
        double singleBestPathScore = previousNodes.get(0).bestScoreToThisTree; //init to arbitrary score of one of em (first entry)
        WindowPathGraphNode bestFinalNode = previousNodes.get(0);

        //previous nodes should be currently holding the last ones
        for(WindowPathGraphNode singleFinalNode : previousNodes){
            if(singleFinalNode.bestScoreToThisTree < singleBestPathScore){
                singleBestPathScore = singleFinalNode.bestScoreToThisTree;
                singleBestPath = singleFinalNode.bestPathToThisTree;
                bestFinalNode = singleFinalNode;
            }
        }

        TPrint("[comp.above] rf total path length "+singleBestPathScore);
        //TPrint();
        System.out.println("rf total path length "+singleBestPathScore);

        ArrayList<String> shortestPath = singleBestPath;


        //print the trees out
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(whereToRunSmoothWinds + "/treefileSmoothWinds"));

            //start at 1 to skip 'start', and end one short to skip 'end'
            //howManyWindows = (shortestPath.size()-2); //already populated
            ArrayList<Integer> stretchLengths = GetStretchLengths(windowSize,windowOffset,genomeLength, doCentering);
            for (int i = 0; i < howManyWindows; i++) {

                //grab tree from file (removing +1 in get cuz i dont have a 'start' node this time
                String bestTreeString = shortestPath.get(i).split(";")[0] + ";"; //strip end num off

                //if Im using offset do the more complicated version where the last stretch is longer
                out.write("[" + stretchLengths.get(i) + "]" + bestTreeString + "\n"); //now just doing this in stretch length getter function
                /*
                if (i == howManyWindows - 1) {
                    out.write("[" + windowSize + "]" + bestTreeString + "\n"); //last tree is longer (in the small offset case)
                } else {
                    out.write("[" + windowOffset + "]" + bestTreeString + "\n");
                }*/

            }

            //close the writers
            out.close();
        }catch(Exception e){
            e.printStackTrace();
        }


        int checkk=0;





        TPrint("end|"+"smoothwindsLessMem"+windowSize+"bp"+windowOffset+"os"); //timing info in case i need it later
    }

    public static void SmoothWindsLessMemNoTPrint(String whereToRun, int windowSize, int windowOffset, int genomeLength, boolean useBL, boolean doCentering){


        //add \ms folder to where we are running all this stuff
        String whereTheWindowsAre = whereToRun+"/windows"+windowSize+"bp"+windowOffset+"os";
        String whereIsRAXML = whereToRun+"/raxml"+windowSize+"bp"+windowOffset+"os";
        String whereToRunSmoothWinds;

        String centeringLabel = "";
        if(doCentering){
            centeringLabel = "-CNTRD";
        }

        if(useBL){
            //TPrint("start|"+"smoothwindsLessMemBL"+centeringLabel+windowSize+"bp"+windowOffset+"os"); //timing info in case i need it later
            whereToRunSmoothWinds = whereToRun + "/smoothwindsLessMemBL"+centeringLabel + windowSize + "bp" + windowOffset + "os";
        }else {
            //TPrint("start|"+"smoothwindsLessMem"+centeringLabel+windowSize+"bp"+windowOffset+"os"); //timing info in case i need it later
            whereToRunSmoothWinds = whereToRun + "/smoothwindsLessMem"+centeringLabel + windowSize + "bp" + windowOffset + "os";
        }

        //this will be populated after reading in all the bootstrap files, to be made into a graph after
        //ArrayList<HashSet<String>> orderedWindowBootstrapTrees = new ArrayList<HashSet<String>>();

        MakeFolder(whereToRunSmoothWinds); //create ms output folder

        //how many 'raxml bootstrap files' are there?
        File dir = new File(whereTheWindowsAre);
        File [] files = dir.listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.startsWith("RAxML_bootstrap");
            }
        });
        int howManyWindows = files.length;

        ArrayList<WindowPathGraphNode> previousNodes = new ArrayList<WindowPathGraphNode>();
        ArrayList<WindowPathGraphNode> theseNodes;

        //loop through bootstrap files, im now going to do everything inside this loop so im only doing one 'layer' at a time
        try {

            //BufferedWriter out = new BufferedWriter(new FileWriter(whereToRunSmoothWinds + "/treefileSmoothWinds"));



            for (int i = 0; i < howManyWindows; i++) {

                HashSet<String> theseTrees = new HashSet<String>();

                //grab trees from file
                BufferedReader bsTreeFile = new BufferedReader(new FileReader(whereTheWindowsAre+"/RAxML_bootstrap.w"+i+".fasta"));
                //String bestTreeString = "";
                String line;
                while ((line = bsTreeFile.readLine()) != null) {
                    theseTrees.add(line);
                    //bestTreeString+=line;
                }

                // THIS IF STATEMENT HAS ALL THE GRAPH PROCESSING STUFF, FOLLOWED BY GRABBING THE BEST PATH
                //process the very first trees slightly differently
                if(i == 0){
                    //init the nodes and add to previous nodes list
                    previousNodes = new ArrayList<WindowPathGraphNode>();
                    for(String singleInitialTree : theseTrees){
                        WindowPathGraphNode makeThisTreeNode = new WindowPathGraphNode(singleInitialTree);
                        makeThisTreeNode.InitFirstLayerNode(); //sets best score to 0 from -1.0
                        makeThisTreeNode.addSelfToOwnBestPath();
                        previousNodes.add(makeThisTreeNode);
                    }
                }else{
                    //create this set of window path graph nodes and compare each one to all the previous nodes
                    theseNodes = new ArrayList<WindowPathGraphNode>();
                    for(String singleCurrentTree : theseTrees){
                        WindowPathGraphNode thisTreeNode = new WindowPathGraphNode(singleCurrentTree);

                        //compare to all previous nodes
                        for(WindowPathGraphNode previousNode : previousNodes){
                            thisTreeNode.ProcessThisPathRF(previousNode);
                        }

                        //after, a best path to this node has been found, add self to own best path
                        thisTreeNode.addSelfToOwnBestPath();
                        theseNodes.add(thisTreeNode);
                    }

                    //this layer has been process, wipe out previous nodes (make previous nodes = these nodes)
                    previousNodes = theseNodes;
                }

                //orderedWindowBootstrapTrees.add(theseTrees);

            }

            //close the writers
            //out.close();


        } catch(Exception e){
            e.printStackTrace();
        }



        //once its done, the nodes should be populated with best paths to each. lets grab the best of the best for the single best path
        ArrayList<String> singleBestPath = previousNodes.get(0).bestPathToThisTree; //init arbitrarily to first entry
        double singleBestPathScore = previousNodes.get(0).bestScoreToThisTree; //init to arbitrary score of one of em (first entry)
        WindowPathGraphNode bestFinalNode = previousNodes.get(0);

        //previous nodes should be currently holding the last ones
        for(WindowPathGraphNode singleFinalNode : previousNodes){
            if(singleFinalNode.bestScoreToThisTree < singleBestPathScore){
                singleBestPathScore = singleFinalNode.bestScoreToThisTree;
                singleBestPath = singleFinalNode.bestPathToThisTree;
                bestFinalNode = singleFinalNode;
            }
        }

        //TPrint("[comp.above] rf total path length "+singleBestPathScore);
        //TPrint();
        //System.out.println("rf total path length "+singleBestPathScore);

        ArrayList<String> shortestPath = singleBestPath;


        //print the trees out
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(whereToRunSmoothWinds + "/treefileSmoothWinds"));

            //start at 1 to skip 'start', and end one short to skip 'end'
            //howManyWindows = (shortestPath.size()-2); //already populated
            ArrayList<Integer> stretchLengths = GetStretchLengths(windowSize,windowOffset,genomeLength, doCentering);
            for (int i = 0; i < howManyWindows; i++) {

                //grab tree from file (removing +1 in get cuz i dont have a 'start' node this time
                String bestTreeString = shortestPath.get(i).split(";")[0] + ";"; //strip end num off

                //if Im using offset do the more complicated version where the last stretch is longer
                out.write("[" + stretchLengths.get(i) + "]" + bestTreeString + "\n"); //now just doing this in stretch length getter function
                /*
                if (i == howManyWindows - 1) {
                    out.write("[" + windowSize + "]" + bestTreeString + "\n"); //last tree is longer (in the small offset case)
                } else {
                    out.write("[" + windowOffset + "]" + bestTreeString + "\n");
                }*/

            }

            //close the writers
            out.close();
        }catch(Exception e){
            e.printStackTrace();
        }


        int checkk=0;





        //TPrint("end|"+"smoothwindsLessMem"+windowSize+"bp"+windowOffset+"os"); //timing info in case i need it later
    }



    public static double CalculateEdgeWeightV1(String fromTree, String toTree) {


        //lets just start by making it the RF distance between the two trees

        SymmetricDifference aaaaa = new SymmetricDifference();
        aaaaa.computeDifference(Trees.readTree(fromTree), Trees.readTree(toTree), false);
        double theRFDistCalculated = aaaaa.getUnweightedAverage();

        return theRFDistCalculated;

    }

    public static double CalculateEdgeWeightBLV1(String fromTree, String toTree) {


        //newer formulation with NRBS stuff to take into account BLs

        //unweighted rf dist
        //SymmetricDifference aaaaa = new SymmetricDifference();
        //aaaaa.computeDifference(Trees.readTree(fromTree), Trees.readTree(toTree), true);
        //double theRFDistCalculated = aaaaa.getUnweightedAverage();

        double nrbsDist = (new RFDistanceBL(Trees.readTree(fromTree), Trees.readTree(toTree))).getDistance();

        return nrbsDist;

    }



    //code/tutorial for running command line commands in java from https://www.linglom.com/programming/java/how-to-run-command-line-or-execute-external-application-from-java/
    public static String RunCommandLineCommand(String command) {

        //have to do this to get piping to work correctly
        String[] cmd = {
                "/bin/sh",
                "-c",
                command
        };

        StringBuilder stringToOutput = new StringBuilder();
        try {
            Runtime rt = Runtime.getRuntime();

            Process pr = rt.exec(cmd);

            BufferedReader input = new BufferedReader(new InputStreamReader(pr.getInputStream()));

            String line=null;

            while((line=input.readLine()) != null) {
                //System.out.println(line);
                stringToOutput.append(line+"\n");
            }

            int exitVal = pr.waitFor();
            System.out.println("Exited with error code " + exitVal);

        } catch(Exception e) {
            System.out.println(e.toString());
            e.printStackTrace();
        }

        return stringToOutput.toString();
    }



    public static String RunCommandLineCommandFromDir(String command,String runDir) {

        //have to do this to get piping to work correctly
        String[] cmd = {
                "/bin/sh",
                "-c",
                command
        };

        StringBuilder stringToOutput = new StringBuilder();
        try {
            Runtime rt = Runtime.getRuntime();

            Process pr = rt.exec(cmd,null,new File(runDir));

            BufferedReader input = new BufferedReader(new InputStreamReader(pr.getInputStream()));

            String line=null;

            while((line=input.readLine()) != null) {
                //System.out.println(line);
                stringToOutput.append(line+"\n");
            }

            int exitVal = pr.waitFor();
            System.out.println("Exited with error code " + exitVal);

        } catch(Exception e) {
            System.out.println(e.toString());
            e.printStackTrace();
        }

        return stringToOutput.toString();
    }




    //helpful code snippet from https://stackoverflow.com/questions/13008526/runtime-getruntime-execcmd-hanging
    //to prevent argweaver from hanging (the reason why it hangs seems to be semi complicated having to do with streams and things running on a new thread)

    //it seems there may be a problem tho with there not being correct waiting for things to finish before, say, trying to access files the process creates
    public static String GetOutputFromProgramNoHangFromDir(String program, String runDir) throws IOException {

        //have to do this to get piping to work correctly
        String[] cmd = {
                "/bin/sh",
                "-c",
                program
        };


        Process proc;
        if(runDir == null){
            proc = Runtime.getRuntime().exec(cmd);
        }else{
            proc = Runtime.getRuntime().exec(cmd, null, new File(runDir));
        }



        return Stream.of(proc.getErrorStream(), proc.getInputStream()).parallel().map((InputStream isForOutput) -> {
            StringBuilder output = new StringBuilder();
            try (BufferedReader br = new BufferedReader(new InputStreamReader(isForOutput))) {
                String line;
                while ((line = br.readLine()) != null) {
                    output.append(line);
                    output.append("\n");
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            //adding this block to do a wait for to not having things going faster than they should when outputs are needed in future steps
            //not doing it tho..
            //try {
            //    proc.waitFor();
            //} catch (InterruptedException e) {
            //    e.printStackTrace();
            //}

            return output;
        }).collect(Collectors.joining());
    }


    //helper function stolen from leo script temp to figure out which arg had highest likelihood to get point estimate of local topologies from
    //might later add offset and every how many samples one was actually written out
    private static String getHighestLikelihoodArg(String folderWithStats) {

        double highestLikelihood = Double.NEGATIVE_INFINITY;
        int highestLikelihoodFileNum = 0;

        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(folderWithStats+"/arg-sample.stats"));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        String line;

        try {

            //get rid of first two unnecessary lines
            line = br.readLine();
            line = br.readLine();

            while ((line = br.readLine()) != null) {

                int sampleNum = Integer.parseInt(line.split("\t")[1]);
                double sampleLikelihood = Double.parseDouble(line.split("\t")[3]);

                //if multiple of ten sample (only ones output), do stuff, otherwise dont
                //this would be the line to change later if i change how many samples are outputted
                if(sampleNum % 10 == 0){

                    if(sampleLikelihood > highestLikelihood){
                        highestLikelihoodFileNum = sampleNum;
                        highestLikelihood = sampleLikelihood;
                    }

                }else{
                    //do nothing
                }

            }
        } catch (IOException e) {
            e.printStackTrace();
        }


        return folderWithStats+"/arg-sample."+highestLikelihoodFileNum+".smc";
    }

    //helper function stolen from leo script temp to turn an arg into an ms style local trees file
    private static void argweaverToMS(String theFileSMC, String whereToWriteTo) {

        int startingIndex = 0;
        String tree = "";
        int index = 0;

        int sanityCheck = 0;

        try {

            String line;

            PrintWriter writer = new PrintWriter(whereToWriteTo, "UTF-8");

            boolean EOF = false;


            BufferedReader br = new BufferedReader(new FileReader(theFileSMC));

            //first line has names to numbers mapping...

            line = br.readLine();
            HashMap<String,String> numsToNames = new HashMap<String,String>();
            int i = -1;
            for(String numMapMe : line.split("\t")){
                numsToNames.put(Integer.toString(i),numMapMe.replace(" ",""));
                i++;
            }


            while ((line = br.readLine()) != null) {

                if(!line.contains("TREE")){
                    continue;
                }

                String[] splitLine = line.split("\t");

                int length = Integer.parseInt(splitLine[2]) - Integer.parseInt(splitLine[1]) + 1;
                String[] treeStringStuff = splitLine[3].split("\\[&&NHX:age=\\w*.\\w*\\]");
                String finalTreeString = String.join("", treeStringStuff);
                //System.out.println(finalTreeString);
                //finalTreeString.replaceAll("4:\\w*.\\w*","");
                //finalTreeString.replaceAll("6","");
                Tree tr = Trees.readTree(finalTreeString);

                for(TNode renameMe : tr.getNodes()){

                    //only messing with leaves
                    if(renameMe.getChildCount()==0) {
                        //rename each node based on name mapping in smc file
                        String currentNodeName = renameMe.getName();
                        String newNodeName = "temp_"+numsToNames.get(currentNodeName);
                        ((STINode) renameMe).setName(newNodeName);
                        //System.out.println("renamed "+currentNodeName+" to "+newNodeName);
                    }else{
                        //they really didnt make this easy on me...duplicate names to my ms output names can exist on internal tree nodes >.< ...this should take care of it
                        ((STINode) renameMe).setName("i_"+renameMe.getName());
                    }
                }


                //then go back and remove temp_ from names since i had to do that to avoid duplicate naming issues
                for(TNode renameMe : tr.getNodes()){
                    //only messing with leaves
                    if(renameMe.getChildCount()==0) {
                        //rename each node based on name mapping in smc file
                        String currentNodeName = renameMe.getName();
                        String newNodeName = currentNodeName.substring(5);
                        ((STINode) renameMe).setName(newNodeName);
                    }
                }


                String finalTreeStringWithReplacement = tr.toNewick();



                //System.out.println("["+length+"]"+tr.toNewick());

                /*
                //print it out after replacing things to be correct (argweaver species names out of order)
                finalTreeString = finalTreeString.replace("(0:","(OrOut:").replace("(1:", "(ChWes:").replace("(2:","(ChEas:").replace("(3:","(ChCen:").replace("(4:","(ChNig:").replace("(5:","(Human:").replace(",0:", ",OrOut:").replace(",1:", ",ChWes:").replace(",2:",",ChEas:").replace(",3:",",ChCen:").replace(",4:",",ChNig:").replace(",5:",",Human:");
                //this line is a total hack for right now, as the species nums go up this will become annoying (its cuz i index 1 to num and argweaver is indexing 0 to num)
                finalTreeString = finalTreeString.replace("6:","99:");
                finalTreeString = finalTreeString.replace("ChEas", "3").replace("ChCen","2").replace("ChWes","5").replace("ChNig","4").replace("Human","1").replace("OrOut","6");
                */

                //have to now do this dumb crap dynamically
                //Integer test1 = namesToNums.get("OrOut");
                //finalTreeString = finalTreeString.replace("("+namesToNums.get("OrOut")+":","(OrOut:").replace("("+namesToNums.get("ChWes")+":", "(ChWes:").replace("("+namesToNums.get("ChEas")+":","(ChEas:").replace("("+namesToNums.get("ChCen")+":","(ChCen:").replace("("+namesToNums.get("ChNig")+":","(ChNig:").replace("("+namesToNums.get("Human")+":","(Human:").replace(","+namesToNums.get("OrOut")+":", ",OrOut:").replace(","+namesToNums.get("ChWes")+":", ",ChWes:").replace(","+namesToNums.get("ChEas")+":",",ChEas:").replace(","+namesToNums.get("ChCen")+":",",ChCen:").replace(","+namesToNums.get("ChNig")+":",",ChNig:").replace(","+namesToNums.get("Human")+":",",Human:");
                //this line is a total hack for right now, as the species nums go up this will become annoying (its cuz i index 1 to num and argweaver is indexing 0 to num)
                //finalTreeString = finalTreeString.replace("6:","99:");
                //finalTreeString = finalTreeString.replace("ChEas", "3").replace("ChCen","2").replace("ChWes","5").replace("ChNig","4").replace("Human","1").replace("OrOut","6");
                //finalTreeString = finalTreeString.replace(";","999;");

                //System.out.println("["+length+"]"+finalTreeString);
                writer.write("[" + length + "]" + finalTreeStringWithReplacement + "\n");

                //tr.rerootTreeAtNode(tr.getNode("6"));



            }



            writer.close();



        } catch (Exception e) {
            System.out.println("san check : " + sanityCheck);
            e.printStackTrace();
        }



    }

    //helper function taken from leo script temp to compare two MS files and get results on accuracy
    public static void compareTwoMSFiles(String trueTrees, String predictedTrees,String whereToWriteResults, int numBP){

        //arg1 - the truth (the actual output from MS)
        //arg2 - the predictions (how well did whatever method do?

        //antiquated code not used for RF calc:
        boolean predictionsRooted = false;


        try {

            int predictionDebugBPInitial = 0;

            BufferedReader brTruth = new BufferedReader(new FileReader(trueTrees));
            BufferedReader brPrediction = new BufferedReader(new FileReader(predictedTrees));

            BufferedWriter writeResults = new BufferedWriter(new FileWriter(whereToWriteResults));
            StringBuilder resultsText = new StringBuilder();

            ArrayList<String> truthText = new ArrayList<>();
            ArrayList<String> predictionText = new ArrayList<>();

            //read in the two files text and create map of the trees in the truth
            String line;
            while ((line = brPrediction.readLine()) != null) {
                predictionText.add(line);
            }
            HashMap<STITree,Integer> topologyTally = new HashMap<STITree,Integer>();

            //while ((line = brTruth.readLine()) != null && line.contains("]")) {
            //    truthText.add(line);
            //}
            //System.out.print(truthText.get(truthText.size()-2));
            //System.out.print(truthText.get(truthText.size()-1));


            while ((line = brTruth.readLine()) != null && line.contains("]")) {
                truthText.add(line);

                //region OLD STUFF FROM LEOSCRIPT_TEMP
                //add entry to topology tally
                //String firstParse = line.split("\\[")[1];
                /*
                int howManyBases = Integer.parseInt((line.split("]")[0]).split("\\[")[1]);

                Tree thisTree = Trees.readTree(line.split("]")[1]);

                //root on outgroup
                thisTree.rerootTreeAtNode(thisTree.getNode("6"));
                //set branch lengths to use RFDistBL
                for (TNode thisNode : thisTree.getNodes()) {
                    ((STINode)thisNode).setParentDistance(1.0);
                }
                boolean foundOne = false;
                for(Map.Entry<STITree,Integer> blah : topologyTally.entrySet()){
                    double similarity = (new RFDistanceBL(blah.getKey(),thisTree)).getDistance();
                    if(similarity == 0.0){
                        blah.setValue(blah.getValue()+howManyBases);
                        foundOne = true;
                    }
                }
                if(!foundOne){
                    topologyTally.put((STITree)thisTree,howManyBases);
                }
                //end add entry to topology Tally
                */
                //endregion FRO

            }

            //region OLD STUFF FROM LEOSCRIPT_TEMP
            /*
            //turn topology tally into ordered list
            //ArrayList<Object> orderedTopologyTally = new ArrayList<>();

            int sumOfSites = 0;
            for(Object printTreeInfo : topologyTally.entrySet()){
                STITree theTreeToPrint = (STITree)((HashMap.Entry)printTreeInfo).getKey();
                System.out.println(theTreeToPrint.toString() + " occurs " + topologyTally.get(theTreeToPrint) + " times");
                resultsText.append(theTreeToPrint.toString() + " occurs " + topologyTally.get(theTreeToPrint) + " times"+"\n");
                sumOfSites += topologyTally.get(theTreeToPrint);
            }
            System.out.println("Total of sites(debug): " + sumOfSites);
            resultsText.append("Total of sites(debug): " + sumOfSites+"\n");
            //int debug = 0;


            //now for each type of tree, calculate TPR and FPR, as well as calculate weighted average TPr and FPR
            double averageTPR = 0.0;
            double averageFPR = 0.0;
            for(Object printTreeInfo : topologyTally.entrySet()){
                STITree theTreeToPrint = (STITree)((HashMap.Entry)printTreeInfo).getKey();
                int numOfSites = topologyTally.get(theTreeToPrint);

                double TPR = 0.0;
                double FPR = 0.0;

                //make two int arrays, 1 if topology, 0 if not topology
                int[] truthInt = new int[sumOfSites];
                int[] predictionInt = new int[sumOfSites];

                //POPULATE TRUTH INT ARRAY
                int whereInArrayAmI = 0;
                for(String oneLine : truthText){
                    int howManyBases = Integer.parseInt((oneLine.split("]")[0]).split("\\[")[1]);
                    int newArrayPosition = whereInArrayAmI + howManyBases;
                    Tree thisTree = Trees.readTree(oneLine.split("]")[1]);

                    int isThisMatch = 0;
                    //check if RFdistBL is the same
                    //set branch lengths to use RFDistBL
                    for (TNode thisNode : thisTree.getNodes()) {
                        ((STINode)thisNode).setParentDistance(1.0);
                    }
                    double similarity = (new RFDistanceBL(theTreeToPrint,thisTree)).getDistance();
                    if(similarity == 0.0){
                        isThisMatch = 1;
                    }

                    while(whereInArrayAmI < newArrayPosition){
                        truthInt[whereInArrayAmI] = isThisMatch;
                        whereInArrayAmI++;
                    }
                }

                //POPULATE PREDICTION INT ARRAY
                whereInArrayAmI = 0;

                for(String oneLine : predictionText){
                    int howManyBases = Integer.parseInt((oneLine.split("]")[0]).split("\\[")[1]);
                    predictionDebugBPInitial += howManyBases;
                    int newArrayPosition = whereInArrayAmI + howManyBases;
                    Tree thisTree = Trees.readTree(oneLine.split("]")[1]);

                    int isThisMatch = 0;
                    //check if RFdistBL is the same
                    //set branch lengths to use RFDistBL and reroot on outgroup
                    if(predictionsRooted == false) {
                        //System.out.println("debug:"+howManyBases);
                        thisTree.rerootTreeAtNode(thisTree.getNode("6"));
                    }
                    for (TNode thisNode : thisTree.getNodes()) {
                        ((STINode)thisNode).setParentDistance(1.0);
                    }
                    double similarity = (new RFDistanceBL(theTreeToPrint,thisTree)).getDistance();
                    if(similarity == 0.0){
                        isThisMatch = 1;
                    }

                    while(whereInArrayAmI < newArrayPosition){
                        predictionInt[whereInArrayAmI] = isThisMatch;
                        whereInArrayAmI++;
                    }
                }

                //calculate TPR and FPR
                double TP = 0;
                double P = 0;
                double FP = 0;
                double N = 0;
                //calculate TP and P, FP and N
                for(int i = 0; i < sumOfSites; i++){
                    if(predictionInt[i] == 1 && truthInt[i] == 1){
                        TP++;
                    }
                    if(predictionInt[i] == 1 && truthInt[i] == 0){
                        FP++;
                    }
                    if(truthInt[i] == 1){
                        P++;
                    }else{
                        N++;
                    }
                }
                TPR = TP / P;
                FPR = FP / N;

                //add to weighted average TPR and FPR
                averageTPR += TPR * numOfSites / sumOfSites;
                averageFPR += FPR * numOfSites / sumOfSites;



                System.out.println(theTreeToPrint.toString() + " occurs " + topologyTally.get(theTreeToPrint) + " times");
                System.out.println("\tTPR: " + TPR);
                System.out.println("\tFPR: " + FPR);
                resultsText.append(theTreeToPrint.toString() + " occurs " + topologyTally.get(theTreeToPrint) + " times"+"\n");
                resultsText.append("\tTPR: " + TPR+"\n");
                resultsText.append("\tFPR: " + FPR+"\n");


                sumOfSites += topologyTally.get(theTreeToPrint);
            }



            //no longer doing this part

            System.out.println("Total Weighted Average TPR and FPR");
            System.out.println("\tTPR: " + averageTPR);
            System.out.println("\tFPR: " + averageFPR);
            resultsText.append("Total Weighted Average TPR and FPR" + "\n");
            resultsText.append("\tTPR: " + averageTPR + "\n");
            resultsText.append("\tFPR: " + averageFPR + "\n");
            */
            //endregion


            //do RF dist. calcs
            int a = 0;
            double totalRFDist = 0.0;
            int totalBP = numBP; //(will divide total by number of sites to get average rf dist)
            double[] eachSiteRF = new double[totalBP];

            int truthIndex = 0;
            Tree truthTree = Trees.readTree(truthText.get(truthIndex).split("]")[1]);
            int truthLength = Integer.parseInt((truthText.get(truthIndex).split("]")[0]).split("\\[")[1]);
            int predIndex = 0;
            Tree predTree = Trees.readTree(predictionText.get(predIndex).split("]")[1]);
            int predLength = Integer.parseInt((predictionText.get(predIndex).split("]")[0]).split("\\[")[1]);

            int debugTruthBP = truthLength;
            int debugPredBP = predLength;



            for(int bpI = 0; bpI < totalBP; bpI++){

                //if either has run out of length, move index forward and grab new trees
                if(truthLength<1){
                    truthIndex++;
                    truthTree = Trees.readTree(truthText.get(truthIndex).split("]")[1]);
                    truthLength = Integer.parseInt((truthText.get(truthIndex).split("]")[0]).split("\\[")[1]);
                    debugTruthBP+=truthLength;
                }
                if(predLength<1){
                    predIndex++;

                    predTree = Trees.readTree(predictionText.get(predIndex).split("]")[1]);
                    predLength = Integer.parseInt((predictionText.get(predIndex).split("]")[0]).split("\\[")[1]);
                    debugPredBP += predLength;
                }

                //calc RF distance between the two trees (withoutBL)
                //double theRFDistCalculated = (new RFDistanceBL(truthTree,predTree)).getDistance();
                SymmetricDifference aaaaa = new SymmetricDifference();
                aaaaa.computeDifference(truthTree, predTree, false); //FALSE - TREES R UNROOTED
                double theRFDistCalculated = aaaaa.getUnweightedAverage();

                //put rf dist in my array for potentially graphing later
                //can also use this to help a little bit with debugging and making sure its working right
                eachSiteRF[bpI] = theRFDistCalculated;

                totalRFDist+=theRFDistCalculated;

                //decrement remaining length for these trees
                truthLength--;
                predLength--;

            }

            //finally, finish RF dist calc with a simple division
            double averageRFDist = totalRFDist / totalBP;

            System.out.println("Average RF Distance");
            System.out.println("\tvalue: " + averageRFDist);

            resultsText.append("Average RF Distance value:" + averageRFDist + "\n");

            //System.out.println("DEBUG");
            //System.out.println("\ttruth total bp: " + debugTruthBP);
            //System.out.println("\tprediction total bp: " + debugPredBP);


            writeResults.write(resultsText.toString());
            writeResults.close();

            brTruth.close();
            brPrediction.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    //helper function taken from leo script temp to compare two MS files and get results on accuracy
    public static void compareMSFileVsOptimal(String trueTrees, String baseFolder, String whereToWriteResults, int numBP, int windowSize, int windowOffset, boolean doCentering){

        //arg1 - the truth (the actual output from MS)
        //arg2 - the predictions (how well did whatever method do?

        String whereTheWindowsAre = baseFolder+"/windows"+windowSize+"bp"+windowOffset+"os";

        //this will be populated after reading in all the bootstrap files, to be made into a graph after
        ArrayList<HashSet<String>> orderedWindowBootstrapTrees = new ArrayList<HashSet<String>>();

        //how many 'raxml bootstrap files' are there?
        File dir = new File(whereTheWindowsAre);
        File [] files = dir.listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.startsWith("RAxML_bootstrap");
            }
        });
        int howManyWindows = files.length;

        try {
            for (int i = 0; i < howManyWindows; i++) {

                HashSet<String> theseTrees = new HashSet<String>();

                //grab trees from file
                BufferedReader bsTreeFile = new BufferedReader(new FileReader(whereTheWindowsAre + "/RAxML_bootstrap.w" + i + ".fasta"));
                //String bestTreeString = "";
                String line;
                while ((line = bsTreeFile.readLine()) != null) {
                    theseTrees.add(line);
                    //bestTreeString+=line;
                }

                orderedWindowBootstrapTrees.add(theseTrees);

            }
        }catch(Exception e){
            e.printStackTrace();
        }

        ArrayList<Integer> stretchLengths = GetStretchLengths(windowSize, windowOffset, numBP, doCentering);

        int[] numsByRF = new int[11];

        //for bebugging rf distances
        String debugTreeString = "";

        try {

            //current window
            int currentWindow = 0;

            int predictionDebugBPInitial = 0;

            BufferedReader brTruth = new BufferedReader(new FileReader(trueTrees));
            //BufferedReader brPrediction = new BufferedReader(new FileReader(predictedTrees));

            BufferedWriter writeResults = new BufferedWriter(new FileWriter(whereToWriteResults));
            StringBuilder resultsText = new StringBuilder();

            ArrayList<String> truthText = new ArrayList<>();
            ArrayList<String> predictionText = new ArrayList<>();

            //read in the two files text and create map of the trees in the truth
            String line;

            while ((line = brTruth.readLine()) != null && line.contains("]")) {
                truthText.add(line);

            }

            //do RF dist. calcs
            int a = 0;
            double totalRFDist = 0.0;
            int totalBP = numBP; //(will divide total by number of sites to get average rf dist)
            double[] eachSiteRF = new double[totalBP];

            int truthIndex = 0;
            Tree truthTree = Trees.readTree(truthText.get(truthIndex).split("]")[1]);
            int truthLength = Integer.parseInt((truthText.get(truthIndex).split("]")[0]).split("\\[")[1]);


            //int predIndex = 0;
            //editting to optimality stuff
            //Tree predTree = Trees.readTree(predictionText.get(predIndex).split("]")[1]);
            HashSet<String> predTrees = orderedWindowBootstrapTrees.get(currentWindow);
            int predLength = stretchLengths.get(currentWindow);
                    //Integer.parseInt((predictionText.get(predIndex).split("]")[0]).split("\\[")[1]);

            int debugTruthBP = truthLength;
            int debugPredBP = predLength;


            for(int bpI = 0; bpI < totalBP; bpI++){

                //if either has run out of length, move index forward and grab new trees

                //manual inspection debugging - crap, was doing rooted RF dist b4, fixing now
                /*
                if((truthLength<1 || predLength<1) && !debugTreeString.equals(truthTree.toNewick())) {
                    System.out.println("TRUE TREE");
                    System.out.println(truthTree.toNewick());
                    //STITree a
                    truthTree.toNewickWD();
                    System.out.println("BS TREES");
                    for(String oneBSTree : predTrees) {
                    System.out.println(oneBSTree);
                        //with rooted=true
                        SymmetricDifference aaaaa = new SymmetricDifference();
                        aaaaa.computeDifference(truthTree, Trees.readTree(oneBSTree), true);
                        System.out.println(aaaaa.getUnweightedAverage());
                        //with rooted=false
                        aaaaa = new SymmetricDifference();
                        aaaaa.computeDifference(truthTree, Trees.readTree(oneBSTree), false);
                        System.out.println(aaaaa.getUnweightedAverage());
                        if(aaaaa.getUnweightedAverage()==0.0){
                            int stopHereSec = 0;
                        }
                    }

                    debugTreeString = truthTree.toNewick();
                    int stopHereSec2 = 0;
                }
                */

                if(truthLength<1){
                    truthIndex++;
                    truthTree = Trees.readTree(truthText.get(truthIndex).split("]")[1]);
                    truthLength = Integer.parseInt((truthText.get(truthIndex).split("]")[0]).split("\\[")[1]);
                    debugTruthBP+=truthLength;
                }
                if(predLength<1){
                    currentWindow++;

                    //predTree = Trees.readTree(predictionText.get(predIndex).split("]")[1]);
                    //predLength = Integer.parseInt((predictionText.get(predIndex).split("]")[0]).split("\\[")[1]);

                    predTrees = orderedWindowBootstrapTrees.get(currentWindow);
                    predLength = stretchLengths.get(currentWindow);


                    debugPredBP += predLength;
                }

                //calc RF distance between the two trees (withoutBL)
                //double theRFDistCalculated = (new RFDistanceBL(truthTree,predTree)).getDistance();
                //CHANGED FOR OPTIMALITY CALC : NOW CALCULATING BEST POSSIBLE RF DISTANCE BETWEEN THE truE
                //TREE AND ALL POSSIBLE WINDOW TREES

                double bestRFScore = 999999999999999.0;
                String debugBestTree = null;
                for(String singleBSTree : predTrees){
                    SymmetricDifference aaaaa = new SymmetricDifference();
                    aaaaa.computeDifference(truthTree, Trees.readTree(singleBSTree), false); //FALSE - TREES R UNROOTED

                    if(aaaaa.getUnweightedAverage() < bestRFScore){
                        bestRFScore = aaaaa.getUnweightedAverage();
                        debugBestTree = singleBSTree;
                    }
                    //double theRFDistCalculated = aaaaa.getUnweightedAverage();
                }

                //now tracking number of occurences of 0 rf distances as best tree and others:
                //(using the x2 this time so true RF dist.)
                int trueRFscore = (int)(bestRFScore*2);
                if(trueRFscore<11){
                    numsByRF[trueRFscore]++;
                }

                //debugging..why cant i get a 0 RF dist..?
                //if (trueRFscore == 3){
                //    int seeWhatsGoingOn = 0;
                //}

                //the old way for a single tree for the normal ms truth vs pred calculation
                /*
                SymmetricDifference aaaaa = new SymmetricDifference();
                aaaaa.computeDifference(truthTree, predTree, true);
                double theRFDistCalculated = aaaaa.getUnweightedAverage();
                */

                //put rf dist in my array for potentially graphing later
                //can also use this to help a little bit with debugging and making sure its working right
                eachSiteRF[bpI] = bestRFScore;

                totalRFDist+=bestRFScore;

                //decrement remaining length for these trees
                truthLength--;
                predLength--;

            }

            //finally, finish RF dist calc with a simple division
            double averageRFDist = totalRFDist / totalBP;

            System.out.println("OPTIMAL Average RF Distance");
            System.out.println("\tvalue: " + averageRFDist);
            resultsText.append("OPTIMAL Average RF Distance value:" + averageRFDist + "\n");

            //System.out.println("DEBUG");
            //System.out.println("\ttruth total bp: " + debugTruthBP);
            //System.out.println("\tprediction total bp: " + debugPredBP);


            writeResults.write(resultsText.toString());
            writeResults.close();

            brTruth.close();
            //brPrediction.close();


            //FINALLY, lets write out the new thing we are tracking (true RF scores)
            BufferedWriter writeResults2 = new BufferedWriter(new FileWriter(whereToWriteResults+"BestRFvsNumSites"));
            writeResults2.write("RF dists (true RF a.k.a. sym dist x2) for best BS trees vs. number of occurences for " + numBP + " bps\n");
            for(int rfI = 0; rfI<11; rfI++){
                writeResults2.write(rfI+" : "+numsByRF[rfI]+"\n");
            }
            writeResults2.close();


        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }



    public static void TestGraphLibrary(){

        // constructs a directed graph with the specified vertices and edges

        org.jgrapht.Graph<String, DefaultWeightedEdge> graph = new SimpleDirectedWeightedGraph<String, DefaultWeightedEdge>(DefaultWeightedEdge.class);
        graph.addVertex("vertex1");
        graph.addVertex("vertex2");
        graph.addVertex("vertex3");
        graph.addVertex("vertex4");
        graph.addVertex("vertex5");


        DefaultWeightedEdge e1 = graph.addEdge("vertex1", "vertex2");
        graph.setEdgeWeight(e1, 1);

        DefaultWeightedEdge e2 = graph.addEdge("vertex2", "vertex3");
        graph.setEdgeWeight(e2, 3);

        DefaultWeightedEdge e3 = graph.addEdge("vertex4", "vertex5");
        graph.setEdgeWeight(e3, 6);

        DefaultWeightedEdge e4 = graph.addEdge("vertex2", "vertex4");
        graph.setEdgeWeight(e4, 2);

        DefaultWeightedEdge e5 = graph.addEdge("vertex5", "vertex4");
        graph.setEdgeWeight(e5, 4);


        DefaultWeightedEdge e6 = graph.addEdge("vertex2", "vertex5");
        graph.setEdgeWeight(e6, 9);

        DefaultWeightedEdge e7 = graph.addEdge("vertex4", "vertex1");
        graph.setEdgeWeight(e7, 7);

        DefaultWeightedEdge e8 = graph.addEdge("vertex3", "vertex2");
        graph.setEdgeWeight(e8, 2);



        DefaultWeightedEdge e10 = graph.addEdge("vertex3", "vertex5");
        graph.setEdgeWeight(e10, 1);

        DefaultWeightedEdge e11 = graph.addEdge("vertex1", "vertex5");
        graph.setEdgeWeight(e11, 7);

        DefaultWeightedEdge e9 = graph.addEdge("vertex1", "vertex3");
        graph.setEdgeWeight(e9, 5);


        System.out.println("Shortest path from vertex1 to vertex5:");

        //List<DefaultWeightedEdge> sp;

        //these two lines are the only out of a ton of internet examples that actually compiles and works 0.0
        DijkstraShortestPath<String, DefaultWeightedEdge> runDijkstra = new DijkstraShortestPath<String, DefaultWeightedEdge>(graph);
        List<String> shortestPath = runDijkstra.getPath("vertex1","vertex5").getVertexList();

        //print the ordered list of vertices on shortest path
        System.out.println(shortestPath);

    }

    public static void MakeFolder(String folderName){
        File directory = new File(String.valueOf(folderName));
        if(!directory.exists()){
            directory.mkdir();
        }
    }



    private static void TestRFCalc() {

        compareTwoMSFiles("/Users/leo/rice/tmp/testRF/treefileConst","/Users/leo/rice/tmp/testRF/treefileConst2","/Users/leo/rice/tmp/testRF/same",10); //should be 0
        compareTwoMSFiles("/Users/leo/rice/tmp/testRF/treefileConst","/Users/leo/rice/tmp/testRF/treefile1Diff","/Users/leo/rice/tmp/testRF/1Diff",10); // should be 0.1 (2*1 / 10 /2)

        compareTwoMSFiles("/Users/leo/rice/tmp/testRF/treefileConst","/Users/leo/rice/tmp/testRF/treefile5Diff","/Users/leo/rice/tmp/testRF/5Diff",10); // should be 0.5 (2*5 / 10 /2)



    }






}
