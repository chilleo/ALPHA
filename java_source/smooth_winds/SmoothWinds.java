package edu.rice.cs.bioinfo.programs.phylonet.util;

import java.util.Date;

import static edu.rice.cs.bioinfo.programs.phylonet.LeoScript_FullPaperCwR.*;
import org.jgrapht.*;
import org.jgrapht.alg.shortestpath.DijkstraShortestPath;
import org.jgrapht.graph.*;

public class SmoothWinds {


    public static void main(String[] args) {



        String sequenceFullPath = args[0];
        String seqFileName = sequenceFullPath.split("/")[sequenceFullPath.split("/").length-1];
        String folderContainingSequence = sequenceFullPath.replace(seqFileName,"");

        String genomeLengthString = args[1];
        int genomeLength = Integer.parseInt(genomeLengthString);

        String windowSizeString = args[2];
        String windowOffsetString = args[3];
        int iSize = Integer.parseInt(windowSizeString);
        int jOff = Integer.parseInt(windowOffsetString);

        //SW stuff

        //params
        //int genomeLength = 10000;


        String baseFolder = folderContainingSequence+"/SmoothWinds"; //base folder where the exp. happens
        MakeFolder(baseFolder); //create this run's folder
        //MakeFolder(baseFolder + "/results"); //create folder for results


        //System.out.println("wind " + iSize + " off "+ jOff);

        //can do this one time
        WindowerWithSeqPath(baseFolder, sequenceFullPath, iSize, jOff);

        //must have raxmlHPC installed and in PATH
        RunWindowedRAXMLNoTPrint(baseFolder, iSize, jOff);

        boolean doCentering = true;

        //do i need to do process raxml? (i think this was just for labelling trees solely based off raxml ML trees
        //ProcessRAXML(baseFolder, iSize, jOff,genomeLength,doCentering);

        SmoothWindsLessMemNoTPrint(baseFolder, iSize, jOff, genomeLength, false, doCentering);

        //System.out.println("less mem done ws "+iSize+" off "+jOff);


        int debugTesting = 0;



    }

}
