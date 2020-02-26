package test_mainfunc;

import test_mainfunc.optimization.BendersIntModelAlter6;
import test_mainfunc.optimization.BendersIntModelAlter6ReversedCut;
import test_mainfunc.optimization.BendersStolletz;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.util.Stopwatch;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;

import static java.lang.Math.ceil;


public class MainExperiment1_2 {

    public static void main(String argv[]) throws Exception {
        String programPath = System.getProperty("user.dir");

        //***   Input files   *********************************************************
        String in_System = programPath + File.separator+"INPUT"+File.separator+"SerialLine_experiment12_highEta.yaml";
        InputStream in_SystemFile = null;
        try {
            in_SystemFile = new FileInputStream(in_System);
        } catch (FileNotFoundException e) {
            System.err.println("Error opening 'System file'");
            e.printStackTrace();
            System.exit(-1);
        }



        //***   Output summary file   *********************************************************
        String out_resFileSummary = programPath +File.separator+"OUTPUT"+File.separator+"BAP_experiment12_higEta_summary.txt";

        OutputStream outRessummary= null;
        try {
            outRessummary = new FileOutputStream(out_resFileSummary);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(-1);
        }

        // Ouput format
        DecimalFormat df;
        df = new DecimalFormat("#.#####");
        df.setRoundingMode(RoundingMode.CEILING);

        // System read from INPUT file
        SystemCombinationsForDOE myDOE=new SystemCombinationsForDOE(in_SystemFile);
        PrintWriter writersum = new PrintWriter(outRessummary, true);

        writersum.write( "nbStage BNEfficiency BNposition Sigma noBN_CT " +
                "Alter6Rev_numit Alter6Rev_TotalTime Alter6Rev_CplexTime Alter6Rev_totalBuffer BAP");
        writersum.println();


        // Range of stage numbers under study
        int stageNumberLowerBound = myDOE.Jfactor[0];
        int stageNumberUpperBound = myDOE.Jfactor[1];

        //here the DoE starts
        for(int r=1;r<=1;r++){
            String[] bnPositions = {"ML"};
            for(String bn: bnPositions){
                for (int etaFac=0; etaFac < myDOE.etaFactor.length;etaFac++){
                    for(int alfac=0; alfac< myDOE.alphafactor.length; alfac++){
                        for(int noBNctfac=0; noBNctfac< myDOE.noBNfactor.length;noBNctfac++){
                            for(int varfac=0; varfac< myDOE.varfactor.length;varfac++){
                                for(int stageNumber = stageNumberLowerBound; stageNumber<=stageNumberUpperBound; stageNumber=stageNumber+=2){
                                    SerialLine mySystem=myDOE.getOneSystemConfiguration(stageNumber, bn, alfac, noBNctfac, varfac);
                                    int BN1=0;
                                    if (bn.equals("MM")){
                                        BN1 = stageNumber/2;
                                    }
                                    else if(bn.equals("ML")){
                                        BN1 = stageNumber/2+1;
                                    }
                                    else if(bn.equals("FF")){
                                        BN1 = 1;
                                    }
                                    double meanBnCt=mySystem.CT[BN1].getMean();

                                    int[] lB=new int[mySystem.nbStage];
                                    int[] uB=new int[mySystem.nbStage];
                                    for(int j=0;j<=mySystem.nbStage-1;j++){
                                        lB[j]=myDOE.Lj;
                                        uB[j]=myDOE.Uj;
                                    }

                                    // sampling
                                    int simLength = myDOE.Njobs * (int) Math.ceil(mySystem.nbStage/10.0);
                                    double[][] tij=new double[simLength+1][mySystem.nbStage+1];
                                    int seed =(int) System.currentTimeMillis();
                                    mySystem.procTimeGeneration(simLength,tij,seed);

                                    writersum.write( mySystem.nbStage + " "+myDOE.etaFactor[etaFac] +" "+bn +" "+myDOE.alphafactor[alfac]+" "+myDOE.noBNfactor[noBNctfac]+" ");
                                    myDOE.tempinstance = "J"+mySystem.nbStage+"_TH_"+myDOE.etaFactor[etaFac] +"_BN_"+bn+"_alpha_"+myDOE.alphafactor[alfac]+"_BNf_"+myDOE.noBNfactor[noBNctfac];

                                    // Start optimization with Alter 6 reversed cut
                                    BendersIntModelAlter6ReversedCut myReversedAlter6=new BendersIntModelAlter6ReversedCut(mySystem, myDOE.etaFactor[etaFac]/meanBnCt, lB, uB, simLength, myDOE.W);
                                    myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream(), true);

                                    Stopwatch totalAlter6RevTime=new Stopwatch();
                                    totalAlter6RevTime.start();
                                    try{
                                        myReversedAlter6.solveBAPWithIntModel(tij);
                                    }catch(Exception exc){exc.printStackTrace();}

                                    totalAlter6RevTime.stop();

                                    int totcap=0;
                                    for(int j=1;j<=mySystem.nbStage-1;j++){
                                        totcap = totcap+ mySystem.buffer[j];
                                    }
                                    writersum.write(myReversedAlter6.numit + " " + df.format(totalAlter6RevTime.elapseTimeSeconds)+ " " +df.format(myReversedAlter6.cplexTimeMeasure.elapseTimeSeconds)+ " "+totcap+" ");

                                    String BAP = String.valueOf(mySystem.buffer[1]);
                                    for(int j=2;j<=mySystem.nbStage-1;j++)
                                        BAP += ","+mySystem.buffer[j];
                                    writersum.println(BAP);
                                    // End Optimization with Alter6 reversed cut

                                    if (myReversedAlter6.cplexTimeMeasure.elapseTimeSeconds>17999){
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        try {
            outRessummary.close();
            in_SystemFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
