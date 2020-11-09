package test_mainfunc;

import test_mainfunc.optimization.BendersIntModelAlter6;
import test_mainfunc.optimization.BendersStolletz;
import test_mainfunc.simulation.Failure;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.util.Stopwatch;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;

public class Main_SimDOE_Failures {
    public static void main(String argv[]) throws Exception {
        String programPath = System.getProperty("user.dir");

        //***   Input files   *********************************************************
        String in_System = programPath + File.separator+"INPUT"+File.separator+"SerialLine_test_DoE_Fail_sim.yaml";
        InputStream in_SystemFile = null;
        try {
            in_SystemFile = new FileInputStream(in_System);
        } catch (FileNotFoundException e) {
            System.err.println("Error opening 'System file'");
            System.exit(-1);
        }


        //***   Output summary file   *********************************************************
        String out_resFileSummary = programPath +File.separator+"OUTPUT"+File.separator+"BAP_DOE_sim_Fail.txt";

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

        writersum.println( "Nstage DiffFailed diffUp Diffdown maxCT min_TH minCT max_TH");


        int[] BNpositions=new int[3];
        int[] BNpositions4={0,1,2,3};
        int[] BNpositions6={0,2,5,8};
        int[] BNpositions2={0,1,2};
        //here the DoE starts
        int totstage= 0;
        for(int Jfac=0; Jfac < myDOE.Jfactor.length; Jfac++){
            if (myDOE.Jfactor[Jfac]==2){
                BNpositions =BNpositions2;
            }
            else if (myDOE.Jfactor[Jfac]==4){
                BNpositions =BNpositions4;
            }
            else if (myDOE.Jfactor[Jfac]==6){
                BNpositions = BNpositions6;
            }
            if(Jfac == 0)
                totstage=2;
            else totstage=3;
            for(int etaFac=0; etaFac < myDOE.etaFactor.length;etaFac++){
                for(int BNfac: BNpositions){
                    for(int alfac=0; alfac< myDOE.alphafactor.length; alfac++){
                        for(int noBNctfac=0; noBNctfac< myDOE.noBNfactor.length;noBNctfac++){
                            for(int varfac=0; varfac< myDOE.varfactor.length;varfac++){
                                for(int failfac=0; failfac < myDOE.diffuprateFactor.length; failfac++) {
                                    for (int fm = 0; fm < totstage; fm++) {
                                        SerialLine mySystem = myDOE.getOneSystemConfiguration(Jfac, BNfac, alfac, noBNctfac, varfac);

                                        int[] lB = new int[mySystem.nbStage];
                                        int[] uB = new int[mySystem.nbStage];
                                        for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                                            lB[j] = myDOE.Lj;
                                            uB[j] = myDOE.Uj;
                                        }

                                        // sampling
                                        double[][] tij = new double[myDOE.Njobs + 1][mySystem.nbStage + 1];
                                        int seed = (int) System.currentTimeMillis();
                                        mySystem.procTimeGeneration(myDOE.Njobs, tij, seed);

                                        //FAILURES: it is assumed that the factor diffuprateFactor and iiduprateFactor have same length, the same for the repair factor
                                        double[] Machinept = new double[myDOE.Njobs+1];
                                        Failure myFailure = new Failure();
                                        for (int j = 1; j <= myDOE.Jfactor[Jfac]; j++) {
                                            for (int i = 1; i <= myDOE.Njobs; i++) {
                                                Machinept[i] = tij[i][j];
                                            }
                                            if (j == myDOE.difffailedstage[fm]) {
                                                myFailure.repairTimeGeneration(Machinept, myDOE.diffuprateFactor[failfac], myDOE.diffdownrateFactor[failfac]);
                                            } else {
                                                myFailure.repairTimeGeneration(Machinept, myDOE.iiduprateFactor[failfac], myDOE.iiddownrateFactor[failfac]);
                                            }
                                            for (int i = 1; i <= myDOE.Njobs; i++) {
                                                tij[i][j] = tij[i][j] + myFailure.repairTimeSamples[i];
                                            }
                                        }
                                        double Availabilitydiff = (myDOE.diffuprateFactor[failfac]) / (myDOE.diffuprateFactor[failfac] + myDOE.diffdownrateFactor[failfac]);
                                        double meanBnCt = mySystem.CT[myDOE.difffailedstage[fm]].getMean() / Availabilitydiff;
                                        mySystem.mySimulation = mySystem.new SimulationBAS(myDOE.Njobs, myDOE.W, tij);

                                        writersum.write( myDOE.Jfactor[Jfac] + " " + myDOE.difffailedstage[fm] + " " + myDOE.diffuprateFactor[failfac] + " " + myDOE.diffdownrateFactor[failfac] + " ");
                                        //simulation with LB of buffer
                                        for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                                            mySystem.buffer[j] = lB[j];
                                        }
                                        mySystem.mySimulation.simBAS(false,0);
                                        writersum.write(mySystem.OverallCT  + " " + meanBnCt / mySystem.OverallCT  + " ");

                                        //simulation with UB of buffer
                                        for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                                            mySystem.buffer[j] = uB[j];
                                        }
                                        mySystem.mySimulation.simBAS(false,0);
                                        writersum.println(mySystem.OverallCT + " " +  meanBnCt / mySystem.OverallCT);

                                    }
                                }
                            }
                        }
                    }//end Mfac
                }//end BNfac
            }//end THfac
        }//end jfac

        try {
            outRessummary.close();
            in_SystemFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }







}
