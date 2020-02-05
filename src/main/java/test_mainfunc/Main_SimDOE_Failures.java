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
        String out_resFileSummary = programPath +File.separator+"OUTPUT"+File.separator+"BAP_DOE_sim.txt";

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

        writersum.println( "nbStage BN1 BN2 Sigma maxCT min_eta minCT max_eta");


        int[] BNpositions=new int[3];
        int[] BNpositions4={0,1,2,3};
        int[] BNpositions6={0,2,5,8};
        int[] BNpositions2={0,1,2};
        //here the DoE starts
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
            for(int etaFac=0; etaFac < myDOE.etaFactor.length;etaFac++){
                for(int BNfac: BNpositions){
                    for(int alfac=0; alfac< myDOE.alphafactor.length; alfac++){
                        for(int noBNctfac=0; noBNctfac< myDOE.noBNfactor.length;noBNctfac++){
                            for(int varfac=0; varfac< myDOE.varfactor.length;varfac++){
                                for(int ttffac=0; ttffac < myDOE.diffuprateFactor.length; ttffac++) {
                                    for (int ttrfac=0; ttrfac < myDOE.diffdownrateFactor.length; ttrfac++) {
                                        SerialLine mySystem = myDOE.getOneSystemConfiguration(Jfac, BNfac, alfac, noBNctfac, varfac);
                                        double meanBnCt = mySystem.CT[myDOE.BN1[BNfac]].getMean();

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

                                        double [] RepairVector = new double[myDOE.Njobs];
                                        double [] Machinept = new double[myDOE.Njobs];

                                        for (int j =0; j < myDOE.Jfactor[Jfac]; j++)
                                        {
                                            for(int row = 0; row < myDOE.Njobs; row++) {
                                            Machinept[row] = tij[row][j];
                                        }
                                            Failure myFailure = new Failure();
                                            myFailure.repairTimeGeneration(Machinept, RepairVector, myDOE.iiduprateFactor[ttffac],myDOE.iiddownrateFactor[ttrfac]);
                                            myFailure.ProctimeUpdateWithRep(Machinept, RepairVector);
                                            for(int row = 0; row < myDOE.Njobs; row++) {
                                                tij[row][j] = Machinept[row];
                                            }
                                        }

                                        mySystem.mySimulation = mySystem.new SimulationBAS(myDOE.Njobs, myDOE.W, tij);

                                        writersum.write(mySystem.nbStage + " " + myDOE.BN1[BNfac] + " " + myDOE.BN2[BNfac] + " " + myDOE.alphafactor[alfac] + " ");
                                        //simulation with LB of buffer
                                        for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                                            mySystem.buffer[j] = lB[j];
                                        }
                                        mySystem.mySimulation.simBAS(false);
                                        writersum.write(mySystem.OverallCT + " " + meanBnCt / mySystem.OverallCT + " ");

                                        //simulation with UB of buffer
                                        for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                                            mySystem.buffer[j] = uB[j];
                                        }
                                        mySystem.mySimulation.simBAS(false);
                                        writersum.println(mySystem.OverallCT + " " + meanBnCt / mySystem.OverallCT);
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
