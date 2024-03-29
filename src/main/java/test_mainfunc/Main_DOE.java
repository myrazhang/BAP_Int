package test_mainfunc;

import test_mainfunc.optimization.BendersIntModelAlter5;
import test_mainfunc.optimization.BendersIntModelAlter6;
import test_mainfunc.optimization.BendersIntModelAlter6ReversedCut;
import test_mainfunc.optimization.BendersStolletz;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.util.Stopwatch;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;

import static java.lang.Math.exp;


public class Main_DOE {

    public static void main(String argv[]) throws Exception {
        String programPath = System.getProperty("user.dir");

        //***   Input files   *********************************************************
        String in_System = programPath + File.separator + "INPUT" + File.separator + "SerialLine_test_DoE_46.yaml";
        InputStream in_SystemFile = null;
        try {
            in_SystemFile = new FileInputStream(in_System);
        } catch (FileNotFoundException e) {
            System.err.println("Error opening 'System file'");
            System.exit(-1);
        }

        //***   Output summary file   *********************************************************
        String out_resFileSummary = programPath + File.separator + "OUTPUT" + File.separator + "SerialLine_test_DoE_46.yaml_decomposition_summary.txt";

        OutputStream outRessummary = null;
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
        SystemCombinationsForDOE myDOE = new SystemCombinationsForDOE(in_SystemFile);
        PrintWriter writersum = new PrintWriter(outRessummary, true);

        writersum.write("nbStage BNEfficiency BN1 BN2 Sigma noBN_CT " +
                //"Alter6_numit Alter6_TotalTime Alter6_CplexTime Alter6_totalBuffer " +
                "Alter6Rev_numit Alter6Rev_TotalTime Alter6Rev_CplexTime Alter6Rev_totalBuffer Buffer_allocation " +
                "Stolletz_numit Stolletz_TotalTime Stolletz_CplexTime Stolletz_totalBuffer");
        writersum.println();

        int[] BNpositions = null;
        int[] BNpositions6 = {6, 7, 8, 9, 10, 11};
        int[] BNpositions4 = {0, 1, 2, 3, 4, 5};

        //here the DoE starts
        for (int r = 1; r <= 1; r++) {
            for (int Jfac = 0; Jfac < myDOE.Jfactor.length; Jfac++) {
                if (myDOE.Jfactor[Jfac] == 4) {
                    BNpositions = BNpositions4;
                } else if (myDOE.Jfactor[Jfac] == 6) {
                    BNpositions = BNpositions6;
                }
                for (int etaFac = 0; etaFac < myDOE.etaFactor.length; etaFac++) {
                    for (int BNfac : BNpositions) {
                        for (int alfac = 0; alfac < myDOE.alphafactor.length; alfac++) {
                            for (int noBNctfac = 0; noBNctfac < myDOE.noBNfactor.length; noBNctfac++) {
                                for (int varfac = 0; varfac < myDOE.varfactor.length; varfac++) {

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

                                    writersum.write(mySystem.nbStage + " " + myDOE.etaFactor[etaFac] + " " + myDOE.BN1[BNfac] + " " + myDOE.BN2[BNfac] + " " + myDOE.alphafactor[alfac] + " " + myDOE.noBNfactor[noBNctfac] + " ");
                                    myDOE.tempinstance = "J" + mySystem.nbStage + "_TH_" + myDOE.etaFactor[etaFac] + "_BN_" + myDOE.BN1[BNfac] + myDOE.BN2[BNfac] + "_alpha_" + myDOE.alphafactor[alfac] + "_BNf_" + myDOE.noBNfactor[noBNctfac];

                                    /*// Start optimization with Alter 6
                                    String out_resFile6 = programPath +File.separator+"OUTPUT"+File.separator+"Out_"+ myDOE.tempinstance + "_Alter6_"+(r)+".txt";
                                    OutputStream outRes6= null;
                                    try {
                                        outRes6 = new FileOutputStream(out_resFile6);
                                    } catch (FileNotFoundException e) {
                                        e.printStackTrace();
                                        System.exit(-1);
                                    }
                                    BendersIntModelAlter6 myAlter6=new BendersIntModelAlter6(mySystem, myDOE.etaFactor[etaFac]/meanBnCt, lB, uB, myDOE.Njobs, myDOE.W);
                                    myAlter6.writer = new PrintWriter(outRes6, true);


                                    Stopwatch totalAlter6Time=new Stopwatch();
                                    totalAlter6Time.start();
                                    try{
                                        myAlter6.solveBAPWithIntModel(tij);
                                    }catch(Exception exc){exc.printStackTrace();}

                                    totalAlter6Time.stop();

                                    int totcap=0;
                                    for(int j=1;j<=mySystem.nbStage-1;j++){
                                        totcap = totcap+ mySystem.buffer[j];
                                    }
                                    writersum.write(myAlter6.numit + " " + df.format(totalAlter6Time.elapseTimeSeconds)+ " " +df.format(myAlter6.cplexTimeMeasure.elapseTimeSeconds)+ " "+totcap+" ");
                                    // End Optimization with Alter6*/

                                    // Start optimization with Alter 6 reversed cut
                                    String out_resFile6Reversed = programPath + File.separator + "OUTPUT" + File.separator + "Out_" + myDOE.tempinstance + "_Alter6RevCut_" + (r) + ".txt";
                                    OutputStream outRes6RevCut = null;
                                    try {
                                        outRes6RevCut = new FileOutputStream(out_resFile6Reversed);
                                    } catch (FileNotFoundException e) {
                                        e.printStackTrace();
                                        System.exit(-1);
                                    }
                                    BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, myDOE.etaFactor[etaFac] / meanBnCt, lB, uB, myDOE.Njobs, myDOE.W);
                                    myReversedAlter6.writer = new PrintWriter(outRes6RevCut, true);

                                    Stopwatch totalAlter6RevTime = new Stopwatch();
                                    totalAlter6RevTime.start();
                                    try {
                                        myReversedAlter6.solveWithDecomposition(tij);
                                    } catch (Exception exc) {
                                        exc.printStackTrace();
                                    }

                                    totalAlter6RevTime.stop();

                                    int totcap = 0;
                                    for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                                        totcap = totcap + mySystem.buffer[j];
                                    }
                                    writersum.write(myReversedAlter6.numit + " " + df.format(totalAlter6RevTime.elapseTimeSeconds) + " " + df.format(myReversedAlter6.cplexTimeMeasure.elapseTimeSeconds) + " " + totcap + " ");
                                    for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                                        writersum.write(mySystem.buffer[j] + ",");
                                    }
                                    writersum.write(" ");

                                    // End Optimization with Alter6 reversed cut


                                    // Start optimization with stolletz
                                    //if (totcap <= 18 || mySystem.nbStage == 4) {
                                        BendersStolletz myStolletz = new BendersStolletz(mySystem, myDOE.etaFactor[etaFac] / meanBnCt, lB, uB, myDOE.Njobs, myDOE.W);

                                        Stopwatch totalStolletzTime = new Stopwatch();
                                        totalStolletzTime.start();
                                        try {
                                            //myStolletz.solveBAPWithStolletz(tij, 0);
                                            myStolletz.solveWithDecomposition(tij);
                                        } catch (Exception exc) {
                                            exc.printStackTrace();
                                        }

                                        totalStolletzTime.stop();

                                        totcap = 0;
                                        for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                                            totcap = totcap + mySystem.buffer[j];
                                        }

                                        writersum.write(myStolletz.numit + " " + df.format(totalStolletzTime.elapseTimeSeconds) + " " + df.format(myStolletz.cplexTimeMeasure.elapseTimeSeconds) + " " + totcap);
                                        writersum.println();
                                    //} else {
                                    //    writersum.println();
                                    //}
                                    //End of optimization with Stolletz


                                    //close single instance file
                                    /*try {
                                        outRes6.close();
                                    } catch (IOException e) {
                                        e.printStackTrace();
                                    }*/
                                }
                            }
                        }//end Mfac
                    }//end BNfac
                }//end THfac
            }//end jfac
        }


        try {
            outRessummary.close();
            in_SystemFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}

