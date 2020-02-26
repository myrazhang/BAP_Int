package test_mainfunc;

import test_mainfunc.optimization.BendersIntModelAlter5;
import test_mainfunc.optimization.BendersIntModelAlter6;
import test_mainfunc.optimization.BendersIntModelAlter6ReversedCut;
import test_mainfunc.optimization.BendersStolletz;
import test_mainfunc.simulation.Failure;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.util.Stopwatch;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import static java.lang.Math.exp;


public class Main_DOE_Failures {

    public static void main(String argv[]) throws Exception {
        String programPath = System.getProperty("user.dir");

        //***   Input files   *********************************************************
        String in_System = programPath + File.separator+"INPUT"+File.separator+"SerialLine_test_DoE_Fail.yaml";
        InputStream in_SystemFile = null;
        try {
            in_SystemFile = new FileInputStream(in_System);
        } catch (FileNotFoundException e) {
            System.err.println("Error opening 'System file'");
            System.exit(-1);
        }


        //***   Output summary file   *********************************************************
        String out_resFileSummary = programPath +File.separator+"OUTPUT"+File.separator+"BAP_DOE_summary.txt";

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

        writersum.write( "nbStage BNEfficiency DiffFailed diffUp Diffdown " +
                //"Alter5_numit Alter5_TotalTime Alter5_CplexTime Alter5_totalBuffer " +
                "Alter6_numit Alter6_TotalTime Alter6_CplexTime Alter6_totalBuffer " +
                "Alter6Rev_numit Alter6Rev_TotalTime Alter6Rev_CplexTime Alter6Rev_totalBuffer " +
                "Stolletz_numit Stolletz_TotalTime Stolletz_CplexTime Stolletz_totalBuffer");
        writersum.println();

        //here the DoE starts
        int BNfac = 1;
        for(int r=1;r<=5;r++){
            for(int Jfac=0; Jfac < myDOE.Jfactor.length; Jfac++){
                for(int etaFac=0; etaFac < myDOE.etaFactor.length;etaFac++){
                        for(int alfac=0; alfac< myDOE.alphafactor.length; alfac++){
                            for(int noBNctfac=0; noBNctfac< myDOE.noBNfactor.length;noBNctfac++) {
                                for (int varfac = 0; varfac < myDOE.varfactor.length; varfac++) {
                                    for (int ttffac = 0; ttffac < myDOE.diffuprateFactor.length; ttffac++) {
                                        for (int ttrfac = 0; ttrfac < myDOE.diffdownrateFactor.length; ttrfac++) {
                                            for (int fm = 0; fm < myDOE.difffailedstage.length; fm++) {
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

                                                //FAILURES
                                                double[] Machinept = new double[myDOE.Njobs];
                                                Failure myFailure = new Failure();
                                                for (int j = 1; j < myDOE.Jfactor[Jfac]; j++) {
                                                    for (int i = 1; i < myDOE.Njobs; i++) {
                                                        Machinept[i] = tij[i][j];
                                                    }
                                                    if (j == myDOE.difffailedstage[fm]) {
                                                        myFailure.repairTimeGeneration(Machinept, myDOE.diffuprateFactor[ttffac], myDOE.diffdownrateFactor[ttrfac]);
                                                    } else {
                                                        myFailure.repairTimeGeneration(Machinept, myDOE.iiduprateFactor[ttffac], myDOE.iiddownrateFactor[ttrfac]);
                                                    }
                                                    for (int i = 1; i < myDOE.Njobs; i++) {
                                                        tij[i][j] = tij[i][j] + myFailure.repairTimeSamples[i];
                                                    }
                                                }
                                                double Availabilitydiff = (1 / myDOE.diffuprateFactor[ttffac]) / ((1 / myDOE.diffuprateFactor[ttffac]) + (1 / myDOE.diffdownrateFactor[ttrfac]));
                                                double meanBnCt = mySystem.CT[myDOE.difffailedstage[fm]].getMean() / Availabilitydiff;

                                                writersum.write(mySystem.nbStage + " " + myDOE.etaFactor[etaFac] + " " + myDOE.difffailedstage[fm] + " " + myDOE.diffuprateFactor[ttffac] + " " + myDOE.diffdownrateFactor[ttrfac] + " ");
                                                myDOE.tempinstance = "TH_" + myDOE.etaFactor[etaFac] + "_MF_" + myDOE.difffailedstage[fm] + "_Fr" + myDOE.diffuprateFactor[ttffac] + "_Rr_" + myDOE.diffdownrateFactor[ttrfac];

                                                // Start optimization with Alter 6
                                                String out_resFile6 = programPath + File.separator + "OUTPUT" + File.separator + "Out_" + myDOE.tempinstance + "_Alter6_" + (r) + ".txt";
                                                OutputStream outRes6 = null;
                                                try {
                                                    outRes6 = new FileOutputStream(out_resFile6);
                                                } catch (FileNotFoundException e) {
                                                    e.printStackTrace();
                                                    System.exit(-1);
                                                }
                                                BendersIntModelAlter6 myAlter6 = new BendersIntModelAlter6(mySystem, myDOE.etaFactor[etaFac] / meanBnCt, lB, uB, myDOE.Njobs, myDOE.W);
                                                myAlter6.writer = new PrintWriter(outRes6, true);


                                                Stopwatch totalAlter6Time = new Stopwatch();
                                                totalAlter6Time.start();
                                                try {
                                                    myAlter6.solveBAPWithIntModel(tij);
                                                } catch (Exception exc) {
                                                    exc.printStackTrace();
                                                }

                                                totalAlter6Time.stop();

                                                int totcap = 0;
                                                for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                                                    totcap = totcap + mySystem.buffer[j];
                                                }
                                                writersum.write(myAlter6.numit + " " + df.format(totalAlter6Time.elapseTimeSeconds) + " " + df.format(myAlter6.cplexTimeMeasure.elapseTimeSeconds) + " " + totcap + " ");
                                                // End Optimization with Alter6

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
                                                    myReversedAlter6.solveBAPWithIntModel(tij);
                                                } catch (Exception exc) {
                                                    exc.printStackTrace();
                                                }

                                                totalAlter6RevTime.stop();

                                                totcap = 0;
                                                for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                                                    totcap = totcap + mySystem.buffer[j];
                                                }
                                                writersum.write(myReversedAlter6.numit + " " + df.format(totalAlter6RevTime.elapseTimeSeconds) + " " + df.format(myReversedAlter6.cplexTimeMeasure.elapseTimeSeconds) + " " + totcap + " ");
                                                // End Optimization with Alter6 reversed cut


                                                // Start optimization with stolletz
                                                BendersStolletz myStolletz = new BendersStolletz(mySystem, myDOE.etaFactor[etaFac] / meanBnCt, lB, uB, myDOE.Njobs, myDOE.W);

                                                Stopwatch totalStolletzTime = new Stopwatch();
                                                totalStolletzTime.start();
                                                try {
                                                    myStolletz.solveBAPWithStolletz(tij);
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
                                                //End of optimization with Stolletz

                                                //close single instance file
                                                try {
                                                    outRes6.close();
                                                } catch (IOException e) {
                                                    e.printStackTrace();
                                                }
                                            }
                                        }
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

////////////////////////////////////////////////////////////////////////
////////CODE TO CHECK THE CORRECTNESS OF THE DESIGN OF EXPERIMENT///////
////////////////////////////////////////////////////////////////////////
/*
writer.println("ist: "+numistances);
        writer.println("Jfac THfac BNfac alfac noBNctfac varfac ARE: ");
        writer.println(Jfac + " " + THfac + " " + BNfac + " " + alfac + " " + noBNctfac + " " + varfac);
        writer.write("J: " + mySystem.NbStage + " and ");
        writer.write("Buffer length: " + mySystem.Buffer.length+ " and ");
        writer.write("TH star: " + mySystem.THstar+ " and ");
        writer.println("tij factors are: ");
        writer.write("alpha: "+ mySystem.alphafactor[alfac] + "and ");
        writer.write("beta: "+ mySystem.betafactor[alfac] + "and ");
        writer.println("var: "+ mySystem.varfactor[varfac] + "and ");
        writer.println("BN positions are: " + mySystem.BN1[BNfac] + " and "+ mySystem.BN2[BNfac]);
        writer.println("tij parameters are: ");
        for (int j = 0; j < mySystem.NbStage; j++)
        {
        writer.println( "j: "+ j + ": " + mySystem.CT[j].distribution + " " + mySystem.CT[j].Para1 + " " +mySystem.CT[j].Para2 + " " +mySystem.CT[j].Para3 + " " +mySystem.CT[j].Para4);
        }

        writer.println("------------------");*/
