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

        writersum.write( "nbStage TH DiffFailed diffUp Diffdown iidUp iidDown " +
                //"Alter5_numit Alter5_TotalTime Alter5_CplexTime Alter5_totalBuffer " +
                "Alter6Rev_numit Alter6Rev_TotalTime Alter6Rev_CplexTime Alter6Rev_totalBuffer Alter6Rev_singleBuffers");
        writersum.println();

        //here the DoE starts
        int totstage=0;
        int BNfac = 1;
        for(int r=1;r<=5;r++){
            for(int Jfac=0; Jfac < myDOE.Jfactor.length; Jfac++){
                if(Jfac == 0)
                    totstage=3;
                else totstage=3;
                for(int etaFac=0; etaFac < myDOE.etaFactor.length;etaFac++){
                        for(int alfac=0; alfac< myDOE.alphafactor.length; alfac++){
                            for(int noBNctfac=0; noBNctfac< myDOE.noBNfactor.length;noBNctfac++) {
                                for (int varfac = 0; varfac < myDOE.varfactor.length; varfac++) {
                                    for (int failfac = 0; failfac < myDOE.diffuprateFactor.length; failfac++) {
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

                                            //FAILURES
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
                                            double Availabilitydiff = (myDOE.diffuprateFactor[failfac]) / (myDOE.diffuprateFactor[failfac] + myDOE.diffdownrateFactor[failfac] );
                                            double meanBnCt = mySystem.CT[myDOE.difffailedstage[fm]].getMean() / Availabilitydiff;

                                                writersum.write(mySystem.nbStage + " " + myDOE.etaFactor[etaFac] + " " + myDOE.difffailedstage[fm] + " " + myDOE.diffuprateFactor[failfac] + " " + myDOE.diffdownrateFactor[failfac] + " "+ myDOE.iiduprateFactor[failfac] + " "+ myDOE.iiddownrateFactor[failfac] + " ");
                                            myDOE.tempinstance = "TH_" + myDOE.etaFactor[etaFac] + "_MF_" + myDOE.difffailedstage[fm] + "_Fr" + myDOE.diffuprateFactor[failfac] + "_Rr_" + myDOE.diffdownrateFactor[failfac];

                                            // Start optimization with Alter 6 reversed cut
                                            String out_resFile6Reversed = programPath + File.separator + "OUTPUT" + File.separator + "Out_" + myDOE.tempinstance + "_Alter6RevCut_" + (r) + ".txt";
                                            OutputStream outRes6RevCut = null;
                                            try {
                                                outRes6RevCut = new FileOutputStream(out_resFile6Reversed);
                                            } catch (FileNotFoundException e) {
                                                e.printStackTrace();
                                                System.exit(-1);
                                            }
                                            BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, myDOE.etaFactor[etaFac]/meanBnCt, lB, uB, myDOE.Njobs, myDOE.W);
                                            myReversedAlter6.writer = new PrintWriter(outRes6RevCut, true);

                                            Stopwatch totalAlter6RevTime = new Stopwatch();
                                            totalAlter6RevTime.start();
                                            try {
                                                myReversedAlter6.solveBAPWithIntModel(tij);
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
                                                writersum.write(mySystem.buffer[j] + ", ");
                                            }
                                            writersum.println();
                                            // End Optimization with Alter6 reversed cut

                                            //close single instance file
                                            try {
                                                outRes6RevCut.close();
                                            } catch (IOException e) {
                                                e.printStackTrace();
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
