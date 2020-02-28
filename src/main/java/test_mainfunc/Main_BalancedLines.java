package test_mainfunc;

import test_mainfunc.optimization.BendersIntModelAlter6;
import test_mainfunc.optimization.BendersIntModelAlter6ReversedCut;
import test_mainfunc.optimization.BendersStolletz;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.simulation.StochNum;
import test_mainfunc.util.Stopwatch;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;


public class Main_BalancedLines {

    public static void main(String argv[]) throws Exception {
        String programPath = System.getProperty("user.dir");

        //***   Output summary file   *********************************************************
        String out_resFileSummary = programPath + File.separator + "OUTPUT" + File.separator + "DoE_balanced_line_summary.txt";

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

        PrintWriter writersum = new PrintWriter(outRessummary, true);

        writersum.write("nbStage BNEfficiency BN1 BN2 Sigma noBN_CT " +
                //"Alter5_numit Alter5_TotalTime Alter5_CplexTime Alter5_totalBuffer " +
                "Alter6_numit Alter6_TotalTime Alter6_CplexTime Alter6_totalBuffer " +
                "Alter6Rev_numit Alter6Rev_TotalTime Alter6Rev_CplexTime Alter6Rev_totalBuffer " +
                "Stolletz_numit Stolletz_TotalTime Stolletz_CplexTime Stolletz_totalBuffer");
        writersum.println();


        // System configuration
        int Njobs = 1000000, W = 1, lb = 1, ub = 20;
        double etaM4 = 0.66, meanCT = 1.5;
        double[] etaM6 = {0.4, 0.5};

        SerialLine myM4System = new SerialLine();
        myM4System.nbStage = 4;
        myM4System.buffer = new int[4];
        myM4System.CT = new StochNum[5];
        for (int j = 1; j <= 4; j++) {
            myM4System.CT[j] = new StochNum();
            myM4System.CT[j].distribution = "LogNorm";
            myM4System.CT[j].para1 = meanCT;
            myM4System.CT[j].para2 = 0.85;
        }

        SerialLine myM6System = new SerialLine();
        myM6System.nbStage = 6;
        myM6System.buffer = new int[6];
        myM6System.CT = new StochNum[7];
        for (int j = 1; j <= 6; j++) {
            myM6System.CT[j] = new StochNum();
            myM6System.CT[j].distribution = "LogNorm";
            myM6System.CT[j].para1 = meanCT;
            myM6System.CT[j].para2 = 0.9;
        }

        ArrayList<SerialLine> experimentDesign = new ArrayList<>();
        //experimentDesign.add(myM4System);
        experimentDesign.add(myM6System);


        for (SerialLine mySystem : experimentDesign) {
            for(int k= 0;k<etaM6.length;k++){
                for (int r = 1; r <= 5; r++) {
                    double meanBnCt = mySystem.CT[1].getMean();
                    int[] lB = new int[mySystem.nbStage];
                    int[] uB = new int[mySystem.nbStage];
                    for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                        lB[j] = lb;
                        uB[j] = ub;
                    }
                    double eta = 0;
                    if (mySystem.nbStage == 4)
                        eta = etaM4;
                    if (mySystem.nbStage == 6)
                        eta = etaM6[k];

                    // sampling
                    double[][] tij = new double[Njobs + 1][mySystem.nbStage + 1];
                    int seed = (int) System.currentTimeMillis();
                    mySystem.procTimeGeneration(Njobs, tij, seed);

                    writersum.write(mySystem.nbStage + " " + eta + " " + 0 + " " + 0 + " " + mySystem.CT[1].para2 + " " + meanCT + " ");


                    // Start optimization with Alter 6
                    BendersIntModelAlter6 myAlter6 = new BendersIntModelAlter6(mySystem, eta / meanBnCt, lB, uB, Njobs, W);
                    myAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());

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
                    BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, eta / meanBnCt, lB, uB, Njobs, W);
                    myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());

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
                    if (totcap <= 18 || mySystem.nbStage == 4) {
                        BendersStolletz myStolletz = new BendersStolletz(mySystem, eta / meanBnCt, lB, uB, Njobs, W);

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
                    } else {
                        writersum.println();
                    }
                    //End of optimization with Stolletz
                }

            }
        }

        try {
            outRessummary.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}


