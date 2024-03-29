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

public class Main_DOE_Templ {
    public static void main(String argv[]) throws Exception {
        String programPath = System.getProperty("user.dir");

        //***   Input files   *********************************************************
        String in_System = programPath + File.separator + "INPUT" + File.separator + "SerialLine_test_DoE_Templ.yaml";
        InputStream in_SystemFile = null;
        try {
            in_SystemFile = new FileInputStream(in_System);
        } catch (FileNotFoundException e) {
            System.err.println("Error opening 'System file'");
            System.exit(-1);
        }


        //***   Output summary file   *********************************************************
        String out_resFileSummary = programPath + File.separator + "OUTPUT" + File.separator + "BAP_DOE_Templ_summary.txt";

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

        writersum.write("nbStage TH* MTTF " +
                //"Alter5_numit Alter5_TotalTime Alter5_CplexTime Alter5_totalBuffer " +
//                "Alter6Rev_numit Alter6Rev_TotalTime Alter6Rev_CplexTime Alter6Rev_totalBuffer Alter6Rev_singleBuffers ");
                "Alter6Rev_numit Alter6Rev_TotalTime Alter6Rev_CplexTime Alter6Rev_totalBuffer Alter6Rev_singleBuffers ");
        writersum.println();

        //here the DoE starts
        int BNfac = 1;
        for (int r = 1; r <= 4; r++) {

            //////////LINE A
/*
            String out_pt = programPath + File.separator + "OUTPUT" + File.separator + "processingtime_A.txt";

            OutputStream outpt = null;
            try {
                outpt = new FileOutputStream(out_pt);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
                System.exit(-1);
            }
            // System read from INPUT file
            PrintWriter outptwrite = new PrintWriter(outpt, true);

            SerialLine mySystem = myDOE.getTemplA();
            int[] lB = new int[mySystem.nbStage];
            int[] uB = new int[mySystem.nbStage];
            for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                lB[j] = myDOE.Lj;
                uB[j] = myDOE.Uj;
            }
            // System.out.println("LB: " + lB[0] + " UB: " + uB[0]);
            // sampling
            double[][] tij = new double[myDOE.Njobs + 1][mySystem.nbStage + 1];
            int seed = (int) System.currentTimeMillis();
            mySystem.procTimeGeneration(myDOE.Njobs, tij, seed);

            //FAILURES
            double[] Machinept = new double[myDOE.Njobs + 1];
            double[] meantr = new double[]{0, 49, 51, 55, 28, 47, 28, 68, 35, 43, 33, 36, 31, 60, 39, 38, 48, 47, 55, 43};
            double[] meantf = new double[]{0, 2178.273, 2009, 3382.5, 3972, 2890.5, 1189.391, 5162.769, 1131.667, 1260.03, 8217, 11964, 7719, 1940, 3211, 1482, 4315.636, 2426.684, 6820, 1610.846};
            Failure myFailure = new Failure();
            for (int j = 1; j <= mySystem.nbStage; j++) {
                for (int i = 1; i <= myDOE.Njobs; i++) {
                    Machinept[i] = tij[i][j];
                }
                myFailure.repairTimeGeneration(Machinept, meantf[j], meantr[j]);
                for (int i = 1; i < myDOE.Njobs; i++) {
                    tij[i][j] = tij[i][j] + myFailure.repairTimeSamples[i];
                }
            }

            for (int i = 1; i <= myDOE.Njobs; i++) {
                for (int j = 1; j <= mySystem.nbStage; j++) {
                    outptwrite.write(tij[i][j] + " ");
                }
                outptwrite.println();
            }

            //find THstar from simulation
            int[] bbuffer = new int[]{0, 14, 9, 9, 16, 28, 23, 27, 8, 12, 28, 6, 9, 12, 8, 24, 29, 10, 13};

            mySystem.mySimulation = mySystem.new SimulationBAS(myDOE.Njobs, myDOE.W, tij);
            for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                mySystem.buffer[j] = bbuffer[j];
            }
            mySystem.mySimulation.simBAS(false);

            //double thstar  = mySystem.TH + 0.0001;
            double thstar = mySystem.TH;
            writersum.write("TH: " + thstar + " ");
            System.out.println("TH star: " + thstar + " ");

            writersum.write("line_A " + thstar + " ");
            myDOE.tempinstance = "lineA_TH_" + thstar;

            // Start optimization with Alter 6 reversed cut
            String out_resFile6Reversed = programPath + File.separator + "OUTPUT" + File.separator + "Out_" + myDOE.tempinstance + "_Alter6rev_" + (r) + ".txt";
            OutputStream outRes6RevCut = null;
            try {
                outRes6RevCut = new FileOutputStream(out_resFile6Reversed);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
                System.exit(-1);
            }

           BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, thstar, lB, uB, myDOE.Njobs, myDOE.W);
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
                writersum.write(mySystem.buffer[j] + ",");
            }
            writersum.println();

            //simulate with optimal buffer capacity
            for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                bbuffer[j] = mySystem.buffer[j];
            }
            mySystem.mySimulation = mySystem.new SimulationBAS(myDOE.Njobs, myDOE.W, tij);
            for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                mySystem.buffer[j] = bbuffer[j];
            }
            mySystem.mySimulation.simBAS(false);
            writersum.write(" TH: " + mySystem.TH + " ");
            System.out.println("TH reached: " + mySystem.TH + " ");
            writersum.println();
            // End Optimization with Alter6 reversed cut


            //close single instance file
            try {
                outRes6RevCut.close();
            } catch (IOException e) {
                e.printStackTrace();
            }

    */
            //////////////////////////////////////////////////////////
            //////////LINE B
/*
            String out_pt = programPath + File.separator + "OUTPUT" + File.separator + "processingtime_B.txt";

            OutputStream outpt = null;
            try {
                outpt = new FileOutputStream(out_pt);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
                System.exit(-1);
            }

            // System read from INPUT file
            PrintWriter outptwrite = new PrintWriter(outpt, true);

            SerialLine mySystem = myDOE.getTemplB();

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
            double[] Machinept = new double[myDOE.Njobs + 1];
            double[] meantr = new double[]{0, 29.69, 35.89, 34.52, 44.4, 27.270, 30.400, 54.970, 33.980, 0, 27.750, 47.530, 45.230, 61.280, 0, 35.430, 33.19, 214.19, 51.07, 51.64, 75.83, 45.39, 229.85, 229.85};
            double[] meantf = new double[]{0, 1716.781, 3553.11, 1346.28, 638.677, 3381.48, 749.087, 909.416, 1027.895, 1000, 442.589, 1272.748, 1921.292, 6066.72, 1000, 2178.945, 2984.083, 26621.56, 2380.835, 2407.408, 639.547, 221.61, 5242.769, 5242.769};
            Failure myFailure = new Failure();
            for (int j = 1; j <= mySystem.nbStage; j++) {
                if(meantr[j]>0){
                    for (int i = 1; i <= myDOE.Njobs; i++) {
                        Machinept[i] = tij[i][j];
                    }
                    myFailure.repairTimeGeneration(Machinept, meantf[j], meantr[j]);
                    for (int i = 1; i < myDOE.Njobs; i++) {
                        tij[i][j] = tij[i][j] + myFailure.repairTimeSamples[i];
                    }
                }
            }

            //print processing time sample path
            for (int i = 1; i <= myDOE.Njobs; i++) {
                for (int j = 1; j <= mySystem.nbStage; j++) {
                    outptwrite.write(tij[i][j] + " ");
                }
                outptwrite.println();
            }

            writersum.write("line_B " );
            myDOE.tempinstance = "lineB_TH_";

            // Start optimization with Alter 6 reversed cut
            String out_resFile6Reversed = programPath + File.separator + "OUTPUT" + File.separator + "Out_" + myDOE.tempinstance + "_Alter6rev_" + (r) + ".txt";
            OutputStream outRes6RevCut = null;
            try {
                outRes6RevCut = new FileOutputStream(out_resFile6Reversed);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
                System.exit(-1);
            }

            //find TH* from simulation
            int[] bbuffer = new int[]{0,  9, 12, 12, 9, 6, 5, 6, 10, 9, 30, 6, 9, 9, 6, 7, 8, 6, 6, 7, 31, 128,128};
            mySystem.mySimulation = mySystem.new SimulationBAS(myDOE.Njobs, myDOE.W, tij);
            for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                mySystem.buffer[j] = bbuffer[j];
            }
            mySystem.mySimulation.simBAS(false);

            double thstar  = mySystem.TH;
            writersum.write("TH: " + thstar + " ");
            System.out.println("TH star: " + thstar+ " ");

            //find optimal solution
            BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, thstar, lB, uB, myDOE.Njobs, myDOE.W);
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
                writersum.write(mySystem.buffer[j] + ",");
            }
            writersum.println();
            // End Optimization with Alter6 reversed cut

            //simulate with optimal buffer capacity
            for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                bbuffer[j] = mySystem.buffer[j];
            }
            mySystem.mySimulation = mySystem.new SimulationBAS(myDOE.Njobs, myDOE.W, tij);
            for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                mySystem.buffer[j] = bbuffer[j];
            }
            mySystem.mySimulation.simBAS(false);
            writersum.write(" TH: " + mySystem.TH + " ");
            writersum.println();
            System.out.println("TH reached: " + mySystem.TH + " ");

            //close single instance file
            try {
                outRes6RevCut.close();
            } catch (IOException e) {
                e.printStackTrace();
            }

*/
            ////////////////////////////////////////////////////
            //////////LINE C
/*
            //write proc time in a file
            String out_pt = programPath + File.separator + "OUTPUT" + File.separator + "processingtime_C_" + r +".txt";

            OutputStream outpt = null;
            try {
                outpt = new FileOutputStream(out_pt);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
                System.exit(-1);
            }

            // System read from INPUT file
            PrintWriter outptwrite = new PrintWriter(outpt, true);


            SerialLine mySystem = myDOE.getTemplC();
            int[] lB = new int[mySystem.nbStage];
            int[]uB = new int[mySystem.nbStage];
            for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                lB[j] = myDOE.Lj;
                uB[j] = myDOE.Uj;
            }
            // sampling
            double [][] tij = new double[myDOE.Njobs + 1][mySystem.nbStage + 1];
            int seed = (int) System.currentTimeMillis();
            mySystem.procTimeGeneration(myDOE.Njobs, tij, seed);

            //FAILURES
            double[] Machinept = new double[myDOE.Njobs + 1];
            double mmeantr = 300;
            double[] meantf = new double[]{0, 2427.273, 2200, 2100, 2007.692, 1200, 2007.762, 7200, 7200};
           Failure myFailure = new Failure();
            for (int j = 1; j <= mySystem.nbStage; j++) {
                for (int i = 1; i <= myDOE.Njobs; i++) {
                    Machinept[i] = tij[i][j];
                }
                myFailure.repairTimeGeneration(Machinept, meantf[j], mmeantr);
                for (int i = 1; i < myDOE.Njobs; i++) {
                    tij[i][j] = tij[i][j] + myFailure.repairTimeSamples[i];
                }
            }
            myDOE.tempinstance = "lineC_";

            for (int i = 1; i <= myDOE.Njobs; i++) {
                for (int j = 1; j <= mySystem.nbStage; j++) {
                    outptwrite.write(tij[i][j] + " ");
                }
                outptwrite.println();
            }

/*
            // open output file
            String out_resFile6Reversed = programPath + File.separator + "OUTPUT" + File.separator + "Out_" + myDOE.tempinstance + "_Alter6rev_" + (r) + ".txt";
            OutputStream outRes6RevCut = null;
            try {
                outRes6RevCut = new FileOutputStream(out_resFile6Reversed);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
                System.exit(-1);
            }

            //find THstar from simulation
            int[] bbuffer = new int[]{0,  21, 31, 20, 24, 19, 4, 12};

            mySystem.mySimulation = mySystem.new SimulationBAS(myDOE.Njobs, myDOE.W, tij);
            for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                mySystem.buffer[j] = bbuffer[j];
            }
            mySystem.mySimulation.simBAS(false);

            //double thstar  = mySystem.TH + 0.0001;
            double thstar  = mySystem.TH;
            writersum.write("TH: " + thstar + " ");
            System.out.println("TH star: " + thstar+ " ");


            // Start optimization with Alter 6 reversed cut
            BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, thstar, lB, uB, myDOE.Njobs, myDOE.W);
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
                writersum.write(mySystem.buffer[j] + ",");
            }
            writersum.println();
            // End Optimization with Alter6 reversed cut

            //simulate with optimal buffer capacity
            for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                bbuffer[j] = mySystem.buffer[j];
            }
            mySystem.mySimulation = mySystem.new SimulationBAS(myDOE.Njobs, myDOE.W, tij);
            for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                mySystem.buffer[j] = bbuffer[j];
            }
            mySystem.mySimulation.simBAS(false);
            writersum.write(" TH: " + mySystem.TH + " ");
            System.out.println("TH reached: " + mySystem.TH + " ");
            //close single instance file
            try {
                outRes6RevCut.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
*/

            ////////////////////////////////////////////////////
            //////////LINE D

            //write proc time in a file
            String out_pt = programPath + File.separator + "OUTPUT" + File.separator + "processingtime_D_" + r + ".txt";

            OutputStream outpt = null;
            try {
                outpt = new FileOutputStream(out_pt);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
                System.exit(-1);
            }

            // System read from INPUT file
            PrintWriter outptwrite = new PrintWriter(outpt, true);

            //design of the line
            SerialLine mySystem = myDOE.getTemplD();
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
            double[] Machinept = new double[myDOE.Njobs + 1];
            int mmeantr = 36;
            double[] meantf = new double[]{0, 1764, 1164, 3564, 1764, 3564, 1164, 1164, 1164, 1164, 1164, 1764, 1164, 1164, 3564};
            Failure myFailure = new Failure();
            for (int j = 1; j <= mySystem.nbStage; j++) {
                for (int i = 1; i <= myDOE.Njobs; i++) {
                    Machinept[i] = tij[i][j];
                }
                myFailure.repairTimeGeneration(Machinept, meantf[j], mmeantr);
                for (int i = 1; i < myDOE.Njobs; i++) {
                    tij[i][j] = tij[i][j] + myFailure.repairTimeSamples[i];
                }

            }

            for (int i = 1; i <= myDOE.Njobs; i++) {
                for (int j = 1; j <= mySystem.nbStage; j++) {
                    outptwrite.write(tij[i][j] + " ");
                }
                outptwrite.println();
            }
/*
            //find THstar from simulation
            int[] bbuffer = new int[]{0,8,25,1,2,16,7,32,1,8,20,9,21,16};

            mySystem.mySimulation = mySystem.new SimulationBAS(myDOE.Njobs, myDOE.W, tij);
            for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                mySystem.buffer[j] = bbuffer[j];
            }
            mySystem.mySimulation.simBAS(false);

            //double thstar  = mySystem.TH + 0.0001;
            double thstar  = mySystem.TH;
            writersum.write("TH: " + thstar + " ");
            System.out.println("TH star: " + thstar+ " ");

            writersum.write("line_D " + thstar + " ");
            myDOE.tempinstance = "lineD_TH_" + thstar;

            // Start optimization with Alter 6 reversed cut
            String out_resFile6Reversed = programPath + File.separator + "OUTPUT" + File.separator + "Out_" + myDOE.tempinstance + "_Alter6rev_" + (r) + ".txt";
            OutputStream outRes6RevCut = null;
            try {
                outRes6RevCut = new FileOutputStream(out_resFile6Reversed);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
                System.exit(-1);
            }

            /*BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, thstar, lB, uB, myDOE.Njobs, myDOE.W);
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
                writersum.write(mySystem.buffer[j] + ",");
            }
            writersum.println();
            // End Optimization with Alter6 reversed cut

            //simulate with optimal buffer capacity
            for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                bbuffer[j] = mySystem.buffer[j];
            }
            mySystem.mySimulation = mySystem.new SimulationBAS(myDOE.Njobs, myDOE.W, tij);
            for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                mySystem.buffer[j] = bbuffer[j];
            }
            mySystem.mySimulation.simBAS(false);
            writersum.write(" TH: " + mySystem.TH + " ");
            System.out.println("TH reached: " + mySystem.TH + " ");

            //close single instance file
            try {
                outRes6RevCut.close();
            } catch (IOException e) {
                e.printStackTrace();
            }

        } //end replicates
        try {
            outRessummary.close();
            in_SystemFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
*/
        }
    }
}