package test_mainfunc;

import test_mainfunc.optimization.BendersIntModelAlter6ReversedCut;
import test_mainfunc.simulation.Failure;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.util.Stopwatch;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.Scanner;

public class Main_DOE_Templ_SimValidation {
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
        String out_resFileSummary = programPath + File.separator + "OUTPUT" + File.separator + "BAP_DOE_Templ_Sim.txt";

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

        writersum.write("System TH singleBuffers ");
        writersum.println();

        //here the DoE starts
        int BNfac = 1;
        double[] thTarget = new double[]{0.027, 0.028, 0.029};
        //for (int thIndex = 0; thIndex <= 2; thIndex++) {
        for (int r = 1; r <= 1; r++) {

            ////////////////////////////////////////////////////
            //////////LINE C
            {
                SerialLine mySystem = myDOE.getTemplC();
                int[] lB = new int[mySystem.nbStage];
                int[] uB = new int[mySystem.nbStage];
                for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                    lB[j] = myDOE.Lj;
                    uB[j] = myDOE.Uj;
                }
                int nJobs = 500000;
                double[][] tij = new double[nJobs + 1][mySystem.nbStage + 1];
                // sampling
            /*String ptC = programPath + File.separator + "INPUT" + File.separator + "processingtime_C.txt";
            InputStream ptStream = null;
            try {
                ptStream = new FileInputStream(ptC);
            } catch (FileNotFoundException e) {
                System.err.println("Error opening 'processing time' file.");
                System.exit(-1);
            }
            Scanner ptReader = new Scanner(ptStream);
            for (int i = 1; i <= nJobs; i++) {
                for (int j = 1; j <= mySystem.nbStage; j++) {
                    tij[i][j] = ptReader.nextDouble();
                }
            }*/

                //FAILURES
                int seed = (int) System.currentTimeMillis();
                mySystem.procTimeGeneration(nJobs, tij, seed);
                double[] Machinept = new double[nJobs + 1];
                double mmeantr = 300;
                double[] meantf = new double[]{0, 2427.273, 2200, 2100, 2007.692, 1200, 2007.762, 7200, 7200};
                Failure myFailure = new Failure();
                for (int j = 1; j <= mySystem.nbStage; j++) {
                    for (int i = 1; i <= nJobs; i++) {
                        Machinept[i] = tij[i][j];
                    }
                    myFailure.repairTimeGeneration(Machinept, meantf[j], mmeantr);
                    for (int i = 1; i < nJobs; i++) {
                        tij[i][j] = tij[i][j] + myFailure.repairTimeSamples[i];
                    }
                }

                mySystem.buffer = new int[]{0, 12,13,11,21,14,6,5};
                mySystem.mySimulation = mySystem.new SimulationBAS(nJobs, myDOE.W, tij);
                mySystem.mySimulation.simBAS(false, 0);
                double thstar = mySystem.TH;
                writersum.print("line_C " + thstar + " ");
                for(int j=1;j<mySystem.nbStage;j++) writersum.print(mySystem.buffer[j]+",");
                writersum.println();
                System.out.println("Th of Line C is " +thstar);
                //myDOE.tempinstance = "lineC_TH_" + thstar;

                // Start optimization with Alter 6 reversed cut
                /*out_resFile6Reversed = programPath + File.separator + "OUTPUT" + File.separator + "Out_" + myDOE.tempinstance + "_Alter6rev_" + (r) + ".txt";
                outRes6RevCut = null;
                try {
                    outRes6RevCut = new FileOutputStream(out_resFile6Reversed);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }*/
            }


            //}
            ////////////////////////////////////////////////////
            //////////LINE D
            SerialLine mySystem = myDOE.getTemplD();
            int[] lB = new int[mySystem.nbStage];
            int[] uB = new int[mySystem.nbStage];
            for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                lB[j] = myDOE.Lj;
                uB[j] = myDOE.Uj;
            }

            // sampling
            int nJobs = 1000000;
            double[][] tij = new double[nJobs + 1][mySystem.nbStage + 1];
            /*String ptD = programPath + File.separator + "INPUT" + File.separator + "processingtime_D.txt";
            try {
                ptStream = new FileInputStream(ptD);
            } catch (FileNotFoundException e) {
                System.err.println("Error opening 'processing time' file.");
                System.exit(-1);
            }
            ptReader = new Scanner(ptStream);

            for (int i = 1; i <= nJobs; i++) {
                for (int j = 1; j <= mySystem.nbStage; j++) {
                    tij[i][j] = ptReader.nextDouble();
                }
            }*/

            int seed = (int) System.currentTimeMillis();
            mySystem.procTimeGeneration(nJobs, tij, seed);
            //FAILURES
                double[] Machinept = new double[nJobs + 1];
                int mmeantr = 36;
                double[] meantf = new double[]{0, 1764, 1164, 3564, 1764, 3564, 1164, 1164, 1164, 1164, 1164, 1764, 1164, 1164, 3564};
                Failure myFailure = new Failure();
                for (int j = 1; j <= mySystem.nbStage; j++) {
                    for (int i = 1; i <= nJobs; i++) {
                        Machinept[i] = tij[i][j];
                    }
                    myFailure.repairTimeGeneration(Machinept, meantf[j], mmeantr);
                    for (int i = 1; i < nJobs; i++) {
                        tij[i][j] = tij[i][j] + myFailure.repairTimeSamples[i];
                    }
                }

            mySystem.buffer = new int[]{0,1,1,2,2,6,7,2,2,1,2,3,1,1};
            mySystem.mySimulation = mySystem.new SimulationBAS(nJobs, myDOE.W, tij);
            mySystem.mySimulation.simBAS(false, 0);
            double thstar = mySystem.TH;
            //thstar = thTarget[0];
            System.out.println("Th of Line D is " +thstar);
            writersum.write("line_D " + thstar + " ");
            for(int j=1;j<mySystem.nbStage;j++) writersum.print(mySystem.buffer[j]+",");
            writersum.println();


            //////////LINE A
            mySystem = myDOE.getTemplA();
            lB = new int[mySystem.nbStage];
            uB = new int[mySystem.nbStage];
            for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                lB[j] = myDOE.Lj;
                uB[j] = myDOE.Uj;
            }
            System.out.println("LB: " + lB[0] + " UB: " + uB[0]);
            // sampling
            nJobs = 500000;
            tij = new double[nJobs + 1][mySystem.nbStage + 1];

            /*String ptA = programPath + File.separator + "INPUT" + File.separator + "processingtime_A.txt";
            FileInputStream ptStream = null;
            try {
                ptStream = new FileInputStream(ptA);
            } catch (FileNotFoundException e) {
                System.err.println("Error opening 'processing time' file.");
                System.exit(-1);
            }
            Scanner ptReader = new Scanner(ptStream);
            for (int i = 1; i <= nJobs; i++) {
                for (int j = 1; j <= mySystem.nbStage; j++) {
                    tij[i][j] = ptReader.nextDouble();
                }
            }*/

            //tij = new double[nJobs + 1][mySystem.nbStage + 1];
            seed = (int) System.currentTimeMillis();
            mySystem.procTimeGeneration(nJobs, tij, seed);

            //FAILURES
                Machinept = new double[nJobs + 1];
                double[] meantr = new double[]{0, 49, 51, 55, 28, 47, 28, 68, 35, 43, 33, 36, 31, 60, 39, 38, 48, 47, 55, 43};
                meantf = new double[]{0, 2178.273, 2009, 3382.5, 3972, 2890.5, 1189.391, 5162.769, 1131.667, 1260.03, 8217, 11964, 7719, 1940, 3211, 1482, 4315.636, 2426.684, 6820, 1610.846};
                myFailure = new Failure();
                for (int j = 1; j <= mySystem.nbStage; j++) {
                    for (int i = 1; i <= nJobs; i++) {
                        Machinept[i] = tij[i][j];
                    }
                    myFailure.repairTimeGeneration(Machinept, meantf[j], meantr[j]);
                    for (int i = 1; i < nJobs; i++) {
                        tij[i][j] = tij[i][j] + myFailure.repairTimeSamples[i];
                    }
                }

            mySystem.buffer = new int[]{0, 1,1,1,1,1,2,1,4,1,1,3,5,2,2,1,1,1,1};
            mySystem.mySimulation = mySystem.new SimulationBAS(nJobs, myDOE.W, tij);
            mySystem.mySimulation.simBAS(false, 0);
            thstar = mySystem.TH;
            //thstar = 1.68984;
            writersum.write("line_A " + thstar + " ");
            for(int j=1;j<mySystem.nbStage;j++) writersum.print(mySystem.buffer[j]+",");
            writersum.println();
            System.out.println("Th of Line A is " +thstar);


            //}


            //////////////////////////////////////////////////////////
            //////////LINE B
            mySystem = myDOE.getTemplB();

            lB = new int[mySystem.nbStage];
            uB = new int[mySystem.nbStage];
            for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                lB[j] = myDOE.Lj;
                uB[j] = myDOE.Uj;
            }

            // sampling
            nJobs = 500000;
            tij = new double[nJobs + 1][mySystem.nbStage + 1];

            /*String ptB = programPath + File.separator + "INPUT" + File.separator + "processingtime_B.txt";
            try {
                ptStream = new FileInputStream(ptB);
            } catch (FileNotFoundException e) {
                System.err.println("Error opening 'processing time' file.");
                System.exit(-1);
            }
            ptReader = new Scanner(ptStream);
            tij = new double[nJobs + 1][mySystem.nbStage + 1];
            for (int i = 1; i <= nJobs; i++) {
                for (int j = 1; j <= mySystem.nbStage; j++) {
                    tij[i][j] = ptReader.nextDouble();
                }
            }*/

            //tij = new double[nJobs + 1][mySystem.nbStage + 1];
            seed = (int) System.currentTimeMillis();
            mySystem.procTimeGeneration(nJobs, tij, seed);

            //FAILURES
                Machinept = new double[nJobs + 1];
                meantr = new double[]{0, 29.69, 35.89, 34.52, 44.4, 27.270, 30.400, 54.970, 33.980, 0, 27.750, 47.530, 45.230, 61.280, 0, 35.430, 33.19, 214.19, 51.07, 51.64, 75.83, 45.39, 229.85, 229.85};
                meantf = new double[]{0, 1716.781, 3553.11, 1346.28, 638.677, 3381.48, 749.087, 909.416, 1027.895, 1000, 442.589, 1272.748, 1921.292, 6066.72, 1000, 2178.945, 2984.083, 26621.56, 2380.835, 2407.408, 639.547, 221.61, 5242.769, 5242.769};
                myFailure = new Failure();
                for (int j = 1; j <= mySystem.nbStage; j++) {
                    if (meantr[j] > 0) {
                        for (int i = 1; i <= nJobs; i++) {
                            Machinept[i] = tij[i][j];
                        }
                        myFailure.repairTimeGeneration(Machinept, meantf[j], meantr[j]);
                        for (int i = 1; i < nJobs; i++) {
                            tij[i][j] = tij[i][j] + myFailure.repairTimeSamples[i];
                        }
                    }
                }

            mySystem.buffer = new int[]{0, 1,1,1,1,1,1,1,1,1,2,1,2,1,1,2,1,1,1,2,32,43,1};
            mySystem.mySimulation = mySystem.new SimulationBAS(nJobs, myDOE.W, tij);
            mySystem.mySimulation.simBAS(false, 0);
            thstar = mySystem.TH;
            //thstar = 0.13474;
            writersum.write("line_B " + thstar + " ");
            for(int j=1;j<mySystem.nbStage;j++) writersum.print(mySystem.buffer[j]+",");
            writersum.println();
            System.out.println("Th of Line B is " +thstar);

        }
        //}//end replicates
        try {
            outRessummary.close();
            in_SystemFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}