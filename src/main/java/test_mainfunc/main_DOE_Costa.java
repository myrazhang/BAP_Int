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


public class main_DOE_Costa {

    public static void main(String argv[]) throws Exception {
        String programPath = System.getProperty("user.dir");

        //***   Input files   *********************************************************
        String inPT = programPath + File.separator + "INPUT" + File.separator + "processingtime_D_2.txt";
        InputStream in_PTFile = null;

        try {
            in_PTFile = new FileInputStream(inPT);
        } catch (FileNotFoundException e) {
            System.err.println("Error opening 'Processing time file'");
            System.exit(-1);
        }


        //***   Output summary file   *********************************************************
        String out_resFileSummary = programPath +File.separator+"OUTPUT"+File.separator+"ResultsSummaryCosta.txt";

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

        //System configuration
        int Njobs = 10000, W = 1, lb = 0, ub = 40;

        double[] THcosta = {0.344827586};
        SerialLine mySystem = new SerialLine();
        mySystem.nbStage = 14;
        mySystem.buffer = new int[14];


        // System read from INPUT file
        PrintWriter writersum = new PrintWriter(outRessummary, true);
        writersum.write( "nbStage TH " +
                //"Alter5_numit Alter5_TotalTime Alter5_CplexTime Alter5_totalBuffer " +
                "Alter6Rev_numit Alter6Rev_TotalTime Alter6Rev_CplexTime Alter6Rev_totalBuffer Alter6Rev_singleBuffers"
                );
        writersum.println();

        // read processing time
        // read pt from file
        double[][] tij = new double[Njobs + 1][mySystem.nbStage + 1];
        int seed = (int) System.currentTimeMillis();

        BufferedReader buffer = new BufferedReader(new InputStreamReader(in_PTFile));
        String line;
        int row = 1;
        while ((line = buffer.readLine()) != null && row<(Njobs+1)) {
            String[] vals =  line.trim().split(" ");
            double[] tmp = new double[vals.length+1];
            for (int i=0;i<vals.length;i++){
                tmp[i] = Double.parseDouble(vals[i]);
            }
            tij[row]=tmp;
            if(row == 1)
                System.out.printf("t11" + tij[row][1]);
            row++;
        }

        for (int i =1; i<= Njobs; i++)
        {
            for(int j=mySystem.nbStage;j> 0; j--)
            {
                tij[i][j] = tij[i][j-1];
            }
            tij[i][0]=0;
        }

        //here the DoE starts

        for(int k= 0;k<THcosta.length;k++){
            for (int r = 1; r <= 1; r++) {
                int[] lB = new int[mySystem.nbStage];
                int[] uB = new int[mySystem.nbStage];
                for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                    lB[j] = lb;
                    uB[j] = ub;
                }
                double TH = THcosta[k];

                writersum.write(mySystem.nbStage + " " + TH  + " ");
                writersum.println();

                //Start Simulation with given buffer
                int[] bbuffer = new int[]{0,5,5,5,5,5,5,5,5,5,5,5,5,5};
                mySystem.mySimulation = mySystem.new SimulationBAS(Njobs, W, tij);
                for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                    mySystem.buffer[j] = bbuffer[j];
                }
                mySystem.mySimulation.simBAS(false);

                for (int i = W+1; i <= Njobs; i++){
                    for (int j = 1; j <= mySystem.nbStage; j++) {
                        writersum.write( mySystem.mySimulation.Dij[i][j] + " ");
                    }
                    writersum.println();
                }

                double thsim  = mySystem.TH;
                System.out.println("TH star: " + thsim+ " ");
/*

                // Start optimization with Alter 6 reversed cut
                BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, TH, lB, uB, Njobs, W);
                myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());

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
*/
                ///////////SIMULATE A SPECIFIC BUFFER DISTRIBUTION:
                /*mySystem.buffer[1]=5;
                mySystem.buffer[2]=4;
                mySystem.buffer[3]=6;
                mySystem.buffer[4]=5;

                mySystem.mySimulation = mySystem.new SimulationBAS(Njobs, W, tij);
                mySystem.mySimulation.simBAS(false);
                writersum.print("CT: " + mySystem.OverallCT+" " + "TH: " + mySystem.TH +" ");
                */
            }

        }


        try {
            outRessummary.close();
            in_PTFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}

