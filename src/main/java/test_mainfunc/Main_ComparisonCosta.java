package test_mainfunc;

import test_mainfunc.optimization.BendersIntModelAlter6ReversedCut;
import test_mainfunc.optimization.BendersStolletz;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.util.Stopwatch;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.Scanner;


public class Main_ComparisonCosta {

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
        String out_resFileSummary = programPath + File.separator + "OUTPUT" + File.separator + "Comparison_Costa_summary.txt";

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

        writersum.write("Id nbStage TargetThp BN1 BN2 Sigma " +
                "Alter6Rev_numit Alter6Rev_TotalTime Alter6Rev_CplexTime Alter6Rev_totalBuffer Buffer_allocation "
        );
        writersum.println();


        int nbJobs = 10000;

        for (int r = 2; r <= 5; r++) {

            // ************************************** System 1 ********************************************
            {
                writersum.write("1 ");
                int Jfac = 0;
                int BNfac = 0;
                int etaFac = 0;
                int alfac = 0;
                String out_pt = programPath + File.separator + "OUTPUT" + File.separator + "pt_1_" + (r) + ".txt";

                SerialLine mySystem = myDOE.getOneSystemConfiguration(Jfac, BNfac, alfac, 0, 0); // variance  factor (always 0)
                double meanBnCt = mySystem.CT[myDOE.BN1[BNfac]].getMean();
                int[] lB = new int[mySystem.nbStage];
                int[] uB = new int[mySystem.nbStage];
                for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                    lB[j] = myDOE.Lj;
                    uB[j] = myDOE.Uj;
                }

                // sampling
                double[][] tij = new double[nbJobs + 1][mySystem.nbStage + 1];
                /*
                InputStream inpt = null;
                try {
                    inpt = new FileInputStream(out_pt);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
                Scanner ptScan = new Scanner(inpt);
                for(int i=1;i<=nbJobs;i++){
                    for(int j=1;j<=mySystem.nbStage;j++) tij[i][j] = ptScan.nextDouble();
                }*/
                int seed = (int) System.currentTimeMillis();
                mySystem.procTimeGeneration(nbJobs, tij, seed);

                writersum.write(mySystem.nbStage + " " + myDOE.etaFactor[etaFac] / meanBnCt + " " + myDOE.BN1[BNfac] + " " + myDOE.BN2[BNfac] + " " + myDOE.alphafactor[alfac] + " ");
                BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, myDOE.etaFactor[etaFac] / meanBnCt, lB, uB, myDOE.Njobs, myDOE.W);
                myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());

                Stopwatch totalAlter6RevTime = new Stopwatch();
                totalAlter6RevTime.start();
                try {
                    myReversedAlter6.solveBAPWithIntModel(tij, false);
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

                OutputStream outpt = null;
                try {
                    outpt = new FileOutputStream(out_pt);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
                PrintWriter ptwriter = new PrintWriter(outpt, true);
                for (int i = 1; i <= nbJobs; i++) {
                    for (int j = 1; j <= mySystem.nbStage; j++) ptwriter.print(tij[i][j] + " ");
                    ptwriter.println();
                }

                try {
                    outpt.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }


            // ************************************** System 2 ********************************************
            {
                writersum.write("2 ");
                int Jfac = 0;
                int BNfac = 1;
                int etaFac = 0;
                int alfac = 0;
                String out_pt = programPath + File.separator + "OUTPUT" + File.separator + "pt_2_" + (r) + ".txt";
                /*InputStream inpt = null;
                try {
                    inpt = new FileInputStream(out_pt);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
                Scanner ptScan = new Scanner(inpt);*/

                SerialLine mySystem = myDOE.getOneSystemConfiguration(Jfac, BNfac, alfac, 0, 0); // variance  factor (always 0)
                double meanBnCt = mySystem.CT[myDOE.BN1[BNfac]].getMean();
                int[] lB = new int[mySystem.nbStage];
                int[] uB = new int[mySystem.nbStage];
                for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                    lB[j] = myDOE.Lj;
                    uB[j] = myDOE.Uj;
                }

                // sampling
                double[][] tij = new double[nbJobs + 1][mySystem.nbStage + 1];
                /*for(int i=1;i<=nbJobs;i++){
                    for(int j=1;j<=mySystem.nbStage;j++) tij[i][j] = ptScan.nextDouble();
                }*/
                int seed = (int) System.currentTimeMillis();
                mySystem.procTimeGeneration(nbJobs, tij, seed);

                writersum.write(mySystem.nbStage + " " + myDOE.etaFactor[etaFac] / meanBnCt + " " + myDOE.BN1[BNfac] + " " + myDOE.BN2[BNfac] + " " + myDOE.alphafactor[alfac] + " ");
                BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, myDOE.etaFactor[etaFac] / meanBnCt, lB, uB, myDOE.Njobs, myDOE.W);
                myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());

                Stopwatch totalAlter6RevTime = new Stopwatch();
                totalAlter6RevTime.start();
                try {
                    myReversedAlter6.solveBAPWithIntModel(tij, false);
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

                OutputStream outpt = null;
                try {
                    outpt = new FileOutputStream(out_pt);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
                PrintWriter ptwriter = new PrintWriter(outpt, true);
                for (int i = 1; i <= nbJobs; i++) {
                    for (int j = 1; j <= mySystem.nbStage; j++) ptwriter.print(tij[i][j] + " ");
                    ptwriter.println();
                }

                try {
                    outpt.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

            // ************************************** System 3 ********************************************
            {
                writersum.write("3 ");
                int Jfac = 0;
                int BNfac = 0;
                int etaFac = 1;
                int alfac = 0;
                String out_pt = programPath + File.separator + "OUTPUT" + File.separator + "pt_3_" + (r) + ".txt";
                /*InputStream inpt = null;
                try {
                    inpt = new FileInputStream(out_pt);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
                Scanner ptScan = new Scanner(inpt);*/

                SerialLine mySystem = myDOE.getOneSystemConfiguration(Jfac, BNfac, alfac, 0, 0); // variance  factor (always 0)
                double meanBnCt = mySystem.CT[myDOE.BN1[BNfac]].getMean();
                int[] lB = new int[mySystem.nbStage];
                int[] uB = new int[mySystem.nbStage];
                for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                    lB[j] = myDOE.Lj;
                    uB[j] = myDOE.Uj;
                }

                // sampling
                double[][] tij = new double[nbJobs + 1][mySystem.nbStage + 1];
                /*for(int i=1;i<=nbJobs;i++){
                    for(int j=1;j<=mySystem.nbStage;j++) tij[i][j] = ptScan.nextDouble();
                }*/
                int seed = (int) System.currentTimeMillis();
                mySystem.procTimeGeneration(nbJobs, tij, seed);

                writersum.write(mySystem.nbStage + " " + myDOE.etaFactor[etaFac] / meanBnCt + " " + myDOE.BN1[BNfac] + " " + myDOE.BN2[BNfac] + " " + myDOE.alphafactor[alfac] + " ");
                BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, myDOE.etaFactor[etaFac] / meanBnCt, lB, uB, myDOE.Njobs, myDOE.W);
                myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());

                Stopwatch totalAlter6RevTime = new Stopwatch();
                totalAlter6RevTime.start();
                try {
                    myReversedAlter6.solveBAPWithIntModel(tij, false);
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

                myReversedAlter6.writer.print(myReversedAlter6.cplex.getModel());

                OutputStream outpt = null;
                try {
                    outpt = new FileOutputStream(out_pt);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
                PrintWriter ptwriter = new PrintWriter(outpt, true);
                for (int i = 1; i <= nbJobs; i++) {
                    for (int j = 1; j <= mySystem.nbStage; j++) ptwriter.print(tij[i][j] + " ");
                    ptwriter.println();
                }
                try {
                    outpt.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

            // ************************************** System 4 ********************************************
            {
                writersum.write("4 ");
                int Jfac = 0;
                int BNfac = 1;
                int etaFac = 1;
                int alfac = 0;
                String out_pt = programPath + File.separator + "OUTPUT" + File.separator + "pt_4_" + (r) + ".txt";
                /*InputStream inpt = null;
                try {
                    inpt = new FileInputStream(out_pt);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
                Scanner ptScan = new Scanner(inpt);*/

                SerialLine mySystem = myDOE.getOneSystemConfiguration(Jfac, BNfac, alfac, 0, 0); // variance  factor (always 0)
                double meanBnCt = mySystem.CT[myDOE.BN1[BNfac]].getMean();
                int[] lB = new int[mySystem.nbStage];
                int[] uB = new int[mySystem.nbStage];
                for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                    lB[j] = myDOE.Lj;
                    uB[j] = myDOE.Uj;
                }

                // sampling
                double[][] tij = new double[nbJobs + 1][mySystem.nbStage + 1];
                /*or(int i=1;i<=nbJobs;i++){
                    for(int j=1;j<=mySystem.nbStage;j++) tij[i][j] = ptScan.nextDouble();
                }*/
                int seed = (int) System.currentTimeMillis();
                mySystem.procTimeGeneration(nbJobs, tij, seed);

                writersum.write(mySystem.nbStage + " " + myDOE.etaFactor[etaFac] / meanBnCt + " " + myDOE.BN1[BNfac] + " " + myDOE.BN2[BNfac] + " " + myDOE.alphafactor[alfac] + " ");
                BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, myDOE.etaFactor[etaFac] / meanBnCt, lB, uB, myDOE.Njobs, myDOE.W);
                myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());

                Stopwatch totalAlter6RevTime = new Stopwatch();
                totalAlter6RevTime.start();
                try {
                    myReversedAlter6.solveBAPWithIntModel(tij, false);
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

                OutputStream outpt = null;
                try {
                    outpt = new FileOutputStream(out_pt);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
                PrintWriter ptwriter = new PrintWriter(outpt, true);
                for (int i = 1; i <= nbJobs; i++) {
                    for (int j = 1; j <= mySystem.nbStage; j++) ptwriter.print(tij[i][j] + " ");
                    ptwriter.println();
                }

                try {
                    outpt.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

            // ************************************** System 5 ********************************************
            {
                writersum.write("5 ");
                int Jfac = 1;
                int BNfac = 6;
                int etaFac = 0;
                int alfac = 0;
                String out_pt = programPath + File.separator + "OUTPUT" + File.separator + "pt_5_" + (r) + ".txt";
                /*InputStream inpt = null;
                try {
                    inpt = new FileInputStream(out_pt);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
                Scanner ptScan = new Scanner(inpt);*/

                SerialLine mySystem = myDOE.getOneSystemConfiguration(Jfac, BNfac, alfac, 0, 0); // variance  factor (always 0)
                double meanBnCt = mySystem.CT[myDOE.BN1[BNfac]].getMean();
                int[] lB = new int[mySystem.nbStage];
                int[] uB = new int[mySystem.nbStage];
                for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                    lB[j] = myDOE.Lj;
                    uB[j] = myDOE.Uj;
                }

                // sampling
                double[][] tij = new double[nbJobs + 1][mySystem.nbStage + 1];
                /*for(int i=1;i<=nbJobs;i++){
                    for(int j=1;j<=mySystem.nbStage;j++) tij[i][j] = ptScan.nextDouble();
                }*/
                int seed = (int) System.currentTimeMillis();
                mySystem.procTimeGeneration(nbJobs, tij, seed);

                writersum.write(mySystem.nbStage + " " + myDOE.etaFactor[etaFac] / meanBnCt + " " + myDOE.BN1[BNfac] + " " + myDOE.BN2[BNfac] + " " + myDOE.alphafactor[alfac] + " ");
                BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, myDOE.etaFactor[etaFac] / meanBnCt, lB, uB, myDOE.Njobs, myDOE.W);
                myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());

                Stopwatch totalAlter6RevTime = new Stopwatch();
                totalAlter6RevTime.start();
                try {
                    myReversedAlter6.solveBAPWithIntModel(tij, false);
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

                OutputStream outpt = null;
                try {
                    outpt = new FileOutputStream(out_pt);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
                PrintWriter ptwriter = new PrintWriter(outpt, true);
                for (int i = 1; i <= nbJobs; i++) {
                    for (int j = 1; j <= mySystem.nbStage; j++) ptwriter.print(tij[i][j] + " ");
                    ptwriter.println();
                }

                try {
                    outpt.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

            // ************************************** System 6 ********************************************
            {
                writersum.write("6 ");
                int Jfac = 1;
                int BNfac = 7;
                int etaFac = 0;
                int alfac = 0;
                String out_pt = programPath + File.separator + "OUTPUT" + File.separator + "pt_6_" + (r) + ".txt";
                /*InputStream inpt = null;
                try {
                    inpt = new FileInputStream(out_pt);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
                Scanner ptScan = new Scanner(inpt);*/

                SerialLine mySystem = myDOE.getOneSystemConfiguration(Jfac, BNfac, alfac, 0, 0); // variance  factor (always 0)
                double meanBnCt = mySystem.CT[myDOE.BN1[BNfac]].getMean();
                int[] lB = new int[mySystem.nbStage];
                int[] uB = new int[mySystem.nbStage];
                for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                    lB[j] = myDOE.Lj;
                    uB[j] = myDOE.Uj;
                }

                // sampling
                double[][] tij = new double[nbJobs + 1][mySystem.nbStage + 1];
                /*for(int i=1;i<=nbJobs;i++){
                    for(int j=1;j<=mySystem.nbStage;j++) tij[i][j] = ptScan.nextDouble();
                }*/
                int seed = (int) System.currentTimeMillis();
                mySystem.procTimeGeneration(nbJobs, tij, seed);

                writersum.write(mySystem.nbStage + " " + myDOE.etaFactor[etaFac] / meanBnCt + " " + myDOE.BN1[BNfac] + " " + myDOE.BN2[BNfac] + " " + myDOE.alphafactor[alfac] + " ");
                BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, myDOE.etaFactor[etaFac] / meanBnCt, lB, uB, myDOE.Njobs, myDOE.W);
                myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());

                Stopwatch totalAlter6RevTime = new Stopwatch();
                totalAlter6RevTime.start();
                try {
                    myReversedAlter6.solveBAPWithIntModel(tij, false);
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

                OutputStream outpt = null;
                try {
                    outpt = new FileOutputStream(out_pt);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
                PrintWriter ptwriter = new PrintWriter(outpt, true);
                for (int i = 1; i <= nbJobs; i++) {
                    for (int j = 1; j <= mySystem.nbStage; j++) ptwriter.print(tij[i][j] + " ");
                    ptwriter.println();
                }
                try {
                    outpt.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

            // ************************************** System 7 ********************************************
            {
                writersum.write("7 ");
                int Jfac = 1;
                int BNfac = 6;
                int etaFac = 1;
                int alfac = 0;
                String out_pt = programPath + File.separator + "OUTPUT" + File.separator + "pt_7_" + (r) + ".txt";
                /*InputStream inpt = null;
                try {
                    inpt = new FileInputStream(out_pt);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
                Scanner ptScan = new Scanner(inpt);*/

                SerialLine mySystem = myDOE.getOneSystemConfiguration(Jfac, BNfac, alfac, 0, 0); // variance  factor (always 0)
                double meanBnCt = mySystem.CT[myDOE.BN1[BNfac]].getMean();
                int[] lB = new int[mySystem.nbStage];
                int[] uB = new int[mySystem.nbStage];
                for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                    lB[j] = myDOE.Lj;
                    uB[j] = myDOE.Uj;
                }

                // sampling
                double[][] tij = new double[nbJobs + 1][mySystem.nbStage + 1];
                /*for(int i=1;i<=nbJobs;i++){
                    for(int j=1;j<=mySystem.nbStage;j++) tij[i][j] = ptScan.nextDouble();
                }*/
                int seed = (int) System.currentTimeMillis();
                mySystem.procTimeGeneration(nbJobs, tij, seed);

                writersum.write(mySystem.nbStage + " " + myDOE.etaFactor[etaFac] / meanBnCt + " " + myDOE.BN1[BNfac] + " " + myDOE.BN2[BNfac] + " " + myDOE.alphafactor[alfac] + " ");
                BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, myDOE.etaFactor[etaFac] / meanBnCt, lB, uB, myDOE.Njobs, myDOE.W);
                myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());

                Stopwatch totalAlter6RevTime = new Stopwatch();
                totalAlter6RevTime.start();
                try {
                    myReversedAlter6.solveBAPWithIntModel(tij, false);
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

                OutputStream outpt = null;
                try {
                    outpt = new FileOutputStream(out_pt);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
                PrintWriter ptwriter = new PrintWriter(outpt, true);
                for (int i = 1; i <= nbJobs; i++) {
                    for (int j = 1; j <= mySystem.nbStage; j++) ptwriter.print(tij[i][j] + " ");
                    ptwriter.println();
                }
                try {
                    outpt.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

            // ************************************** System 8 ********************************************
            {
                writersum.write("8 ");
                int Jfac = 1;
                int BNfac = 7;
                int etaFac = 1;
                int alfac = 0;
                String out_pt = programPath + File.separator + "OUTPUT" + File.separator + "pt_8_" + (r) + ".txt";
                /*InputStream inpt = null;
                try {
                    inpt = new FileInputStream(out_pt);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
                Scanner ptScan = new Scanner(inpt);*/

                SerialLine mySystem = myDOE.getOneSystemConfiguration(Jfac, BNfac, alfac, 0, 0); // variance  factor (always 0)
                double meanBnCt = mySystem.CT[myDOE.BN1[BNfac]].getMean();
                int[] lB = new int[mySystem.nbStage];
                int[] uB = new int[mySystem.nbStage];
                for (int j = 0; j <= mySystem.nbStage - 1; j++) {
                    lB[j] = myDOE.Lj;
                    uB[j] = myDOE.Uj;
                }

                // sampling
                double[][] tij = new double[nbJobs + 1][mySystem.nbStage + 1];
                /*for(int i=1;i<=nbJobs;i++){
                    for(int j=1;j<=mySystem.nbStage;j++) tij[i][j] = ptScan.nextDouble();
                }*/
                int seed = (int) System.currentTimeMillis();
                mySystem.procTimeGeneration(nbJobs, tij, seed);

                writersum.write(mySystem.nbStage + " " + myDOE.etaFactor[etaFac] / meanBnCt + " " + myDOE.BN1[BNfac] + " " + myDOE.BN2[BNfac] + " " + myDOE.alphafactor[alfac] + " ");
                BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, myDOE.etaFactor[etaFac] / meanBnCt, lB, uB, myDOE.Njobs, myDOE.W);
                myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());

                Stopwatch totalAlter6RevTime = new Stopwatch();
                totalAlter6RevTime.start();
                try {
                    myReversedAlter6.solveBAPWithIntModel(tij, false);
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

                OutputStream outpt = null;
                try {
                    outpt = new FileOutputStream(out_pt);
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                    System.exit(-1);
                }
                PrintWriter ptwriter = new PrintWriter(outpt, true);
                for (int i = 1; i <= nbJobs; i++) {
                    for (int j = 1; j <= mySystem.nbStage; j++) ptwriter.print(tij[i][j] + " ");
                    ptwriter.println();
                }
                try {
                    outpt.close();
                } catch (IOException e) {
                    e.printStackTrace();
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

