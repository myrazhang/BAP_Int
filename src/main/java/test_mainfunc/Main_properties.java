package test_mainfunc;

import test_mainfunc.optimization.BendersIntModelAlter6ReversedCut;
import test_mainfunc.optimization.BendersStolletz;
import test_mainfunc.simulation.Failure;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.simulation.StochNum;
import test_mainfunc.util.Stopwatch;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import static java.lang.Math.max;


public class Main_properties {

    public static void main(String argv[]) throws Exception {
        String programPath = System.getProperty("user.dir");

        //***   Output summary file   *********************************************************
        String out_resFileSummary = programPath + File.separator + "OUTPUT" + File.separator + "Properties_4.txt";

        OutputStream outRessummary = null;
        try {
            outRessummary = new FileOutputStream(out_resFileSummary);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(-1);
        }
        PrintWriter writersum = new PrintWriter(outRessummary, true);
        writersum.println("system" + "numit TotalTime totalBuffer BufferAllocation ");


        // Ouput format
        DecimalFormat df;
        df = new DecimalFormat("#.#####");
        df.setRoundingMode(RoundingMode.CEILING);


        // Systems
        int nbJobs = 1000000;
        int nbStages = 6;
        int[] lB = new int[]{1,1,1,1,1,1};
        int[] uB = new int[]{20,20,20,20,20,20};
        double unrelMachineEfficiency = 0.8;
        int[] referBuffer = new int[]{0,4,4,4,4,4};

        for(int r= 1;r<=5;r++){


            // Reliable lines
            /*HashMap<String, SerialLine> systems = new HashMap<>();
            systems.put("Balance", new SerialLine());
            {
                SerialLine theSystem = systems.get("Balance");
                theSystem.nbStage = nbStages;
                theSystem.buffer = new int[theSystem.nbStage];
                //save distribution information
                theSystem.CT = new StochNum[theSystem.nbStage+1];
                for (int j = 1; j <= theSystem.nbStage; j++)
                {
                    theSystem.CT[j]=new StochNum();
                    theSystem.CT[j].distribution= "LogNorm";
                    theSystem.CT[j].para1 = 2;
                    theSystem.CT[j].para2 = 0.85;
                }

            }

            systems.put("Bowl", new SerialLine());
            {
                SerialLine theSystem = systems.get("Bowl");
                theSystem.nbStage = 6;
                theSystem.buffer = new int[theSystem.nbStage];
                //save distribution information
                theSystem.CT = new StochNum[theSystem.nbStage+1];
                for (int j = 1; j <= theSystem.nbStage; j++)
                {
                    theSystem.CT[j]=new StochNum();
                    theSystem.CT[j].distribution= "LogNorm";
                    theSystem.CT[j].para2 = 0.85;
                }
                theSystem.CT[1].para1 = 1.5;
                theSystem.CT[2].para1 = 2;
                theSystem.CT[3].para1 = 2.5;
                theSystem.CT[4].para1 = 2.5;
                theSystem.CT[5].para1 = 2;
                theSystem.CT[6].para1 = 1.5;
            }

            systems.put("InvertedBowl", new SerialLine());
            {
                SerialLine theSystem = systems.get("InvertedBowl");
                theSystem.nbStage = 6;
                theSystem.buffer = new int[theSystem.nbStage];
                //save distribution information
                theSystem.CT = new StochNum[theSystem.nbStage+1];
                for (int j = 1; j <= theSystem.nbStage; j++)
                {
                    theSystem.CT[j]=new StochNum();
                    theSystem.CT[j].distribution= "LogNorm";
                    theSystem.CT[j].para2 = 0.85;
                }
                theSystem.CT[1].para1 = 2.5;
                theSystem.CT[2].para1 = 2;
                theSystem.CT[3].para1 = 1.5;
                theSystem.CT[4].para1 = 1.5;
                theSystem.CT[5].para1 = 2;
                theSystem.CT[6].para1 = 2.5;
            }

            systems.put("Increase", new SerialLine());
            {
                SerialLine theSystem = systems.get("Increase");
                theSystem.nbStage = 6;
                theSystem.buffer = new int[theSystem.nbStage];
                //save distribution information
                theSystem.CT = new StochNum[theSystem.nbStage+1];
                for (int j = 1; j <= theSystem.nbStage; j++)
                {
                    theSystem.CT[j]=new StochNum();
                    theSystem.CT[j].distribution= "LogNorm";
                    theSystem.CT[j].para2 = 0.85;
                }
                theSystem.CT[1].para1 = 1.5;
                theSystem.CT[2].para1 = 2;
                theSystem.CT[3].para1 = 2.5;
                theSystem.CT[4].para1 = 3;
                theSystem.CT[5].para1 = 3.5;
                theSystem.CT[6].para1 = 4;
            }

            systems.put("Decrease",new SerialLine());
            {
                SerialLine theSystem = systems.get("Decrease");
                theSystem.nbStage = 6;
                theSystem.buffer = new int[theSystem.nbStage];
                //save distribution information
                theSystem.CT = new StochNum[theSystem.nbStage+1];
                for (int j = 1; j <= theSystem.nbStage; j++)
                {
                    theSystem.CT[j]=new StochNum();
                    theSystem.CT[j].distribution= "LogNorm";
                    theSystem.CT[j].para2 = 0.85;
                }
                theSystem.CT[1].para1 = 4;
                theSystem.CT[2].para1 = 3.5;
                theSystem.CT[3].para1 = 3;
                theSystem.CT[4].para1 = 2.5;
                theSystem.CT[5].para1 = 2;
                theSystem.CT[6].para1 = 1.5;
            }

            for(String str: systems.keySet()){
                SerialLine mySystem = systems.get(str);
                double[][] tij = new double[nbJobs + 1][mySystem.nbStage + 1];
                int seed = (int) System.currentTimeMillis();
                mySystem.procTimeGeneration(nbJobs, tij, seed);
                writersum.write(str+" ");
                mySystem.buffer=referBuffer.clone();
                mySystem.mySimulation = mySystem.new SimulationBAS(nbJobs, 1, tij);
                mySystem.mySimulation.simBAS(false, 0);
                double thstar = mySystem.TH;

                BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(mySystem, thstar, lB, uB, nbJobs, 1);
                myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());

                Stopwatch totalAlter6RevTime = new Stopwatch();
                totalAlter6RevTime.start();
                try {
                    myReversedAlter6.solveBAPWithIntModel(tij,true);
                } catch (Exception exc) {
                    exc.printStackTrace();
                }
                totalAlter6RevTime.stop();

                int totcap = 0;
                for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                    totcap = totcap + mySystem.buffer[j];
                }
                writersum.write(myReversedAlter6.numit + " " + df.format(totalAlter6RevTime.elapseTimeSeconds) + " " + totcap + " ");
                for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                    writersum.write(mySystem.buffer[j] + ",");
                }
                writersum.println();
            }*/

            // mixed lines
            /*HashMap<String, int[]> unreMachinePosition = new HashMap<>();
            unreMachinePosition.put("mix1",new int[]{1,3,4,6});
            unreMachinePosition.put("mix2",new int[]{3,4});
            unreMachinePosition.put("mix3",new int[]{1,3,5});
            unreMachinePosition.put("mix4",new int[]{1,3,6});

            for(String str: unreMachinePosition.keySet()){
                SerialLine theSystem = new SerialLine();
                theSystem.nbStage = 6;
                theSystem.buffer = new int[theSystem.nbStage];
                theSystem.CT = new StochNum[theSystem.nbStage+1];

                HashSet<Integer> unre = new HashSet<>();
                for(int j: unreMachinePosition.get(str)) unre.add(j);
                StochNum sn=new StochNum();
                sn.distribution= "LogNorm";
                sn.para1 = 1.5;
                sn.para2 = 0.85;
                double pt = sn.getMean()*unrelMachineEfficiency;

                for(int j=1;j<=nbStages;j++){
                    if (unre.contains(j)){
                        theSystem.CT[j]=new StochNum();
                        theSystem.CT[j].distribution= "Deterministic";
                        theSystem.CT[j].para1 = pt;
                    }else{
                        theSystem.CT[j]=new StochNum();
                        theSystem.CT[j].distribution= "LogNorm";
                        theSystem.CT[j].para1 = 1.5;
                        theSystem.CT[j].para2 = 0.85;
                    }
                }

                double[][] tij = new double[nbJobs + 1][nbStages + 1];
                int seed = (int) System.currentTimeMillis();
                theSystem.procTimeGeneration(nbJobs,tij,seed);
                for(int j: unre){
                    double[] Machinept = new double[nbJobs + 1];
                    Failure myFailure = new Failure();
                    for (int i = 1; i <= nbJobs; i++) {
                        Machinept[i] = tij[i][j];
                    }
                    myFailure.repairTimeGeneration(Machinept, 80*pt, 20*pt);
                    for (int i = 1; i <= nbJobs; i++) {
                        tij[i][j] = tij[i][j] + myFailure.repairTimeSamples[i];
                    }
                }

                writersum.write(str+" ");
                theSystem.buffer=referBuffer.clone();
                theSystem.mySimulation = theSystem.new SimulationBAS(nbJobs, 1, tij);
                theSystem.mySimulation.simBAS(false, 0);
                double thstar = theSystem.TH;

                BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(theSystem, thstar, lB, uB, nbJobs, 1);
                myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());

                Stopwatch totalAlter6RevTime = new Stopwatch();
                totalAlter6RevTime.start();
                try {
                    myReversedAlter6.solveBAPWithIntModel(tij,false);
                } catch (Exception exc) {
                    exc.printStackTrace();
                }
                totalAlter6RevTime.stop();

                int totcap = 0;
                for (int j = 1; j <= theSystem.nbStage - 1; j++) {
                    totcap = totcap + theSystem.buffer[j];
                }
                writersum.write(myReversedAlter6.numit + " " + df.format(totalAlter6RevTime.elapseTimeSeconds) + " " + totcap + " ");
                for (int j = 1; j <= theSystem.nbStage - 1; j++) {
                    writersum.write(theSystem.buffer[j] + ",");
                }
                writersum.println();
            }*/

            // unreliable lines
            HashMap<String, int[]> MTTRs = new HashMap<>();
            MTTRs.put("UnreBalance3", new int[]{0,1,1,1,1,1,1});
            MTTRs.put("UnreBalance4", new int[]{0,100,100,100,100,100,100});
            /*MTTRs.put("UnreBalance1", new int[]{0,4,4,4,4,4,4});
            MTTRs.put("UnreBalance2", new int[]{0,40,40,40,40,40,40});
            MTTRs.put("UnreDecrease", new int[]{0,60,50,40,30,20,10});
            MTTRs.put("UnreIncrease", new int[]{0,10,20,30,40,50,60});
            MTTRs.put("UnreV-shape", new int[]{0,10,20,30,30,20,10});
            MTTRs.put("UnreInvertedV-shape", new int[]{0,30,20,10,10,20,30});
            MTTRs.put("UnreStep1", new int[]{0,10,10,10,30,30,30});
            MTTRs.put("UnreStep2", new int[]{0,60,60,60,30,30,30});*/

            for(String str: MTTRs.keySet()){
                SerialLine theSystem = new SerialLine();
                theSystem.nbStage = 6;
                theSystem.buffer = new int[theSystem.nbStage];
                theSystem.CT = new StochNum[theSystem.nbStage+1];

                double pt = 1;
                for (int j =1;j<=nbStages;j++){
                    theSystem.CT[j]=new StochNum();
                    theSystem.CT[j].distribution= "Deterministic";
                    theSystem.CT[j].para1 = pt;
                }
                double[][] tij = new double[nbJobs + 1][nbStages + 1];
                int seed = (int) System.currentTimeMillis();
                theSystem.procTimeGeneration(nbJobs,tij,seed);
                for(int j =1;j<=nbStages;j++){
                    double[] Machinept = new double[nbJobs + 1];
                    Failure myFailure = new Failure();
                    for (int i = 1; i <= nbJobs; i++) {
                        Machinept[i] = tij[i][j];
                    }
                    myFailure.repairTimeGeneration(Machinept, MTTRs.get(str)[j]*pt*4, MTTRs.get(str)[j]*pt);
                    for (int i = 1; i <= nbJobs; i++) {
                        tij[i][j] = tij[i][j] + myFailure.repairTimeSamples[i];
                    }
                }

                writersum.write(str+" ");
                theSystem.buffer=referBuffer.clone();
                theSystem.mySimulation = theSystem.new SimulationBAS(nbJobs, 1, tij);
                theSystem.mySimulation.simBAS(false, 0);
                double thstar = theSystem.TH;

                BendersIntModelAlter6ReversedCut myReversedAlter6 = new BendersIntModelAlter6ReversedCut(theSystem, thstar, lB, uB, nbJobs, 1);
                myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());

                Stopwatch totalAlter6RevTime = new Stopwatch();
                totalAlter6RevTime.start();
                try {
                    myReversedAlter6.solveBAPWithIntModel(tij,false);
                } catch (Exception exc) {
                    exc.printStackTrace();
                }
                totalAlter6RevTime.stop();

                int totcap = 0;
                for (int j = 1; j <= theSystem.nbStage - 1; j++) {
                    totcap = totcap + theSystem.buffer[j];
                }
                writersum.write(myReversedAlter6.numit + " " + df.format(totalAlter6RevTime.elapseTimeSeconds) + " " + totcap + " ");
                for (int j = 1; j <= theSystem.nbStage - 1; j++) {
                    writersum.write(theSystem.buffer[j] + ",");
                }
                writersum.println();
            }
        }

        try {
            outRessummary.close();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

}

