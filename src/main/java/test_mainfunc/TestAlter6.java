package test_mainfunc;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;

import test_mainfunc.optimization.BendersIntModelAlter5;
import test_mainfunc.optimization.BendersIntModelAlter6;
import test_mainfunc.optimization.BendersIntModelAlter6ReversedCut;
import test_mainfunc.optimization.BendersStolletz;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.util.Stopwatch;

public class TestAlter6 {

    public static void main(String argv[]){
        String programPath = System.getProperty("user.dir");

        // 1. an object of serial line is constructed by reading the input file /INPUT/SerialLine_test_6stage.txt
        //***   Input system files   *********************************************************
        String in_System = programPath + "\\INPUT\\SerialLine_test_6stage.txt";
        InputStream in_SystemFile = null;
        try {
            in_SystemFile = new FileInputStream(in_System);
        } catch (FileNotFoundException e) {
            System.err.println("Error opening 'System file'");
            System.exit(-1);
        }
        SerialLine mySystem=new SerialLine(in_SystemFile);


        //  2. parameters of BAP, such as buffer lowerbound, buffer upperbound, target throughput and simulation length, are defined
        int[] myLB=new int[mySystem.nbStage+1]; // buffer lowerbound
        int[] myUB=new int[mySystem.nbStage+1]; // buffer upperbound
        for(int j=1;j<=mySystem.nbStage;j++){
            myLB[j]=1;
            myUB[j]=40;
        }
        double myTHstar=0.45; //target throughput
        int N=100000; //simulation length


        // 3. the samples of processing time tij are generated.
        double[][] tij=new double[N+1][mySystem.nbStage+1];
        int seed =(int) System.currentTimeMillis();
        mySystem.procTimeGeneration(N,tij,seed);


        // 4. an object myReversedAlter6 of class BendersIntModelAlter6ReversedCut is constructed by providing an object of SerialLine, target throughput, buffer lowerbound, buffer upperbound and simulation length.
        BendersIntModelAlter6ReversedCut myReversedAlter6=new BendersIntModelAlter6ReversedCut(mySystem, myTHstar, myLB, myUB, N);
        myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());
        Stopwatch totalReversedAlter6Time=new Stopwatch();
        totalReversedAlter6Time.start();
        try{
            // 5. the method solveBAPWithIntModel(tij,false) is called to solve the BAP.
            myReversedAlter6.solveBAPWithIntModel(tij,false);
        }catch(Exception exc){exc.printStackTrace();}
        totalReversedAlter6Time.stop();
        int totalBuffer=0;
        for(int j=1;j<mySystem.nbStage;j++)
            totalBuffer+=mySystem.buffer[j];


        // 6. the solution is then printed in the output file /OUTPUT/Alter6_6stage.txt.
        //***   Output summary file   *********************************************************
        String out_resFileSummary = programPath + "\\OUTPUT\\Alter6_6stage.txt";
        OutputStream outRessummary= null;
        try {
            outRessummary = new FileOutputStream(out_resFileSummary);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(-1);
        }
        PrintWriter writersum = new PrintWriter(outRessummary, true);
        DecimalFormat df;
        df = new DecimalFormat("#.#####");
        df.setRoundingMode(RoundingMode.CEILING);
        writersum.println("Method totaltime iterations totalcap bj");
        writersum.print("ReversedAlter6 "+df.format(totalReversedAlter6Time.elapseTimeSeconds)+" seconds "+ myReversedAlter6.numit+" iterations "+totalBuffer+" ");
        for(int j=1;j<mySystem.nbStage;j++)
            writersum.print(mySystem.buffer[j]+" ");
        writersum.println();

    }
}
