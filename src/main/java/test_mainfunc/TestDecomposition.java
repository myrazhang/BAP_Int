package test_mainfunc;

import test_mainfunc.optimization.BendersIntModelAlter6;
import test_mainfunc.optimization.BendersIntModelAlter6ReversedCut;
import test_mainfunc.optimization.BendersStolletz;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.util.Stopwatch;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;

public class TestDecomposition {

    public static void main(String argv[]){
        String programPath = System.getProperty("user.dir");

        //***   Input files   *********************************************************
        String in_System = programPath + "\\INPUT\\SerialLine_test_6stage.txt";
        InputStream in_SystemFile = null;
        try {
            in_SystemFile = new FileInputStream(in_System);
        } catch (FileNotFoundException e) {
            System.err.println("Error opening 'System file'");
            System.exit(-1);
        }


        //***   Output summary file   *********************************************************
        String out_resFileSummary = programPath + "\\OUTPUT\\Decomposition_6stage.txt";

        OutputStream outRessummary= null;
        try {
            outRessummary = new FileOutputStream(out_resFileSummary);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(-1);
        }
        PrintWriter writersum = new PrintWriter(outRessummary, true);

        // Ouput format
        DecimalFormat df;
        df = new DecimalFormat("#.#####");
        df.setRoundingMode(RoundingMode.CEILING);


        SerialLine mySystem=new SerialLine(in_SystemFile);
        int[] myLB=new int[mySystem.nbStage+1];
        int[] myUB=new int[mySystem.nbStage+1];
        for(int j=1;j<=mySystem.nbStage;j++){
            myLB[j]=1;
            myUB[j]=40;
        }
        double myTHstar=0.45;
        int N=100000;

        // sampling
        double[][] tij=new double[N+1][mySystem.nbStage+1];
        int seed = 2019;
        mySystem.procTimeGeneration(N,tij,seed);

        // output file
        String out_resFile = programPath +"\\OUTPUT\\Out.txt";
        OutputStream outRes= null;
        try {
            outRes = new FileOutputStream(out_resFile);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(-1);
        }

        //output
        writersum.println("Method totaltime iterations totalcap bj");
        int totalBuffer;

        //BAP with reversed_Alter6 + decomposition
        BendersIntModelAlter6ReversedCut myReversedAlter6Decomp=new BendersIntModelAlter6ReversedCut(mySystem, myTHstar, myLB, myUB, N, 1);
        myReversedAlter6Decomp.writer = new PrintWriter(outRes, true);
        Stopwatch totalReversedAlter6DecompTime=new Stopwatch();
        totalReversedAlter6DecompTime.start();
        try{
            myReversedAlter6Decomp.solveWithDecomposition(tij);
        }catch(Exception exc){exc.printStackTrace();}
        totalReversedAlter6DecompTime.stop();
        totalBuffer=0;
        for(int j=1;j<mySystem.nbStage;j++)
            totalBuffer+=mySystem.buffer[j];
        writersum.print("ReversedAlter6Decomp "+df.format(totalReversedAlter6DecompTime.elapseTimeSeconds)+" seconds "+ myReversedAlter6Decomp.numit+" iterations "+totalBuffer+" ");
        for(int j=1;j<mySystem.nbStage;j++)
            writersum.print(mySystem.buffer[j]+" ");
        writersum.println();


        //BAP with Stolletz + decomposition
        BendersStolletz myStolletzDecomposition=new BendersStolletz(mySystem, myTHstar, myLB, myUB, N, 1);
        Stopwatch totalStolletzDecompositionTime=new Stopwatch();
        totalStolletzDecompositionTime.start();
        try{
            myStolletzDecomposition.solveWithDecomposition(tij);
        }catch(Exception exc){exc.printStackTrace();}
        totalStolletzDecompositionTime.stop();
        totalBuffer=0;
        for(int j=1;j<mySystem.nbStage;j++)
            totalBuffer+=mySystem.buffer[j];
        writersum.print("StolletzDecomp "+df.format(totalStolletzDecompositionTime.elapseTimeSeconds)+"seconds "+ myStolletzDecomposition.numit+"iterations "+totalBuffer+" ");
        for(int j=1;j<mySystem.nbStage;j++)
            writersum.print(mySystem.buffer[j]+" ");
        writersum.println();


        //BAP with reversed_Alter6
        BendersIntModelAlter6ReversedCut myReversedAlter6=new BendersIntModelAlter6ReversedCut(mySystem, myTHstar, myLB, myUB, N, 1);
        myReversedAlter6.writer = new PrintWriter(outRes, true);
        Stopwatch totalReversedAlter6Time=new Stopwatch();
        totalReversedAlter6Time.start();
        try{
            myReversedAlter6.solveBAPWithIntModel(tij,true);
        }catch(Exception exc){exc.printStackTrace();}
        totalReversedAlter6Time.stop();
        totalBuffer=0;
        for(int j=1;j<mySystem.nbStage;j++)
            totalBuffer+=mySystem.buffer[j];
        writersum.print("ReversedAlter6 "+df.format(totalReversedAlter6Time.elapseTimeSeconds)+" seconds "+ myReversedAlter6.numit+" iterations "+totalBuffer+" ");
        for(int j=1;j<mySystem.nbStage;j++)
            writersum.print(mySystem.buffer[j]+" ");
        writersum.println();
        //BAP with Stolletz
        /*BendersStolletz myStolletz=new BendersStolletz(mySystem, myTHstar, myLB, myUB, N, 1);
        Stopwatch totalStolletzTime=new Stopwatch();
        totalStolletzTime.start();
        try{
            myStolletz.solveBAPWithStolletz(tij,0);
        }catch(Exception exc){exc.printStackTrace();}
        totalStolletzTime.stop();
        totalBuffer=0;
        for(int j=1;j<mySystem.nbStage;j++)
            totalBuffer+=mySystem.buffer[j];
        writersum.print("Stolletz "+df.format(totalStolletzTime.elapseTimeSeconds)+"seconds "+ myStolletz.numit+"iterations "+totalBuffer+" ");
        for(int j=1;j<mySystem.nbStage;j++)
            writersum.print(mySystem.buffer[j]+" ");
        writersum.println();*/

    }

}
