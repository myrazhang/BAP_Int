package test_mainfunc;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;

import test_mainfunc.optimization.BendersIntModelAlter5;
import test_mainfunc.optimization.BendersIntModelAlter6;
import test_mainfunc.optimization.BendersStolletz;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.util.Stopwatch;

public class TestAlter6 {

    public static void main(String argv[]){
        String programPath = System.getProperty("user.dir");

        //***   Input files   *********************************************************
        String in_System = programPath + "\\INPUT\\SerialLine_test_10stage.txt";
        InputStream in_SystemFile = null;
        try {
            in_SystemFile = new FileInputStream(in_System);
        } catch (FileNotFoundException e) {
            System.err.println("Error opening 'System file'");
            System.exit(-1);
        }


        //***   Output summary file   *********************************************************
        String out_resFileSummary = programPath + "\\OUTPUT\\Alter6_test.txt";

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
        int seed =(int) System.currentTimeMillis();
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
        // BAP with Alter 5
        BendersIntModelAlter5 myAlter5=new BendersIntModelAlter5(mySystem, myTHstar, myLB, myUB, N, 1);
        myAlter5.writer = new PrintWriter(outRes, true);
        Stopwatch totalAlter5Time=new Stopwatch();
        totalAlter5Time.start();
        try{
            myAlter5.solveBAPWithIntModel(tij);
        }catch(Exception exc){exc.printStackTrace();}
        totalAlter5Time.stop();
        int totalBuffer=0;
        for(int j=1;j<mySystem.nbStage;j++)
            totalBuffer+=mySystem.buffer[j];
        writersum.print("Alter5 "+df.format(totalAlter5Time.elapseTimeSeconds)+"seconds "+ myAlter5.numit+"iterations "+totalBuffer+" ");
        for(int j=1;j<mySystem.nbStage;j++)
            writersum.print(mySystem.buffer[j]+" ");
        writersum.println();


        //BAP with Alter6
        BendersIntModelAlter6 myAlter6=new BendersIntModelAlter6(mySystem, myTHstar, myLB, myUB, N, 1);
        myAlter6.writer = new PrintWriter(outRes, true);
        Stopwatch totalAlter6Time=new Stopwatch();
        totalAlter6Time.start();
        try{
            myAlter6.solveBAPWithIntModel(tij);
        }catch(Exception exc){exc.printStackTrace();}
        totalAlter6Time.stop();
        totalBuffer=0;
        for(int j=1;j<mySystem.nbStage;j++)
            totalBuffer+=mySystem.buffer[j];
        writersum.print("Alter6 "+df.format(totalAlter6Time.elapseTimeSeconds)+"seconds "+ myAlter6.numit+"iterations "+totalBuffer+" ");
        for(int j=1;j<mySystem.nbStage;j++)
            writersum.print(mySystem.buffer[j]+" ");
        writersum.println();


        //BAP with Stolletz
        /*BendersStolletz myStolletz=new BendersStolletz(mySystem, myTHstar, myLB, myUB, N, 1);
        Stopwatch totalStolletzTime=new Stopwatch();
        totalStolletzTime.start();
        try{
            myStolletz.solveBAPWithStolletz(tij);
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
