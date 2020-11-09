package test_mainfunc;

import test_mainfunc.optimization.BendersIntModelAlter6;
import test_mainfunc.optimization.BendersIntModelAlter6ReversedCut;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.util.Stopwatch;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;

public class Main_VerifyAdditive {
    public static void main(String argv[]) {
        String programPath = System.getProperty("user.dir");

        int[] myLB=new int[7];
        int[] myUB=new int[7];
        for(int j=1;j<=5;j++){
            myLB[j]=0;
            myUB[j]=40;
        }
        double myTHstar;
        int N=1000000;


        // Ouput format
        DecimalFormat df;
        df = new DecimalFormat("#.#####");
        df.setRoundingMode(RoundingMode.CEILING);

        //***   Output summary file   *********************************************************
        String out_resFileSummary = programPath + File.separator+ "OUTPUT"+File.separator+"AddVerify.txt";

        OutputStream outRessummary= null;
        try {
            outRessummary = new FileOutputStream(out_resFileSummary);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(-1);
        }
        PrintWriter writersum = new PrintWriter(outRessummary, true);

        writersum.println("System totalcap bj");

        for(int r=1;r<=3;r++){

            int seed =(int) System.currentTimeMillis();

            // ********************* 1BN_1 ******************************************************* //
            //***   Input files   *********************************************************
            String in_System = programPath + File.separator+"INPUT"+File.separator+"Add_1BN_1.txt";
            InputStream in_SystemFile = null;
            try {
                in_SystemFile = new FileInputStream(in_System);
            } catch (FileNotFoundException e) {
                System.err.println("Error opening 'System file'");
                System.exit(-1);
            }
            SerialLine mySystem=new SerialLine(in_SystemFile);

            // sampling
            double[][] tij=new double[N+1][mySystem.nbStage+1];

            mySystem.procTimeGeneration(N,tij,seed);
            myTHstar =0.72/ mySystem.CT[1].getMean();

            BendersIntModelAlter6ReversedCut myReversedAlter6=new BendersIntModelAlter6ReversedCut(mySystem, myTHstar, myLB, myUB, N, 1);
            myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());
            Stopwatch totalReversedAlter6Time=new Stopwatch();
            totalReversedAlter6Time.start();
            try{
                myReversedAlter6.solveBAPWithIntModel(tij);
            }catch(Exception exc){exc.printStackTrace();}
            totalReversedAlter6Time.stop();
            int totalBuffer=0;
            for(int j=1;j<mySystem.nbStage;j++)
                totalBuffer+=mySystem.buffer[j];
            writersum.print("1BN_1 "+totalBuffer+" ");
            for(int j=1;j<mySystem.nbStage;j++)
                writersum.print(mySystem.buffer[j]+" ");
            writersum.println();

            try{
                in_SystemFile.close();
            }catch (Exception exc){
                exc.printStackTrace();
            }


            // ********************* 1BN_1 END ******************************************************* //


            // ********************* 1BN_6 ******************************************************* //
            //***   Input files   *********************************************************
            in_System = programPath + File.separator+"INPUT"+File.separator+"Add_1BN_6.txt";
            in_SystemFile = null;
            try {
                in_SystemFile = new FileInputStream(in_System);
            } catch (FileNotFoundException e) {
                System.err.println("Error opening 'System file'");
                System.exit(-1);
            }
            mySystem=new SerialLine(in_SystemFile);
            myTHstar =0.72/ mySystem.CT[6].getMean();
            mySystem.procTimeGeneration(N,tij,seed);

            myReversedAlter6=new BendersIntModelAlter6ReversedCut(mySystem, myTHstar, myLB, myUB, N, 1);
            myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());
            totalReversedAlter6Time.start();
            try{
                myReversedAlter6.solveBAPWithIntModel(tij);
            }catch(Exception exc){exc.printStackTrace();}
            totalReversedAlter6Time.stop();
            totalBuffer=0;
            for(int j=1;j<mySystem.nbStage;j++)
                totalBuffer+=mySystem.buffer[j];
            writersum.print("1BN_6 "+totalBuffer+" ");
            for(int j=1;j<mySystem.nbStage;j++)
                writersum.print(mySystem.buffer[j]+" ");
            writersum.println();
            try{
                in_SystemFile.close();
            }catch (Exception exc){
                exc.printStackTrace();
            }
            // ********************* 1BN_6 END ******************************************************* //


            // ********************* 2BN_16 ******************************************************* //
            //***   Input files   *********************************************************
            in_System = programPath + File.separator+"INPUT"+File.separator+"Add_2BN_16.txt";
            in_SystemFile = null;
            try {
                in_SystemFile = new FileInputStream(in_System);
            } catch (FileNotFoundException e) {
                System.err.println("Error opening 'System file'");
                System.exit(-1);
            }
            mySystem=new SerialLine(in_SystemFile);
            myTHstar = 0.72/mySystem.CT[1].getMean();
            mySystem.procTimeGeneration(N,tij,seed);

            myReversedAlter6=new BendersIntModelAlter6ReversedCut(mySystem, myTHstar, myLB, myUB, N, 1);
            myReversedAlter6.writer = new PrintWriter(OutputStream.nullOutputStream());
            totalReversedAlter6Time.start();
            try{
                myReversedAlter6.solveBAPWithIntModel(tij);
            }catch(Exception exc){exc.printStackTrace();}
            totalReversedAlter6Time.stop();
            totalBuffer=0;
            for(int j=1;j<mySystem.nbStage;j++)
                totalBuffer+=mySystem.buffer[j];
            writersum.print("2BN_16 "+totalBuffer+" ");
            for(int j=1;j<mySystem.nbStage;j++)
                writersum.print(mySystem.buffer[j]+" ");
            writersum.println();
            try{
                in_SystemFile.close();
            }catch (Exception exc){
                exc.printStackTrace();
            }
            // ********************* 2BN_16 END ******************************************************* //
        }






    }

}
