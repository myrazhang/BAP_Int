package test_mainfunc;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;


import static java.lang.Math.max;

public class Test_SerialLineSimulation {


    public static void main(String argv[]) {
        String programPath = System.getProperty("user.dir");

        //***   Input files   *********************************************************
        String in_System = programPath + "\\INPUT\\SerialLine_test.txt";
        InputStream in_SystemFile = null;
        try {
            in_SystemFile = new FileInputStream(in_System);
        } catch (FileNotFoundException e) {
            System.err.println("Error opening 'System file'");
            System.exit(-1);
        }



        //***   Output files   *********************************************************
        String out_resFile = programPath + "\\OUTPUT\\SimOutput_test.txt";
        OutputStream outRes= null;
        try {
            outRes = new FileOutputStream(out_resFile);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(-1);
        }

        // Simulation parameters
        int N=10000;
        int W=1000;
        // System read from INPUT file
        SerialLine mySystem= new SerialLine(in_SystemFile);
        double[][] tij=new double[N][mySystem.NbStage];
        int[][] uij=new int[N][mySystem.NbStage];
        int[][] vij=new int[N][mySystem.NbStage];
        int[][] wij=new int[N][mySystem.NbStage];
        int[][] sij=new int[N][mySystem.NbStage];
        double[] bar_Dij=new double[mySystem.NbStage];
        double thetapar=0.0;
        int B = 0;
        for (int j = 0; j < mySystem.NbStage - 1; j++) {
            B = max(B, mySystem.Buffer[j]);
        }
        double[][] bar_Sij = new double [B][mySystem.NbStage];

        mySystem.ProcTimeGeneration(N,tij);
        mySystem.SIM_Serial_BAS(N,W,tij,uij,vij,wij,sij,thetapar,bar_Sij,bar_Dij,0,false, false);

        PrintWriter writer = new PrintWriter(outRes, true);
        // cycle time is calculated between W and N, whatever value SteadyState will take.
        writer.print("Cycle time with simulation from 0: ");
        writer.println(mySystem.OverallCT);

        mySystem.SIM_Serial_BAS(N,W,tij,uij,vij,wij,sij,thetapar,bar_Sij,bar_Dij,0,false,false );
        writer.print("Cycle time with simulation from steady state: ");
        writer.println(mySystem.OverallCT);


        //optimization part
        long StartOpt = System.currentTimeMillis();
        mySystem.solveBender(N,W,tij,uij,vij,wij,sij,bar_Sij,bar_Dij,true);
        long elapsedTimeMillisC = System.currentTimeMillis()- StartOpt;
        float elapsedTimeSecC = elapsedTimeMillisC/1000F;
        double OptTime = elapsedTimeSecC;

        writer.print("Cycle time with Benders Decomposition: ");
        double tempCT=1/mySystem.TH;
        writer.println(tempCT);

        DecimalFormat df = new DecimalFormat("#.#####");
        df.setRoundingMode(RoundingMode.CEILING);
        writer.write("total time CG: ");
        writer.write(df.format(OptTime));

        writer.print("Optimal Buffer: ");
        for(int j=0;j<mySystem.NbStage-1;j++){
            writer.println(mySystem.BJsol[j][mySystem.totit]);
            writer.println(" ");
        }

        writer.print("Total number of iterations: ");
        writer.println(mySystem.totit);

        try {
            outRes.close();
            in_SystemFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
