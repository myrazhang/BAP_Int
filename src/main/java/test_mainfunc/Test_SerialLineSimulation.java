package test_mainfunc;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;


import static java.lang.Math.max;

public class Test_SerialLineSimulation {


    public static void main(String argv[]) {
        String programPath = System.getProperty("user.dir");

        //***   Input files   *********************************************************
        String in_System = programPath + "\\INPUT\\SerialLine_test_10stage_exp.txt";
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
        int N=50000;
        int W=1;

        // System read from INPUT file
        SerialLine mySystem= new SerialLine(in_SystemFile);
        double[][] tij=new double[N][mySystem.NbStage];
        mySystem.resc = new double [mySystem.Maxit];
        mySystem.ProcTimeGeneration(N,tij);

        PrintWriter writer = new PrintWriter(outRes, true);

        //optimization part with our cut
        long StartOpt;
        long elapsedTimeMillisC;
        float elapsedTimeSecC;
        double OptTime;
        double tempCT;
        DecimalFormat df;

        //Alter 5
        writer.print("---------");
        writer.println("Alter 5");
        StartOpt = System.currentTimeMillis();
        mySystem.solveBender(N,W,tij,false,5);
        elapsedTimeMillisC = System.currentTimeMillis()- StartOpt;
        elapsedTimeSecC = elapsedTimeMillisC/1000F;
        OptTime = elapsedTimeSecC;

        writer.print("TH: ");
        tempCT=mySystem.TH;
        writer.println(tempCT);

        df = new DecimalFormat("#.#####");
        df.setRoundingMode(RoundingMode.CEILING);
        writer.write("total time CG: ");
        writer.write(df.format(OptTime));
        writer.println();

        writer.write("cplex time CG: ");
        writer.write(df.format(mySystem.our_cplextime[2]));
        writer.println();

        writer.print("Optimal Buffer: ");
        for(int j=0;j<mySystem.NbStage-1;j++){
            writer.print(mySystem.Buffer[j]);
            writer.print(" ");
        }
        writer.println();

        writer.print("Total number of iterations: ");
        writer.println(mySystem.totit);
        /*writer.println("BAP of each iteration: ");
        for(int k=0;k<mySystem.numit;k++){
            for(int j=0;j<mySystem.NbStage-1;j++){
                writer.print(mySystem.BJsol[j][k]);
                writer.print(' ');
            }
            writer.println();
        }*/
        //Alter 4
        writer.print("---------");
        writer.println("Alter 4");
        StartOpt = System.currentTimeMillis();
        mySystem.solveBender(N,W,tij,false,4);
        elapsedTimeMillisC = System.currentTimeMillis()- StartOpt;
        elapsedTimeSecC = elapsedTimeMillisC/1000F;
        OptTime = elapsedTimeSecC;

        writer.print("TH: ");
        tempCT=mySystem.TH;
        writer.println(tempCT);

        df = new DecimalFormat("#.#####");
        df.setRoundingMode(RoundingMode.CEILING);
        writer.write("total time CG: ");
        writer.write(df.format(OptTime));
        writer.println();

        writer.write("cplex time CG: ");
        writer.write(df.format(mySystem.our_cplextime[1]));
        writer.println();

        writer.print("Optimal Buffer: ");
        for(int j=0;j<mySystem.NbStage-1;j++){
            writer.print(mySystem.Buffer[j]);
            writer.print(" ");
        }
        writer.println();

        writer.print("Total number of iterations: ");
        writer.println(mySystem.totit);

        /*writer.println("BAP of each iteration: ");
        for(int k=0;k<mySystem.numit;k++){
            for(int j=0;j<mySystem.NbStage-1;j++){
                writer.print(mySystem.BJsol[j][k]);
                writer.print(' ');
            }
            writer.println();
        }*/

        //Alter 3
        writer.print("---------");
        writer.println("Alter 3");
        StartOpt = System.currentTimeMillis();
        mySystem.solveBender(N,W,tij,false,3);
        elapsedTimeMillisC = System.currentTimeMillis()- StartOpt;
        elapsedTimeSecC = elapsedTimeMillisC/1000F;
        OptTime = elapsedTimeSecC;

        writer.print("TH: ");
        tempCT=mySystem.TH;
        writer.println(tempCT);

        df = new DecimalFormat("#.#####");
        df.setRoundingMode(RoundingMode.CEILING);
        writer.write("total time CG: ");
        writer.write(df.format(OptTime));
        writer.println();

        writer.write("cplex time CG: ");
        writer.write(df.format(mySystem.our_cplextime[0]));
        writer.println();

        writer.print("Optimal Buffer: ");
        for(int j=0;j<mySystem.NbStage-1;j++){
            writer.print(mySystem.Buffer[j]);
            writer.print(" ");
        }
        writer.println();

        writer.print("Total number of iterations: ");
        writer.println(mySystem.totit);

        /*writer.println("BAP of each iteration: ");
        for(int k=0;k<mySystem.numit;k++){
            for(int j=0;j<mySystem.NbStage-1;j++){
                writer.print(mySystem.BJsol[j][k]);
                writer.print(' ');
            }
            writer.println();
        }*/






        //optimization part with stolletz
        writer.print("---------");
        writer.println("Stolletz optimization");
        StartOpt = System.currentTimeMillis();
        mySystem.solveBender(N,W,tij,true,0);
        elapsedTimeMillisC = System.currentTimeMillis()- StartOpt;
        elapsedTimeSecC = elapsedTimeMillisC/1000F;
        OptTime = elapsedTimeSecC;

        writer.print("TH: ");
        tempCT=mySystem.TH;
        writer.println(tempCT);

        df = new DecimalFormat("#.#####");
        df.setRoundingMode(RoundingMode.CEILING);
        writer.write("total time CG: ");
        writer.write(df.format(OptTime));
        writer.println();
        writer.write("cplex time CG: ");
        writer.write(df.format(mySystem.stolletz_cplextime));
        writer.println();

        writer.print("Optimal Buffer: ");
        for(int j=0;j<mySystem.NbStage-1;j++){
            writer.print(mySystem.Buffer[j]);
            writer.print(" ");
        }
        writer.println();

        writer.print("Total number of iterations: ");
        writer.println(mySystem.totit);

        /*writer.println("BAP of each iteration: ");
        for(int k=0;k<mySystem.numit;k++){
            for(int j=0;j<mySystem.NbStage-1;j++){
                writer.print(mySystem.BJsol[j][k]);
                writer.print(' ');
            }
            writer.println();
        }*/



        try {
            outRes.close();
            in_SystemFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
