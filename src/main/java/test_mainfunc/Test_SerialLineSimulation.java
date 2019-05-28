package test_mainfunc;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;


import static java.lang.Math.max;

public class Test_SerialLineSimulation {


    public static void main(String argv[]) {
        String programPath = System.getProperty("user.dir");

        //***   Input files   *********************************************************
        String in_System = programPath + "\\INPUT\\SerialLine_test_DoE_prova.txt";
        InputStream in_SystemFile = null;
        try {
            in_SystemFile = new FileInputStream(in_System);
        } catch (FileNotFoundException e) {
            System.err.println("Error opening 'System file'");
            System.exit(-1);
        }



        //***   Output summary file   *********************************************************
        String out_resFileSummary = programPath + "\\OUTPUT\\SimOutput_summary.txt";

        OutputStream outRessummary= null;
        try {
            outRessummary = new FileOutputStream(out_resFileSummary);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(-1);
        }

        // Simulation parameters
        //int N=100000;
        //int W=1;

        // System read from INPUT file
        SerialLine mySystem= new SerialLine(in_SystemFile);
        PrintWriter writersum = new PrintWriter(outRessummary, true);

        writersum.println("J TH BN1 BN2 alpha CTnoBN var A5_numiter A5_comptime A5_cplext ST_numiter ST_comptime ST_cplext Optsol");

        int BNpositions=0;
        //here the DoE starts
        for(int Jfac=0; Jfac < mySystem.Jfactor.length; Jfac++){
            if (mySystem.Jfactor[Jfac]==4){
                BNpositions =4;
            }
            else if (mySystem.Jfactor[Jfac]==7){
                BNpositions = 12;
            }
            for(int THfac=0; THfac < mySystem.THfactor.length;THfac++){
                for(int BNfac=0; BNfac < BNpositions; BNfac++){
                    for(int alfac=0; alfac< mySystem.alphafactor.length; alfac++){
                        for(int noBNctfac=0; noBNctfac< mySystem.noBNfactor.length;noBNctfac++){
                            for(int varfac=0; varfac< mySystem.varfactor.length;varfac++){

                                //configuration of a single instance
                                mySystem.SystemConfiguration(Jfac, THfac, BNfac, alfac, noBNctfac, varfac);

                                //open output file
                                mySystem.tempinstance = "J"+mySystem.NbStage+"_TH_"+mySystem.THstar +"_BN_"+mySystem.BN1[BNfac]+mySystem.BN2[BNfac]+"_alpha_"+mySystem.alphafactor[alfac]+"_BNf_"+mySystem.noBNfactor[noBNctfac]+"_var_"+mySystem.varfactor[varfac];
                                String out_resFile = programPath +"\\OUTPUT\\Out_"+ mySystem.tempinstance + ".txt";
                                OutputStream outRes= null;
                                try {
                                    outRes = new FileOutputStream(out_resFile);
                                } catch (FileNotFoundException e) {
                                    e.printStackTrace();
                                    System.exit(-1);
                                }
                                mySystem.writer = new PrintWriter(outRes, true);

                                //start optimization
                                double[][] tij=new double[mySystem.Njobs][mySystem.NbStage];
                                mySystem.resc = new double [mySystem.Maxit];
                                mySystem.ProcTimeGeneration(mySystem.Njobs,tij);

                                //optimization part with our cut
                                long StartOpt;
                                long elapsedTimeMillisC;
                                float elapsedTimeSecC;
                                double OptTime;
                                DecimalFormat df;

                                //Alter 5
                                StartOpt = System.currentTimeMillis();
                                mySystem.solveBender(mySystem.Njobs,mySystem.W,tij,false,5);
                                elapsedTimeMillisC = System.currentTimeMillis()- StartOpt;
                                elapsedTimeSecC = elapsedTimeMillisC/1000F;
                                OptTime = elapsedTimeSecC;
                                df = new DecimalFormat("#.#####");
                                df.setRoundingMode(RoundingMode.CEILING);

                                //write results on summary file
                                writersum.write(" " + mySystem.NbStage + " "+mySystem.THstar +" "+mySystem.BN1[BNfac]+" " + mySystem.BN2[BNfac] +" "+mySystem.alphafactor[alfac]+" "+mySystem.noBNfactor[noBNctfac]+" "+mySystem.varfactor[varfac] + " ");
                                writersum.write(mySystem.totit + " " + df.format(OptTime)+ " " +df.format(mySystem.our_cplextime[2])+ " ");

                                //optimization part with stolletz
                                StartOpt = System.currentTimeMillis();
                                mySystem.solveBender(mySystem.Njobs,mySystem.W,tij,true,0);
                                elapsedTimeMillisC = System.currentTimeMillis()- StartOpt;
                                elapsedTimeSecC = elapsedTimeMillisC/1000F;
                                OptTime = elapsedTimeSecC;
                                df = new DecimalFormat("#.#####");
                                df.setRoundingMode(RoundingMode.CEILING);
                                int totcap=0;
                                for(int j=0;j<mySystem.NbStage-1;j++){
                                    totcap = totcap+ mySystem.Buffer[j];
                                }
                                //write on summary file
                                writersum.write(mySystem.totit + " " + df.format(OptTime)+ " " +df.format(mySystem.stolletz_cplextime)+ " " + totcap);
                                writersum.println();

                                //close single instance file
                                try {
                                    outRes.close();
                                } catch (IOException e) {
                                    e.printStackTrace();
                                }
                            }
                        }

                    }//end Mfac
                }//end BNfac
            }//end THfac
        }//end jfac

        try {
            outRessummary.close();
            in_SystemFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

          /*writer.println("BAP of each iteration: ");
        for(int k=0;k<mySystem.numit;k++){
            for(int j=0;j<mySystem.NbStage-1;j++){
                writer.print(mySystem.BJsol[j][k]);
                writer.print(' ');
            }
            writer.println();
        }*/


        /*writer.println("BAP of each iteration: ");
        for(int k=0;k<mySystem.numit;k++){
            for(int j=0;j<mySystem.NbStage-1;j++){
                writer.print(mySystem.BJsol[j][k]);
                writer.print(' ');
            }
            writer.println();
        }*/
        //Alter 4
        /*writer.print("---------");
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
        writer.println(mySystem.totit);*/

        /*writer.println("BAP of each iteration: ");
        for(int k=0;k<mySystem.numit;k++){
            for(int j=0;j<mySystem.NbStage-1;j++){
                writer.print(mySystem.BJsol[j][k]);
                writer.print(' ');
            }
            writer.println();
        }*/

        //Alter 3
       /* writer.print("---------");
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
        writer.println(mySystem.totit);*/

        /*writer.println("BAP of each iteration: ");
        for(int k=0;k<mySystem.numit;k++){
            for(int j=0;j<mySystem.NbStage-1;j++){
                writer.print(mySystem.BJsol[j][k]);
                writer.print(' ');
            }
            writer.println();
        }*/


    }

}

////////////////////////////////////////////////////////////////////////
////////CODE TO CHECK THE CORRECTNESS OF THE DESIGN OF EXPERIMENT///////
////////////////////////////////////////////////////////////////////////
/*
writer.println("ist: "+numistances);
        writer.println("Jfac THfac BNfac alfac noBNctfac varfac ARE: ");
        writer.println(Jfac + " " + THfac + " " + BNfac + " " + alfac + " " + noBNctfac + " " + varfac);
        writer.write("J: " + mySystem.NbStage + " and ");
        writer.write("Buffer length: " + mySystem.Buffer.length+ " and ");
        writer.write("TH star: " + mySystem.THstar+ " and ");
        writer.println("tij factors are: ");
        writer.write("alpha: "+ mySystem.alphafactor[alfac] + "and ");
        writer.write("beta: "+ mySystem.betafactor[alfac] + "and ");
        writer.println("var: "+ mySystem.varfactor[varfac] + "and ");
        writer.println("BN positions are: " + mySystem.BN1[BNfac] + " and "+ mySystem.BN2[BNfac]);
        writer.println("tij parameters are: ");
        for (int j = 0; j < mySystem.NbStage; j++)
        {
        writer.println( "j: "+ j + ": " + mySystem.CT[j].distribution + " " + mySystem.CT[j].Para1 + " " +mySystem.CT[j].Para2 + " " +mySystem.CT[j].Para3 + " " +mySystem.CT[j].Para4);
        }

        writer.println("------------------");*/
