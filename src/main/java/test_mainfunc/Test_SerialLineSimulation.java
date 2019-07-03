package test_mainfunc;

import test_mainfunc.optimization.BendersIntModelAlter5;
import test_mainfunc.optimization.BendersStolletz;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.util.Stopwatch;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;

import static java.lang.Math.exp;


public class Test_SerialLineSimulation {

    public static void main(String argv[]) {
        String programPath = System.getProperty("user.dir");

        //***   Input files   *********************************************************
        String in_System = programPath + "\\INPUT\\SerialLine_test_DoE.txt";
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

        // Ouput format
        DecimalFormat df;
        df = new DecimalFormat("#.#####");
        df.setRoundingMode(RoundingMode.CEILING);

        // System read from INPUT file
        SystemCombinationsForDOE myDOE=new SystemCombinationsForDOE(in_SystemFile);
        PrintWriter writersum = new PrintWriter(outRessummary, true);

       //todo: for test
        //writersum.println("J TH BN1 BN2 alpha CTnoBN var A5_numiter A5_comptime A5_cplext ST_numiter ST_comptime ST_cplext Optsol");
        writersum.println("J BN1 BN2 alpha CTnoBN var TH_B1 BN_eta_B1 TH_B20 BN_eta_B20");
        //todo: end of test

        int BNpositions=0;
        //here the DoE starts
        for(int Jfac=0; Jfac < myDOE.Jfactor.length; Jfac++){
            if (myDOE.Jfactor[Jfac]==4){
                BNpositions =4;
            }
            else if (myDOE.Jfactor[Jfac]==7){
                BNpositions = 12;
            }
            for(int etaFac=0; etaFac < myDOE.etaFactor.length;etaFac++){
                for(int BNfac=0; BNfac < BNpositions; BNfac++){
                    for(int alfac=0; alfac< myDOE.alphafactor.length; alfac++){
                        for(int noBNctfac=0; noBNctfac< myDOE.noBNfactor.length;noBNctfac++){
                            for(int varfac=0; varfac< myDOE.varfactor.length;varfac++){

                                SerialLine mySystem=myDOE.getOneSystemConfiguration(Jfac, BNfac, alfac, noBNctfac, varfac);
                                double meanBnCt=mySystem.CT[myDOE.BN1[BNfac]].getMean();

                                int[] lB=new int[mySystem.nbStage];
                                int[] uB=new int[mySystem.nbStage];
                                for(int j=0;j<=mySystem.nbStage-1;j++){
                                    lB[j]=myDOE.Lj;
                                    uB[j]=myDOE.Uj;
                                }

                                // sampling
                                double[][] tij=new double[myDOE.Njobs+1][mySystem.nbStage+1];
                                int seed =(int) System.currentTimeMillis();
                                mySystem.procTimeGeneration(myDOE.Njobs,tij,seed);

                                // todo: simulation to check throughput
                                mySystem.mySimulation=mySystem.new SimulationBAS(myDOE.Njobs,myDOE.W,tij);
                                for(int j=1;j<=mySystem.nbStage-1;j++)
                                    mySystem.buffer[j]=lB[j];
                                mySystem.mySimulation.simDualBAS(false);


                                for(int j=1;j<=mySystem.nbStage-1;j++)
                                    mySystem.buffer[j]=lB[j];
                                mySystem.mySimulation.simDualBAS(false);


                                double BN_eta=mySystem.TH*meanBnCt;

                                writersum.write(" " + mySystem.nbStage + " "+myDOE.BN1[BNfac]+" " + myDOE.BN2[BNfac] +" "+myDOE.alphafactor[alfac]+" "+ +myDOE.noBNfactor[noBNctfac]+" "+myDOE.varfactor[varfac]+" " + mySystem.TH+" "+BN_eta);
                                for(int j=1;j<=mySystem.nbStage-1;j++)
                                    mySystem.buffer[j]=uB[j];
                                mySystem.mySimulation.simDualBAS(false);
                                BN_eta=mySystem.TH*meanBnCt;
                                writersum.write(" "+mySystem.TH+" "+BN_eta);
                                writersum.println();

                                // todo: end of test


                                //Start optimization with Alter 5
                                /*myDOE.tempinstance = "J"+mySystem.nbStage+"_TH_"+myDOE.etaFactor[etaFac] +"_BN_"+myDOE.BN1[BNfac]+myDOE.BN2[BNfac]+"_alpha_"+myDOE.alphafactor[alfac]+"_BNf_"+myDOE.noBNfactor[noBNctfac]+"_var_"+myDOE.varfactor[varfac];
                                String out_resFile = programPath +"\\OUTPUT\\Out_"+ myDOE.tempinstance + ".txt";
                                OutputStream outRes= null;
                                try {
                                    outRes = new FileOutputStream(out_resFile);
                                } catch (FileNotFoundException e) {
                                    e.printStackTrace();
                                    System.exit(-1);
                                }
                                BendersIntModelAlter5 myAlter5=new BendersIntModelAlter5(mySystem, myDOE.etaFactor[etaFac]/meanBnCt, lB, uB, myDOE.Njobs, myDOE.W);
                                myAlter5.writer = new PrintWriter(outRes, true);


                                Stopwatch totalAlter5Time=new Stopwatch();
                                totalAlter5Time.start();
                                try{
                                    myAlter5.solveBAPWithIntModel(tij);
                                }catch(Exception exc){exc.printStackTrace();}

                                totalAlter5Time.stop();


                                writersum.write(" " + mySystem.nbStage + " "+myAlter5.THstar +" "+myDOE.BN1[BNfac]+" " + myDOE.BN2[BNfac] +" "+myDOE.alphafactor[alfac]+" "+myDOE.noBNfactor[noBNctfac]+" "+myDOE.varfactor[varfac] + " ");
                                writersum.write(myAlter5.numit + " " + df.format(totalAlter5Time.elapseTimeSeconds)+ " " +df.format(myAlter5.cplexTimeMeasure.elapseTimeSeconds)+ " ");
                                // End Optimization with Alter5


                                // Start optimization with stolletz
                                BendersStolletz myStolletz=new BendersStolletz(mySystem, myDOE.etaFactor[etaFac]/meanBnCt, lB, uB, myDOE.Njobs, myDOE.W);

                                Stopwatch totalStolletzTime=new Stopwatch();
                                totalStolletzTime.start();
                                try{
                                    myStolletz.solveBAPWithStolletz(tij);
                                }catch(Exception exc){exc.printStackTrace();}

                                totalStolletzTime.stop();

                                int totcap=0;
                                for(int j=1;j<=mySystem.nbStage-1;j++){
                                    totcap = totcap+ mySystem.buffer[j];
                                }

                                writersum.write(myStolletz.numit + " " + df.format(totalStolletzTime.elapseTimeSeconds)+ " " +df.format(myStolletz.cplexTimeMeasure.elapseTimeSeconds)+ " " + totcap);
                                writersum.println();
                                //End of optimization with Stolletz

                                //close single instance file
                                try {
                                    outRes.close();
                                } catch (IOException e) {
                                    e.printStackTrace();
                                }*/
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
