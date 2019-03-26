package test_mainfunc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.io.InputStream;
import java.util.Date;
import java.util.Random;
import java.util.Scanner;

import ilog.concert.*;
import ilog.cplex.*;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
// import org.apache.commons.math3.distribution.UniformRealDistribution;
// import org.apache.commons.math3.distribution.WeibullDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;

import static java.lang.Math.max;
import static java.lang.Math.min;

public class SerialLine {
    public int NbStage;
    public int[] Buffer;
    public StochNum[] CT;
    public int[] Uj;
    public int[] Lj;
    //public int DeltaUj;
    //public int DeltaLj;
    public double THstar;

    //Sim_Output
    public double OverallCT;
    public double TH;
    public int Maxit = 20000;
    public int numit;

    public double[][] bar_Sij;
    public double[] bar_Dij;
    public int[][] ubarij;
    public int[][] wbarij;
    public double stolletz_cplextime;
    public double[] our_cplextime=new double[3];
    public boolean[][] Alter5_Delta;


    //cplex public variables
    public IloCplex cplex;

    //ERICA: CHECK IF NEEDED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    public IloNumVar[] BJ;
    public IloNumVar[][] DeltaBJrp;
    public IloNumVar[][] DeltaBJrm;
    public IloNumVar[][] DeltaJKp;
    public IloNumVar[][] DeltaJKm;

    public IloLinearNumExpr objective;
    public IloNumVar[][] YJK;
    public IloNumVar[][] YJR;
    public IloNumVar[][] ZJK;

    // Solutions from all iterations
    public int totit;
    public int BJsol[][];
    public double[][] DeltapPar;
    public double[][] DeltamPar;
    public double[] resc;


    // Create a new test_mainfunc.SerialLine from txt file
    public SerialLine(InputStream system) {

        Scanner scanner=new Scanner (system);
        scanner.useDelimiter("\\s+");


        // scan J size
        scanner.next();
        this.NbStage=scanner.nextInt();

        //scan Lowerbound Lj
        scanner.next();
        this.Lj = new int [this.NbStage - 1];
        int templj = scanner.nextInt();
        //this.DeltaLj=templj;
        for(int j=0;j< this.NbStage -1;j++)
        {
            this.Lj[j]= templj;
        }

        //scan upperbound Uj
        scanner.next();
        this.Uj = new int [this.NbStage - 1];
        int tempuj = scanner.nextInt();
        //this.DeltaUj=tempuj;
        for(int j=0;j< this.NbStage -1;j++)
        {
            this.Uj[j]= tempuj;
        }

        //scan optimization TH*
        scanner.next();
        this.THstar=scanner.nextDouble();

        //scan simulation buffer configuration
        this.Buffer = new int[this.NbStage - 1];
        scanner.next();
        for (int j = 0; j < this.NbStage - 1; j++)
        {
            this.Buffer[j]=scanner.nextInt();
        }

        //scan distributions of machine processing time
        scanner.next();
        this.CT = new StochNum[this.NbStage];
        scanner.next();
        for (int j = 0; j < this.NbStage; j++)
        {
            this.CT[j]=new StochNum();
            this.CT[j].distribution=scanner.next(); // ???????????? this.CT[j].distribution 是字符串
        }
        scanner.next();
        for (int j = 0; j < this.NbStage; j++)
        {
            this.CT[j].Para1=scanner.nextDouble();
        }
        scanner.next();
        for (int j = 0; j < this.NbStage; j++)
        {
            this.CT[j].Para2=scanner.nextDouble();
        }
        scanner.next();
        for (int j = 0; j < this.NbStage; j++)
        {
            this.CT[j].Para3=scanner.nextDouble();
        }
        scanner.next();
        for (int j = 0; j < this.NbStage; j++)
        {
            this.CT[j].Para4=scanner.nextDouble();
        }
    }

    // Processing time generation
    public void ProcTimeGeneration(int N, double[][] pij){
        RandomGenerator generator = RandomGeneratorFactory.createRandomGenerator(new Random());

        //Seed
        generator.setSeed(new Date().getTime());

        for (int j = 0; j < NbStage; j++) {
            //Processing time
            //Normal distribution
            //CT[j].Para1 -  mean
            //CT[j].Para2 -  coefficient of variance
            if (CT[j].distribution.equals("Norm")) {
                NormalDistribution pt = new NormalDistribution(generator, CT[j].Para1, CT[j].Para1 * CT[j].Para2);
                double p = pt.sample();
                for (int i = 0; i < N; i++) {
                    while (p < 0 || p > 2 * CT[j].Para1) {
                        p = pt.sample();
                    }
                    pij[i][j] = p;
                    p = pt.sample();
                }
            }
            else if(CT[j].distribution.equals("Beta")){
                BetaDistribution pt=new BetaDistribution(generator, CT[j].Para1, CT[j].Para2);

                double p = pt.sample();
                for (int i = 0; i < N; i++) {
                    pij[i][j] =CT[j].Para3+ p * (CT[j].Para4 - CT[j].Para3);
                    p = pt.sample();
                }
            }
            // else if(CT[j].distribution.equals("Exp")){...}
        }
    }


    //*****************************************************************************************************************
    // Simulation
    // Dual variables: uij, vij, wij, sij
    /////////////////////////// OUTPUT //////////////////////////////////////////////////////
    /////uij[i][j] : D[i][j] trigerred by S[i][j]                ////////////////////////////
    /////wij[i][j] : D[i][j] trigerred by S[i-b_j][j+1] BLOCKAGE ////////////////////////////
    /////vij[i][j] : S[i][j] trigerred by D[i-1][j]              ////////////////////////////
    /////sij[i][j] : S[i][j] trigerred by D[i][j-1]   STARVATION ////////////////////////////

    /////////////////////////// INPUT ///////////////////////////////////////////////////////
    ///// N : Simulation length /////////////////////////////////////////////////////////////
    ///// W : Warmup length (in terms of part number) ///////////////////////////////////////
    ///// tij : processing time of a part at a station //////////////////////////////////////
    ///// bar_Dij, bar_Sij : event time in warmup period ////////////////////////////////////
    ///// SteadyState: simulation from W or from 0:
    /////           if SteadyState=false, simulate from 0; bar_Dij and bar_Sij will be output
    /////           if SteadyState=true, simulate from W; bar_Dij and bar_Sij will be used as input
    /////////////////////////////////////////////////////////////////////////////////////////

    public void SIM_Serial_BAS(
            int N,
            int W,
            double[][] tij,
            int[][] uij,
            int[][] vij,
            int[][] wij,
            int[][] sij,
            double[][] bar_Sij,
            double[] bar_Dij,
            boolean SteadyState) {
        //*****************************************************************************
        //*****************************************************************************
        //** The notations of u,v,w,s is consistent with the journal paper ************
        //** on throughput improvement (2018).                             ************
        //*****************************************************************************
        //*****************************************************************************

        //***************************************
        //********** EVENT TIME *****************
        //***************************************
        double[][] Sij = new double[N][NbStage];
        double[][] Dij = new double[N][NbStage];
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < NbStage; j++)
            {
                Sij[i][j] = 0;
                Dij[i][j] = 0;
            }
        }


        //***************************************
        //********** SIM ERG variables   ********
        //***************************************
        int[][] bsij = new int[N][this.NbStage];
        int[][] buij = new int[N][this.NbStage];
        int[][] bvij = new int[N][this.NbStage];
        int[][] bwij = new int[N][this.NbStage];


        for (int i = 0; i <N; i++)
        {

            for (int j = 0; j < NbStage; j++)
            {
                bsij[i][j] = 0;
                buij[i][j] = 0;
                bvij[i][j] = 0;
                bwij[i][j] = 0;

                // initialize dual variables
                sij[i][j] = 0;
                uij[i][j] = 0;
                vij[i][j] = 0;
                wij[i][j] = 0;
            }
        }

        if (!SteadyState) {

            /////////////////////////////////////////////////////////////////////////////////////////////
            ///// Simulate from i=0 /////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////
            /////buij[i][j] = 1: D[i][j] trigerred by S[i][j]               /////////////////////////////
            /////bwij[i][j] = 1: D[i][j] trigerred by S[i-b_j][j+1] BLOCKAGE ////////////////////////////
            /////bvij[i][j] = 1: S[i][j] trigerred by D[i-1][j]             /////////////////////////////
            /////bsij[i][j] = 1: S[i][j] trigerred by D[i][j-1]   STARVATION ////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////


            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < NbStage; j++)
                {
                    ////////Starting event ////////////////////////////
                    //first arrival
                    if (i == 0 && j == 0)
                    {
                        Sij[i][j] = 0;
                        bsij[i][j] = 1;
                    }
                    //the first part facing the empty line
                    else if (i == 0 && j > 0)
                    {
                        Sij[i][j] = Dij[i][j - 1];
                        bsij[i][j] = 1;
                    }
                    //first stage: no starvation
                    else if (i > 0 && j == 0)
                    {
                        Sij[i][j] = Dij[i - 1][j];
                        bvij[i][j] = 1;
                    }
                    //others: max{D[i-1][j],D[i][j-1]}
                    else
                    {
                        if (Dij[i - 1][j] > Dij[i][j - 1])
                        {
                            Sij[i][j] = Dij[i - 1][j];
                            bvij[i][j] = 1;
                        }
                        else
                        {
                            Sij[i][j] = Dij[i][j - 1];
                            bsij[i][j] = 1;
                        }
                    }


                    //////// Departure event ////////////////////////////
                    // First parts: no blockage
                    if (j < NbStage - 1 && i < Buffer[j])
                    {
                        Dij[i][j] = Sij[i][j] + tij[i][j];
                        buij[i][j] = 1;
                    }
                    //last stage: no blockage
                    else if (j == NbStage - 1)
                    {
                        Dij[i][j] = Sij[i][j] + tij[i][j];
                        buij[i][j] = 1;
                    }
                    //others: max{S[i][j]+tij,S[i-b_j][j]}
                    else
                    {
                        if (Sij[i][j] + tij[i][j] > Sij[i - Buffer[j]][j + 1])
                        {
                            Dij[i][j] = Sij[i][j] + tij[i][j];
                            buij[i][j] = 1;
                        }
                        else
                        {
                            Dij[i][j] = Sij[i - Buffer[j]][j + 1];
                            bwij[i][j] = 1;
                        }
                    }
                }
            }

            /////////////////////////////////////////////////////////////////////////
            /////   Dual Optimal Solution               /////////////////////////////
            /////////////////////////////////////////////////////////////////////////
            for (int i = N - 1; i >= 0; i--)
            {
                for (int j = NbStage - 1; j >=0 ; j--)
                {
                    /////////Departure event////////////////////
                    if (i < N - 1 && j < NbStage - 1)
                    {
                        uij[i][j] = (vij[i + 1][j] + sij[i][j + 1])*buij[i][j];
                        wij[i][j] = (vij[i + 1][j] + sij[i][j + 1])*bwij[i][j];
                    }
                    else if (i == N - 1 && j < NbStage - 1)
                    {
                        uij[i][j] = (sij[i][j + 1])*buij[i][j];
                        wij[i][j] = (sij[i][j + 1])*bwij[i][j];
                    }
                    else if (i < N - 1 && j == NbStage - 1)
                    {
                        uij[i][j] = (vij[i + 1][j])*buij[i][j];
                    }
                    else if (i == N - 1 && j == NbStage - 1)
                    {
                        uij[i][j] = 1;
                    }

                    /////////Starting event///////////////////////
                    if (j > 0 && i < N - Buffer[j - 1])
                    {
                        sij[i][j] = (uij[i][j] + wij[i + Buffer[j - 1]][j - 1])*bsij[i][j];
                        vij[i][j] = (uij[i][j] + wij[i + Buffer[j - 1]][j - 1])*bvij[i][j];
                    }
                    else
                    {
                        sij[i][j] = (uij[i][j])*bsij[i][j];
                        vij[i][j] = (uij[i][j])*bvij[i][j];
                    }

                }
            }

            /*for(int j=0;j<NbStage-1;j++)
            {
                System.out.println("buffer at " + j + " is: " + Buffer[j]);
            }*/

            this.OverallCT = (Dij[N - 1][NbStage - 1] - Dij[W - 1][NbStage - 1]) / (double)(N - W);

            //***************   SAVE bar_Sij, bar_Dij    ***************
            /*for (int j = 0; j < NbStage; j++) {
                bar_Dij[j] = Dij[W - 1][j];
                if (j != NbStage -1) {
                    for (int i = 0 ; i < this.Uj[j]; i++) {
                        bar_Sij[i][j+1] = Sij[W- Uj[j]+i][j+1];
                    }
                }
            }*/


        }
        else {
            //***************   Fix event times in transient period   ***************
            for (int j = 0; j < NbStage; j++) {
                Dij[W - 1][j] = bar_Dij[j];
                if (j != NbStage - 1) {
                    for (int i = 0; i < this.Uj[j]; i++) {
                        Sij[W- Uj[j]+i][j+1] = bar_Sij[i][j+1];
                    }
                }
            }

            /////////////////////////////////////////////////////////////////////////////////////////////
            /////bzij[i][j] = 1: D[i][j] trigerred by S[i][j]                ////////////////////////////
            /////bwij[i][j] = 1: D[i][j] trigerred by S[i-b_j][j+1] BLOCKAGE ////////////////////////////
            /////buij[i][j] = 1: S[i][j] trigerred by D[i-1][j]              ////////////////////////////
            /////bvij[i][j] = 1: S[i][j] trigerred by D[i][j-1]   STARVATION ////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////
            for (int i = W; i < N; i++)
            {
                for (int j = 0; j < NbStage; j++)
                {
                    ////////Starting event ////////////////////////////

                    //first stage: no starvation
                    if (j == 0)
                    {
                        Sij[i][j] = Dij[i - 1][j];
                        bvij[i][j] = 1;
                    }
                    //others: max{D[i-1][j],D[i][j-1]}
                    else
                    {
                        if (Dij[i - 1][j] > Dij[i][j - 1])
                        {
                            Sij[i][j] = Dij[i - 1][j];
                            bvij[i][j] = 1;
                        }
                        else
                        {
                            Sij[i][j] = Dij[i][j - 1];
                            bsij[i][j] = 1;
                        }
                    }


                    ////////Departure event ////////////////////////////
                    //last stage: no blockage
                    if (j == NbStage - 1)
                    {
                        Dij[i][j] = Sij[i][j] + tij[i][j];
                        buij[i][j] = 1;
                    }
                    //others: max{S[i][j]+tij,S[i-b_j][j]}
                    else
                    {
                        if (Sij[i][j] + tij[i][j] > Sij[i - Buffer[j]][j + 1])
                        {
                            Dij[i][j] = Sij[i][j] + tij[i][j];
                            buij[i][j] = 1;
                        }
                        else
                        {
                            Dij[i][j] = Sij[i - Buffer[j]][j + 1];
                            bwij[i][j] = 1;
                        }
                    }
                }
            }

            for (int i = N - 1; i >= W; i--)
            {
                for (int j = NbStage - 1; j >= 0; j--)
                {
                    /////////Departure event////////////////////
                    if (i < N - 1 && j < NbStage - 1)
                    {
                        uij[i][j] = (vij[i + 1][j] + sij[i][j + 1])*buij[i][j];
                        wij[i][j] = (vij[i + 1][j] + sij[i][j + 1])*bwij[i][j];
                    }
                    else if (i == N - 1 && j < NbStage - 1)
                    {
                        uij[i][j] = (sij[i][j + 1])*buij[i][j];
                        wij[i][j] = (sij[i][j + 1])*bwij[i][j];
                    }
                    else if (i < N - 1 && j == NbStage - 1)
                    {
                        uij[i][j] = (vij[i + 1][j])*buij[i][j];
                    }
                    else if (i == N - 1 && j == NbStage - 1)
                    {
                        uij[i][j] = 1;
                    }

                    /////////Starting event///////////////////////
                    if (j > 0 && i < N - Buffer[j - 1])
                    {
                        sij[i][j] = (uij[i][j] + wij[i + Buffer[j - 1]][j - 1])*bsij[i][j];
                        vij[i][j] = (uij[i][j] + wij[i + Buffer[j - 1]][j - 1])*bvij[i][j];
                    }
                    else
                    {
                        sij[i][j] = (uij[i][j])*bsij[i][j];
                        vij[i][j] = (uij[i][j])*bvij[i][j];
                    }

                }
            }


        }
        this.OverallCT = (Dij[N-1][NbStage-1]- Dij[W-1][NbStage-1])/ (double)(N-W);

        /*for (int j=0; j<this.NbStage;j++)
        {
            for(int i=0; i<N;i++)
            {
                this.ubarij[i][j][iter]=uij[i][j];
                this.wbarij[i][j][iter]=wij[i][j];
            }

        }*/
    }

    public void solveBender(
            int jobs,
            int W,
            double[][] tij,
            boolean Stolletz,
            int AlterID) {


        try {
            for(int i=0;i<3;i++)
                our_cplextime[i]=0.0;
            stolletz_cplextime=0.0;
            double theta = jobs - W;
            //Initialization properties of the class
            int[][] vij=new int[jobs][this.NbStage];
            int[][] sij=new int[jobs][this.NbStage];
            DeltapPar = new double[this.NbStage-1][this.Maxit];
            DeltamPar = new double[this.NbStage-1][this.Maxit];
            ubarij=new int[jobs][this.NbStage];
            wbarij=new int[jobs][this.NbStage];
            this.bar_Dij=new double[this.NbStage];
            int B = 0;
            for (int j = 0; j < this.NbStage - 1; j++) {
                B = max(B, this.Uj[j]);
            }
            this.bar_Sij = new double [B][this.NbStage];
            BJsol=new int[this.NbStage-1][Maxit];

            //run master problem
            this.MasterBenders(Stolletz,AlterID);

            //run simulation (as subproblem)
            numit = 0;
            this.SIM_Serial_BAS(jobs,W,tij,this.ubarij,vij,this.wbarij,sij,this.bar_Sij,this.bar_Dij, false);
            for(int j=0;j<this.NbStage-1;j++){
                this.BJsol[j][0]=this.Buffer[j];
            }
            this.TH = (double)1/this.OverallCT;
            mM_value mm = new test_mainfunc.mM_value(this.NbStage, jobs, this.Uj, this.Lj);
            mm.get_mM_value(tij);

            while((this.THstar-this.TH > 0.0001) && (numit < this.Maxit-1))
            {
                numit++;
                long StartT = System.currentTimeMillis();
                this.AddFeasibilityCut(jobs,numit, tij, theta, Stolletz, mm,AlterID);
                long elapsedT = System.currentTimeMillis()- StartT;
                float elapsedTime = elapsedT/1000F;
                if(Stolletz)
                    stolletz_cplextime+=elapsedTime;
                else if(AlterID==3)
                    our_cplextime[0]+=elapsedTime;
                else if(AlterID==4)
                    our_cplextime[1]+=elapsedTime;
                else if(AlterID==5)
                    our_cplextime[2]+=elapsedTime;

                this.SIM_Serial_BAS(jobs,W,tij,this.ubarij,vij,this.wbarij,sij,this.bar_Sij,this.bar_Dij,false);
                this.TH = 1/this.OverallCT;
            }
            this.totit = numit;

            cplex.end();
        }


        catch (Exception exc) {
            exc.printStackTrace();
        }

    }


    //Initial master problem.
    //No feasibility cut is included.
    public void MasterBenders(
            boolean Stolletz,
            int AlterID
    ) {

        try {
            //Environment definition
            cplex = new IloCplex();

            //Variable definition
            BJ = cplex.intVarArray(this.NbStage-1,this.Lj, this.Uj);


            if(!Stolletz)
            {
                String label;
                for (int j=0; j<this.NbStage-1; j++)
                {
                    label = "b_"+(j+1);
                    BJ[j] = cplex.intVar(this.Lj[j],this.Uj[j],label);
                }

                if(AlterID==3 || AlterID==4){
                    DeltaBJrp = new IloNumVar[this.NbStage-1][];
                    DeltaBJrm = new IloNumVar[this.NbStage-1][];
                    for (int j=0; j<this.NbStage-1; j++)
                    {
                        DeltaBJrp[j] = cplex.intVarArray(this.Maxit,0,this.Uj[j]-this.Lj[j]);
                        DeltaBJrm[j] = cplex.intVarArray(this.Maxit,0,this.Uj[j]-this.Lj[j]);
                        for (int r=0; r<this.Maxit; r++)
                        {
                            label = "Deltap_" + (j+1) +'^'+ (r+1);
                            DeltaBJrp[j][r] = cplex.intVar(0,this.Uj[j]-this.Lj[j],label);
                            label = "Deltam_" + (j+1) +'^'+ (r+1);
                            DeltaBJrm[j][r] = cplex.intVar(0, this.Uj[j]-this.Lj[j],label);
                        }
                    }
                }
                else{
                    DeltaJKp = new IloNumVar[this.NbStage-1][];
                    DeltaJKm = new IloNumVar[this.NbStage-1][];
                    this.Alter5_Delta=new boolean[this.NbStage-1][];
                    for (int j=0; j<this.NbStage-1; j++)
                    {
                        DeltaJKp[j] = cplex.intVarArray(this.Uj[j],0,this.Uj[j]-this.Lj[j]);
                        DeltaJKm[j] = cplex.intVarArray(this.Uj[j],0,this.Uj[j]-this.Lj[j]);
                        for (int r=this.Lj[j]; r<this.Uj[j]; r++)
                        {
                            label = "Deltap_" + (j+1) +'^'+ (r+1);
                            DeltaJKp[j][r] = cplex.intVar(0, this.Uj[j]-r,label);
                            label = "Deltam_" + (j+1) +'^'+ (r+1);
                            DeltaJKm[j][r] = cplex.intVar(0, r-this.Lj[j],label);
                        }
                        this.Alter5_Delta[j]=new boolean[this.Uj[j]];
                        for(int k=0;k<this.Uj[j];k++)
                            this.Alter5_Delta[j][k]=false;
                    }
                }

                if(AlterID==3){
                    YJR=new IloNumVar[this.NbStage-1][];
                    for (int j=0; j<this.NbStage-1; j++)
                    {
                        YJR[j] = cplex.intVarArray(this.Maxit,0,1);
                        for (int r=0; r<this.Maxit; r++)
                        {
                            label = "Yrj_" + (j+1) +'^'+ (r+1);
                            YJR[j][r] = cplex.intVar(0, 1,label);
                        }
                    }
                }
                else{
                    YJK = new IloNumVar[this.NbStage-1][];
                    for (int j=0; j<this.NbStage-1; j++)
                    {
                        YJK[j] = cplex.intVarArray(this.Uj[j],0,1);
                        for (int r=0; r<this.Uj[j]; r++)
                        {
                            label = "Yjk_" + (j+1) +'^'+ (r+1);
                            YJK[j][r] = cplex.intVar(0, 1,label);
                        }
                    }


                    for (int j=0; j<this.NbStage-1; j++)
                    {
                        IloLinearNumExpr sumyjk_expr = cplex.linearNumExpr();
                        IloRange rng;
                        for (int k=1; k<this.Uj[0]; k++)
                        {
                            sumyjk_expr.addTerm(1,YJK[j][k]);
                        }
                        sumyjk_expr.addTerm(-1,BJ[j]);
                        rng = cplex.addEq(sumyjk_expr,0);
                        rng.setName("sumyjk" + j);
                    }

                    for (int j=0; j<this.NbStage-1; j++)
                    {
                        for (int k=1; k<this.Uj[0]; k++)
                        {
                            IloLinearNumExpr yjk_expr = cplex.linearNumExpr();
                            IloRange rng;
                            yjk_expr.addTerm(1,YJK[j][k]);
                            yjk_expr.addTerm(-1,YJK[j][k-1]);
                            rng = cplex.addLe(yjk_expr,0);
                            rng.setName("yjk" + j);
                        }
                    }
                    for (int j=0; j<this.NbStage-1; j++)
                    {
                        for(int k=0;k<=this.Lj[j];k++){
                            IloLinearNumExpr yjk_expr = cplex.linearNumExpr();
                            IloRange rng;
                            yjk_expr.addTerm(1,YJK[j][k]);
                            rng = cplex.addEq(yjk_expr,1);
                            rng.setName("yj0" + j);
                        }
                    }
                }


                /*if(AlterID==5){
                    for (int j=0; j<this.NbStage-1; j++){
                        for(int k=1;k<this.Uj[j];k++){
                            IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                            IloRange rng;
                            singlebj_expr.addTerm(1,DeltaJKp[j][k]);
                            singlebj_expr.addTerm(-1,BJ[j]);
                            for(int l=1;l<=k;l++)
                                singlebj_expr.addTerm(1,YJK[j][l]);
                            rng = cplex.addEq(singlebj_expr, 0) ;
                            rng.setName("def: Deltap_" + (j+1) +'^'+ (k+1));
                        }
                    }

                    for (int j=0; j<this.NbStage-1; j++){
                        for(int k=1;k<this.Uj[j];k++){
                            IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                            IloRange rng;
                            singlebj_expr.addTerm(1,DeltaJKm[j][k]);
                            for(int l=1;l<=k;l++)
                                singlebj_expr.addTerm(1,YJK[j][l]);
                            rng = cplex.addEq(singlebj_expr, k) ;
                            rng.setName("def: Deltam_" + (j+1) +'^'+ (k+1));
                        }
                    }
                }*/
            }

            if(Stolletz)
            {
                ZJK = new IloNumVar[this.NbStage-1][];
                for (int j=0; j<this.NbStage-1; j++)
                {
                    ZJK[j] = cplex.boolVarArray(this.Uj[0]);
                }
            }

            //objective function
            this.objective = cplex.linearNumExpr();
            for (int j=0; j<this.NbStage-1; j++)
                this.objective.addTerm(1,BJ[j]);


            if(Stolletz)
            {
                //constraint (8) Stolletz
                for (int j=0; j<this.NbStage-1; j++)
                {
                    IloLinearNumExpr sumzjk_expr = cplex.linearNumExpr();
                    IloRange rng;
                    for (int k=0; k<this.Uj[0]; k++)
                    {
                        sumzjk_expr.addTerm(1,ZJK[j][k]);
                    }
                    rng = cplex.addEq(sumzjk_expr,1.0);
                    rng.setName("sumk" + j);
                }

                //constraint (9) Stolletz
                for (int j=0; j<this.NbStage-1; j++)
                {
                    IloLinearNumExpr sumzjk_expr = cplex.linearNumExpr();
                    IloRange rng;
                    sumzjk_expr.addTerm(1,BJ[j]);
                    for (int k=0; k<this.Uj[0]; k++)
                    {
                        sumzjk_expr.addTerm(-k,ZJK[j][k]);
                    }
                    rng = cplex.addEq(sumzjk_expr,0);
                    rng.setName("cap" + j);
                }
            }

            cplex.addMinimize(this.objective);


            //E: RISOLUZIONE DEL MODELLO
            if (cplex.solve())
            {
                //save current solution
                for(int j=0;j<this.NbStage-1;j++)
                {
                    //this.BJsol[j][0] = (int) (cplex.getValue(this.BJ[j]));
                    this.BJsol[j][0] = this.Lj[j];
                    this.Buffer[j]=this.Lj[j];
                }
            }
            else
            {
                System.out.println("problem not solved");
            }

        }
        catch (Exception exc) {
            exc.printStackTrace();
        }

    }

    public void AddFeasibilityCut(
            int n,
            int numint,
            double[][] tij,
            double thetabar,
            boolean Stolletz,
            mM_value mm,
            int AlterID) {

        //i define the last already explored iteration as nnint
        int nnint = numint -1;

        try {


            if(!Stolletz)
            {

                if(AlterID==3||AlterID==4){
                    for (int j=0; j<this.NbStage-1; j++)
                    {
                        IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                        IloRange rng;
                        singlebj_expr .addTerm(1,BJ[j]);
                        singlebj_expr .addTerm(-1,DeltaBJrp[j][nnint]);
                        singlebj_expr .addTerm(+1,DeltaBJrm[j][nnint]);
                        rng = cplex.addEq(singlebj_expr, BJsol[j][nnint]) ;
                        rng.setName("def: bj_" + (j+1) +'^'+ (nnint+1));
                    }

                    if(AlterID==3){
                        for (int j=0; j<this.NbStage-1; j++)
                        {
                            IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                            IloRange rng;
                            singlebj_expr .addTerm(this.Uj[j]-BJsol[j][nnint],YJR[j][nnint]);
                            singlebj_expr .addTerm(-1,DeltaBJrp[j][nnint]);
                            rng = cplex.addGe(singlebj_expr, 0) ;
                            rng.setName("def: Deltap_" + (j+1) +'^'+ (nnint+1));
                        }

                        for (int j=0; j<this.NbStage-1; j++)
                        {
                            IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                            IloRange rng;
                            singlebj_expr .addTerm(this.Lj[j]-this.BJsol[j][nnint],YJR[j][nnint]);
                            singlebj_expr .addTerm(-1,DeltaBJrm[j][nnint]);
                            rng = cplex.addGe(singlebj_expr, this.Lj[j]-this.BJsol[j][nnint]) ;
                            rng.setName("def: Deltam_" + (j+1) +'^'+ (nnint+1));
                        }
                    }
                    else{
                        for (int j=0; j<this.NbStage-1; j++)
                        {
                            IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                            IloRange rng;
                            singlebj_expr .addTerm(this.Uj[j]-BJsol[j][nnint],YJK[j][min(BJsol[j][nnint]+1,this.Uj[j]-1)]);
                            singlebj_expr .addTerm(-1,DeltaBJrp[j][nnint]);
                            rng = cplex.addGe(singlebj_expr, 0) ;
                            rng.setName("def: Deltap_" + (j+1) +'^'+ (nnint+1));
                        }

                        for (int j=0; j<this.NbStage-1; j++)
                        {
                            IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                            IloRange rng;
                            singlebj_expr .addTerm(this.Lj[j]-this.BJsol[j][nnint],YJK[j][BJsol[j][nnint]]);
                            singlebj_expr .addTerm(-1,DeltaBJrm[j][nnint]);
                            rng = cplex.addGe(singlebj_expr, this.Lj[j]-this.BJsol[j][nnint]) ;
                            rng.setName("def: Deltam_" + (j+1) +'^'+ (nnint+1));
                        }
                    }
                }
                else{
                    for(int j=0;j<this.NbStage-1;j++){
                        if(!this.Alter5_Delta[j][this.Buffer[j]]){

                                    IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                                    IloRange rng;
                                    singlebj_expr.addTerm(1,DeltaJKp[j][this.Buffer[j]]);
                                    singlebj_expr.addTerm(-1,BJ[j]);
                                    for(int l=1;l<=this.Buffer[j];l++)
                                        singlebj_expr.addTerm(1,YJK[j][l]);
                                    rng = cplex.addEq(singlebj_expr, 0) ;
                                    rng.setName("def: Deltap_" + (j+1) +'^'+ (this.Buffer[j]));


                                    singlebj_expr = cplex.linearNumExpr();
                                    singlebj_expr.addTerm(1,DeltaJKm[j][this.Buffer[j]]);
                                    for(int l=1;l<=this.Buffer[j];l++)
                                        singlebj_expr.addTerm(1,YJK[j][l]);
                                    rng = cplex.addEq(singlebj_expr, this.Buffer[j]) ;
                                    rng.setName("def: Deltam_" + (j+1) +'^'+ (this.Buffer[j]));

                            this.Alter5_Delta[j][this.Buffer[j]]=true;
                        }
                    }



                }

                ///************************* feasibility cut ****************************************///////////////////////////////////
                for(int j=0;j<this.NbStage -1;j++)                                                                              ////////
                {                                                                                                               ////////
                    for (int i=0;i<n;i++)                                                                                       ////////
                        {                                                                                                       ////////
                       DeltapPar[j][nnint]= DeltapPar[j][nnint] + (double)this.wbarij[i][j]*mm.Mijk[i][j][BJsol[j][nnint]];     ////////
                       DeltamPar[j][nnint]= DeltamPar[j][nnint] + (double)this.wbarij[i][j]*mm.mijk[i][j][BJsol[j][nnint]];     ////////
                    }                                                                                                           ////////
                }                                                                                                               ////////                                                                                  ///////////////////////////////
                                                                                                                                ////////
                                                                                                                                ////////
                double tijpar=0.0;                                                                                              ////////
                for(int j=0;j<this.NbStage;j++)                                                                                 ////////
                {                                                                                                               ////////
                    for (int i=0;i<n;i++)                                                                                       ////////
                    {                                                                                                           ////////
                        tijpar = tijpar + tij[i][j]*(double)this.ubarij[i][j];                                                  ////////
                    }                                                                                                           ////////
                }                                                                                                               ////////
                                                                                                                                ////////
                                                                                                                                ////////
                double thetapar = thetabar/this.THstar;//E: if theta is an array, then change this formulation                  ////////
                IloLinearNumExpr sumBJcut_expr = cplex.linearNumExpr();                                                         ////////
                IloRange rng;                                                                                                   ////////
                ///*************************MP-alter3/4 ****************************************////////////////////////////////////////
                if(AlterID==3 || AlterID==4){
                    for (int j=0; j<this.NbStage-1; j++)
                    {                                                                                                               ////////
                        sumBJcut_expr.addTerm(-DeltapPar[j][nnint],DeltaBJrp[j][nnint]);                                            ////////
                        sumBJcut_expr.addTerm(+DeltamPar[j][nnint],DeltaBJrm[j][nnint]);                                            ////////
                    }
                }
                else{
                    ///*************************MP-alter5 ****************************************//////////////////////////////////////////
                    for (int j=0; j<this.NbStage-1; j++)                                                                            ////////
                    {                                                                                                               ////////
                        sumBJcut_expr.addTerm(-DeltapPar[j][nnint],DeltaJKp[j][BJsol[j][nnint]]);                                   ////////
                        sumBJcut_expr.addTerm(+DeltamPar[j][nnint],DeltaJKm[j][BJsol[j][nnint]]);                                   ////////
                    }
                }
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                                                                                                                ////////
                this.resc[nnint] = thetapar - tijpar;                                                                           ////////
                rng = cplex.addLe(sumBJcut_expr,resc[nnint]);                                                                   ////////
                rng.setName("feascut of iter: "+ (nnint+1));                                                                    ////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            }

            //if we don't want to use (25) but combinatorial cut:
            if(Stolletz)
            {
                IloLinearNumExpr sumBJcut_expr = cplex.linearNumExpr();
                IloRange rng;
                for (int j=0; j<this.NbStage-1; j++)
                {
                    for (int k= (int) this.Buffer[j] +1; k<this.Uj[0]; k++)
                    {
                        sumBJcut_expr.addTerm(1,ZJK[j][k]);
                    }
                }
                rng = cplex.addGe(sumBJcut_expr,1.0);
                rng.setName("combcut of iter"+nnint);
            }

            //Model solution
            if (cplex.solve())
            {


                //finished=true;
                String program = System.getProperty("user.dir");
                String prova = program + "\\OUTPUT\\model.lp";
                cplex.exportModel(prova);
                //save current solution
                for(int j=0;j<this.NbStage-1;j++)
                {
                    //MP-Alter3
                    this.BJsol[j][nnint+1] =  (int)  (cplex.getValue(this.BJ[j])+0.1);
                    Buffer[j] =  (int) (cplex.getValue(this.BJ[j])+0.1);

                    //MP-Alter0
                    /*this.BJsol[j][nnint+1] =  (int) (this.BJsol[j][nnint] +  cplex.getValue(this.DeltaBJrp[j][nnint]) - cplex.getValue(this.DeltaBJrm[j][nnint]));
                    Buffer[j]=(int) (this.BJsol[j][nnint] +  cplex.getValue(this.DeltaBJrp[j][nnint]) - cplex.getValue(this.DeltaBJrm[j][nnint]));*/

                }

            }
            else
            {
                System.out.println("problem not solved");
            }

        }
        catch (Exception exc) {
            exc.printStackTrace();
        }
    }


}
