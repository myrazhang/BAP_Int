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

import org.apache.commons.math3.distribution.NormalDistribution;
// import org.apache.commons.math3.distribution.UniformRealDistribution;
// import org.apache.commons.math3.distribution.WeibullDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;

public class SerialLine {
    public int NbStage;
    public int[] Buffer;
    public StochNum[] CT;
    public int Uj;
    public double THstar;

    //Sim_Output
    public double OverallCT;
    public double TH;
    public double[] Pblock;
    public double[] Pstarve;
    public int Maxit = 100;

    //cplex public variables
    public IloCplex cplex;
    //ERICA: CHECK IF NEEDED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    public IloNumVar[] BJ;
    public IloNumVar[][] DeltaBJrp;
    public IloNumVar[][] DeltaBJrm;
    public IloLinearNumExpr objective;
    public int BJsol[][];
    public double[][][] ubarij;
    public double[][][] wbarij;
    public double thetabar;
    public int totit;


    // Create a new test_mainfunc.SerialLine from txt file
    public SerialLine(InputStream system) {

        Scanner scanner=new Scanner (system);
        scanner.useDelimiter("\\s+");

        System.out.println();

        scanner.next();
        this.NbStage=scanner.nextInt();

        scanner.next();
        this.Uj=scanner.nextInt();

        this.Buffer = new int[this.NbStage - 1];
        scanner.next();
        for (int j = 0; j < this.NbStage - 1; j++)
        {
            this.Buffer[j]=scanner.nextInt();
        }

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
            // else if(CT[j].distribution.equals("Exp")){...}
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
            double thetabar,
            double[][] bar_Sij,
            double[] bar_Dij,
            int iter,
            boolean SteadyState,
            boolean bendersuse) {
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

        if(bendersuse)
        {
            for(int j=0;j< NbStage -1;j++){
                Buffer[j]=this.BJsol[j][iter];
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

            for(int j=0;j<NbStage-1;j++)
            {
                System.out.println("buffer at " + j + " is: " + Buffer[j]);
            }

            this.OverallCT = (Dij[N - 1][NbStage - 1] - Dij[W - 1][NbStage - 1]) / (double)(N - W);
            //***************   SAVE bar_Sij, bar_Dij    ***************
            for (int j = 0; j < NbStage; j++) {
                bar_Dij[j] = Dij[W - 1][j];
                if (j != NbStage -1) {
                    for (int i = 0 ; i < Buffer[j]; i++) {
                        bar_Sij[i][j+1] = Sij[W- Buffer[j]+i][j+1];
                    }
                }
            }


        }
        else {
            //***************   Fix event times in transient period   ***************
            for (int j = 0; j < NbStage; j++) {
                Dij[W - 1][j] = bar_Dij[j];
                if (j != NbStage - 1) {
                    for (int i = 0; i < Buffer[j]; i++) {
                        Sij[W - Buffer[j] + i][j + 1] = bar_Sij[i][j+1];
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
            this.OverallCT = (Dij[N-1][NbStage-1]- Dij[W-1][NbStage-1])/ (double)(N-W);
            this.ubarij = new double [N][NbStage][this.Maxit];
            this.wbarij = new double [N][NbStage][this.Maxit];
            for (int j=0; j<this.NbStage;j++)
            {
                for(int i=0; i<N;i++)
                {
                    this.ubarij[i][j][iter]=uij[i][j];
                    this.wbarij[i][j][iter]=wij[i][j];
                }
            }
            /////////////////DEFINE HOW MUCH IS THETA FROM THE SIMULATION
            thetabar=1;
        }


    }

    public void solveBender(
            int jobs,
            int W,
            double[][] tij,
            int[][] uij,
            int[][] vij,
            int[][] wij,
            int[][] sij,
            double[][] bar_Sij,
            double[] bar_Dij,
            boolean Stolletz) {

        double theta = 0;

        try {

            //run master problem
            this.MasterBenders();
            //run simulation (as subproblem)
            int numit = 0;
            this.SIM_Serial_BAS(jobs,W,tij,uij,vij,wij,sij,theta,bar_Sij,bar_Dij, numit,false,true);
            System.out.println("CT is:" + this.OverallCT );
            this.TH = 1/this.OverallCT;
            while((this.TH < this.THstar) && (numit < this.Maxit-1))
            {
                numit++;
                this.AddFeasibilityCut(jobs,numit, tij, theta, Stolletz);
                this.SIM_Serial_BAS(jobs,W,tij,uij,vij,wij,sij,theta,bar_Sij,bar_Dij,numit,false,true);
                this.TH = 1/this.OverallCT;
                System.out.println("simulated TH" + this.TH );
            }
            this.totit = numit;

            cplex.end();

            System.out.println("num. iteraz: " + numit);
        }

        catch (Exception exc) {
            exc.printStackTrace();
        }

    }

    public void MasterBenders() {

        int J = this.NbStage;
        int Uj = this.Uj;

        try {
            //E: Model definition
            cplex = new IloCplex();

            //E: Variable definition

            BJ = cplex.intVarArray(J-1,0, Uj);
            DeltaBJrp = new IloNumVar[J-1][];

            for (int j=0; j<J-1; j++)
            {
                DeltaBJrp[j] = cplex.intVarArray(this.Maxit,0,Uj);
            }
            DeltaBJrm = new IloNumVar[J-1][];
            for (int j=0; j<J-1; j++)
            {
                DeltaBJrm[j] = cplex.intVarArray(this.Maxit,0,Uj);
            }

            //E: objective function
            this.objective = cplex.linearNumExpr();
            for (int j=0; j<J-1; j++)
            {
                this.objective.addTerm(1,BJ[j]);
                for(int r = 0;r<this.Maxit;r++){
                    this.objective.addTerm(1, DeltaBJrp[j][r]);
                    this.objective.addTerm(1, DeltaBJrm[j][r]);
                }
            }
            cplex.addMinimize(this.objective);

            //E: Delta definition constraints (22)
//            for (int j=0; j<J-1; j++)
////            {
////                    IloLinearNumExpr sumBJ_expr = cplex.linearNumExpr();
////                    IloRange rng;
////                    sumBJ_expr.addTerm(1,DeltaBJrp[j][0]);
////                    sumBJ_expr.addTerm(-1,BJ[j]);
////                    rng = cplex.addGe(sumBJ_expr, -1*BJsol[j][0]);
////                    rng.setName("deltaBJp" + j + "0");
////
////            }
////
////            //E: Delta definition constraints (23)
////            for (int j=0; j<J-1; j++)
////            {
////                    IloLinearNumExpr sumDeltaBJ_expr = cplex.linearNumExpr();
////                    IloRange rng;
////                    sumDeltaBJ_expr.addTerm(1,DeltaBJrm[j][0]);
////                    sumDeltaBJ_expr.addTerm(1,BJ[j]);
////                    rng = cplex.addGe(sumDeltaBJ_expr, BJsol[j][0]);
////                    rng.setName("deltaBJm" + j + "0");
////            }

            //E: RISOLUZIONE DEL MODELLO
            if (cplex.solve())
            {
                System.out.println("obj = "+cplex.getObjValue());
                //finished=true;

                //save current solution
                this.BJsol = new int[J-1][this.Maxit];
                for(int j=0;j<J-1;j++)
                {
                    //System.out.println(cplex.getValue(this.BJ[j]));
                    this.BJsol[j][0] = (int) cplex.getValue(this.BJ[j]);
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
            int jobs,
            int numint,
            double[][] tij,
            double thetabar,
            boolean Stolletz) {

        int J = this.NbStage;
        int n = jobs;
        //i define the last already explored iteration as nnint
        int nnint = numint -1;

        try {
            ////E: insert new nnint-delta in objective function
            //this.objective = cplex.linearNumExpr();
//            for (int j=0; j<J-1; j++)
//            {
//                this.objective.addTerm(1, DeltaBJrp[j][nnint]);
//                this.objective.addTerm(1, DeltaBJrm[j][nnint]);
//            }
//            cplex.addMinimize(this.objective);


            //E: Delta definition constraints (22)
            for (int j=0; j<J-1; j++)
            {
                IloLinearNumExpr sumBJ_expr = cplex.linearNumExpr();
                IloRange rng;
                sumBJ_expr.addTerm(1,DeltaBJrp[j][nnint]);
                sumBJ_expr.addTerm(-1,BJ[j]);
                /////////////////CHECK THE CORRECT ITERATION OF BJSOL
                rng = cplex.addGe(sumBJ_expr, -1*BJsol[j][nnint]);
                rng.setName("deltaBJp" + j + nnint);

            }

            //E: Delta definition constraints (23)
            for (int j=0; j<J-1; j++)
            {
                IloLinearNumExpr sumDeltaBJ_expr = cplex.linearNumExpr();
                IloRange rng;
                sumDeltaBJ_expr.addTerm(1,DeltaBJrm[j][nnint]);
                sumDeltaBJ_expr.addTerm(1,BJ[j]);
                /////////////////CHECK THE CORRECT ITERATION OF BJSOL
                rng = cplex.addGe(sumDeltaBJ_expr, BJsol[j][nnint]);
                rng.setName("deltaBJm" + j + nnint);
            }

            if(!Stolletz)
            {
                double[] DeltapPar;
                double[] DeltamPar;
                DeltapPar = new double[J-1];
                DeltamPar = new double[J-1];

                mM_value mm = new mM_value(J, n, this.Uj, 0);

                for(int j=0;j<J -1;j++)
                {
                    for (int i=0;i<n;i++)
                    {
                        DeltapPar[j]= DeltapPar[j] + this.wbarij[i][j][nnint]*mm.Mijk[i][j][BJsol[j][nnint]];
                        DeltamPar[j]= DeltamPar[j] + this.wbarij[i][j][nnint]*mm.mijk[i][j][BJsol[j][nnint]];
                    }
                }

                double tijpar=0.0;
                for(int j=0;j<J -1;j++)
                {
                    for (int i=0;i<n;i++)
                    {
                        tijpar = tijpar + tij[i][j]*this.ubarij[i][j][nnint];
                    }
                }

                double thetapar = thetabar/this.THstar;


                IloLinearNumExpr sumBJcut_expr = cplex.linearNumExpr();
                IloRange rng;
                for (int j=0; j<J-1; j++)
                {
                    sumBJcut_expr.addTerm(-DeltapPar[j],DeltaBJrp[j][nnint]);
                    sumBJcut_expr.addTerm(-DeltamPar[j],DeltaBJrm[j][nnint]);
                }
                double resc = thetapar - tijpar;
                rng = cplex.addLe(sumBJcut_expr,resc);
                rng.setName("feascut of iter"+nnint);

            }
            if(Stolletz)
            {
                IloLinearNumExpr sumBJcut_expr = cplex.linearNumExpr();
                IloRange rng;
                int lastcapsum = 0;
                for(int j=0;j< J-1;j++)
                {
                    sumBJcut_expr.addTerm(1,BJ[j]);
                    lastcapsum = lastcapsum + this.BJsol[j][nnint];
                }
                lastcapsum=lastcapsum+1;
                rng = cplex.addGe(sumBJcut_expr,lastcapsum);
                rng.setName("combcut of iter"+nnint);
            }

            //E: Model solution
            if (cplex.solve())
            {
                System.out.println("obj = "+cplex.getObjValue());
                //finished=true;
                System.out.println("obj = "+cplex.getObjValue());

                //save current solution
                for(int j=0;j<this.NbStage-1;j++)
                {
                    this.BJsol[j][nnint+1] = (int) cplex.getValue(this.BJ[j]);
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