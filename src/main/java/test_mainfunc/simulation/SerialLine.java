package test_mainfunc.simulation;

import java.io.InputStream;
import java.util.Scanner;

public class SerialLine {
    public int nbStage;
    public int[] buffer;
    public StochNum[] CT;
    public SimulationBAS mySimulation;


    //Sim_Output
    public double OverallCT;
    public double TH;


    // constructing a SINGLE serial line from a input file is not needed currently. For example:
    /*public SerialLine(InputStream system) {}*/

    public SerialLine(){}

    public SerialLine(InputStream inputFile){
        Scanner scanner=new Scanner (inputFile);
        scanner.useDelimiter("\\s+");

        scanner.next();
        this.nbStage=scanner.nextInt();
        //System.out.println(this.nbStage);

        buffer=new int[nbStage];
        CT=new StochNum[nbStage+1];

        //for(int j=1;j<=7;j++)
            scanner.next();
        for(int j=1;j<=nbStage-1;j++)
        { buffer[j]=scanner.nextInt();}
        scanner.next();
        scanner.next();
        for(int j=1;j<=nbStage;j++){
            CT[j]=new StochNum();
            CT[j].distribution=scanner.next();
        }

        scanner.next();
        for(int j=1;j<=nbStage;j++){
            CT[j].para1=scanner.nextDouble();
        }

        scanner.next();
        for(int j=1;j<=nbStage;j++){
            CT[j].para2=scanner.nextDouble();
        }

        scanner.next();
        for(int j=1;j<=nbStage;j++){
            CT[j].para3=scanner.nextDouble();
        }

        scanner.next();
        for(int j=1;j<=nbStage;j++){
            CT[j].para4=scanner.nextDouble();
        }
    }

    // Processing time generation
    public void procTimeGeneration(int N, double[][] pij,int seed){
        for(int j=1;j<=this.nbStage;j++){
            double[] pi=new double[N+1];
                this.CT[j].iidGeneration(N,pi,seed+j);
            for(int i=1;i<=N;i++){
                double a=pi[i];
                pij[i][j]=a;
            }
        }
    }

    public class SimulationBAS{

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


        double[][] tij;
        public int[][] uij;
        public int[][] wij;
        int[][] vij;
        int[][] sij;
        int N;
        int W;
        private double[][] bar_Sij;
        private double[][] bar_Dij;
        private double[][] Sij;
        public double[][] Dij;
        private int[][] bsij;
        private int[][] buij;
        private int[][] bvij;
        private int[][] bwij;

        public SimulationBAS(int N,int W,double[][] tij){
            this.N=N;
            this.W=W;
            this.tij=tij;
            this.Sij=new double[N+1][nbStage+1];
            this.Dij=new double[N+1][nbStage+1];
            this.bsij = new int[N+1][nbStage+1];
            this.buij = new int[N+1][nbStage+1];
            this.bvij = new int[N+1][nbStage+1];
            this.bwij = new int[N+1][nbStage+1];
            this.sij = new int[N+1][nbStage+1];
            this.uij = new int[N+1][nbStage+1];
            this.vij = new int[N+1][nbStage+1];
            this.wij = new int[N+1][nbStage+1];
            this.bar_Dij=new double[W+1][nbStage+1];
            this.bar_Sij=new double[W+1][nbStage+1];
        }

        private void initialization(){
            for (int i = 1; i <= this.N; i++){
                for (int j = 1; j <= nbStage; j++){
                    this.Sij[i][j] = 0;
                    this.Dij[i][j] = 0;
                }
            }
            for (int i = 1; i <= this.N; i++){
                for (int j = 1; j <= nbStage; j++){
                    bsij[i][j] = 0;
                    buij[i][j] = 0;
                    bvij[i][j] = 0;
                    bwij[i][j] = 0;
                }
            }
        }

        public void simBAS(boolean steadyState){
            initialization();

            if(!steadyState){
                for (int i = 1; i <= this.N; i++) {
                    for (int j = 1; j <= nbStage; j++){
                        this.start_BAS(i,j);
                        this.departure_BAS(i,j);
                    }
                }
                for (int j = 1; j <= nbStage; j++) {
                    for (int i = 1; i <= this.W; i++) {
                        this.bar_Dij[i][j]=this.Dij[i][j];
                        this.bar_Sij[i][j]=this.Sij[i][j];
                    }
                }
                OverallCT = (this.Dij[this.N][nbStage] - this.Dij[this.W][nbStage]) / (double)(this.N - this.W);
                TH=(double) 1/OverallCT;
            }
            else{
                for (int j = 1; j <= nbStage; j++) {
                    for (int i = 1; i <= this.W; i++) {
                        this.Dij[i][j] = this.bar_Dij[i][j];
                        this.Sij[i][j] = this.bar_Sij[i][j];
                    }
                }

                for (int i = this.W+1; i <= this.N; i++){
                    for (int j = 1; j <= nbStage; j++){
                        this.start_BAS(i,j);
                        this.departure_BAS(i,j);
                    }
                }

                OverallCT = (this.Dij[N][nbStage]- Dij[W][nbStage])/ (double)(this.N-this.W);
                TH=(double) 1/OverallCT;

            }
        }

        public void simDualBAS(boolean steadyState) {



            //*****************************************************************************
            //*****************************************************************************
            //** The notations of u,v,w,s is consistent with the journal paper ************
            //** on throughput improvement (2019).                             ************
            //*****************************************************************************
            //*****************************************************************************

            for (int i = 1; i <= this.N; i++){
                for (int j = 1; j <= nbStage; j++){
                    this.sij[i][j] = 0;
                    this.uij[i][j] = 0;
                    this.vij[i][j] = 0;
                    this.wij[i][j] = 0;
                }
            }

            this.simBAS(steadyState);
            if (!steadyState) {

                /////////////////////////////////////////////////////////////////////////////////////////////
                ///// Simulate from i=0 /////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////
                /////buij[i][j] = 1: D[i][j] trigerred by S[i][j]               /////////////////////////////
                /////bwij[i][j] = 1: D[i][j] trigerred by S[i-b_j][j+1] BLOCKAGE ////////////////////////////
                /////bvij[i][j] = 1: S[i][j] trigerred by D[i-1][j]             /////////////////////////////
                /////bsij[i][j] = 1: S[i][j] trigerred by D[i][j-1]   STARVATION ////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////
                for (int i = this.N; i >= 1; i--) {
                    for (int j = nbStage; j >=1 ; j--) {
                        this.departure_dual(i,j);
                        this.start_dual(i,j);
                    }
                }
            }
            else {
                for (int i = this.N; i >= this.W+1; i--){
                    for (int j = nbStage; j >= 1; j--){
                        this.departure_dual(i,j);
                        this.start_dual(i,j);
                    }
                }
            }
        }

        private void start_BAS(int i0,int j0){

            if(i0==1&& j0==1){
                this.Sij[i0][j0]=0;
                this.bsij[i0][j0] = 1;
            }
            else if(i0==1){
                this.Sij[i0][j0]=this.Dij[i0][j0-1];
                this.bsij[i0][j0] = 1;
            }
            else if(j0==1){
                this.Sij[i0][j0]=this.Dij[i0-1][j0];
                this.bvij[i0][j0] = 1;
            }
            else{
                if (this.Dij[i0 - 1][j0] >= this.Dij[i0][j0 - 1]){
                    this.Sij[i0][j0] = this.Dij[i0 - 1][j0];
                    this.bvij[i0][j0] = 1;
                }
                else{
                    this.Sij[i0][j0] = this.Dij[i0][j0 - 1];
                    this.bsij[i0][j0] = 1;
                }
            }

        }

        private void departure_BAS(int i0,int j0){

            if (j0 == nbStage) {
                this.Dij[i0][j0] = this.Sij[i0][j0] + this.tij[i0][j0];
                this.buij[i0][j0] = 1;
            }
            else if (i0 <= buffer[j0] + 1) {
                this.Dij[i0][j0] = this.Sij[i0][j0] + this.tij[i0][j0];
                this.buij[i0][j0] = 1;
            }

            else {
                if (this.Sij[i0][j0] + this.tij[i0][j0] >= this.Dij[i0 - buffer[j0]-1][j0 + 1]) {
                    this.Dij[i0][j0] = this.Sij[i0][j0] + this.tij[i0][j0];
                    this.buij[i0][j0] = 1;
                }
                else{
                    this.Dij[i0][j0] = this.Dij[i0 - buffer[j0]-1][j0 + 1];
                    bwij[i0][j0] = 1;
                }
            }
        }

        private void departure_dual(int i0, int j0 ){

            int flowOut = 0;
            if(i0 < this.N){
                flowOut += this.vij[i0 + 1][j0];
            }
            if( j0 < nbStage){
                flowOut += this.sij[i0][j0 + 1];
            }
            if(j0 >= 2 && i0 <= this.N-buffer[j0-1]-1){
                flowOut += this.wij[i0+buffer[j0-1]+1][j0-1];
            }

            if (i0 == this.N && j0 == nbStage){
                this.uij[i0][j0] = 1;
            }
            else{
                this.uij[i0][j0] = flowOut * this.buij[i0][j0];
                this.wij[i0][j0] = flowOut * this.bwij[i0][j0];
            }

        }

        private void start_dual(int i0,int j0){
            this.sij[i0][j0] = (this.uij[i0][j0])*this.bsij[i0][j0];
            this.vij[i0][j0] = (this.uij[i0][j0])*this.bvij[i0][j0];
        }

    }



















    //Initial master problem.
    //No feasibility cut is included.
    /*public void MasterBenders(
            boolean Stolletz,
            int AlterID
    ) {

        try {

            if(!Stolletz)
            {
                String label;
                for (int j=0; j<this.nbStage-1; j++)
                {
                    label = "b_"+(j+1);
                    BJ[j] = cplex.intVar(this.Lj,this.Uj,label);// BendersBAP
                }

                if(AlterID==3 || AlterID==4){ //Alter3 & Alter4
                    DeltaBJrp = new IloNumVar[this.nbStage-1][];//Alter3 & Alter4
                    DeltaBJrm = new IloNumVar[this.nbStage-1][];//Alter3 & Alter4
                    for (int j=0; j<this.nbStage-1; j++)//Alter3 & Alter4
                    {//Alter3 & Alter4
                        DeltaBJrp[j] = cplex.intVarArray(this.Maxit,0,this.Uj-this.Lj);//Alter3 & Alter4
                        DeltaBJrm[j] = cplex.intVarArray(this.Maxit,0,this.Uj-this.Lj);//Alter3 & Alter4
                        for (int r=0; r<this.Maxit; r++)//Alter3 & Alter4
                        {//Alter3 & Alter4
                            label = "Deltap_" + (j+1) +'^'+ (r+1);//Alter3 & Alter4
                            DeltaBJrp[j][r] = cplex.intVar(0,this.Uj-this.Lj,label);//Alter3 & Alter4
                            label = "Deltam_" + (j+1) +'^'+ (r+1);//Alter3 & Alter4
                            DeltaBJrm[j][r] = cplex.intVar(0, this.Uj-this.Lj,label);//Alter3 & Alter4
                        }//Alter3 & Alter4
                    }//Alter3 & Alter4
                }//Alter3 & Alter4
                else{//Alter5
                    DeltaJKp = new IloNumVar[this.nbStage-1][];//Alter5
                    DeltaJKm = new IloNumVar[this.nbStage-1][];//Alter5
                    this.Alter5_Delta=new boolean[this.nbStage-1][];
                    for (int j=0; j<this.nbStage-1; j++)//
                    {//
                        DeltaJKp[j] = cplex.intVarArray(this.Uj,0,this.Uj-this.Lj);//Alter5
                        DeltaJKm[j] = cplex.intVarArray(this.Uj,0,this.Uj-this.Lj);//Alter5
                        for (int r=this.Lj; r<this.Uj; r++)//Alter5
                        {//Alter5
                            label = "Deltap_" + (j+1) +'^'+ (r+1);//Alter5
                            DeltaJKp[j][r] = cplex.intVar(0, this.Uj-r,label);//Alter5
                            label = "Deltam_" + (j+1) +'^'+ (r+1);//Alter5
                            DeltaJKm[j][r] = cplex.intVar(0, r-this.Lj,label);//Alter5
                        }//Alter5
                        this.Alter5_Delta[j]=new boolean[this.Uj];
                        for(int k=0;k<this.Uj;k++)
                            this.Alter5_Delta[j][k]=false;
                    }
                }

                if(AlterID==3){
                    YJR=new IloNumVar[this.nbStage-1][];////Alter3
                    for (int j=0; j<this.nbStage-1; j++)//Alter3
                    {//Alter3
                        YJR[j] = cplex.intVarArray(this.Maxit,0,1);//Alter3
                        for (int r=0; r<this.Maxit; r++)//Alter3
                        {//Alter3
                            label = "Yrj_" + (j+1) +'^'+ (r+1);//Alter3
                            YJR[j][r] = cplex.intVar(0, 1,label);//Alter3
                        }//Alter3
                    }//Alter3
                }//Alter3
                else{
                    YJK = new IloNumVar[this.nbStage-1][];//Alter4 & 5
                    for (int j=0; j<this.nbStage-1; j++)//Alter4 & 5
                    {
                        YJK[j] = cplex.intVarArray(this.Uj,0,1);//Alter4 & 5
                        for (int r=0; r<this.Uj; r++)//Alter4 & 5
                        {
                            label = "Yjk_" + (j+1) +'^'+ (r+1);//Alter4 & 5
                            YJK[j][r] = cplex.intVar(0, 1,label);//Alter4 & 5
                        }
                    }


                    for (int j=0; j<this.nbStage-1; j++)//Alter4 & 5
                    {
                        IloLinearNumExpr sumyjk_expr = cplex.linearNumExpr();//Alter4 & 5
                        IloRange rng;//Alter4 & 5
                        for (int k=1; k<this.Uj; k++)//Alter4 & 5
                        {
                            sumyjk_expr.addTerm(1,YJK[j][k]);//Alter4 & 5
                        }
                        sumyjk_expr.addTerm(-1,BJ[j]);//Alter4 & 5
                        rng = cplex.addEq(sumyjk_expr,0);//Alter4 & 5
                        rng.setName("sumyjk" + j);//Alter4 & 5
                    }

                    for (int j=0; j<this.nbStage-1; j++)//Alter4 & 5
                    {
                        for (int k=1; k<this.Uj; k++)//Alter4 & 5
                        {
                            IloLinearNumExpr yjk_expr = cplex.linearNumExpr();//Alter4 & 5
                            IloRange rng;//Alter4 & 5
                            yjk_expr.addTerm(1,YJK[j][k]);//Alter4 & 5
                            yjk_expr.addTerm(-1,YJK[j][k-1]);//Alter4 & 5
                            rng = cplex.addLe(yjk_expr,0);//Alter4 & 5
                            rng.setName("yjk" + j);//Alter4 & 5
                        }
                    }
                    for (int j=0; j<this.nbStage-1; j++)//Alter4 & 5
                    {
                        for(int k=0;k<=this.Lj;k++){//Alter4 & 5
                            IloLinearNumExpr yjk_expr = cplex.linearNumExpr();//Alter4 & 5
                            IloRange rng;//Alter4 & 5
                            yjk_expr.addTerm(1,YJK[j][k]);//Alter4 & 5
                            rng = cplex.addEq(yjk_expr,1);//Alter4 & 5
                            rng.setName("yj0" + j);//Alter4 & 5
                        }
                    }
                }
/////////////////////////////////////////以下代码已被上面的某一部分替换//////////////////////////////////////////////////

                if(AlterID==5){
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
                }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            }

            if(Stolletz)
            {
                ZJK = new IloNumVar[this.nbStage-1][];// Stolletz
                for (int j=0; j<this.nbStage-1; j++)
                {
                    ZJK[j] = cplex.boolVarArray(this.Uj);// Stolletz
                }
            }

            //objective function
            this.objective = cplex.linearNumExpr();// BendersBAP
            for (int j=0; j<this.nbStage-1; j++)// BendersBAP
                this.objective.addTerm(1,BJ[j]);// BendersBAP


            if(Stolletz)
            {
                //constraint (8) Stolletz
                for (int j=0; j<this.nbStage-1; j++) // Stolletz
                {
                    IloLinearNumExpr sumzjk_expr = cplex.linearNumExpr();// Stolletz
                    IloRange rng;// Stolletz
                    for (int k=0; k<this.Uj; k++)// Stolletz
                    {
                        sumzjk_expr.addTerm(1,ZJK[j][k]);// Stolletz
                    }
                    rng = cplex.addEq(sumzjk_expr,1.0);// Stolletz
                    rng.setName("sumk" + j);// Stolletz
                }

                //constraint (9) Stolletz
                for (int j=0; j<this.nbStage-1; j++)// Stolletz
                {
                    IloLinearNumExpr sumzjk_expr = cplex.linearNumExpr();// Stolletz
                    IloRange rng;// Stolletz
                    sumzjk_expr.addTerm(1,BJ[j]);// Stolletz
                    for (int k=0; k<this.Uj; k++)// Stolletz
                    {
                        sumzjk_expr.addTerm(-k,ZJK[j][k]);// Stolletz
                    }
                    rng = cplex.addEq(sumzjk_expr,0);// Stolletz
                    rng.setName("cap" + j);// Stolletz
                }
            }

            cplex.addMinimize(this.objective);// BendersBAP


            //E: RISOLUZIONE DEL MODELLO
            if (cplex.solve())
            {
                if(!Stolletz) {
                    String program = System.getProperty("user.dir");
                    String prova = program + "\\OUTPUT\\" + this.tempinstance + "model.lp";
                    cplex.exportModel(prova);
                }
                //save current solution
                int totcap = 0;
                for(int j=0;j<this.nbStage-1;j++)
                {
                    //this.BJsol[j][0] = (int) (cplex.getValue(this.BJ[j]));
                    this.BJsol[j][0] = this.Lj;
                    this.buffer[j]=this.Lj;
                    totcap = this.buffer[j] + totcap;
                }
                this.writer.write("it 0 OF: "+ totcap + "\r\n");
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


    /*public void AddFeasibilityCut(
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
                    for (int j=0; j<this.nbStage-1; j++)
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
                        for (int j=0; j<this.nbStage-1; j++)
                        {
                            IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                            IloRange rng;
                            singlebj_expr .addTerm(this.Uj-BJsol[j][nnint],YJR[j][nnint]);
                            singlebj_expr .addTerm(-1,DeltaBJrp[j][nnint]);
                            rng = cplex.addGe(singlebj_expr, 0) ;
                            rng.setName("def: Deltap_" + (j+1) +'^'+ (nnint+1));
                        }

                        for (int j=0; j<this.nbStage-1; j++)
                        {
                            IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                            IloRange rng;
                            singlebj_expr .addTerm(this.Lj-this.BJsol[j][nnint],YJR[j][nnint]);
                            singlebj_expr .addTerm(-1,DeltaBJrm[j][nnint]);
                            rng = cplex.addGe(singlebj_expr, this.Lj-this.BJsol[j][nnint]) ;
                            rng.setName("def: Deltam_" + (j+1) +'^'+ (nnint+1));
                        }
                    }
                    else{
                        for (int j=0; j<this.NbStage-1; j++)
                        {
                            IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                            IloRange rng;
                            singlebj_expr .addTerm(this.Uj-BJsol[j][nnint],YJK[j][min(BJsol[j][nnint]+1,this.Uj-1)]);
                            singlebj_expr .addTerm(-1,DeltaBJrp[j][nnint]);
                            rng = cplex.addGe(singlebj_expr, 0) ;
                            rng.setName("def: Deltap_" + (j+1) +'^'+ (nnint+1));
                        }

                        for (int j=0; j<this.NbStage-1; j++)
                        {
                            IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                            IloRange rng;
                            singlebj_expr .addTerm(this.Lj-this.BJsol[j][nnint],YJK[j][BJsol[j][nnint]]);
                            singlebj_expr .addTerm(-1,DeltaBJrm[j][nnint]);
                            rng = cplex.addGe(singlebj_expr, this.Lj-this.BJsol[j][nnint]) ;
                            rng.setName("def: Deltam_" + (j+1) +'^'+ (nnint+1));
                        }
                    }
                }
                else{
                    for(int j=0;j<this.nbStage-1;j++){
                        if(!this.Alter5_Delta[j][this.buffer[j]]){

                                    IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                                    IloRange rng;
                                    singlebj_expr.addTerm(1,DeltaJKp[j][this.buffer[j]]);
                                    singlebj_expr.addTerm(-1,BJ[j]);
                                    for(int l=1;l<=this.buffer[j];l++)
                                        singlebj_expr.addTerm(1,YJK[j][l]);
                                    rng = cplex.addEq(singlebj_expr, 0) ;
                                    rng.setName("def: Deltap_" + (j+1) +'^'+ (this.buffer[j]));


                                    singlebj_expr = cplex.linearNumExpr();
                                    singlebj_expr.addTerm(1,DeltaJKm[j][this.buffer[j]]);
                                    for(int l=1;l<=this.buffer[j];l++)
                                        singlebj_expr.addTerm(1,YJK[j][l]);
                                    rng = cplex.addEq(singlebj_expr, this.buffer[j]) ;
                                    rng.setName("def: Deltam_" + (j+1) +'^'+ (this.buffer[j]));

                            this.Alter5_Delta[j][this.buffer[j]]=true;
                        }
                    }



                }

                ///////////////////////////////// feasibility cut ///////////////////////////////////
                for(int j=0;j<this.NbStage -1;j++)                                                                              ////////
                {                                                                                                               ////////
                    for (int i=0;i<n;i++)                                                                                       ////////
                        {                                                                                                       ////////
                       DeltapPar[j][nnint]= DeltapPar[j][nnint] + (double)this.wbarij[i][j]*mm.Mijk[i][j][BJsol[j][nnint]]; //done    ////////
                       DeltamPar[j][nnint]= DeltamPar[j][nnint] + (double)this.wbarij[i][j]*mm.mijk[i][j][BJsol[j][nnint]]; //done     ////////
                    }                                                                                                           ////////
                }                                                                                                               ////////                                                                                  ///////////////////////////////
                                                                                                                                ////////
                                                                                                                                ////////
                double tijpar=0.0;   //done                                                                                           ////////
                for(int j=0;j<this.NbStage;j++)                                                                                 ////////
                {                                                                                                               ////////
                    for (int i=0;i<n;i++)                                                                                       ////////
                    {                                                                                                           ////////
                        tijpar = tijpar + tij[i][j]*(double)this.ubarij[i][j]; //done                                                  ////////
                    }                                                                                                           ////////
                }                                                                                                               ////////
                                                                                                                                ////////
                                                                                                                                ////////
                double thetapar = thetabar/this.THstar;//done //E: if theta is an array, then change this formulation                  ////////
                IloLinearNumExpr sumBJcut_expr = cplex.linearNumExpr();                                                         ////////
                IloRange rng;                                                                                                   ////////
                /////////////////////////////////MP-alter3/4 ////////////////////////////////////////
                if(AlterID==3 || AlterID==4){
                    for (int j=0; j<this.nbStage-1; j++)
                    {                                                                                                               ////////
                        sumBJcut_expr.addTerm(-DeltapPar[j][nnint],DeltaBJrp[j][nnint]);                                            ////////
                        sumBJcut_expr.addTerm(+DeltamPar[j][nnint],DeltaBJrm[j][nnint]);                                            ////////
                    }
                }
                else{
                    /////////////////////////////////MP-alter5 //////////////////////////////////////////
                    for (int j=0; j<this.nbStage-1; j++)                                                                            ////////
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
                    for (int k= (int) this.Buffer[j] +1; k<this.Uj; k++)
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
                String prova = program + "\\OUTPUT\\"+this.tempinstance +"model.lp";
                cplex.exportModel(prova);
                //save current solution
                for(int j=0;j<this.NbStage-1;j++)
                {
                    //MP-Alter3
                    this.BJsol[j][nnint+1] =  (int)  (cplex.getValue(this.BJ[j])+0.1);
                    Buffer[j] =  (int) (cplex.getValue(this.BJ[j])+0.1);
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
    }*/


}
