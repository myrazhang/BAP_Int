package test_mainfunc.optimization;

import ilog.concert.*;
import ilog.cplex.*;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.util.Stopwatch;

import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;


public abstract class BendersBAP {

    //Input
    public double THstar;
    public SerialLine mySystem;
    int upperBoundj[];
    int lowerBoundj[];
    int simulationLength;
    int warmupLength;
    int MAX_ITE=20000;
    double MAX_TIME=18000;
    double TOLERANCE_RATIO = 0.000;
    // End of input

    //Output
    public int numit; // Iteration number
    public Stopwatch cplexTimeMeasure; // Total cplex time
    public Stopwatch totalTimeMeasure; // Total time
    public PrintWriter writer; // Output file: this file
    private List<IterationSolution> bjsol;// Master problem solution of each iteration
    boolean solvability; // whether the problem is feasible or not.
    // End of output

    // Cplex model elements
    public IloCplex cplex;
    IloNumVar[] bj; //decision variable: buffer spaces
    List<IloLinearNumExpr> InitialEqConstraints; // Initial constraints of "="
    List<IloLinearNumExpr> InitialLeConstraints; // Initial constraints of "<="
    IloLinearNumExpr objective;
    // End of Cplex model elements


    // Constructor
    BendersBAP(SerialLine system, double THstar,int[] lB, int[] uB, int N, int W){
        // Input
        this.mySystem=system;
        this.THstar=THstar;
        this.lowerBoundj=new int[this.mySystem.nbStage];
        this.upperBoundj=new int[this.mySystem.nbStage];
        for(int j=1 ; j<= this.mySystem.nbStage-1;j++){
            this.lowerBoundj[j]=lB[j];
            this.upperBoundj[j]=uB[j];
        }
        this.simulationLength=N;
        this.warmupLength=W;

        // Output
        this.bjsol=new ArrayList<>();
        this.numit=0;
        this.cplexTimeMeasure=new Stopwatch();
        this.totalTimeMeasure=new Stopwatch();
        this.writer = new PrintWriter(OutputStream.nullOutputStream());

        // Cplex
        try {

            this.cplex=new IloCplex();

            InitialEqConstraints=new ArrayList<>();
            InitialLeConstraints=new ArrayList<>();

            this.setupVariables();
            this.initialConstraints();
            this.initializeMasterProb();

        }catch(Exception exc){exc.printStackTrace();}
    }

    BendersBAP(SerialLine system, double THstar, int[] lB, int[] uB, int N, int W,
               // **************************************************************
               int j1, // index of first machine
               int j2, // index of last machine
               HashMap<String,Integer> initialBounds
               // ***************************************************************
    ){
        // Input
        this.mySystem=system.subsystem(j1,j2);
        this.THstar=THstar;
        this.lowerBoundj=new int[this.mySystem.nbStage];
        this.upperBoundj=new int[this.mySystem.nbStage];
        for(int j=1 ; j<= this.mySystem.nbStage-1;j++){
            this.lowerBoundj[j]=lB[j+j1-1];
            this.upperBoundj[j]=uB[j+j1-1];
        }
        this.simulationLength=N;
        this.warmupLength=W;

        // Output
        this.bjsol=new ArrayList<>();
        this.numit=0;
        this.cplexTimeMeasure=new Stopwatch();
        this.totalTimeMeasure=new Stopwatch();
        this.writer = new PrintWriter(OutputStream.nullOutputStream());


        // Cplex
        try {

            this.cplex=new IloCplex();

            InitialEqConstraints=new ArrayList<>();
            InitialLeConstraints=new ArrayList<>();

            this.setupVariables();
            this.initialConstraints();
            this.initializeMasterProb();

            if(initialBounds != null){
                for(int j=j1;j<=j2-1;j++){
                    for(int j0=j+1;j0<=j2;j0++){
                        if(initialBounds.containsKey(j+"_"+j0)){
                            IloLinearNumExpr constBound = cplex.linearNumExpr();
                            for(int k=j-j1+1;k<=j0-j1;k++){
                                constBound.addTerm(1,bj[k]);
                            }
                            cplex.addGe(constBound,initialBounds.get(j+"_"+j0),"bound_"+j+"_"+j0);
                        }
                    }
                }
            }

        }catch(Exception exc){exc.printStackTrace();}
    }

    public void addInitialBounds(HashMap<String,Integer> initialBounds,int j1,int j2,int bound){
        initialBounds.put(j1+"_"+j2,bound);
    }


    // Public methods
    public void setupVariables() throws IloException{
        String label;
        this.bj=new IloNumVar[this.mySystem.nbStage];
        for (int j = 0; j <= this.mySystem.nbStage - 1; j++) {
            label = "b_" + (j);
            this.bj[j] = cplex.intVar(this.lowerBoundj[j], this.upperBoundj[j], label);
        }
    }
    public void initialConstraints(){}
    public void solveMasterProb(double timeLimit) throws IloException{
        try{
            //double timeLimit=MAX_CPLEX_TIME-cplexTimeMeasure.elapseTimeSeconds;
            if(timeLimit<=0) {
                this.solvability = false;
                System.out.println("Time is used up!");
                return;
            }
            cplex.setParam(IloCplex.DoubleParam.TiLim,timeLimit);
            cplex.setParam(IloCplex.Param.Threads, 16);

            this.cplexTimeMeasure.restart();
            this.solvability=this.cplex.solve();
            this.cplexTimeMeasure.pause();

            if(!this.solvability)
                System.out.println("No solution was found by Cplex!");
            else if(cplex.getCplexStatus()== IloCplex.CplexStatus.AbortTimeLim){
                this.solvability=false;
            }

        }catch (Exception exc) {exc.printStackTrace();}
    }
    public class IterationSolution{
        public int[] solution;

        public IterationSolution(){
            this.solution=new int [mySystem.nbStage];
        }
    }
    public void setMaxIte(int i){MAX_ITE = i;}
    public void setMaxTime(double t){MAX_TIME=t;}

    // Private & package-private methods
    void endMasterProb() throws IloException{
        this.cplex.end();
    }
    void saveIterationSolution(){
        IterationSolution sol=new IterationSolution();
        for(int j=1;j<=mySystem.nbStage-1;j++)
            sol.solution[j]=this.mySystem.buffer[j];
        this.bjsol.add(sol);
    }
    void updateBufferSpaces(){
        try {
            if (this.solvability) {
                for (int j = 1; j <= this.mySystem.nbStage - 1; j++)
                    this.mySystem.buffer[j] = (int) (this.cplex.getValue(this.bj[j]) + 0.1);
            }
        }catch(Exception exc){exc.printStackTrace();}
    }
    private void initializeMasterProb(){
        try{
            this.objective = this.cplex.linearNumExpr();
            for (int j=1; j<= this.mySystem.nbStage-1; j++)
                this.objective.addTerm(1,this.bj[j]);
            this.cplex.addMinimize(this.objective);

            IloRange rng;
            for(IloLinearNumExpr i: InitialEqConstraints){
                rng = this.cplex.addEq(i,0);
                rng.setName("Initial constraint");
            }
            for(IloLinearNumExpr i: InitialLeConstraints){
                rng = this.cplex.addLe(i,0);
                rng.setName("Initial constraint");
            }
        }catch(Exception exc) {exc.printStackTrace();}
    }
}

