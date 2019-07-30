package test_mainfunc.optimization;

import ilog.concert.*;
import ilog.cplex.*;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.util.Stopwatch;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;


public abstract class BendersBAP {

    //Input
    public double THstar;
    SerialLine mySystem;
    int upperBoundj[];
    int lowerBoundj[];
    int simulationLength;
    int warmupLength;
    final int MAX_ITE=10000;
    final double MAX_CPLEX_TIME=86400;
    // End of input

    //Output
    public int numit; // Iteration number
    public Stopwatch cplexTimeMeasure; // Total cplex time
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
        this.mySystem=new SerialLine();
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
    public void solveMasterProb() throws IloException{
        try{
            double timeLimit=MAX_CPLEX_TIME-cplexTimeMeasure.elapseTimeSeconds;
            cplex.setParam(IloCplex.DoubleParam.TiLim,timeLimit);



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

