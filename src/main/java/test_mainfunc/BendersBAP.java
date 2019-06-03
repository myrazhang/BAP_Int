package test_mainfunc;

import ilog.concert.*;
import ilog.cplex.*;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;


public class BendersBAP {

    public SerialLine mySystem;
    public double THstar;
    public PrintWriter writer;
    public int upperBoundj[];
    public int lowerBoundj[];
    public List<IterationSolution> bjsol;
    public int simulationLength;
    public int warmupLength;
    public int numit;
    public StopWatch cplexTimeMeasure;
    public final int MAX_ITE=10000;


    public IloNumVar[] bj; //buffer spaces
    public IloCplex cplex;
    public IloLinearNumExpr objective;
    public List<IloLinearNumExpr> InitialEqConstraints;
    public List<IloLinearNumExpr> InitialLeConstraints;

    public boolean solvability;

    public BendersBAP(SerialLine system,double THstar,int[] lB, int[] uB, int N, int W){
        String label;
        this.simulationLength=N;
        this.warmupLength=W;
        this.mySystem=new SerialLine();
        this.mySystem=system;
        this.THstar=THstar;
        this.cplexTimeMeasure=new StopWatch();

        this.lowerBoundj=new int[this.mySystem.nbStage-1];
        this.upperBoundj=new int[this.mySystem.nbStage-1];
        for(int j=0;j<this.mySystem.nbStage-1;j++){
            this.lowerBoundj[j]=lB[j];
            this.upperBoundj[j]=uB[j];
        }

        try {
            this.cplex=new IloCplex();
            this.bj=new IloNumVar[this.mySystem.nbStage-1];
            for (int j = 0; j < this.mySystem.nbStage - 1; j++) {
                label = "b_" + (j + 1);
                this.bj[j] = cplex.intVar(this.lowerBoundj[j], this.upperBoundj[j], label);
            }
            InitialEqConstraints=new ArrayList<>();
            InitialLeConstraints=new ArrayList<>();

            this.setupVariables();
            this.initialConstraints();
            this.initializeMasterProb();

        }catch(Exception exc){exc.printStackTrace();}
    }
    public void setupVariables(){}
    public void initialConstraints(){}
    public void initializeMasterProb(){
        try{
            this.objective = this.cplex.linearNumExpr();
            for (int j=0; j<this.mySystem.nbStage-1; j++)
                this.objective.addTerm(1,bj[j]);
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

    public void solveMasterProb(){
        try{
            this.cplexTimeMeasure.restart();
            this.solvability=this.cplex.solve();
            this.cplexTimeMeasure.pause();

            if(this.solvability){

                int totcap = 0;
                for(int j=0;j<this.mySystem.nbStage-1;j++)
                {
                    totcap += this.cplex.getValue(this.bj[j]);
                }
                this.writer.println("it " + numit + " OF: "+ totcap);
            }
            else
                System.out.println("No solution was found by Cplex!");

        }catch (Exception exc) {exc.printStackTrace();}
    }

    public void endMasterProb(){
        this.cplex.end();
    }

    public class IterationSolution{
        public int[] solution;
    }

    public void saveIterationSolution(){
        IterationSolution sol=new IterationSolution();
        sol.solution=new int[mySystem.nbStage-1];
        for(int j=0;j<mySystem.nbStage-1;j++)
            sol.solution[j]=this.mySystem.buffer[j];
        this.bjsol.add(sol);
    }

    public void updateBufferSpaces(){
        try {
            if (this.solvability) {
                for (int j = 0; j < this.mySystem.nbStage - 1; j++)
                    this.mySystem.buffer[j] = (int) (this.cplex.getValue(this.bj[j]) + 0.1);
            }
        }catch(Exception exc){exc.printStackTrace();}
    }
}
