package test_mainfunc.optimization;

import ilog.concert.*;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.util.Stopwatch;

import java.io.*;
import java.util.HashMap;

public class BendersStolletz extends BendersBAP {
    private IloNumVar[][] zjk;

    // Master problem is initialized inside the constructor.
    public BendersStolletz(SerialLine system, double THstar, int[] lB, int[] uB, int N, int W){
        super(system,THstar, lB, uB,N,W);
    }


    public BendersStolletz(SerialLine system, double THstar, int[] lB, int[] uB, int N, int W, int j1, int j2,
                           HashMap<String,Integer> initialBounds){
        super(system, THstar,  lB, uB, N, W, j1, j2, initialBounds);
    }


    public void solveBAPWithStolletz(double[][] tij,int offset) throws IloException{
        totalTimeMeasure.start();
        this.mySystem.mySimulation =this.mySystem.new SimulationBAS(this.simulationLength,this.warmupLength,tij);
        for(int j=1;j<=mySystem.nbStage-1;j++)
            this.mySystem.buffer[j]=lowerBoundj[j];

        this.mySystem.mySimulation.simBAS(false,offset);
        this.saveIterationSolution();

        solvability=true;
        while((this.THstar-this.mySystem.TH > 0)&&(numit < this.MAX_ITE)&&solvability){
            numit++;
            this.addFeasibilityCut();
            totalTimeMeasure.pause();
            double timeLimit = MAX_TIME-totalTimeMeasure.elapseTimeSeconds;
            totalTimeMeasure.restart();
            this.solveMasterProb(timeLimit);
            this.updateBufferSpaces();
            this.saveIterationSolution();
            this.mySystem.mySimulation.simBAS(false,offset);
        }
        this.endMasterProb();

    }

    public void solveWithDecomposition(double[][] tij) throws IloException {
        for(int j=1;j<=mySystem.nbStage-1;j++){
            mySystem.buffer[j] = lowerBoundj[j];
        }

        HashMap<String, Integer> initialBounds = new HashMap<>();
        numit=0;
        double remainTime = MAX_TIME;
        boolean stopLoop = false;
        for(int m=2;m<=mySystem.nbStage && !stopLoop;m++){
            for(int startMachine = 1; startMachine<=mySystem.nbStage-1 && !stopLoop;startMachine++){
                int endMachine = startMachine + m-1;
                if(endMachine<=mySystem.nbStage){
                    Stopwatch sw = new Stopwatch();
                    sw.start();
                    BendersStolletz subSystemBap = new BendersStolletz(mySystem, THstar, lowerBoundj, upperBoundj, simulationLength, warmupLength,
                            startMachine, endMachine, initialBounds);
                    subSystemBap.writer = new PrintWriter(OutputStream.nullOutputStream());
                    subSystemBap.setMaxTime(remainTime);
                    subSystemBap.solveBAPWithStolletz(tij,startMachine-1);
                    int solSubProb = 0;
                    for(int j=1;j<=subSystemBap.mySystem.nbStage-1;j++){
                        solSubProb+=subSystemBap.mySystem.buffer[j];
                    }
                    addInitialBounds(initialBounds,startMachine,endMachine,solSubProb);
                    numit+=subSystemBap.numit;
                    sw.stop();
                    remainTime -= sw.elapseTimeSeconds;

                    if(remainTime<0){
                        if(m<mySystem.nbStage)
                            subSystemBap = new BendersStolletz(mySystem, THstar, lowerBoundj, upperBoundj, simulationLength, warmupLength,
                                1, mySystem.nbStage, initialBounds);
                        subSystemBap.setMaxIte(1);
                        subSystemBap.setMaxTime(MAX_TIME);
                        subSystemBap.solveBAPWithStolletz(tij,0);
                        if (mySystem.nbStage - 1 >= 0)
                            System.arraycopy(subSystemBap.mySystem.buffer, 1, mySystem.buffer, 1, mySystem.nbStage - 1);
                        stopLoop = true;
                    }
                    else if(m==mySystem.nbStage){
                        if (mySystem.nbStage - 1 >= 0)
                            System.arraycopy(subSystemBap.mySystem.buffer, 1, mySystem.buffer, 1, mySystem.nbStage - 1);
                        stopLoop = true;
                    }
                }
            }
        }
    }

    @Override
    public void setupVariables() throws IloException{
        super.setupVariables();
        this.zjk=new IloNumVar[this.mySystem.nbStage][];
        for(int j=0; j <= this.mySystem.nbStage-1;j++){
            this.zjk[j]=new IloNumVar[this.upperBoundj[j]+1];
            for(int k=0;k<=this.upperBoundj[j];k++){
                String label = "z_" + (j)+"_"+(k);
                this.zjk[j][k]=cplex.boolVar(label);
            }
        }
    }

    @Override
    public void initialConstraints(){
        try{
            for (int j=1; j <= this.mySystem.nbStage-1; j++)
            {
                IloLinearNumExpr sumzjk_expr = cplex.linearNumExpr();
                for (int k=0; k<=this.upperBoundj[j]; k++)
                    sumzjk_expr.addTerm(1,zjk[j][k]);
                sumzjk_expr.setConstant(-1.0);
                this.InitialEqConstraints.add(sumzjk_expr);
            }

            for (int j=1; j<=this.mySystem.nbStage-1; j++)
            {
                IloLinearNumExpr bj_expr = cplex.linearNumExpr();
                bj_expr.addTerm(1,bj[j]);
                for (int k=0; k<=this.upperBoundj[j]; k++)
                {
                    bj_expr.addTerm(-k,zjk[j][k]);
                }
                this.InitialEqConstraints.add(bj_expr);
            }
        }catch(Exception exc){exc.printStackTrace();}

    }

    private void addFeasibilityCut() throws IloException{

        IloLinearNumExpr sumBJcut_expr = cplex.linearNumExpr();
        IloRange rng;
        for (int j=1; j<=this.mySystem.nbStage-1; j++)
        {
            for (int k = this.mySystem.buffer[j]+1; k<= this.upperBoundj[j]; k++)
            {
                sumBJcut_expr.addTerm(1,this.zjk[j][k]);
            }
        }
        rng = cplex.addGe(sumBJcut_expr,1.0);
        rng.setName("combcut of iter"+numit);

    }

}
