package test_mainfunc.optimization;

import ilog.concert.*;
import test_mainfunc.simulation.SerialLine;

public class BendersStolletz extends BendersBAP {
    private IloNumVar[][] zjk;

    // Master problem is initialized inside the constructor.
    public BendersStolletz(SerialLine system, double THstar, int[] lB, int[] uB, int N, int W){
        super(system,THstar, lB, uB,N,W);
    }

    public void solveBAPWithStolletz(double[][] tij) throws IloException{

        this.mySystem.mySimulation =this.mySystem.new SimulationBAS(this.simulationLength,this.warmupLength,tij);
        for(int j=1;j<=mySystem.nbStage-1;j++)
            this.mySystem.buffer[j]=lowerBoundj[j];

        this.mySystem.mySimulation.simBAS(false);
        this.saveIterationSolution();

        solvability=true;
        while((this.THstar-this.mySystem.TH > 0.0001)&&(numit < this.MAX_ITE)&&solvability){
            numit++;

            this.addFeasibilityCut();
            this.solveMasterProb();
            this.updateBufferSpaces();
            this.saveIterationSolution();
            this.mySystem.mySimulation.simBAS(false);
        }
        this.endMasterProb();

    }

    @Override
    public void setupVariables() throws IloException{
        super.setupVariables();
        this.zjk=new IloNumVar[this.mySystem.nbStage][];
        for(int j=0; j <= this.mySystem.nbStage-1;j++){
            this.zjk[j]=new IloNumVar[this.upperBoundj[j]+1];
            for(int k=0;k<=this.upperBoundj[j];k++){
                this.zjk[j][k]=cplex.boolVar();
            }
        }
    }

    @Override
    public void initialConstraints(){
        try{
            for (int j=1; j <= this.mySystem.nbStage-1; j++)
            {
                IloLinearNumExpr sumzjk_expr = cplex.linearNumExpr();
                for (int k=1; k<=this.upperBoundj[j]; k++)
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
