package test_mainfunc.optimization;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import test_mainfunc.simulation.SerialLine;

import static java.lang.Math.min;

public class BendersIntModelAlter4 extends BendersIntModel{

    private IloNumVar[][] yjk;

    public BendersIntModelAlter4(SerialLine system, double THstar, int[] lB, int[] uB, int N, int W){
        super(system, THstar,  lB, uB, N,  W);
    }

    @Override
    public void solveBAPWithIntModel(double[][] tij)throws IloException {
        this.writer.println("Alter 4:");
        super.solveBAPWithIntModel(tij);
    }

    @Override
    public void addFeasibilityCut(){
        String label;
        try{
            IloNumVar[] newdeltaPbjr=new IloNumVar[this.mySystem.nbStage];
            IloNumVar[] newdeltaMbjr=new IloNumVar[this.mySystem.nbStage];

            for(int j=0;j<=this.mySystem.nbStage-1;j++){
                label = "Deltap_" + (j) +'^'+ (numit);
                newdeltaPbjr[j]=cplex.intVar(0,this.upperBoundj[j]-this.mySystem.buffer[j],label);
                label = "Deltam_" + (j) +'^'+ (numit);
                newdeltaMbjr[j]=cplex.intVar(0,this.mySystem.buffer[j]-this.lowerBoundj[j],label);
            }

            for (int j=1; j<=this.mySystem.nbStage-1; j++)
            {
                IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                IloRange rng;
                singlebj_expr .addTerm(1,bj[j]);
                singlebj_expr .addTerm(-1,newdeltaPbjr[j]);
                singlebj_expr .addTerm(+1,newdeltaMbjr[j]);
                rng = cplex.addEq(singlebj_expr, this.mySystem.buffer[j]) ;
                rng.setName("def: bj_" + (j) +'^'+ (this.numit));
            }

            for (int j=1; j<= this.mySystem.nbStage-1; j++)
            {
                IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                IloRange rng;
                singlebj_expr .addTerm(this.upperBoundj[j]-this.mySystem.buffer[j],this.yjk[j][min(this.mySystem.buffer[j]+1,this.upperBoundj[j])]);
                singlebj_expr .addTerm(-1,newdeltaPbjr[j]);
                rng = cplex.addGe(singlebj_expr, 0) ;
                rng.setName("def: Deltap_" + (j) +'^'+ (numit));
            }

            for (int j=1; j<=this.mySystem.nbStage-1; j++)
            {
                IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                IloRange rng;
                singlebj_expr .addTerm(this.lowerBoundj[j]-this.mySystem.buffer[j],this.yjk[j][this.mySystem.buffer[j]]);
                singlebj_expr .addTerm(-1,newdeltaMbjr[j]);
                rng = cplex.addGe(singlebj_expr, this.lowerBoundj[j]-this.mySystem.buffer[j]) ;
                rng.setName("def: Deltam_" + (j) +'^'+ (numit));
            }

            //feasibility cut
            IloLinearNumExpr sumBJcut_expr = cplex.linearNumExpr();
            sumBJcut_expr.setConstant(this.newCut.constantTerm);
            IloRange rng;
            for (int j=1; j<=this.mySystem.nbStage-1; j++)
            {
                sumBJcut_expr.addTerm(+this.newCut.coefdeltaP[j],newdeltaPbjr[j]);
                sumBJcut_expr.addTerm(+this.newCut.coefdeltaM[j],newdeltaMbjr[j]);
            }
            rng = cplex.addLe(sumBJcut_expr,0);
            rng.setName("feascut of iter: "+ (numit));


        }catch(Exception exc){exc.printStackTrace();}



    }

    @Override
    public void setupVariables(){
        try{
            this.yjk=new IloNumVar[this.mySystem.nbStage][];
            for(int j=0;j<=this.mySystem.nbStage;j++){
                this.yjk[j]=new IloNumVar[this.upperBoundj[j]+1];
                for(int k=0;k<=this.upperBoundj[j];k++)
                    this.yjk[j][k]=cplex.boolVar();
            }
        }catch(Exception exc){exc.printStackTrace();}
    }

    @Override
    public void initialConstraints(){
        try{
            for (int j=1; j<=this.mySystem.nbStage-1; j++){
                IloLinearNumExpr sumyjk_expr = cplex.linearNumExpr();
                for (int k=1; k<=this.upperBoundj[j]; k++)
                    sumyjk_expr.addTerm(1,yjk[j][k]);
                sumyjk_expr.addTerm(-1,bj[j]);
                this.InitialEqConstraints.add(sumyjk_expr);
            }

            for (int j=1; j<=this.mySystem.nbStage-1; j++){
                for (int k=1; k<=this.upperBoundj[j]; k++){
                    IloLinearNumExpr yjk_expr = cplex.linearNumExpr();
                    yjk_expr.addTerm(1,yjk[j][k]);
                    yjk_expr.addTerm(-1,yjk[j][k-1]);
                    this.InitialLeConstraints.add(yjk_expr);
                }
            }
            for (int j=1; j<=this.mySystem.nbStage-1; j++)
            {
                for(int k=1;k<=this.lowerBoundj[j];k++){
                    IloLinearNumExpr yjk_expr = cplex.linearNumExpr();
                    yjk_expr.addTerm(1,yjk[j][k]);
                    yjk_expr.setConstant(-1);
                    this.InitialEqConstraints.add(yjk_expr);
                }
            }

        }catch(Exception exc){exc.printStackTrace();}
    }


}
