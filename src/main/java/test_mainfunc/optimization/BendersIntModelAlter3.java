package test_mainfunc.optimization;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import test_mainfunc.simulation.SerialLine;

public class BendersIntModelAlter3 extends BendersIntModel {

    public BendersIntModelAlter3(SerialLine system, double THstar, int[] lB, int[] uB, int N){
        super(system, THstar,  lB, uB, N);
    }

    /*@Override
    public void solveBAPWithIntModel(double[][] tij,int offset, boolean resetMm)throws IloException {
        this.writer.println("Alter 3:");
        super.solveBAPWithIntModel(tij, offset,resetMm);
    }*/

    @Override
    public void addFeasibilityCut(){
        String label;
        try{
            IloNumVar[] newdeltaPbjr=new IloNumVar[this.mySystem.nbStage];
            IloNumVar[] newdeltaMbjr=new IloNumVar[this.mySystem.nbStage];
            IloNumVar[] newyjr=new IloNumVar[this.mySystem.nbStage];

            for(int j=0;j<=this.mySystem.nbStage-1;j++){
                label = "Deltap_" + (j) +'^'+ (numit);
                newdeltaPbjr[j]=cplex.intVar(0,this.upperBoundj[j]-this.mySystem.buffer[j],label);
                label = "Deltam_" + (j) +'^'+ (numit);
                newdeltaMbjr[j]=cplex.intVar(0,this.mySystem.buffer[j]-this.lowerBoundj[j],label);
                label = "Yrj_" + (j) +'^'+ (numit);
                newyjr[j] = cplex.boolVar(label);
            }

            for (int j=1; j <= this.mySystem.nbStage-1; j++)
            {
                IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                IloRange rng;
                singlebj_expr .addTerm(1,bj[j]);
                singlebj_expr .addTerm(-1,newdeltaPbjr[j]);
                singlebj_expr .addTerm(+1,newdeltaMbjr[j]);
                rng = cplex.addEq(singlebj_expr, this.mySystem.buffer[j]) ;
                rng.setName("def: bj_" + (j) +'^'+ (this.numit));
            }

            for (int j=1; j<=this.mySystem.nbStage-1; j++)
            {
                IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                IloRange rng;
                singlebj_expr .addTerm(this.upperBoundj[j]-this.mySystem.buffer[j],newyjr[j]);
                singlebj_expr .addTerm(-1,newdeltaPbjr[j]);
                rng = cplex.addGe(singlebj_expr, 0) ;
                rng.setName("def: Deltap_" + (j) +'^'+ (numit));
            }

            for (int j=1; j<=this.mySystem.nbStage-1; j++)
            {
                IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                IloRange rng;
                singlebj_expr .addTerm(this.lowerBoundj[j]-this.mySystem.buffer[j],newyjr[j]);
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
}
