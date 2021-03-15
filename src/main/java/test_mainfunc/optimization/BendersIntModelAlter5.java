package test_mainfunc.optimization;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import test_mainfunc.simulation.SerialLine;

import static java.lang.Math.max;

public class BendersIntModelAlter5 extends BendersIntModel {

    private IloNumVar[][] deltaPjk;
    private IloNumVar[][] deltaMjk;
    private IloNumVar[][] yjk;
    private boolean[][] delta;

    public BendersIntModelAlter5(SerialLine system, double THstar, int[] lB, int[] uB, int N){
        super(system, THstar,  lB, uB, N);
        this.delta=new boolean [this.mySystem.nbStage][];
        for(int j=0;j<=this.mySystem.nbStage-1;j++){
            this.delta[j]=new boolean[this.upperBoundj[j]+1];
            for(int k=0;k<=this.upperBoundj[j];k++)
                delta[j][k]=false;
        }
    }


    /*@Override
    public void solveBAPWithIntModel(double[][] tij,int offset,boolean resetMm)throws IloException {
        this.writer.println("Alter 5:");
        super.solveBAPWithIntModel(tij,offset,resetMm);
    }*/

    @Override
    public void addFeasibilityCut(){
        try{

            for(int j=1;j<=this.mySystem.nbStage-1;j++){
                if(!this.delta[j][this.mySystem.buffer[j]]){

                    IloLinearNumExpr singlebj_expr = cplex.linearNumExpr();
                    IloRange rng;
                    singlebj_expr.addTerm(1,this.deltaPjk[j][this.mySystem.buffer[j]]);
                    singlebj_expr.addTerm(-1,this.bj[j]);
                    for(int l=1;l<=this.mySystem.buffer[j];l++)
                        singlebj_expr.addTerm(1,this.yjk[j][l]);
                    rng = cplex.addEq(singlebj_expr, 0) ;
                    rng.setName("def: Deltap_" + (j) +'^'+ (this.mySystem.buffer[j]));


                    singlebj_expr = cplex.linearNumExpr();
                    singlebj_expr.addTerm(1,this.deltaMjk[j][this.mySystem.buffer[j]]);
                    for(int l=1;l<=this.mySystem.buffer[j];l++)
                        singlebj_expr.addTerm(1,this.yjk[j][l]);
                    rng = cplex.addEq(singlebj_expr, this.mySystem.buffer[j]) ;
                    rng.setName("def: Deltam_" + (j) +'^'+ (this.mySystem.buffer[j]));

                    this.delta[j][this.mySystem.buffer[j]]=true;
                }
            }

            //feasibility cut
            IloLinearNumExpr sumBJcut_expr = cplex.linearNumExpr();
            sumBJcut_expr.setConstant(this.newCut.constantTerm);
            IloRange rng;
            for (int j=1; j<=this.mySystem.nbStage-1; j++)
            {
                sumBJcut_expr.addTerm(+this.newCut.coefdeltaP[j],this.deltaPjk[j][mySystem.buffer[j]]);
                sumBJcut_expr.addTerm(+this.newCut.coefdeltaM[j],this.deltaMjk[j][mySystem.buffer[j]]);
            }
            rng = cplex.addLe(sumBJcut_expr,0);
            rng.setName("feascut of iter: "+ (numit));


        }catch(Exception exc){exc.printStackTrace();}



    }

    @Override
    public void setupVariables() throws IloException {
        super.setupVariables();
        try{
            this.yjk=new IloNumVar[this.mySystem.nbStage][];
            for(int j=0;j<=this.mySystem.nbStage-1;j++){

                this.yjk[j]=new IloNumVar[this.upperBoundj[j]+1];
                for(int k=0;k<=this.upperBoundj[j];k++){
                    String label="y_"+j+"_"+k;
                    this.yjk[j][k]=cplex.boolVar(label);
                }

            }

            this.deltaPjk = new IloNumVar[this.mySystem.nbStage][];
            this.deltaMjk = new IloNumVar[this.mySystem.nbStage][];
            String label;
            for (int j=1; j<=this.mySystem.nbStage-1; j++)
            {
                this.deltaPjk[j] = new IloNumVar[this.upperBoundj[j]+1];
                this.deltaMjk[j] = new IloNumVar[this.upperBoundj[j]+1];

                for(int k=0;k<=this.upperBoundj[j];k++){
                    label = "Deltap_" + (j) +'^'+ (k);
                    deltaPjk[j][k] = cplex.intVar(0, this.upperBoundj[j]-k,label);
                    label = "Deltam_" + (j) +'^'+ (k);
                    deltaMjk[j][k] = cplex.intVar(0, max(0,k-this.lowerBoundj[j]),label);
                }
            }
        }catch(Exception exc){exc.printStackTrace();}
    }

    @Override
    public void initialConstraints(){
        try{
            for (int j=1; j<=this.mySystem.nbStage-1; j++){
                IloLinearNumExpr sumyjk_expr = cplex.linearNumExpr();
                for (int k=1; k<= this.upperBoundj[j]; k++)
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

    @Override
    void endMasterProb()throws IloException{
        this.writer.println(this.cplex.getModel());
        super.endMasterProb();
    }





}
