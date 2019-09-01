package test_mainfunc.optimization;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import test_mainfunc.simulation.SerialLine;

import static java.lang.Math.max;
import static java.lang.Math.min;

public class BendersIntModelAlter6 extends BendersIntModel {

    private IloNumVar[][] yjk;

    public BendersIntModelAlter6(SerialLine system, double THstar, int[] lB, int[] uB, int N, int W){
        super(system, THstar,  lB, uB, N,  W);
    }


    @Override
    public void solveBAPWithIntModel(double[][] tij)throws IloException {
        this.writer.println("Alter 6:");
        super.solveBAPWithIntModel(tij);
    }

    @Override
    public void getMmValue(double[][] tij){
        // initialize
        this.Mijk=new double[this.simulationLength+1][this.mySystem.nbStage][];
        this.mijk=new double[this.simulationLength+1][this.mySystem.nbStage][];
        for (int j = 0; j <= this.mySystem.nbStage-1; j++){
            for (int i = 0; i <= this.simulationLength; i++){
                this.Mijk[i][j]=new double[this.upperBoundj[j]+1];
                this.mijk[i][j]=new double[this.upperBoundj[j]+1];
                for(int k = 1; k <= this.upperBoundj[j]; k++){
                    this.mijk[i][j][k]=0;
                    this.Mijk[i][j][k]=100000;
                }
            }
        }

        // get value of small mijk
        for(int j = 1;j <= this.mySystem.nbStage-1; j++){
            for(int k = this.lowerBoundj[j]; k <= this.upperBoundj[j]; k++){
                for (int i = k+1; i <= this.simulationLength; i++) {
                    this.mijk[i][j][k]=tij[i-k][j+1];
                }
            }
        }

        //get value of big Mijk
        double[][] sij= new double [this.simulationLength+1][this.mySystem.nbStage+1];
        int[][] sumUj=new int[this.mySystem.nbStage+1][this.mySystem.nbStage+1];
        int[][] sumLj=new int[this.mySystem.nbStage+1][this.mySystem.nbStage+1];
        for(int j=1;j<this.mySystem.nbStage;j++){
            for(int j1=j+1;j1<=this.mySystem.nbStage;j1++){
                sumUj[j][j1]=0;
                sumLj[j][j1]=0;
                for(int l=j;l<=j1-1;l++){
                    sumUj[j][j1]+=this.upperBoundj[l];
                    sumLj[j][j1]+=this.lowerBoundj[l];
                }
            }
        }


        for(int i=1;i<=this.simulationLength;i++){
            for(int j=1;j<=this.mySystem.nbStage;j++){
                sij[i][j]=tij[i][j];
                if(j<this.mySystem.nbStage){
                    for(int j1=j+1;j1<=this.mySystem.nbStage;j1++){
                        for(int i1=max(1,i-j1+j-sumUj[j][j1]);i1<=i-j1+j-sumLj[j][j1];i1++){
                            sij[i][j]=max(sij[i][j],tij[i1][j1]);
                        }
                    }
                }
            }
        }

        for(int j=1;j <= this.mySystem.nbStage-1;j++){
            for(int k=this.lowerBoundj[j]; k <= this.upperBoundj[j];k++){
                for (int i = k+1 ; i <= this.simulationLength; i++){
                    this.Mijk[i][j][k]=sij[i-k][j+1];
                }
            }
        }

    }

    @Override
    public void generateFeasibilityCut(double[][] tij){
        double tijpar=0.0;
        for(int j=1;j<=this.mySystem.nbStage;j++)
        {
            for (int i=1;i<=this.simulationLength;i++)
            {
                if(this.mySystem.mySimulation.uij[i][j]>0)
                    tijpar += tij[i][j]*(double)this.mySystem.mySimulation.uij[i][j];
            }
        }
        this.newCut.constantTerm= - this.theta/this.THstar + tijpar;
    }

    @Override
    public void addFeasibilityCut(){
        try{
            double constantTerm=this.newCut.constantTerm;
            IloLinearNumExpr feacut = cplex.linearNumExpr();
            IloRange rng;
            double[][] mcoefficient=new double[mySystem.nbStage][];
            for(int j=1;j<=this.mySystem.nbStage-1;j++) {
                mcoefficient[j] = new double[this.upperBoundj[j] + 1];
                for (int k = this.lowerBoundj[j] + 1; k <= this.upperBoundj[j]; k++)
                    mcoefficient[j][k] = 0;
            }

            // Calculate the coefficient of each term in standard cuts
            for(int j=1;j<=this.mySystem.nbStage-1;j++){
                for(int i=1;i<=this.simulationLength;i++){
                    if(this.mySystem.mySimulation.wij[i][j]>0){
                        for(int k=this.lowerBoundj[j]+1;k<=this.mySystem.buffer[j];k++){
                            mcoefficient[j][k]+=this.mijk[i][j][k]*this.mySystem.mySimulation.wij[i][j];
                        }

                        for(int k=this.mySystem.buffer[j]+1;k<=this.upperBoundj[j];k++){
                            mcoefficient[j][k]+=this.Mijk[i][j][k]*this.mySystem.mySimulation.wij[i][j];
                        }
                    }
                }

                for(int k=this.lowerBoundj[j]+1;k<=this.mySystem.buffer[j];k++){
                    constantTerm+=mcoefficient[j][k];
                }

            }

            // Is combinatorial cut or standard cut tighter?
            boolean combCutIsTighter=true;
            for(int j=1;j<=this.mySystem.nbStage-1;j++){
                double lhs=0;
                for(int k=this.lowerBoundj[j]+1;k<=min(this.mySystem.buffer[j]+1,this.upperBoundj[j]);k++){
                    lhs+=mcoefficient[j][k];
                }
                if(lhs<constantTerm)
                    combCutIsTighter=false;
            }

            // Add the correct cut
            if(combCutIsTighter){
                for(int j=1;j<=this.mySystem.nbStage-1;j++){
                    if(this.mySystem.buffer[j]<this.upperBoundj[j])
                        feacut.addTerm(1,this.yjk[j][this.mySystem.buffer[j]+1]);
                }
                rng=cplex.addGe(feacut,1.0);
                rng.setName("combcut of iter "+ (numit)+": ");
            }
            else{
                for(int j=1;j<=this.mySystem.nbStage-1;j++){
                    double cumulateMjk=0;
                    for(int k=this.lowerBoundj[j]+1;k<=this.upperBoundj[j];k++){
                        cumulateMjk+=mcoefficient[j][k];
                        feacut.addTerm(-mcoefficient[j][k],this.yjk[j][k]);
                        if(cumulateMjk > constantTerm)
                            break;
                    }
                }
                feacut.setConstant(constantTerm);
                rng = cplex.addLe(feacut,0);
                rng.setName("feascut of iter "+ (numit)+": ");
            }

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

            /*for (int j=1; j<=this.mySystem.nbStage-1; j++){
                IloLinearNumExpr sumyjk_expr = cplex.linearNumExpr();
                for (int k=1; k<= this.upperBoundj[j]; k++)
                    sumyjk_expr.addTerm(k,yjk[j][k]);
                sumyjk_expr.addTerm(-1,bj[j]);
                this.InitialLeConstraints.add(sumyjk_expr);
            }*/

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
