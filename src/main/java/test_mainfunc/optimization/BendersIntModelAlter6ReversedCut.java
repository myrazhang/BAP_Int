package test_mainfunc.optimization;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloRange;
import test_mainfunc.simulation.SerialLine;

import static java.lang.Math.*;

public class BendersIntModelAlter6ReversedCut extends BendersIntModelAlter6{

    private BendersIntModelAlter6 myReversedBAPModel;
    private double constantC =1;

    public BendersIntModelAlter6ReversedCut(SerialLine system, double THstar, int[] lB, int[] uB, int N, int W) {
        super(system, THstar, lB, uB, N, W);
        int[] reversedLB = new int[lB.length];
        int[] reversedUB = new int[uB.length];
        for (int j=0;j<lB.length;j++){
            reversedLB[j] = lB[lB.length-1-j];
            reversedUB[j] = uB[uB.length-1-j];
        }
        this.myReversedBAPModel = new BendersIntModelAlter6(system, THstar, reversedLB, reversedUB, N, W);
    }

    @Override
    public void solveBAPWithIntModel(double[][] tij)throws IloException {
        this.writer.println("Alter 6 with reversed cuts:");

        double[][] reversedTij = new double[this.simulationLength+1][this.mySystem.nbStage+1];
        for (int i=1;i<=this.simulationLength;i++){
            for (int j=1;j<=this.mySystem.nbStage;j++){
                reversedTij[i][j] = tij[this.simulationLength+1-i][this.mySystem.nbStage+1-j];
            }
        }
        myReversedBAPModel.getMmValue(reversedTij);

        super.solveBAPWithIntModel(tij);
    }

    @Override
    public void addFeasibilityCut(){
        try{

            // Calculate the coefficient of each term in the reversed cut
            double reversed_constantTerm=this.newCut.constantTerm;
            double[][] reversed_mcoefficient=new double[mySystem.nbStage][];
            for(int j=1;j<=this.mySystem.nbStage-1;j++) {
                reversed_mcoefficient[j] = new double[this.upperBoundj[j] + 1];
                for (int k = this.lowerBoundj[j] + 1; k <= this.upperBoundj[j]; k++)
                    reversed_mcoefficient[j][k] = 0;
            }

            for(int j=1;j<=this.mySystem.nbStage-1;j++){
                for(int i=1;i<=this.simulationLength;i++){
                    if(this.mySystem.mySimulation.wij[i][j]>0){
                        for(int k=this.lowerBoundj[j]+1;k<=this.mySystem.buffer[j];k++){
                            reversed_mcoefficient[j][k]+=this.myReversedBAPModel.mijk[this.simulationLength+this.mySystem.buffer[j]+1-i][this.mySystem.nbStage-j][k] * this.mySystem.mySimulation.wij[i][j];
                        }

                        for(int k=this.mySystem.buffer[j]+1;k<=this.upperBoundj[j];k++){
                            reversed_mcoefficient[j][k]+=this.myReversedBAPModel.Mijk[this.simulationLength+this.mySystem.buffer[j]+1-i][this.mySystem.nbStage-j][k]*this.mySystem.mySimulation.wij[i][j];
                        }
                    }
                }

                for(int k=this.lowerBoundj[j]+1;k<=this.mySystem.buffer[j];k++){
                    reversed_constantTerm+=reversed_mcoefficient[j][k];
                }
                
                double cumulateMjk=0;
                for(int k=this.lowerBoundj[j]+1;k<=this.upperBoundj[j];k++){
                    if(cumulateMjk > reversed_constantTerm)
                        reversed_mcoefficient[j][k] = 0;
                    else
                        cumulateMjk+=reversed_mcoefficient[j][k];                     
                }
            }

            // Calculate the coefficient of each term in the standard cut
            double constantTerm=this.newCut.constantTerm;
            double[][] mcoefficient=new double[mySystem.nbStage][];
            for(int j=1;j<=this.mySystem.nbStage-1;j++) {
                mcoefficient[j] = new double[this.upperBoundj[j] + 1];
                for (int k = this.lowerBoundj[j] + 1; k <= this.upperBoundj[j]; k++)
                    mcoefficient[j][k] = 0;
            }

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

                double cumulateMjk=0;
                for(int k=this.lowerBoundj[j]+1;k<=this.upperBoundj[j];k++){
                    if(cumulateMjk > constantTerm)
                        mcoefficient[j][k] = 0;
                    else
                        cumulateMjk+=mcoefficient[j][k];
                }
            }

            // Is combinatorial cut or standard cut tighter?
            boolean combCutIsTighterThanStandardCut=true;
            for(int j=1;j<=this.mySystem.nbStage-1;j++){
                double lhs=0;
                for(int k=this.lowerBoundj[j]+1;k<=min(this.mySystem.buffer[j]+1,this.upperBoundj[j]);k++){
                    lhs+=mcoefficient[j][k];
                }
                if(lhs<constantTerm)
                    combCutIsTighterThanStandardCut=false;
            }

            // Is combinatorial cut or reversed cut tighter?
            boolean combCutIsTighterThanReversedCut=true;
            for(int j=1;j<=this.mySystem.nbStage-1;j++){
                double lhs=0;
                for(int k=this.lowerBoundj[j]+1;k<=min(this.mySystem.buffer[j]+1,this.upperBoundj[j]);k++){
                    lhs+=reversed_mcoefficient[j][k];
                }
                if(lhs<constantTerm)
                    combCutIsTighterThanReversedCut=false;
            }
            
            // Which cut(s) to be added?
            boolean addCombCut = false;
            boolean addStandardCut = false;
            boolean addReversedCut = false;
            if (combCutIsTighterThanStandardCut && combCutIsTighterThanReversedCut)
                addCombCut = true;
            else if(combCutIsTighterThanStandardCut)
                addReversedCut = true;
            else if(combCutIsTighterThanReversedCut)
                addStandardCut = true;
            else{
                boolean similarCuts = true;
                for (int j=1;j<=this.mySystem.nbStage-1;j++){
                    if (!similarCuts)
                        break;
                    for(int k=this.lowerBoundj[j]+1;k<=min(this.mySystem.buffer[j]+1,this.upperBoundj[j]);k++){
                        if (abs(mcoefficient[j][k]-reversed_mcoefficient[j][k]) >= this.newCut.constantTerm * this.constantC / (this.mySystem.nbStage-1)){
                            similarCuts = false;
                            break;
                        }                         
                    }
                }
                if (similarCuts){
                    if (random()>0.5)
                        addStandardCut = true;
                    else
                        addReversedCut = true;                    
                }
                else{
                    addStandardCut = true;
                    addReversedCut = true;
                }
                    
            }

            // Add the correct cut
            if (addCombCut){
                IloLinearNumExpr feacut = cplex.linearNumExpr();
                IloRange rng;
                for(int j=1;j<=this.mySystem.nbStage-1;j++){
                    if(this.mySystem.buffer[j]<this.upperBoundj[j])
                        feacut.addTerm(1,this.yjk[j][this.mySystem.buffer[j]+1]);
                }
                rng=cplex.addGe(feacut,1.0);
                rng.setName("combcut of iter "+ (numit)+": ");
            }
            if(addStandardCut){
                IloLinearNumExpr feacut = cplex.linearNumExpr();
                IloRange rng;
                for(int j=1;j<=this.mySystem.nbStage-1;j++){
                    for(int k=this.lowerBoundj[j]+1;k<=this.upperBoundj[j];k++)
                        feacut.addTerm(-mcoefficient[j][k],this.yjk[j][k]);
                }
                feacut.setConstant(constantTerm);
                rng = cplex.addLe(feacut,0);
                rng.setName("standard feascut of iter "+ (numit)+": ");
            }
            if(addReversedCut){
                IloLinearNumExpr feacut = cplex.linearNumExpr();
                IloRange rng;
                for(int j=1;j<=this.mySystem.nbStage-1;j++){
                    for(int k=this.lowerBoundj[j]+1;k<=this.upperBoundj[j];k++)
                        feacut.addTerm(-reversed_mcoefficient[j][k],this.yjk[j][k]);
                }
                feacut.setConstant(reversed_constantTerm);
                rng = cplex.addLe(feacut,0);
                rng.setName("reversed feascut of iter "+ (numit)+": ");
            }
        }catch(Exception exc){exc.printStackTrace();}


    }

}
