package test_mainfunc.optimization;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.util.Stopwatch;

import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.HashMap;

import static java.lang.Math.max;
import static java.lang.Math.min;

public class BendersIntModelAlter6 extends BendersIntModel {

    protected IloNumVar[][] yjk;

    public BendersIntModelAlter6(SerialLine system, double THstar, int[] lB, int[] uB, int N, int W){
        super(system, THstar,  lB, uB, N,  W);
    }

    public BendersIntModelAlter6(SerialLine mySystem, double tHstar, int[] lowerBoundj, int[] upperBoundj, int simulationLength, int warmupLength, int startMachine, int endMachine, HashMap<String, Integer> initialBounds) {
        super(mySystem, tHstar,  lowerBoundj, upperBoundj, simulationLength, warmupLength, startMachine, endMachine, initialBounds);
    }

    /*@Override
    public void solveBAPWithIntModel(double[][] tij,int offset,boolean resetMm)throws IloException {
        this.writer.println("Alter 6:");
        super.solveBAPWithIntModel(tij,offset,resetMm);
    }*/



    @Override
    public boolean solveBAPWithIntModel(double[][] tij,boolean redefineBounds)throws IloException {
        this.writer.println("Alter 6:");
        return super.solveBAPWithIntModel(tij,redefineBounds);
    }


    public boolean solveBAPWithIntModel(double[][] tij,HashMap<String,Integer> initialLowerBounds)throws IloException {
        totalTimeMeasure.start();
        this.getMmValue(tij,initialLowerBounds);
        this.mySystem.mySimulation = this.mySystem.new SimulationBAS(this.simulationLength,this.warmupLength,tij);

        /*for(int j=1;j<=mySystem.nbStage-1;j++)
            this.mySystem.buffer[j]=lowerBoundj[j];*/
        solvability=true;
        this.solveMasterProb(MAX_TIME);
        if(!solvability){
            System.out.println("Initial master problem is not solvable!");
            return false;
        }
        this.updateBufferSpaces();
        this.mySystem.mySimulation.simDualBAS(false);
        this.saveIterationSolution();

        double tolerance = TOLERANCE_RATIO*this.THstar;
        while((this.THstar-this.mySystem.TH > tolerance)&&(numit < this.MAX_ITE)&&solvability){
            numit++;
            this.generateFeasibilityCut(tij);
            this.addFeasibilityCut();
            totalTimeMeasure.pause();
            double timeLimit = MAX_TIME-totalTimeMeasure.elapseTimeSeconds;
            totalTimeMeasure.restart();
            this.solveMasterProb(timeLimit);
            this.updateBufferSpaces();
            this.saveIterationSolution();
            this.mySystem.mySimulation.simDualBAS(false);
        }

        writer.print(cplex.getModel());

        //this.endMasterProb();
        return solvability;
    }

    @Override
    public void solveWithDecomposition(double[][] tij) throws IloException {
        /*for(int j=1;j<=mySystem.nbStage-1;j++){
            mySystem.buffer[j] = lowerBoundj[j];
        }
        HashMap<String, Integer> initialBounds = new HashMap<>();
        numit = 0;
        double remainTime = MAX_CPLEX_TIME;
        for(int m=2;m<=mySystem.nbStage;m++){
            for(int startMachine = 1; startMachine<=mySystem.nbStage-1;startMachine++){
                int endMachine = startMachine + m-1;
                if(endMachine<=mySystem.nbStage){
                    Stopwatch sw = new Stopwatch();
                    sw.start();
                    BendersIntModelAlter6 subSystemBap = new BendersIntModelAlter6(mySystem, THstar, lowerBoundj, upperBoundj, simulationLength, warmupLength,
                            startMachine, endMachine, initialBounds);
                    subSystemBap.writer = new PrintWriter(OutputStream.nullOutputStream());
                    subSystemBap.setMaxCplexTime(remainTime);
                    subSystemBap.solveBAPWithIntModel(tij,startMachine-1);
                    int solSubProb = 0;
                    for(int j=1;j<=subSystemBap.mySystem.nbStage-1;j++){
                        solSubProb+=subSystemBap.mySystem.buffer[j];
                    }
                    addInitialBounds(initialBounds,startMachine,endMachine,solSubProb);
                    numit+=subSystemBap.numit;
                    if(m==mySystem.nbStage){
                        for(int j=1;j<=m-1;j++){
                            mySystem.buffer[j] = subSystemBap.mySystem.buffer[j];
                        }
                        break;
                    }
                    sw.stop();
                    remainTime -= sw.elapseTimeSeconds;
                    if(remainTime<0){
                        subSystemBap = new BendersIntModelAlter6(mySystem, THstar, lowerBoundj, upperBoundj, simulationLength, warmupLength,
                                1, mySystem.nbStage, initialBounds);
                        subSystemBap.setMaxIte(1);
                        subSystemBap.solveBAPWithIntModel(tij,0);
                        for(int j=1;j<=m-1;j++){
                            mySystem.buffer[j] = subSystemBap.mySystem.buffer[j];
                        }
                        break;
                    }

                }
            }
        }*/
    }


    @Override
    public void getMmValue(
            double[][] tij){
        // initialize
        this.Mijk=new double[this.simulationLength+1][this.mySystem.nbStage][];
        this.mijk=new double[this.simulationLength+1][this.mySystem.nbStage][];
        for (int j = 0; j <= this.mySystem.nbStage-1; j++){
            for (int i = 0; i <= this.simulationLength; i++){
                this.Mijk[i][j]=new double[this.upperBoundj[j]+1];
                this.mijk[i][j]=new double[this.upperBoundj[j]+1];
                for(int k = 0; k <= this.upperBoundj[j]; k++){
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

    public void getMmValue(
            double[][] tij, HashMap<String,Integer> initialLowerBounds//, int startMachine, int endMachine
            //,HashMap<String,Integer> initialUpperBounds
    ){
        // initialize
        this.Mijk=new double[this.simulationLength+1][this.mySystem.nbStage][];
        this.mijk=new double[this.simulationLength+1][this.mySystem.nbStage][];
        for (int j = 0; j <= this.mySystem.nbStage-1; j++){
            for (int i = 0; i <= this.simulationLength; i++){
                this.Mijk[i][j]=new double[this.upperBoundj[j]+1];
                this.mijk[i][j]=new double[this.upperBoundj[j]+1];
                for(int k = 0; k <= this.upperBoundj[j]; k++){
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
        int incumbent=0;
        int startMachine=-1, endMachine=-1;
        for(int j=1;j<this.mySystem.nbStage;j++){
            if(this.lowerBoundj[j]!=this.upperBoundj[j] && startMachine==-1) startMachine = j;
            if(this.lowerBoundj[j]!=this.upperBoundj[j]) endMachine = j+1;
        }
        for(int j=1;j<this.mySystem.nbStage;j++){
            if(j<startMachine || j>=endMachine) incumbent+=this.upperBoundj[j];
        }
        if(endMachine-startMachine==1){
            incumbent += this.upperBoundj[startMachine];
        }else{
            incumbent += min(initialLowerBounds.get((startMachine+1)+"_"+endMachine)+this.upperBoundj[startMachine],
                    initialLowerBounds.get(startMachine+"_"+(endMachine-1))+this.upperBoundj[endMachine-1]);
        }

        double[][] sij= new double [this.simulationLength+1][this.mySystem.nbStage+1];
        int[][] sumUj=new int[this.mySystem.nbStage+1][this.mySystem.nbStage+1];
        int[][] sumLj=new int[this.mySystem.nbStage+1][this.mySystem.nbStage+1];
        for(int j=startMachine;j<endMachine;j++) {
            for (int j1=j+1;j1<=endMachine;j1++) {
                sumLj[j][j1]=0;
                if(endMachine-startMachine==1) sumLj[j][j1]=lowerBoundj[j];
                else if(initialLowerBounds.containsKey(j+"_"+j1)) sumLj[j][j1]=initialLowerBounds.get(j+"_"+j1);
                else {
                    for(int k = startMachine+1;k<=endMachine-1;k++) {
                        sumLj[j][j1] = max(sumLj[j][j1],initialLowerBounds.get(startMachine+"_"+k)+initialLowerBounds.get(k+"_"+endMachine));
                    }
                }
            }
        }
        for(int j=1;j<this.mySystem.nbStage;j++) {
            for (int j1 = j + 1; j1 <= this.mySystem.nbStage; j1++) {
                if (j1 <= startMachine || j1 > endMachine) {
                    if (j1==2) sumLj[j][j1] = this.lowerBoundj[1];
                    else sumLj[j][j1] = sumLj[j][j1 - 1] + this.lowerBoundj[j1 - 1];
                }
                else if(j<startMachine) sumLj[j][j1] = sumLj[j][startMachine]+sumLj[startMachine][j1];
            }
        }

        for(int j=1;j<this.mySystem.nbStage;j++){
            for(int j1=j+1;j1<=this.mySystem.nbStage;j1++){
                sumUj[j][j1] = incumbent;
                if(j>1) sumUj[j][j1] -= sumLj[1][j];
                if(j1<this.mySystem.nbStage) sumUj[j][j1] -= sumLj[j1][this.mySystem.nbStage];
                int a = 0;
                for(int k=j;k<j1;k++) a += this.upperBoundj[k];
                sumUj[j][j1] = min(sumUj[j][j1],a);
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
    public void generateFeasibilityCut(double[][] tij,int offset){
        double tijpar=0.0;
        for(int j=1;j<=this.mySystem.nbStage;j++)
        {
            for (int i=1;i<=this.simulationLength;i++)
            {
                if(this.mySystem.mySimulation.uij[i][j]>0)
                    tijpar += tij[i][j+offset]*(double)this.mySystem.mySimulation.uij[i][j];
            }
        }
        this.newCut.constantTerm= - this.theta/this.THstar + tijpar;
    }

    @Override
    public void generateFeasibilityCut(double[][] tij){
        /*double tijpar=0.0;
        for(int j=1;j<=this.mySystem.nbStage;j++)
        {
            for (int i=1;i<=this.simulationLength;i++)
            {
                if(this.mySystem.mySimulation.uij[i][j]>0)
                    tijpar += tij[i][j]*(double)this.mySystem.mySimulation.uij[i][j];
            }
        }
        this.newCut.constantTerm= - this.theta/this.THstar + tijpar;*/
        this.newCut.constantTerm= this.theta*(1/this.mySystem.TH - 1/this.THstar);
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
                for(int k=0;k<=this.lowerBoundj[j];k++){
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
