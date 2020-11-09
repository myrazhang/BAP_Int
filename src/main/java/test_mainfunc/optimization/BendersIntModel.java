package test_mainfunc.optimization;

import ilog.concert.IloException;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.util.Stopwatch;

import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.HashMap;

import static java.lang.Math.max;
import static java.lang.Math.min;

public abstract class BendersIntModel extends BendersBAP {

    FeasibilityCutCoef newCut;
    public double[][][] mijk;
    public double[][][] Mijk;
    protected double theta;

    // Constructor
    BendersIntModel(SerialLine system, double THstar, int[] lB, int[] uB, int N, int W){
        super(system,THstar, lB, uB,N,W);

        this.newCut=new FeasibilityCutCoef();
        this.newCut.coefdeltaM=new double[this.mySystem.nbStage];
        this.newCut.coefdeltaP=new double[this.mySystem.nbStage];
        this.theta=(double) (N-W);
    }

    BendersIntModel(SerialLine system, double THstar, int[] lB, int[] uB, int N, int W, int j1, int j2,
                           HashMap<String,Integer> initialBounds){
        super(system, THstar,  lB, uB, N, W, j1, j2, initialBounds);
        this.newCut=new FeasibilityCutCoef();
        this.newCut.coefdeltaM=new double[this.mySystem.nbStage];
        this.newCut.coefdeltaP=new double[this.mySystem.nbStage];
        this.theta=(double) (N-W);
    }

    // Public methods
    /*public void solveBAPWithIntModel(double[][] tij, int offset, boolean resetMm) throws IloException {
        totalTimeMeasure.start();
        if(resetMm){
            this.getMmValue(tij,offset);
        }
        this.mySystem.mySimulation = this.mySystem.new SimulationBAS(this.simulationLength,this.warmupLength,tij);

        this.solveMasterProb(MAX_TIME);
        this.updateBufferSpaces();
        this.mySystem.mySimulation.simDualBAS(false,offset);
        this.saveIterationSolution();


        solvability=true;
        while((this.THstar-this.mySystem.TH > 0)&&(numit < this.MAX_ITE)&&solvability){
            numit++;
            this.generateFeasibilityCut(tij,offset);
            this.addFeasibilityCut();
            totalTimeMeasure.pause();
            double timeLimit = MAX_TIME - totalTimeMeasure.elapseTimeSeconds;
            totalTimeMeasure.restart();
            this.solveMasterProb(timeLimit);
            this.updateBufferSpaces();
            this.saveIterationSolution();
            this.mySystem.mySimulation.simDualBAS(false,offset);
       }

       //this.endMasterProb();
    }*/

    public boolean solveBAPWithIntModel(double[][] tij,boolean redefineBounds) throws IloException {
        totalTimeMeasure.start();
        this.getMmValue(tij);
        this.mySystem.mySimulation = this.mySystem.new SimulationBAS(this.simulationLength,this.warmupLength,tij);

        /*for(int j=1;j<=mySystem.nbStage-1;j++)
            this.mySystem.buffer[j]=lowerBoundj[j];*/
        this.solveMasterProb(MAX_TIME);
        this.updateBufferSpaces();
        this.mySystem.mySimulation.simDualBAS(false);
        this.saveIterationSolution();

        solvability=true;
        double tolerance = TOLERANCE_RATIO*this.THstar;
        while((this.THstar-this.mySystem.TH >= tolerance)&&(numit < this.MAX_ITE)&&solvability){
            numit++;
            this.generateFeasibilityCut(tij);
            this.addFeasibilityCut();
            totalTimeMeasure.pause();
            double timeLimit = MAX_TIME - totalTimeMeasure.elapseTimeSeconds;
            totalTimeMeasure.restart();
            this.solveMasterProb(timeLimit);
            this.updateBufferSpaces();
            this.saveIterationSolution();
            this.mySystem.mySimulation.simDualBAS(false);
        }
        totalTimeMeasure.stop();
        return solvability;
        //this.endMasterProb();
    }

    public abstract void addFeasibilityCut();
    public class FeasibilityCutCoef{
        public double[] coefdeltaP;
        public double[] coefdeltaM;
        public double constantTerm;
    }

    public void solveWithDecomposition(double[][] tij) throws IloException{} /*{
        HashMap<String, Integer> initialBounds = new HashMap<>();

        for(int m=2;m<=mySystem.nbStage;m++){
            for(int startMachine = 1; startMachine<=mySystem.nbStage-1;startMachine++){
                int endMachine = startMachine + m-1;
                if(endMachine<=mySystem.nbStage){
                    BendersIntModel subSystemBap = new BendersIntModel(mySystem, THstar, lowerBoundj, upperBoundj, simulationLength, warmupLength,
                            startMachine, endMachine, initialBounds);
                    subSystemBap.writer = new PrintWriter(OutputStream.nullOutputStream());
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
                    }
                }
            }
        }
    }*/

    // Private & package-private methods
    public void getMmValue(
            double[][] tij,
            int offset){

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
                    this.mijk[i][j][k]=200000;
                    for(int i0 = i - k + 1; i0 <= i-this.lowerBoundj[j] ; i0++){
                        this.mijk[i][j][k] = min( this.mijk[i][j][k], tij[i0][j+1+offset]);
                    }
                }
            }
        }

        //get value of big Mijk
        int[][] U_j_deltaj=new int[this.mySystem.nbStage][this.mySystem.nbStage];
        int[][] L_j_deltaj=new int[this.mySystem.nbStage][this.mySystem.nbStage];

        for (int j = 1; j <= this.mySystem.nbStage-1; j++){
            for(int delta_j = 1; delta_j <= this.mySystem.nbStage-j; delta_j++){
                U_j_deltaj[j][delta_j]=delta_j+1;
                L_j_deltaj[j][delta_j]=delta_j+1;
                for (int l=j; l <= j+delta_j-1 ; l++){
                    U_j_deltaj[j][delta_j]=U_j_deltaj[j][delta_j]+this.upperBoundj[j];
                    L_j_deltaj[j][delta_j]=L_j_deltaj[j][delta_j]+this.lowerBoundj[j];
                }
            }
        }

        double[][] sij= new double [this.simulationLength+1][this.mySystem.nbStage+1];
        for (int i = 2; i <= this.simulationLength; i++){
            for (int j  = 1; j <= this.mySystem.nbStage-1; j++){
                sij[i][j]=tij[i-1][j+offset];
                for (int delta_j = 1; delta_j <= this.mySystem.nbStage-j ; delta_j++){
                    for (int i0= max(i - U_j_deltaj[j][delta_j],1);i0<= min(i - L_j_deltaj[j][delta_j],this.simulationLength);i0++){
                        sij[i][j]=max(sij[i][j], tij[i0][j+delta_j+offset]);
                    }
                }
            }
            sij[i][this.mySystem.nbStage]=tij[i-1][this.mySystem.nbStage+offset];
        }

        for(int j=1;j <= this.mySystem.nbStage-1;j++){
            for(int k=this.lowerBoundj[j]; k <= this.upperBoundj[j];k++){
                for (int i = k+1 ; i <= this.simulationLength; i++){
                    this.Mijk[i][j][k]=0;
                    for (int i0=max(i-this.upperBoundj[j]+1,1);i0<= min(i-k,this.simulationLength);i0++){
                        this.Mijk[i][j][k]=max(this.Mijk[i][j][k],sij[i0][j+1]);
                    }
                }
            }
        }
    }
    public void getMmValue(
            double[][] tij
            ){

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
                    this.mijk[i][j][k]=200000;
                    for(int i0 = i - k + 1; i0 <= i-this.lowerBoundj[j] ; i0++){
                        this.mijk[i][j][k] = min( this.mijk[i][j][k], tij[i0][j+1]);
                    }
                }
            }
        }

        //get value of big Mijk
        int[][] U_j_deltaj=new int[this.mySystem.nbStage][this.mySystem.nbStage];
        int[][] L_j_deltaj=new int[this.mySystem.nbStage][this.mySystem.nbStage];

        for (int j = 1; j <= this.mySystem.nbStage-1; j++){
            for(int delta_j = 1; delta_j <= this.mySystem.nbStage-j; delta_j++){
                U_j_deltaj[j][delta_j]=delta_j+1;
                L_j_deltaj[j][delta_j]=delta_j+1;
                for (int l=j; l <= j+delta_j-1 ; l++){
                    U_j_deltaj[j][delta_j]=U_j_deltaj[j][delta_j]+this.upperBoundj[j];
                    L_j_deltaj[j][delta_j]=L_j_deltaj[j][delta_j]+this.lowerBoundj[j];
                }
            }
        }

        double[][] sij= new double [this.simulationLength+1][this.mySystem.nbStage+1];
        for (int i = 2; i <= this.simulationLength; i++){
            for (int j  = 1; j <= this.mySystem.nbStage-1; j++){
                sij[i][j]=tij[i-1][j];
                for (int delta_j = 1; delta_j <= this.mySystem.nbStage-j ; delta_j++){
                    for (int i0= max(i - U_j_deltaj[j][delta_j],1);i0<= min(i - L_j_deltaj[j][delta_j],this.simulationLength);i0++){
                        sij[i][j]=max(sij[i][j], tij[i0][j+delta_j]);
                    }
                }
            }
            sij[i][this.mySystem.nbStage]=tij[i-1][this.mySystem.nbStage];
        }

        for(int j=1;j <= this.mySystem.nbStage-1;j++){
            for(int k=this.lowerBoundj[j]; k <= this.upperBoundj[j];k++){
                for (int i = k+1 ; i <= this.simulationLength; i++){
                    this.Mijk[i][j][k]=0;
                    for (int i0=max(i-this.upperBoundj[j]+1,1);i0<= min(i-k,this.simulationLength);i0++){
                        this.Mijk[i][j][k]=max(this.Mijk[i][j][k],sij[i0][j+1]);
                    }
                }
            }
        }
    }
    public void generateFeasibilityCut(
            double[][] tij,
            int offset
    ){

        for(int j=1;j<=this.mySystem.nbStage -1;j++)
        {
            this.newCut.coefdeltaP[j]=0.0;
            this.newCut.coefdeltaM[j]=0.0;
            for (int i= this.mySystem.buffer[j]+1;i<=this.simulationLength;i++)
            {
                if(this.mySystem.mySimulation.wij[i][j]>0){
                    this.newCut.coefdeltaP[j]+=  (double)this.mySystem.mySimulation.wij[i][j]*this.Mijk[i][j][this.mySystem.buffer[j]];
                    this.newCut.coefdeltaM[j]+=  (double)this.mySystem.mySimulation.wij[i][j]*this.mijk[i][j][this.mySystem.buffer[j]];
                }
            }
            this.newCut.coefdeltaP[j]=-this.newCut.coefdeltaP[j];
        }

        double tijpar=0.0;
        for(int j=1;j<=this.mySystem.nbStage;j++)
        {
            for (int i=1;i<=this.simulationLength;i++)
            {
                tijpar = tijpar + tij[i][j+offset]*(double)this.mySystem.mySimulation.uij[i][j];
            }
        }
        this.newCut.constantTerm= - this.theta/this.THstar + tijpar;

    }

    public void generateFeasibilityCut(
            double[][] tij
    ){

        for(int j=1;j<=this.mySystem.nbStage -1;j++)
        {
            this.newCut.coefdeltaP[j]=0.0;
            this.newCut.coefdeltaM[j]=0.0;
            for (int i= this.mySystem.buffer[j]+1;i<=this.simulationLength;i++)
            {
                if(this.mySystem.mySimulation.wij[i][j]>0){
                    this.newCut.coefdeltaP[j]+=  (double)this.mySystem.mySimulation.wij[i][j]*this.Mijk[i][j][this.mySystem.buffer[j]];
                    this.newCut.coefdeltaM[j]+=  (double)this.mySystem.mySimulation.wij[i][j]*this.mijk[i][j][this.mySystem.buffer[j]];
                }
            }
            this.newCut.coefdeltaP[j]=-this.newCut.coefdeltaP[j];
        }

        double tijpar=0.0;
        for(int j=1;j<=this.mySystem.nbStage;j++)
        {
            for (int i=1;i<=this.simulationLength;i++)
            {
                tijpar = tijpar + tij[i][j]*(double)this.mySystem.mySimulation.uij[i][j];
            }
        }
        this.newCut.constantTerm= - this.theta/this.THstar + tijpar;

    }


    @Override
    public void solveMasterProb(double timeLimit) throws IloException{
        super.solveMasterProb(timeLimit);
        if(this.solvability){

            int totcap = 0;
            for(int j=1;j <= this.mySystem.nbStage-1;j++)
            {
                totcap += (int) (this.cplex.getValue(this.bj[j])+0.1);
                this.writer.print(Double.toString(this.cplex.getValue(this.bj[j]))+',');
            }
            this.writer.println("it " + numit + " OF: "+ totcap);
        }
    }

}
