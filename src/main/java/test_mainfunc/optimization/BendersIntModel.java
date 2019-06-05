package test_mainfunc.optimization;

import ilog.concert.IloException;
import test_mainfunc.simulation.SerialLine;

import static java.lang.Math.max;
import static java.lang.Math.min;

public abstract class BendersIntModel extends BendersBAP {

    FeasibilityCutCoef newCut;
    private double[][][] mijk;
    private double[][][] Mijk;
    private double theta;

    // Constructor
    BendersIntModel(SerialLine system, double THstar, int[] lB, int[] uB, int N, int W){
        super(system,THstar, lB, uB,N,W);

        this.newCut=new FeasibilityCutCoef();
        this.newCut.coefdeltaM=new double[this.mySystem.nbStage];
        this.newCut.coefdeltaP=new double[this.mySystem.nbStage];
        this.theta=(double) (N-W);
    }

    // Public methods
    public void solveBAPWithIntModel(double[][] tij) throws IloException {

        this.getMmValue(tij);
        this.mySystem.mySimulation = this.mySystem.new SimulationBAS(this.simulationLength,this.warmupLength,tij);

        for(int j=1;j<=mySystem.nbStage-1;j++)
            this.mySystem.buffer[j]=lowerBoundj[j];

        this.mySystem.mySimulation.simDualBAS(false);
        this.saveIterationSolution();

        while((this.THstar-this.mySystem.TH > 0.0001)&&(numit < this.MAX_ITE)){
            numit++;
            this.generateFeasibilityCut(tij);
            this.addFeasibilityCut();
            this.solveMasterProb();
            this.updateBufferSpaces();
            this.saveIterationSolution();
            this.mySystem.mySimulation.simDualBAS(false);
       }


       this.endMasterProb();
    }
    public abstract void addFeasibilityCut();
    public class FeasibilityCutCoef{
        public double[] coefdeltaP;
        public double[] coefdeltaM;
        public double constantTerm;
    }

    // Private & package-private methods
    private void getMmValue(double[][] tij){

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
                /*for(int k = this.lowerBoundj[j]; k <= this.upperBoundj[j]; k++){
                    this.mijk[i][j][k]=200000;
                    this.Mijk[i][j][k]=0;
                }*/
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
    private void generateFeasibilityCut(double[][] tij){

        for(int j=1;j<=this.mySystem.nbStage -1;j++)
        {
            this.newCut.coefdeltaP[j]=0.0;
            this.newCut.coefdeltaM[j]=0.0;
            for (int i= this.mySystem.buffer[j]+1;i<=this.simulationLength;i++)
            {
                this.newCut.coefdeltaP[j]= this.newCut.coefdeltaP[j] + (double)this.mySystem.mySimulation.wij[i][j]*this.Mijk[i][j][this.mySystem.buffer[j]];
                this.newCut.coefdeltaM[j]= this.newCut.coefdeltaM[j] + (double)this.mySystem.mySimulation.wij[i][j]*this.mijk[i][j][this.mySystem.buffer[j]];
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
    public void solveMasterProb() throws IloException{
        super.solveMasterProb();
        if(this.solvability){

            int totcap = 0;
            for(int j=1;j <= this.mySystem.nbStage-1;j++)
            {
                totcap += this.cplex.getValue(this.bj[j]);
                this.writer.print(Double.toString(this.cplex.getValue(this.bj[j]))+',');
            }
            this.writer.println("it " + numit + " OF: "+ totcap);
        }
    }

}
