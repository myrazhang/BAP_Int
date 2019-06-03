package test_mainfunc;

import ilog.concert.IloNumVar;

import static java.lang.Math.max;
import static java.lang.Math.min;

public class BendersIntModel extends BendersBAP{

    public double[][][] mijk;
    public double[][][] Mijk;

    public FeasibilityCutCoef newCut;
    public int[][] wbarij;
    public int[][] ubarij;
    public double theta;


    public BendersIntModel(SerialLine system, double THstar, int[] lB, int[] uB,int N, int W){
        super(system,THstar, lB, uB,N,W);

        this.newCut=new FeasibilityCutCoef();
        this.newCut.coefdeltaM=new double[this.mySystem.nbStage-1];
        this.newCut.coefdeltaP=new double[this.mySystem.nbStage-1];
        this.wbarij=new int[this.simulationLength][this.mySystem.nbStage];
        this.ubarij=new int[this.simulationLength][this.mySystem.nbStage];
        this.theta=(double) (N-W);
    }


    public void solveBAPWithIntModel(double[][] tij){

        this.getMmValue(tij);
        this.mySystem.mySimulation = this.mySystem.new SimulationBAS(this.simulationLength,this.warmupLength,tij);

        for(int j=0;j<mySystem.nbStage-1;j++)
            this.mySystem.buffer[j]=lowerBoundj[j];

        this.mySystem.mySimulation.simDualBAS(false);
        this.saveIterationSolution();

        this.numit=0;
        while((this.THstar-this.mySystem.TH > 0.0001)&&(numit < this.MAX_ITE)){
            numit++;
            this.generateFeasibilityCut(tij);
            this.addFeasibilityCut();
            this.solveMasterProb();
            this.updateBufferSpaces();
            this.saveIterationSolution();
            this.mySystem.mySimulation.simDualBAS(false);
       }

       try{
           this.writer.println(this.cplex.getModel());
       }catch(Exception exc){exc.printStackTrace();}

       this.endMasterProb();

    }


    public void getMmValue(double[][] tij){

        // initialize
        this.Mijk=new double[this.simulationLength][this.mySystem.nbStage-1][];
        this.mijk=new double[this.simulationLength][this.mySystem.nbStage-1][];
        for (int j = 1; j <= this.mySystem.nbStage-1; j++){
            for (int i = 1; i <= this.simulationLength; i++){
                this.Mijk=new double[i][j][this.upperBoundj[j]-1];
                this.mijk=new double[i][j][this.upperBoundj[j]-1];
                for(int k = 1; k <= this.upperBoundj[j]; k++){
                    this.mijk[i-1][j-1][k-1]=0;
                    this.Mijk[i-1][j-1][k-1]=100000;
                }
                for(int k = this.lowerBoundj[j]; k <= this.upperBoundj[j]; k++){
                    this.mijk[i-1][j-1][k-1]=200000;
                    this.Mijk[i-1][j-1][k-1]=0;
                }
            }
        }

        // get value of small mijk
        for(int j = 1;j <= this.mySystem.nbStage-1; j++){
            for(int k = this.lowerBoundj[j]; k <= this.upperBoundj[j]; k++){
                for (int i = k+1; i <= this.simulationLength; i++) {
                    for(int i0 = i - k + 1; i0 <= i-this.lowerBoundj[j] ; i0++){
                        this.mijk[i-1][j-1][k-1] = min( this.mijk[i-1][j-1][k-1], tij[i0-1][j]);
                    }
                }
            }
        }

        //get value of big Mijk
        int[][] U_j_deltaj=new int[this.mySystem.nbStage-1][this.mySystem.nbStage-1];
        int[][] L_j_deltaj=new int[this.mySystem.nbStage-1][this.mySystem.nbStage-1];

        for (int j = 1; j <= this.mySystem.nbStage-1; j++){
            for(int delta_j = 1; delta_j <= this.mySystem.nbStage-j; delta_j++){
                U_j_deltaj[j-1][delta_j-1]=delta_j+1;
                L_j_deltaj[j-1][delta_j-1]=delta_j+1;
                for (int l=j; l <= j+delta_j-1 ; l++){
                    U_j_deltaj[j-1][delta_j-1]=U_j_deltaj[j-1][delta_j-1]+this.upperBoundj[j];
                    L_j_deltaj[j-1][delta_j-1]=L_j_deltaj[j-1][delta_j-1]+this.upperBoundj[j];
                }
            }
        }

        double[][] sij= new double [this.simulationLength][this.mySystem.nbStage];
        for (int i = 2; i <= this.simulationLength; i++){
            for (int j  = 1; j <= this.mySystem.nbStage-1; j++){
                sij[i-1][j-1]=tij[i-2][j-1];
                for (int delta_j = 1; delta_j <= this.mySystem.nbStage-j ; delta_j++){
                    for (int i0= max(i - U_j_deltaj[j-1][delta_j-1],1);i0<= min(i - L_j_deltaj[j-1][delta_j-1],this.simulationLength);i0++){
                        sij[i-1][j-1]=max(sij[i-1][j-1], tij[i0-1][j+delta_j-1]);
                    }
                }
            }
            sij[i-1][this.mySystem.nbStage-1]=tij[i-2][this.mySystem.nbStage-1];
        }

        for(int j=1;j <= this.mySystem.nbStage-1;j++){
            for(int k=this.lowerBoundj[j]; k <= this.upperBoundj[j];k++){
                for (int i = k+1 ; i <= this.simulationLength; i++){
                    for (int i0=max(i-this.upperBoundj[j]+1,1);i0<= min(i-k,this.simulationLength);i0++){
                        this.Mijk[i-1][j-1][k-1]=max(this.Mijk[i-1][j-1][k-1],sij[i0-1][j]);
                    }
                }
            }
        }
    }

    public void generateFeasibilityCut(double[][] tij){

        for(int j=0;j<this.mySystem.nbStage -1;j++)
        {
            this.newCut.coefdeltaP[j]=0.0;
            this.newCut.coefdeltaM[j]=0.0;
            for (int i=0;i<this.simulationLength;i++)
            {
                this.newCut.coefdeltaP[j]= this.newCut.coefdeltaP[j] + (double)this.wbarij[i][j]*this.Mijk[i][j][this.mySystem.buffer[j]];
                this.newCut.coefdeltaM[j]= this.newCut.coefdeltaM[j] + (double)this.wbarij[i][j]*this.mijk[i][j][this.mySystem.buffer[j]];
            }
            this.newCut.coefdeltaP[j]=-this.newCut.coefdeltaP[j];
        }

        double tijpar=0.0;
        for(int j=0;j<this.mySystem.nbStage;j++)
        {
            for (int i=0;i<this.simulationLength;i++)
            {
                tijpar = tijpar + tij[i][j]*(double)this.ubarij[i][j];
            }
        }

        this.newCut.constantTerm= - this.theta/this.THstar + tijpar;

    }

    public void addFeasibilityCut(){ }

    public class FeasibilityCutCoef{
        public double[] coefdeltaP;
        public double[] coefdeltaM;
        public double constantTerm;
    }

    public class IterativelyAddedVariables{
        IloNumVar[] listOfVariables;

        public IterativelyAddedVariables(int J){
            listOfVariables=new IloNumVar[J];
        }
    }
}
