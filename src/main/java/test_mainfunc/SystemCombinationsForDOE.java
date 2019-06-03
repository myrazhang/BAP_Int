package test_mainfunc;

import java.io.InputStream;
import java.util.Scanner;

import static java.lang.Math.sqrt;

public class SystemCombinationsForDOE {
    public int[] Jfactor;
    public double[] THfactor;
    public String distr;
    public double BNct;
    public double[] alphafactor;
    public double[] betafactor;
    public double[] noBNfactor;
    public double[] varfactor;
    public int[] BN1;
    public int[] BN2;
    public int Njobs;
    public int W;
    public int Lj;
    public int Uj;
    public String tempinstance;

    public SystemCombinationsForDOE(InputStream system){
        Scanner scanner=new Scanner (system);
        scanner.useDelimiter("\\s+");

        // scan number of jobs
        scanner.next();
        this.Njobs=scanner.nextInt();

        //scan warmup period
        scanner.next();
        this.W = scanner.nextInt();

        //scan Lowerbound Lj
        scanner.next();
        this.Lj = scanner.nextInt();

        //scan upperbound Uj
        scanner.next();
        this.Uj =  scanner.nextInt();

        //scan number of values for each factor
        //scan factor Number of stages
        scanner.next();
        int tempLength = scanner.nextInt();
        this.Jfactor = new int[tempLength];
        for (int j =0;j<this.Jfactor.length;j++){
            this.Jfactor[j] = scanner.nextInt();
        }

        //scan factor TH
        scanner.next();
        tempLength = scanner.nextInt();
        this.THfactor = new double[tempLength];
        for (int i=0; i<this.THfactor.length;i++){
            this.THfactor[i]=scanner.nextDouble();
        }


        //scan factor Distributions
        scanner.next();
        this.distr = scanner.next();
        //scan parameters size of factors and values
        scanner.next();
        //BN cycle time
        this.BNct = scanner.nextDouble();

        //no BN cycle time
        scanner.next();
        tempLength = scanner.nextInt();
        this.noBNfactor = new double[tempLength];
        for (int i=0; i<this.noBNfactor.length;i++){
            this.noBNfactor[i]=scanner.nextDouble();
        }
        //alpha parameter
        scanner.next();
        tempLength = scanner.nextInt();
        this.alphafactor = new double[tempLength];
        for (int i=0; i<this.alphafactor.length;i++){
            this.alphafactor[i]=scanner.nextDouble();
        }
        //beta parameter
        scanner.next();
        tempLength = scanner.nextInt();
        this.betafactor = new double[tempLength];
        for (int i=0; i<this.betafactor.length;i++){
            this.betafactor[i]=scanner.nextDouble();
        }

        scanner.next();
        tempLength = scanner.nextInt();
        this.varfactor = new double[tempLength];
        for (int i=0; i<this.varfactor.length;i++){
            this.varfactor[i]=scanner.nextDouble();
        }

        //scan bottlenecks information
        scanner.next();
        tempLength = scanner.nextInt();
        this.BN1 = new int[tempLength];
        this.BN2 = new int[tempLength];
        scanner.next();
        for (int i=0; i<this.BN1.length;i++){
            this.BN1[i]=scanner.nextInt();
        }
        scanner.next();
        for (int i=0; i<this.BN2.length;i++){
            this.BN2[i]=scanner.nextInt();
        }
    }

    public SerialLine getOneSystemConfiguration(int Jindex, int BNindex,int alphaindex, int noBNctindex, int varindex){

        SerialLine theSystem=new SerialLine();
        //save number of stages
        theSystem.nbStage = this.Jfactor[Jindex];
        theSystem.buffer = new int[theSystem.nbStage - 1];
        //save distribution information
        theSystem.CT = new StochNum[theSystem.nbStage];
        for (int j = 0; j < theSystem.nbStage; j++)
        {
            theSystem.CT[j]=new StochNum();
            theSystem.CT[j].distribution= this.distr;
            if (theSystem.CT[j].distribution.equals("Beta") ){
                this.CalculateBetaPar(j,BNindex,alphaindex,noBNctindex,varindex,theSystem);
            }//end if beta
            else if(theSystem.CT[j].distribution.equals("Exp")){
                theSystem.CT[j].para1 = this.alphafactor[alphaindex];
            }//end if exponential
        }
        return theSystem;
    }

    public void CalculateBetaPar(int j,int BNindex,int alphaindex, int noBNctindex,int varindex, SerialLine theSystem){
        double alpha = this.alphafactor[alphaindex] ;
        double beta = this.betafactor[alphaindex] ;
        double boolvariance = alpha*beta/((alpha+beta)*(alpha+beta))*(alpha+beta+1);
        double boolmean= alpha/(alpha+beta);
        double lowera=0;
        if(j==this.BN1[BNindex] || j == this.BN2[BNindex]){
            lowera = this.BNct - boolmean* sqrt(this.varfactor[varindex]/boolvariance);

        }//end if BN stage
        else{
            lowera = this.noBNfactor[noBNctindex] - boolmean* sqrt(this.varfactor[varindex]/boolvariance);
        }//end else no BN stage
        double upperb = lowera + sqrt(this.varfactor[varindex]/boolvariance);

        theSystem.CT[j].para1= alpha;
        theSystem.CT[j].para2= beta;
        theSystem.CT[j].para3= lowera;
        theSystem.CT[j].para4= upperb;
    }//end CalculateBetaPar function

}
