package test_mainfunc;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.dataformat.yaml.YAMLMapper;
import test_mainfunc.simulation.Failure;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.simulation.StochNum;
import test_mainfunc.input.DoeInputV1;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.Scanner;
import java.util.stream.StreamSupport;

import static java.lang.Math.sqrt;

public class SystemCombinationsForDOE {
    public int Lj;
    public int Uj;
    int[] Jfactor;
    double[] etaFactor;
    double[] alphafactor;
    double[] noBNfactor;
    double[] varfactor;
    int[] BN1;
    int[] BN2;
    int Njobs;
    int W;
    String tempinstance;
    private String distr;
    double BNct;
    int[] difffailedstage;
    double[] diffuprateFactor;
    double[] diffdownrateFactor;
    double[] iiduprateFactor;
    double[] iiddownrateFactor;
    private double[] betafactor;

    SystemCombinationsForDOE(InputStream system) throws IOException {
        ObjectMapper mapper = new YAMLMapper();
        DoeInputV1 input = mapper.readValue(system, DoeInputV1.class);
        this.Njobs=input.Jobs;
        this.W = input.Warmup;
        this.Lj = input.Lj;
        this.Uj = input.Uj;
        this.Jfactor = input.NbStage_factors;
        this.etaFactor = input.Eta;
        this.distr = input.Distribution;
        this.BNct = input.BN_CT_m;
        this.noBNfactor = input.noBN_CT_m;
        this.alphafactor = input.sigma;
        this.betafactor = input.beta;
        this.varfactor = input.variance;
        this.BN1 = Arrays.stream(input.BN1_position.split(" "))
                .mapToInt(Integer::valueOf).toArray();
        this.BN2 = Arrays.stream(input.BN2_position.split(" "))
                .mapToInt(Integer::valueOf).toArray();
        this.difffailedstage = input.differentfailedmachine;
        this.diffuprateFactor = input.different_uptime_rate;
        this.diffdownrateFactor = input.different_downtime_rate;
        this.iiduprateFactor = input.iid_uptime_rate;
        this.iiddownrateFactor = input.iid_downtime_rate;
    }
    SerialLine getOneSystemConfiguration(int Jindex, int BNindex, int alphaindex, int noBNctindex, int varindex){

        SerialLine theSystem=new SerialLine();
        //save number of stages
        theSystem.nbStage = this.Jfactor[Jindex];
        theSystem.buffer = new int[theSystem.nbStage];
        //save distribution information
        theSystem.CT = new StochNum[theSystem.nbStage+1];
        for (int j = 1; j <= theSystem.nbStage; j++)
        {
            theSystem.CT[j]=new StochNum();
            theSystem.CT[j].distribution= this.distr;
            if (theSystem.CT[j].distribution.equals("Beta") ){
                this.CalculateBetaPar(j,BNindex,alphaindex,noBNctindex,varindex,theSystem);
            }//end if beta
            else if(theSystem.CT[j].distribution.equals("Exp")){
                theSystem.CT[j].para1 = this.alphafactor[alphaindex];
            }//end if exponential
            else if(theSystem.CT[j].distribution.equals("Deterministic")){
                theSystem.CT[j].para1 = this.alphafactor[alphaindex];
            }//end if Deterministic
            else if(theSystem.CT[j].distribution.equals("Erlang")){
                theSystem.CT[j].para1 = this.alphafactor[alphaindex];
                theSystem.CT[j].para2 = this.betafactor[alphaindex];
            }//end if Erlang
            else if(theSystem.CT[j].distribution.equals("LogNorm")){
                if(j==this.BN1[BNindex] || j == this.BN2[BNindex]){
                    theSystem.CT[j].para1 = this.BNct;
                }//end if BN stage
                else{
                    theSystem.CT[j].para1 = this.noBNfactor[noBNctindex];
                }//end else no BN stage
                theSystem.CT[j].para2 = this.alphafactor[alphaindex];
            }//end if LogNorm
            else
                throw new UnsupportedOperationException("Cycle time distribution is not supported!");
        }
        return theSystem;
    }

    SerialLine getStLongLineConfiguration(int Jindex){

        SerialLine theSystem=new SerialLine();
        //save number of stages
        theSystem.nbStage = this.Jfactor[Jindex];
        theSystem.buffer = new int[theSystem.nbStage];
        //save distribution information
        theSystem.CT = new StochNum[theSystem.nbStage+1];
        for (int j = 1; j <= theSystem.nbStage; j++)
        {
            theSystem.CT[j]=new StochNum();
            if(((j == 1) || (j == 3) || (j == 5) || (j == 7) || (j == 9) || (j == 11) || (j == 15) || (j == 17) || (j == 19) || (j==21) || (j==23)))
            {
                theSystem.CT[j].distribution= "Erlang";
                theSystem.CT[j].para1 = this.alphafactor[0];
                theSystem.CT[j].para2 = this.betafactor[0];
                // System.out.println("para 1: " + theSystem.CT[j].para1 + " and para 2: "+ theSystem.CT[j].para2);
            }
            else
            {
                theSystem.CT[j].distribution= "Deterministic";
                theSystem.CT[j].para1 = this.BNct;
            }
        }
        return theSystem;
    }

    SerialLine getTemplA(){

        SerialLine theSystem=new SerialLine();
        //save number of stages
        theSystem.nbStage = 19;
        double[] meanpt = new double[] {0, 0.47, 0.47, 0.493, 0.488, 0.463, 0.463, 0.498, 0.467, 0.473, 0.487, 0.498, 0.478, 0.488, 0.488, 0.482, 0.43, 0.473, 0.5, 0.467};
        theSystem.buffer = new int[theSystem.nbStage];
        //save distribution information
        theSystem.CT = new StochNum[theSystem.nbStage+1];
        for (int j = 1; j <= theSystem.nbStage; j++)
        {
            theSystem.CT[j]=new StochNum();
            theSystem.CT[j].distribution= "Deterministic";
            theSystem.CT[j].para1 = meanpt[j];
        }
        return theSystem;
    }

    SerialLine getTemplB(){

        SerialLine theSystem=new SerialLine();
        //save number of stages
        theSystem.nbStage = 23;
        double[] meanpt = new double[] {0, 5.1, 5.1, 5.3, 5.3, 3.8, 5.2, 5.2, 5, 3.5, 5.1, 5.5, 5.3, 4.7, 3.5, 5, 5.1, 4.1, 4.9, 4.9, 4.9, 6.1, 5.1, 5.1};
        theSystem.buffer = new int[theSystem.nbStage];
        //save distribution information
        theSystem.CT = new StochNum[theSystem.nbStage+1];
        for (int j = 1; j <= theSystem.nbStage; j++)
        {
            theSystem.CT[j]=new StochNum();
            theSystem.CT[j].distribution= "Deterministic";
            theSystem.CT[j].para1 = meanpt[j];
        }
        return theSystem;
    }

    SerialLine getTemplC(){

        SerialLine theSystem=new SerialLine();
        //save number of stages
        theSystem.nbStage = 8;
        double[] meanpt = new double[] {0, 233, 233, 234, 216, 212, 220, 255, 257};
        theSystem.buffer = new int[theSystem.nbStage];
        //save distribution information
        theSystem.CT = new StochNum[theSystem.nbStage+1];
        for (int j = 1; j <= theSystem.nbStage; j++)
        {
            theSystem.CT[j]=new StochNum();
            theSystem.CT[j].distribution= "Deterministic";
            theSystem.CT[j].para1 = meanpt[j];
        }
        return theSystem;
    }

    SerialLine getTemplD(){

        SerialLine theSystem=new SerialLine();
        //save number of stages
        theSystem.nbStage = 14;
        double[] meanpt = new double[] {0, 25, 34, 33.5, 33.5, 22, 36, 22, 34, 26, 30, 26.5, 33, 29.5, 35};
        theSystem.buffer = new int[theSystem.nbStage];
        //save distribution information
        theSystem.CT = new StochNum[theSystem.nbStage+1];
        for (int j = 1; j <= theSystem.nbStage; j++)
        {
            theSystem.CT[j]=new StochNum();
            theSystem.CT[j].distribution= "Deterministic";
            theSystem.CT[j].para1 = meanpt[j];
        }
        return theSystem;
    }

    private void CalculateBetaPar(int j,int BNindex,int alphaindex, int noBNctindex,int varindex, SerialLine theSystem){
        double alpha = this.alphafactor[alphaindex] ;
        double beta = this.betafactor[alphaindex] ;
        double boolvariance = alpha*beta/((alpha+beta)*(alpha+beta)*(alpha+beta+1));
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
