package test_mainfunc;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.dataformat.yaml.YAMLMapper;
import test_mainfunc.input.DoeInputV1;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.simulation.StochNum;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

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
    int difffailedstage;
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
