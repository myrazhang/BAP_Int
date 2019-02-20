package test_mainfunc;
import static java.lang.Math.min;
import static java.lang.Math.max;

public class mM_value {

    public double[][][] mijk;
    public double[][][] Mijk;

    public int nbstage;
    public int N;
    public int[] Uj;
    public int[] Lj;


//    public mM_value(int m, int nbparts, int[] UB, int[] LB){
    public mM_value(int m, int nbparts, int UB, int LB){
        this.nbstage=m;
        this.N=nbparts;
        this.Uj= new int[m-1];
        this.Lj= new int[m-1];

        for(int j=1;j<=this.nbstage-1;j++){
            //this.Uj[j-1]=UB[j-1];
            //this.Lj[j-1]=LB[j-1];
            this.Uj[j-1]=UB;
            this.Lj[j-1]=LB;
        }

        this.Mijk=new double[this.N][this.nbstage-1][];
        this.mijk=new double[this.N][this.nbstage-1][];

        for(int i=1;i<=this.N;i++){
            for(int j=1;j<=this.nbstage-1;j++){
                this.Mijk[i-1][j-1]=new double[this.Uj[j-1]];
                this.mijk[i-1][j-1]=new double[this.Uj[j-1]];
            }
        }

        for (int j = 1; j <= this.nbstage-1; j++){
            for(int k = 1; k <= this.Uj[j-1]; k++){
                for (int i = 1; i <= this.N; i++){
                    this.mijk[i-1][j-1][k-1]=0;
                    this.Mijk[i-1][j-1][k-1]=100000;
                }
            }
        }
    }

    public void get_mM_value(
            double[][] tij
    )
    {
        for (int j = 1; j <= this.nbstage-1; j++){
            for(int k = this.Lj[j-1]; k <= this.Uj[j-1]; k++){
                for (int i = k+1; i <= this.N; i++){
                    this.mijk[i-1][j-1][k-1]=100000;
                    this.Mijk[i-1][j-1][k-1]=0;
                }
            }
        }

        // mijk
        for(int j = 1;j <= this.nbstage-1; j++){
            for(int k = this.Lj[j-1]; k <= this.Uj[j-1]; k++){
                for (int i = k+1; i <= this.N; i++) {
                    for(int i0 = i - k + 1; i0 <= i-this.Lj[j-1] ; i0++){
                        this.mijk[i-1][j-1][k-1] = min( this.mijk[i-1][j-1][k-1], tij[i0-1][j]);
                    }
                }
            }
        }



        // Mijk
        int[][] U_j_deltaj=new int[this.nbstage-1][this.nbstage-1];
        int[][] L_j_deltaj=new int[this.nbstage-1][this.nbstage-1];

        for (int j = 1; j <= this.nbstage-1; j++){
            for(int delta_j = 1; delta_j <= this.nbstage-j; delta_j++){

                U_j_deltaj[j-1][delta_j-1]=delta_j+1;
                L_j_deltaj[j-1][delta_j-1]=delta_j+1;

                for (int l=j; l <= j+delta_j-1 ; l++){
                    U_j_deltaj[j-1][delta_j-1]=U_j_deltaj[j-1][delta_j-1]+this.Uj[l-1];
                    L_j_deltaj[j-1][delta_j-1]=L_j_deltaj[j-1][delta_j-1]+this.Lj[l-1];
                }
            }
        }


        double[][] sij= new double [this.N][this.nbstage];
        for (int i = 2; i <= this.N; i++){
            for (int j  = 1; j <= this.nbstage-1; j++){
                sij[i-1][j-1]=tij[i-2][j-1];
                for (int delta_j = 1; delta_j <= this.nbstage-j ; delta_j++){
                    for (int i0= i - U_j_deltaj[j-1][delta_j-1];i0<= i - L_j_deltaj[j-1][delta_j-1];i0++){
                        sij[i-1][j-1]=max(sij[i-1][j-1], tij[i0-1][j+delta_j-1]);
                    }
                }
            }
            sij[i-1][this.nbstage-1]=tij[i-2][this.nbstage-1];
        }


        for(int j=1;j <= this.nbstage-1;j++){
            for(int k=this.Lj[j-1]; k <= this.Uj[j-1];k++){
                for (int i = k+1 ; i <= this.N; i++){
                    for (int i0=i-this.Uj[j]+1;i0<= i-k;i0++){
                        this.Mijk[i-1][j-1][k-1]=max(this.Mijk[i-1][j-1][k-1],sij[i0-1][j]);
                    }
                }
            }
        }

    }


}
