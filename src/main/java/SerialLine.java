import java.io.InputStream;
import java.util.Scanner;

public class SerialLine {
    public int NbStage;
    public int[] Buffer;
    public StochNum[] CT;

    //Sim_Output
    public double OverallCT;
    public double[] Pblock;
    public double[] Pstarve;
    public double[] IPA_Sensitivity;
    public double[] Tup;
    public double[] Tdown;

    public SerialLine(InputStream system) {

        Scanner scanner=new Scanner (system);
        scanner.useDelimiter("\\s+");

        System.out.println();

        scanner.next();
        this.NbStage=scanner.nextInt();

        this.Buffer = new int[this.NbStage - 1];
        scanner.next();
        for (int j = 0; j < this.NbStage - 1; j++)
        {
            this.Buffer[j]=scanner.nextInt();
        }

        scanner.next();
        this.CT = new StochNum[this.NbStage];
        scanner.next();
        for (int j = 0; j < this.NbStage; j++)
        {
            this.CT[j]=new StochNum();
            this.CT[j].distribution=scanner.next(); // ???????????? this.CT[j].distribution 是字符串
        }
        scanner.next();
        for (int j = 0; j < this.NbStage; j++)
        {
            this.CT[j].Para1=scanner.nextDouble();
        }
        scanner.next();
        for (int j = 0; j < this.NbStage; j++)
        {
            this.CT[j].Para2=scanner.nextDouble();
        }


    }

    public void SIM_Serial_BAS(
            int N,
            int W,
            double[][] tij,
            int[][] uij,
            int[][] vij,
            int[][] wij,
            int[][] sij,
            double[][] bar_Sij,
            double[] bar_Dij,
            boolean SteadyState) {
        //*****************************************************************************
        //*****************************************************************************
        //** The notations of u,v,w,s is consistent with the journal paper ************
        //** on throughput improvement (2018).                             ************
        //*****************************************************************************
        //*****************************************************************************

        //***************************************
        //********** EVENT TIME *****************
        //***************************************
        double[][] Sij = new double[N][];
        double[][] Dij = new double[N][];
        for (int i = 0; i < N; i++)
        {
            Sij[i] = new double[NbStage];
            Dij[i] = new double[NbStage];
            for (int j = 0; j < NbStage; j++)
            {
                Sij[i][j] = 0;
                Dij[i][j] = 0;
            }
        }



        //***************************************
        //********** SIM ERG variables   ********
        //***************************************
        int[][] bsij = new int[N][];
        int[][] buij = new int[N][];
        int[][] bvij = new int[N][];
        int[][] bwij = new int[N][];


        for (int i = 0; i <N; i++)
        {
            bsij[i] = new int[this.NbStage];
            buij[i] = new int[this.NbStage];
            bvij[i] = new int[this.NbStage];
            bwij[i] = new int[this.NbStage];

            for (int j = 0; j < NbStage; j++)
            {
                bsij[i][j] = 0;
                buij[i][j] = 0;
                bvij[i][j] = 0;
                bwij[i][j] = 0;

                // initialize dual variables
                sij[i][j] = 0;
                uij[i][j] = 0;
                vij[i][j] = 0;
                wij[i][j] = 0;
            }
        }



        if (!SteadyState) {

            /////////////////////////////////////////////////////////////////////////////////////////////
            ///// Simulate from i=0 /////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////
            /////buij[i][j] = 1: D[i][j] trigerred by S[i][j]               /////////////////////////////
            /////bwij[i][j] = 1: D[i][j] trigerred by S[i-b_j][j+1] BLOCKAGE ////////////////////////////
            /////bvij[i][j] = 1: S[i][j] trigerred by D[i-1][j]             /////////////////////////////
            /////bsij[i][j] = 1: S[i][j] trigerred by D[i][j-1]   STARVATION ////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////


            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < NbStage; j++)
                {
                    ////////Starting event ////////////////////////////
                    //first arrival
                    if (i == 0 && j == 0)
                    {
                        Sij[i][j] = 0;
                        bsij[i][j] = 1;
                    }
                    //the first part facing the empty line
                    else if (i == 0 && j > 0)
                    {
                        Sij[i][j] = Dij[i][j - 1];
                        bsij[i][j] = 1;
                    }
                    //first stage: no starvation
                    else if (i > 0 && j == 0)
                    {
                        Sij[i][j] = Dij[i - 1][j];
                        bvij[i][j] = 1;
                    }
                    //others: max{D[i-1][j],D[i][j-1]}
                    else
                    {
                        if (Dij[i - 1][j] > Dij[i][j - 1])
                        {
                            Sij[i][j] = Dij[i - 1][j];
                            bvij[i][j] = 1;
                        }
                        else
                        {
                            Sij[i][j] = Dij[i][j - 1];
                            bsij[i][j] = 1;
                        }
                    }


                    //////// Departure event ////////////////////////////
                    // First parts: no blockage
                    if (j < NbStage - 1 && i < Buffer[j])
                    {
                        Dij[i][j] = Sij[i][j] + tij[i][j];
                        buij[i][j] = 1;
                    }
                    //last stage: no blockage
                    else if (j == NbStage - 1)
                    {
                        Dij[i][j] = Sij[i][j] + tij[i][j];
                        buij[i][j] = 1;
                    }
                    //others: max{S[i][j]+tij,S[i-b_j][j]}
                    else
                    {
                        if (Sij[i][j] + tij[i][j] > Sij[i - Buffer[j]][j + 1])
                        {
                            Dij[i][j] = Sij[i][j] + tij[i][j];
                            buij[i][j] = 1;
                        }
                        else
                        {
                            Dij[i][j] = Sij[i - Buffer[j]][j + 1];
                            bwij[i][j] = 1;
                        }
                    }
                }
            }

            /////////////////////////////////////////////////////////////////////////
            /////   Dual Optimal Solution               /////////////////////////////
            /////////////////////////////////////////////////////////////////////////
            for (int i = N - 1; i >= 0; i--)
            {
                for (int j = NbStage - 1; j >=0 ; j--)
                {
                    /////////Departure event////////////////////
                    if (i < N - 1 && j < NbStage - 1)
                    {
                        uij[i][j] = (vij[i + 1][j] + sij[i][j + 1])*buij[i][j];
                        wij[i][j] = (vij[i + 1][j] + sij[i][j + 1])*bwij[i][j];
                    }
                    else if (i == N - 1 && j < NbStage - 1)
                    {
                        uij[i][j] = (sij[i][j + 1])*buij[i][j];
                        wij[i][j] = (sij[i][j + 1])*bwij[i][j];
                    }
                    else if (i < N - 1 && j == NbStage - 1)
                    {
                        uij[i][j] = (vij[i + 1][j])*buij[i][j];
                    }
                    else if (i == N - 1 && j == NbStage - 1)
                    {
                        uij[i][j] = 1;
                    }

                    /////////Starting event///////////////////////
                    if (j > 0 && i < N - Buffer[j - 1])
                    {
                        sij[i][j] = (uij[i][j] + wij[i + Buffer[j - 1]][j - 1])*bsij[i][j];
                        vij[i][j] = (uij[i][j] + wij[i + Buffer[j - 1]][j - 1])*bvij[i][j];
                    }
                    else
                    {
                        sij[i][j] = (uij[i][j])*bsij[i][j];
                        vij[i][j] = (uij[i][j])*bvij[i][j];
                    }

                }
            }
            this.OverallCT = (Dij[N - 1][NbStage - 1] - Dij[W - 1][NbStage - 1]) / (double)(N - W);
            //***************   SAVE bar_Sij, bar_Dij    ***************
            for (int j = 0; j < NbStage; j++) {
                bar_Dij[j] = Dij[W - 1][j];
                if (j != NbStage -1) {
                    for (int i = 0 ; i < Buffer[j]; i++) {
                        bar_Sij[i][j+1] = Sij[W- Buffer[j]+i][j+1];
                    }
                }
            }


        }
        else {
            //***************   Fix event times in transient period   ***************
            for (int j = 0; j < NbStage; j++) {
                Dij[W - 1][j] = bar_Dij[j];
                if (j != NbStage - 1) {
                    for (int i = 0; i < Buffer[j]; i++) {
                        Sij[W - Buffer[j] + i][j + 1] = bar_Sij[i][j+1];
                    }
                }
            }

            /////////////////////////////////////////////////////////////////////////////////////////////
            /////bzij[i][j] = 1: D[i][j] trigerred by S[i][j]                ////////////////////////////
            /////bwij[i][j] = 1: D[i][j] trigerred by S[i-b_j][j+1] BLOCKAGE ////////////////////////////
            /////buij[i][j] = 1: S[i][j] trigerred by D[i-1][j]              ////////////////////////////
            /////bvij[i][j] = 1: S[i][j] trigerred by D[i][j-1]   STARVATION ////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////////
            for (int i = W; i < N; i++)
            {
                for (int j = 0; j < NbStage; j++)
                {
                    ////////Starting event ////////////////////////////

                    //first stage: no starvation
                    if (j == 0)
                    {
                        Sij[i][j] = Dij[i - 1][j];
                        bvij[i][j] = 1;
                    }
                    //others: max{D[i-1][j],D[i][j-1]}
                    else
                    {
                        if (Dij[i - 1][j] > Dij[i][j - 1])
                        {
                            Sij[i][j] = Dij[i - 1][j];
                            bvij[i][j] = 1;
                        }
                        else
                        {
                            Sij[i][j] = Dij[i][j - 1];
                            bsij[i][j] = 1;
                        }
                    }


                    ////////Departure event ////////////////////////////
                    //last stage: no blockage
                    if (j == NbStage - 1)
                    {
                        Dij[i][j] = Sij[i][j] + tij[i][j];
                        buij[i][j] = 1;
                    }
                    //others: max{S[i][j]+tij,S[i-b_j][j]}
                    else
                    {
                        if (Sij[i][j] + tij[i][j] > Sij[i - Buffer[j]][j + 1])
                        {
                            Dij[i][j] = Sij[i][j] + tij[i][j];
                            buij[i][j] = 1;
                        }
                        else
                        {
                            Dij[i][j] = Sij[i - Buffer[j]][j + 1];
                            bwij[i][j] = 1;
                        }
                    }
                }
            }

            for (int i = N - 1; i >= W; i--)
            {
                for (int j = NbStage - 1; j >= 0; j--)
                {
                    /////////Departure event////////////////////
                    if (i < N - 1 && j < NbStage - 1)
                    {
                        uij[i][j] = (vij[i + 1][j] + sij[i][j + 1])*buij[i][j];
                        wij[i][j] = (vij[i + 1][j] + sij[i][j + 1])*bwij[i][j];
                    }
                    else if (i == N - 1 && j < NbStage - 1)
                    {
                        uij[i][j] = (sij[i][j + 1])*buij[i][j];
                        wij[i][j] = (sij[i][j + 1])*bwij[i][j];
                    }
                    else if (i < N - 1 && j == NbStage - 1)
                    {
                        uij[i][j] = (vij[i + 1][j])*buij[i][j];
                    }
                    else if (i == N - 1 && j == NbStage - 1)
                    {
                        uij[i][j] = 1;
                    }

                    /////////Starting event///////////////////////
                    if (j > 0 && i < N - Buffer[j - 1])
                    {
                        sij[i][j] = (uij[i][j] + wij[i + Buffer[j - 1]][j - 1])*bsij[i][j];
                        vij[i][j] = (uij[i][j] + wij[i + Buffer[j - 1]][j - 1])*bvij[i][j];
                    }
                    else
                    {
                        sij[i][j] = (uij[i][j])*bsij[i][j];
                        vij[i][j] = (uij[i][j])*bvij[i][j];
                    }

                }
            }
            this.OverallCT = (Dij[N-1][NbStage-1]- Dij[W-1][NbStage-1])/ (double)(N-W);
        }


    }

}
