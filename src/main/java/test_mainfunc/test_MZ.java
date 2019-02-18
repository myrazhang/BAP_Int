package test_mainfunc;

import java.io.*;

public class test_MZ {
    public static void main(String argv[]){

        int m=5;
        int N=100;
        int[] U=new int[m-1];
        int[] L=new int[m-1];
        for (int j=1;j<=m-1;j++){
            U[j-1]=5;
            L[j-1]=2;
        }
        mM_value A=new mM_value(m, N,U,L);
        double[][] t=new double[N][m];
        for(int i=1;i<=N;i++){
            t[i-1][0]=1;
            t[i-1][1]=2;
            t[i-1][2]=3;
            t[i-1][3]=4;
            t[i-1][4]=1;
        }

        A.get_mM_value(t);

        String programPath = System.getProperty("user.dir");
        String Mijk_File = programPath + "\\OUTPUT\\Mijk.txt";
        OutputStream outResM= null;
        try {
            outResM = new FileOutputStream(Mijk_File);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(-1);
        }
        PrintWriter file_M=new PrintWriter(outResM, true);


        String mijk_File = programPath + "\\OUTPUT\\mmijk.txt";
        OutputStream outResm= null;
        try {
            outResm = new FileOutputStream(mijk_File);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(-1);
        }
        PrintWriter file_m=new PrintWriter(outResm, true);



        for(int i=1;i<=N;i++){
            for(int j=1;j<=m-1;j++){
                for(int k=1;k<=U[j-1];k++){
                    file_m.print(A.mijk[i-1][j-1][k-1]);
                    file_m.print(", ");
                    file_M.print(A.Mijk[i-1][j-1][k-1]);
                    file_M.print(", ");
                }
                file_m.print("; ");
                file_M.print("; ");
            }
            file_m.println();
            file_M.println();
        }

        try {
            outResm.close();
            outResM.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
