package test_mainfunc.optimization;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloRange;
import test_mainfunc.simulation.SerialLine;
import test_mainfunc.simulation.StochNum;
import test_mainfunc.util.Stopwatch;

import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.HashMap;

import static java.lang.Math.*;

public class BendersIntModelAlter6ReversedCut extends BendersIntModelAlter6 {

    public BendersIntModelAlter6 myReversedBAPModel;
    private double constantC = 1;

    public BendersIntModelAlter6ReversedCut(SerialLine system, double THstar, int[] lB, int[] uB, int N, int W) {
        super(system, THstar, lB, uB, N, W);
        int[] reversedLB = new int[lB.length];
        int[] reversedUB = new int[uB.length];
        for (int j = 1; j < lB.length; j++) {
            reversedLB[j] = lB[lB.length - j];
            reversedUB[j] = uB[uB.length - j];
        }
        this.myReversedBAPModel = new BendersIntModelAlter6(system, THstar, reversedLB, reversedUB, N, W);
    }

    public BendersIntModelAlter6ReversedCut(SerialLine system, double THstar, int[] lB, int[] uB, int N, int W, int startMachine, int endMachine, HashMap<String, Integer> initialBounds) {
        super(system, THstar, lB, uB, N, W, startMachine, endMachine, initialBounds);
        int[] reversedLB = new int[this.mySystem.nbStage];
        int[] reversedUB = new int[this.mySystem.nbStage];
        for (int j = 1; j < this.mySystem.nbStage; j++) {
            reversedLB[j] = this.lowerBoundj[this.mySystem.nbStage - j];
            reversedUB[j] = this.upperBoundj[this.mySystem.nbStage - j];
        }
        this.myReversedBAPModel = new BendersIntModelAlter6(this.mySystem, THstar, reversedLB, reversedUB, N, W);
    }


    public boolean solveBAPWithIntModel(double[][] tij
            , boolean redefineBounds
    ) throws IloException {
        this.writer.println("Alter 6 with reversed cuts:");

        if (redefineBounds) {
            // define lowerbound
            for (int j = 1; j <= mySystem.nbStage - 1; j++) {
                int[] tempUb = new int[mySystem.nbStage];
                int[] tempLb = new int[mySystem.nbStage];
                for (int j1 = 1; j1 <= mySystem.nbStage - 1; j1++) {
                    if (j1 != j) {
                        tempLb[j1] = this.upperBoundj[j1];
                        tempUb[j1] = this.upperBoundj[j1];
                    } else {
                        tempLb[j1] = this.lowerBoundj[j1];
                        tempUb[j1] = this.upperBoundj[j1];
                    }
                }
                BendersIntModelAlter6ReversedCut upboundSearchSystem = new BendersIntModelAlter6ReversedCut(mySystem, this.THstar, tempLb, tempUb, this.simulationLength, this.warmupLength);
                upboundSearchSystem.solveBAPWithIntModel(tij, false);
                this.lowerBoundj[j] = upboundSearchSystem.mySystem.buffer[j];
                this.myReversedBAPModel.lowerBoundj[mySystem.nbStage - j] = upboundSearchSystem.mySystem.buffer[j];
            }

            for(int j=1;j<mySystem.nbStage;j++){

                IloLinearNumExpr initialLb = cplex.linearNumExpr();
                IloRange rng;
                initialLb.addTerm(1, this.bj[j]);
                rng = cplex.addGe(initialLb, this.lowerBoundj[j]);
                rng.setName("iniitla lower bound of " + (j) + ": ");
            }

        }

        double[][] reversedTij = new double[this.simulationLength + 1][this.mySystem.nbStage + 1];
        for (int i = 1; i <= this.simulationLength; i++) {
            for (int j = 1; j <= this.mySystem.nbStage; j++) {
                reversedTij[i][j] = tij[this.simulationLength + 1 - i][this.mySystem.nbStage + 1 - j];
            }
        }
        myReversedBAPModel.getMmValue(reversedTij);

        return super.solveBAPWithIntModel(tij,false);
    }


    public boolean solveBAPWithIntModel(double[][] tij, HashMap<String, Integer> initialLowerBounds) throws IloException {
        this.writer.println("Alter 6 with reversed cuts:");

        double[][] reversedTij = new double[this.simulationLength + 1][this.mySystem.nbStage + 1];
        for (int i = 1; i <= this.simulationLength; i++) {
            for (int j = 1; j <= this.mySystem.nbStage; j++) {
                reversedTij[i][j] = tij[this.simulationLength + 1 - i][this.mySystem.nbStage + 1 - j];
            }
        }
        HashMap<String, Integer> reversedLowerBounds = new HashMap<>();
        int m = tij[1].length - 1;
        for (int j = 1; j <= m; j++) {
            for (int j1 = j + 1; j1 <= m; j1++) {
                if (initialLowerBounds.containsKey(j + "_" + j1)) {
                    reversedLowerBounds.put((m + 1 - j1) + "_" + (m + 1 - j), initialLowerBounds.get(j + "_" + j1));
                }
            }
        }
        myReversedBAPModel.getMmValue(reversedTij, reversedLowerBounds);
        return super.solveBAPWithIntModel(tij, initialLowerBounds);
    }

    @Override
    public void solveWithDecomposition(double[][] tij) throws IloException {
        //System.out.println("The method solveWithDecomposition(double[][] tij) of class BendersIntModelAlter6ReversedCut is not defined yet!");

        /*double[][] reversedTij = new double[this.simulationLength + 1][this.mySystem.nbStage + 1];
        for (int i = 1; i <= this.simulationLength; i++) {
            for (int j1 = 1; j1 <= this.mySystem.nbStage; j1++) {
                reversedTij[i][j1] = tij[this.simulationLength + 1 - i][this.mySystem.nbStage + 1 - j1];
            }
        }
        //this.getMmValue(tij,0);
        //myReversedBAPModel.getMmValue(reversedTij,0);
        this.mySystem.mySimulation = this.mySystem.new SimulationBAS(this.simulationLength,this.warmupLength,tij);

        // ======================================== Decomposition with cuts ======================================
        int stepSize =  max(mySystem.nbStage/5,1);
        for (int m = max(stepSize,2) ; m <= mySystem.nbStage; m+=stepSize) {
            for (int startMachine = 1; startMachine <= mySystem.nbStage - 1; startMachine+=2) {
                int endMachine = startMachine + m - 1;
                if (endMachine <= mySystem.nbStage) {
                    for (int j = 1; j <= endMachine - 1; j++) {
                        if (j < startMachine || j >= endMachine) {
                            bj[j].setLB(upperBoundj[j]);
                        }
                    }

                    this.getMmValue(tij, 0);
                    myReversedBAPModel.getMmValue(reversedTij, 0);

                    this.solveMasterProb();
                    this.updateBufferSpaces();
                    this.mySystem.mySimulation.simDualBAS(false, 0);
                    this.saveIterationSolution();

                    solvability = true;
                    while ((this.THstar - this.mySystem.TH > 0) && (numit < this.MAX_ITE) && solvability) {
                        numit++;
                        this.generateFeasibilityCut(tij, 0);
                        this.addFeasibilityCut(startMachine, endMachine);
                        this.solveMasterProb();
                        this.updateBufferSpaces();
                        this.saveIterationSolution();
                        this.mySystem.mySimulation.simDualBAS(false, 0);
                    }

                    int totBufferSubsystem = 0;
                    IloLinearNumExpr constBound = cplex.linearNumExpr();
                    for(int j=startMachine;j<=endMachine-1;j++){
                        constBound.addTerm(1,bj[j]);
                        totBufferSubsystem += this.mySystem.buffer[j];
                    }
                    cplex.addGe(constBound,totBufferSubsystem,"bound_"+startMachine+"_"+endMachine);

                    for (int j = 1; j <= endMachine - 1; j++) {
                        if (j < startMachine || j >= endMachine) {
                            bj[j].setLB(lowerBoundj[j]);
                        }
                    }
                }
            }
            if(m<mySystem.nbStage && m>-stepSize+mySystem.nbStage)
                stepSize=1;
        }*/


        // ======================================== Decomposition with bounds ====================================
        for(int j=1;j<=mySystem.nbStage-1;j++){
            mySystem.buffer[j] = lowerBoundj[j];
        }
        HashMap<String, Integer> initialBounds = new HashMap<>();
        numit=0;
        double remainTime = MAX_TIME;
        boolean stopLoop = false;
        for(int m=2;m<=mySystem.nbStage&& !stopLoop;m++){
            for(int startMachine = 1; startMachine<=mySystem.nbStage-1;startMachine++){
                int endMachine = startMachine + m-1;
                if(endMachine<=mySystem.nbStage){
                    Stopwatch sw = new Stopwatch();
                    sw.start();

                    int[] tempUb = new int[mySystem.nbStage];
                    int[] tempLb = new int[mySystem.nbStage];
                    for (int j1 = 1; j1 <= mySystem.nbStage - 1; j1++) {
                        if (j1 < startMachine || j1>=endMachine) {
                            tempLb[j1] = this.upperBoundj[j1];
                            tempUb[j1] = this.upperBoundj[j1];
                        } else {
                            tempLb[j1] = this.lowerBoundj[j1];
                            tempUb[j1] = this.upperBoundj[j1];
                        }
                    }


                    HashMap<String, Integer> subInitialBounds = new HashMap<>();
                    for(int j=startMachine;j<=endMachine-1;j++){
                        for(int j0=j+1;j0<=endMachine;j0++){
                            if(initialBounds.containsKey(j+"_"+j0)){
                                subInitialBounds.put(j+"_"+j0,initialBounds.get(j+"_"+j0));
                            }
                        }
                    }
                    BendersIntModelAlter6ReversedCut subSystemBap = new BendersIntModelAlter6ReversedCut(mySystem,
                            this.THstar, tempLb, tempUb, this.simulationLength, this.warmupLength,1,mySystem.nbStage,subInitialBounds);
                    /*BendersIntModelAlter6ReversedCut subSystemBap = new BendersIntModelAlter6ReversedCut(mySystem, THstar, lowerBoundj, upperBoundj, simulationLength, warmupLength,
                            startMachine, endMachine, initialBounds);*/
                    subSystemBap.setMaxTime(remainTime);
                    subSystemBap.solveBAPWithIntModel(tij,subInitialBounds);

                    if(subSystemBap.solvability && m==mySystem.nbStage){
                        for(int j=1;j<=mySystem.nbStage-1;j++){
                            mySystem.buffer[j] = subSystemBap.mySystem.buffer[j];
                        }
                        stopLoop = true;
                        break;
                    }else if(subSystemBap.solvability && m<mySystem.nbStage){
                        int solSubProb = 0;
                        for(int j=startMachine;j<=endMachine-1;j++){
                            solSubProb+=subSystemBap.mySystem.buffer[j];
                        }
                        addInitialBounds(initialBounds,startMachine,endMachine,solSubProb);
                    }
                    numit+=subSystemBap.numit;
                    sw.stop();
                    remainTime -= sw.elapseTimeSeconds;

                    if(remainTime<0 ){
                        if(m <mySystem.nbStage)
                            subSystemBap = new BendersIntModelAlter6ReversedCut(mySystem, THstar, lowerBoundj, upperBoundj, simulationLength, warmupLength,
                                1, mySystem.nbStage, initialBounds);
                        subSystemBap.setMaxIte(0);
                        subSystemBap.numit=0;
                        subSystemBap.setMaxTime(MAX_TIME);
                        subSystemBap.solveBAPWithIntModel(tij,false);
                        if (mySystem.nbStage - 1 >= 0)
                            System.arraycopy(subSystemBap.mySystem.buffer, 1, mySystem.buffer, 1, mySystem.nbStage - 1);
                        stopLoop = true;
                        break;
                    }
                }
            }
        }

    }

    @Override
    public void addFeasibilityCut() {
        try {

            // Calculate the coefficient of each term in the reversed cut
            double reversed_constantTerm = this.newCut.constantTerm;
            double[][] reversed_mcoefficient = new double[mySystem.nbStage][];
            for (int j = 1; j <= this.mySystem.nbStage - 1; j++) {
                reversed_mcoefficient[j] = new double[this.upperBoundj[j] + 1];
                for (int k = this.lowerBoundj[j] + 1; k <= this.upperBoundj[j]; k++)
                    reversed_mcoefficient[j][k] = 0;
            }

            for (int j = 1; j <= this.mySystem.nbStage - 1; j++) {
                //if(j>=startMachine && j<endMachine){
                for (int i = 1; i <= this.simulationLength; i++) {
                    if (this.mySystem.mySimulation.wij[i][j] > 0) {
                        for (int k = this.lowerBoundj[j] + 1; k <= this.mySystem.buffer[j]; k++) {
                            reversed_mcoefficient[j][k] += this.myReversedBAPModel.mijk[this.simulationLength + this.mySystem.buffer[j] + 1 - i][this.mySystem.nbStage - j][k] * this.mySystem.mySimulation.wij[i][j];
                        }

                        for (int k = this.mySystem.buffer[j] + 1; k <= this.upperBoundj[j]; k++) {
                            reversed_mcoefficient[j][k] += this.myReversedBAPModel.Mijk[this.simulationLength + this.mySystem.buffer[j] + 1 - i][this.mySystem.nbStage - j][k] * this.mySystem.mySimulation.wij[i][j];
                        }
                    }
                }

                for (int k = this.lowerBoundj[j] + 1; k <= this.mySystem.buffer[j]; k++) {
                    reversed_constantTerm += reversed_mcoefficient[j][k];
                }
            }

            for(int j=1;j<mySystem.nbStage;j++){
                double cumulateMjk = 0;
                for (int k = this.lowerBoundj[j] + 1; k <= this.upperBoundj[j]; k++) {
                    if (cumulateMjk > reversed_constantTerm)
                        reversed_mcoefficient[j][k] = 0;
                    else
                        cumulateMjk += reversed_mcoefficient[j][k];
                }
            }
            //}

            // Calculate the coefficient of each term in the standard cut
            double constantTerm = this.newCut.constantTerm;
            double[][] mcoefficient = new double[mySystem.nbStage][];
            for (int j = 1; j <= this.mySystem.nbStage - 1; j++) {
                mcoefficient[j] = new double[this.upperBoundj[j] + 1];
                for (int k = this.lowerBoundj[j] + 1; k <= this.upperBoundj[j]; k++)
                    mcoefficient[j][k] = 0;
            }

            for (int j = 1; j <= this.mySystem.nbStage - 1; j++) {
                //if(j>=startMachine && j<endMachine){
                for (int i = 1; i <= this.simulationLength; i++) {
                    if (this.mySystem.mySimulation.wij[i][j] > 0) {
                        for (int k = this.lowerBoundj[j] + 1; k <= this.mySystem.buffer[j]; k++) {
                            mcoefficient[j][k] += this.mijk[i][j][k] * this.mySystem.mySimulation.wij[i][j];
                        }

                        for (int k = this.mySystem.buffer[j] + 1; k <= this.upperBoundj[j]; k++) {
                            mcoefficient[j][k] += this.Mijk[i][j][k] * this.mySystem.mySimulation.wij[i][j];
                        }
                    }
                }

                for (int k = this.lowerBoundj[j] + 1; k <= this.mySystem.buffer[j]; k++) {
                    constantTerm += mcoefficient[j][k];
                }
            }
            //}

            for(int j=1;j<mySystem.nbStage;j++){
                double cumulateMjk = 0;
                for (int k = this.lowerBoundj[j] + 1; k <= this.upperBoundj[j]; k++) {
                    if (cumulateMjk > constantTerm)
                        mcoefficient[j][k] = 0;
                    else
                        cumulateMjk += mcoefficient[j][k];
                }
            }

            // Is combinatorial cut or standard cut tighter?
            boolean combCutIsTighterThanStandardCut = true;
            for (int j = 1; j <= this.mySystem.nbStage - 1; j++) {
                double lhs = 0;
                for (int k = this.lowerBoundj[j] + 1; k <= min(this.mySystem.buffer[j] + 1, this.upperBoundj[j]); k++) {
                    lhs += mcoefficient[j][k];
                }
                if (lhs < constantTerm)
                    combCutIsTighterThanStandardCut = false;
            }

            // Is combinatorial cut or reversed cut tighter?
            boolean combCutIsTighterThanReversedCut = true;
            for (int j = 1; j <= this.mySystem.nbStage - 1; j++) {
                double lhs = 0;
                for (int k = this.lowerBoundj[j] + 1; k <= min(this.mySystem.buffer[j] + 1, this.upperBoundj[j]); k++) {
                    lhs += reversed_mcoefficient[j][k];
                }
                if (lhs < constantTerm)
                    combCutIsTighterThanReversedCut = false;
            }

            // Which cut(s) to be added?
            boolean addCombCut = false;
            boolean addStandardCut = false;
            boolean addReversedCut = false;
            if (combCutIsTighterThanStandardCut && combCutIsTighterThanReversedCut)
                addCombCut = true;
            else if (combCutIsTighterThanStandardCut)
                addReversedCut = true;
            else if (combCutIsTighterThanReversedCut)
                addStandardCut = true;
            else {
                boolean similarCuts = true;
                for (int j = 1; j <= this.mySystem.nbStage - 1; j++) {
                    if (!similarCuts)
                        break;
                    for (int k = this.lowerBoundj[j] + 1; k <= min(this.mySystem.buffer[j] + 1, this.upperBoundj[j]); k++) {
                        if (abs(mcoefficient[j][k] - reversed_mcoefficient[j][k]) >= this.newCut.constantTerm * this.constantC / (this.mySystem.nbStage - 1)) {
                            similarCuts = false;
                            break;
                        }
                    }
                }
                if (similarCuts) {
                    if (random() > 0.5)
                        addStandardCut = true;
                    else
                        addReversedCut = true;
                } else {
                    addStandardCut = true;
                    addReversedCut = true;
                }

            }

            // Add the correct cut
            if (addCombCut) {
                IloLinearNumExpr feacut = cplex.linearNumExpr();
                IloRange rng;
                for (int j = 1; j <= this.mySystem.nbStage - 1; j++) {
                    if (this.mySystem.buffer[j] < this.upperBoundj[j])
                        feacut.addTerm(1, this.yjk[j][this.mySystem.buffer[j] + 1]);
                }
                rng = cplex.addGe(feacut, 1.0);
                rng.setName("combcut of iter " + (numit) + ": ");
                writer.println(rng);
            }
            if (addStandardCut) {
                IloLinearNumExpr feacut = cplex.linearNumExpr();
                IloRange rng;
                for (int j = 1; j <= this.mySystem.nbStage - 1; j++) {
                    for (int k = this.lowerBoundj[j] + 1; k <= this.upperBoundj[j]; k++)
                        feacut.addTerm(-mcoefficient[j][k], this.yjk[j][k]);
                }
                feacut.setConstant(constantTerm);
                rng = cplex.addLe(feacut, 0);
                rng.setName("standard feascut of iter " + (numit) + ": ");
                writer.println(rng);
            }
            if (addReversedCut) {
                IloLinearNumExpr feacut = cplex.linearNumExpr();
                IloRange rng;
                for (int j = 1; j <= this.mySystem.nbStage - 1; j++) {
                    for (int k = this.lowerBoundj[j] + 1; k <= this.upperBoundj[j]; k++)
                        feacut.addTerm(-reversed_mcoefficient[j][k], this.yjk[j][k]);
                }
                feacut.setConstant(reversed_constantTerm);
                rng = cplex.addLe(feacut, 0);
                rng.setName("reversed feascut of iter " + (numit) + ": ");
                writer.println(rng);
            }
        } catch (Exception exc) {
            exc.printStackTrace();
        }


    }


}
