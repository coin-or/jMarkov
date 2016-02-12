/**
 * MtjLogRedSolver.java
 * Created: 2/09/2006
 */
package jmarkov.solvers;

import jmarkov.GeomProcess;
import jmarkov.basic.exceptions.NotUnichainException;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;

/**
 * This class implements the Logarithmic reduction algorithm to 
 * find the steady state probabilities of a QBD process
 * @author German Riano. Universidad de los Andes. (C) 2006
 */
public class MtjLogRedSolver extends GeometricSolver {

    private Matrix R = null;
    private GeomProcess mp = null;

    /**
     * @param mp The given QBD
     */
    public MtjLogRedSolver(GeomProcess mp) {
        super(mp);
        this.mp = mp;
    }

    @Override
    public double[][] getRmatrix() throws NotUnichainException {

    	long startTimeR = System.currentTimeMillis();
    	
        double epsilon = 1e-8;
        Matrix A[] = mp.getAMatrices();
        Matrix MA0 = A[0], MA1 = A[1], MA2 = A[2];
        

        int dimen = MA1.numRows();

        Matrix I = new DenseMatrix(Matrices.identity(dimen));

        DenseMatrix mA1 = new DenseMatrix(dimen, dimen);
        DenseMatrix mA1I = new DenseMatrix(dimen, dimen);
        DenseMatrix H = new DenseMatrix(dimen, dimen);
        DenseMatrix L = new DenseMatrix(dimen, dimen);
        
        
        long startTimeS0 = System.currentTimeMillis();
        //MA1.mult(-1, I, mA1);
        mA1.set(MA1);
        mA1.scale(-1);
        long startTimeSolve = System.currentTimeMillis();
        mA1.solve(I, mA1I);
        
        mA1.solve(MA0, H);
        mA1.solve(MA2, L);
        
        long stopTimeS0 = System.currentTimeMillis();
        long elapsedTimeS0 = stopTimeS0 - startTimeS0;
        long elapsedTimeSolve = stopTimeS0 - startTimeSolve;
        System.out.println("Time init mult and solve: "+elapsedTimeS0+" ms");
        
        
        DenseMatrix G = new DenseMatrix(dimen, dimen);
        DenseMatrix T = new DenseMatrix(dimen, dimen);
        DenseMatrix U = new DenseMatrix(dimen, dimen);
        DenseMatrix UA = new DenseMatrix(dimen, dimen);
        DenseMatrix M = new DenseMatrix(dimen, dimen);
        DenseMatrix TA = new DenseMatrix(dimen, dimen);
        
        DenseMatrix one = new DenseMatrix(G.numRows(), 1);
        
        long startTimeE = System.currentTimeMillis();
        for(MatrixEntry e: one){
        	e.set(1.0);
        }
        long stopTimeE = System.currentTimeMillis();
        long elapsedTimeE = stopTimeE - startTimeE;
        System.out.println("Time creating unit vector e: "+elapsedTimeE+" ms");
        
        DenseMatrix check = new DenseMatrix(G.numRows(), 1);
        
        R = new DenseMatrix(dimen, dimen);
        
        /*
        long startTimeInitMult = System.currentTimeMillis();
        mA1I.mult(MA0, H);// H=(-A_1)^{-1}A_0
        mA1I.mult(MA2, L);// L=(-A_1)^{-1}A_2
        long stopTimeInitMult = System.currentTimeMillis();
        long elapsedTimeInitMult = stopTimeInitMult - startTimeInitMult;
        System.out.println("\nTime initial mults: "+elapsedTimeInitMult+" ms\n");
        */
        
        /*
        long startTimeLHinit = System.currentTimeMillis();
        DenseMatrix LHtemp = new DenseMatrix(dimen, dimen);
        DenseMatrix H2 = new DenseMatrix(dimen, dimen);
        DenseMatrix L2 = new DenseMatrix(dimen, dimen);
        
        LHtemp.set(MA1);
        LHtemp.scale(-1);
        LHtemp.solve(MA0, H2);
        //L2.set(MA1);
        //L2.scale(-1);
        LHtemp.solve(MA2, L2);
        long stopTimeLHinit = System.currentTimeMillis();
        long elapsedTimeLHinit = stopTimeLHinit - startTimeLHinit;
        System.out.println("\nTime NEW initial mults: "+elapsedTimeLHinit+" ms\n");
        double err = (H2.add(-1, H)).norm(Matrix.Norm.Maxvalue);
        System.out.println("\nerror new H: "+err+"\n");
        err = (L2.add(-1, L)).norm(Matrix.Norm.Maxvalue);
        System.out.println("\nerror new L: "+err+"\n");*/
        
        G.set(L);
        T.set(H);
        double compare = 1;
        long stopTimeSolve = 0;
        long stopTimeInit = System.currentTimeMillis();
        long elapsedTimeInit = stopTimeInit - startTimeR;
        System.out.println("Time comp R initial section: "+elapsedTimeInit+" ms");
        long elapsedTime = 0;
        long elapsedTime2 = 0;
        while (compare > epsilon) {

        	long startTime = System.currentTimeMillis();
	    	
	    	
            H.mult(L, U);// U=HL
            L.multAdd(H, U);// U=HL+LH
            H.mult(H, M);// M=(H)^2
            H = new DenseMatrix(dimen, dimen);
            //H.zero();
            mA1I = new DenseMatrix(dimen, dimen);
            //mA1I.zero();
            // U.multAdd(-1, I, I, H);// H=(I-U)
            // TODO: check this!!
            //H = (DenseMatrix) I.copy().add(-1, U);// H=(I-U)
            H.set(U);
            H.scale(-1);
            H.add(I);
            
            
            startTimeSolve = System.currentTimeMillis();
            
            long solveTimeS = System.currentTimeMillis();
            H.solve(I, mA1I);// H1I=(I-U)^{-1}
            long solveTimeF = System.currentTimeMillis() - solveTimeS;
            System.out.println("Solve time: "+solveTimeF+" ms");
            
            stopTimeSolve = System.currentTimeMillis();
            elapsedTimeSolve += stopTimeSolve - startTimeSolve;
            //H.zero();
            H = new DenseMatrix(dimen, dimen);
            
            
            solveTimeS = System.currentTimeMillis();
            mA1I.mult(M, H);// H=(I-U)^{-1}*M
            solveTimeF = System.currentTimeMillis() - solveTimeS;
            System.out.println("Mult time: "+solveTimeF+" ms");
            
            //M.zero();
            M = new DenseMatrix(dimen, dimen);
            L.mult(L, M);// M=(L)^2
            
            //L.zero();
            L = new DenseMatrix(dimen, dimen);
            mA1I.mult(M, L);// L=(I-U)^{-1}*M
            T.multAdd(L, G);// G=G+TL
            TA.set(T);
            //T.zero();
            T = new DenseMatrix(dimen, dimen);
            TA.mult(H, T);// T=TH
            //U.zero();
            U = new DenseMatrix(dimen, dimen);

            
            long stopTime = System.currentTimeMillis();
	        elapsedTime += stopTime - startTime;
	        
	        
            // Procedimiento para el calculo de la norma infinito
            // ||1-G1||_{\infty}
	        long startTime2 = System.currentTimeMillis();
	        
	        check.zero();
	        G.mult(one, check);
            //check.add(one);
	        
            compare = Math.abs(1 - check.get(0, 0));

            for (int i = 1; i < check.numRows(); i++) {
                compare = Math.max(compare, Math.abs(1 - check.get(i, 0)));
            }
            
            long stopTime2 = System.currentTimeMillis();
	        elapsedTime2 += stopTime2 - startTime2;

        }
        
        System.out.println("Time matrix multplications: "+elapsedTime+" ms");
        System.out.println("Time norm computation: "+elapsedTime2+" ms");

        long startTimeEnd = System.currentTimeMillis();
        
        U.zero();
        //U = new DenseMatrix(dimen, dimen);
        // MA0.multAdd(G, MA1, U);
        // TODO: check this!!!
        //U = (DenseMatrix) MA0.multAdd(G, MA1.copy());
        
        U.set(MA0);
        U.multAdd(G, MA1.copy());
        
        
        /*        
        mA1I.zero();
        //mA1I = new DenseMatrix(dimen, dimen);
        //U.mult(-1, I, UA);// U=A_1+A_0*G
        UA.set(U); // U=A_1+A_0*G

        UA.scale(-1); // U=A_1+A_0*G
        startTimeSolve = System.currentTimeMillis();
        UA.solve(I, mA1I);// (I-U)^{-1}
        stopTimeSolve = System.currentTimeMillis();
        elapsedTimeSolve += stopTimeSolve - startTimeSolve;
        System.out.println("Time comp R solve: "+elapsedTimeSolve+" ms");
        MA0.mult(mA1I, R);// R=A_0*H1I.
        
                */
        
        /*long startTimeR2 = System.currentTimeMillis();
        DenseMatrix R2 = new DenseMatrix(dimen, dimen);*/
        //U.transpose();
        U.scale(-1);
        U.transSolve(MA0.transpose(), R);
        R.transpose(); 
        /*
        long stopTimeR2 = System.currentTimeMillis();
        long elapsedTimeR2 = stopTimeR2 - startTimeR2;
        System.out.println("\nTime NEW final section: "+elapsedTimeR2+" ms\n");
        double err = (R2.add(-1, R)).norm(Matrix.Norm.Maxvalue);
        System.out.println("\nerror new R: "+err+"\n");*/
        
        long stopTimeEnd = System.currentTimeMillis();
        long elapsedTimeEnd = stopTimeEnd - startTimeEnd;
        System.out.println("Time comp R final section: "+elapsedTimeEnd+" ms");

        
        long stopTimeR = System.currentTimeMillis();
        long elapsedTimeR = stopTimeR - startTimeR;
        System.out.println("Time comp R inside solver: "+elapsedTimeR+" ms");

        
        

        return Matrices.getArray(R);
    }

    /**
     * @see jmarkov.solvers.Solver#label()
     */
    @Override
    public String label() {
        return "MTJ Logarithmic reduction solver";
    }

    public String description() {
        return "MTJ Logarithmic reduction solver. This solver uses the MTJ p"
                + "ackage to handle matrices.";
    }

}
