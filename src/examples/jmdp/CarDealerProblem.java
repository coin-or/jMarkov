package examples.jmdp;

import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.DTMDP;
import jmarkov.jmdp.solvers.PolicyIterationSolver;
import jmarkov.jmdp.solvers.RelativeValueIterationSolver;
import jmarkov.jmdp.solvers.ValueIterationSolver;
import Jama.Matrix;

/**
 * This example is intended to illustrate jMDP's capacity for modeling Markov Decision
 * Processes using events. It uses the classes InvLevel, Order and DemandEvent. 
 * This example models a car dealer's ordering it is a single item, periodic review, 
 * stochastic demand inventory problem. It is modeled like a long-run average cost, 
 * infinite horizon, Markov Decision Problem. Demand is assumed to be random according to 
 * a Poisson Process distribution with given rate per period. The system is a periodic
 * review problem in which an entity periodically checks the inventory level and
 * takes decisions according to the states he finds. There is a price of selling
 * each item and a cost for buying it. Besides, there is a holding cost incurred
 * when holding one item in stock from one period to another. There is also a
 * truck cost for ordering as units must be placed in trucks and only whole trucks can 
 * be sent, even if the 'last' truck is not full. We assume that backorders are not allowed.
 * The objective is to minimize the expected long-run average cost. We include a switch that 
 * allows modelling the problem with the objective of minimizing total discounted cost.
 * 
 * @author Daniel F. Silva
 */
public class CarDealerProblem extends DTMDPEv<InvLevel, Order, DemandEvent> {
    // Problem parameters:
    private int maxInv, truckSize;
    // Cost and demand parameters:
    private double truckCost, holdingCost, expDemand, price, cost, discountFactor;
    // Demand distribution information:
    private double[] demPMF, demCDF, demandLoss1;
    private boolean isdisc = false;


    // Constructor
    /**
     * @param maxInv
     *            maximum physical capacity in inventory, warehouse size.
     * @param truckSize
     *            maximum items in each fixed cost order. Orders can be greater
     *            than this value, but will be charged more than one fixed cost.
     * @param truckCost
     *            fixed cost per order
     * @param price
     *            unit price
     * @param cost
     *            unit aquistion cost
     * @param holdingCost
     *            non-finantial holding cost (it does NOT include finantial cost)
     * @param discountFactor
     *            interest per period
     * @param expDemand
     *            demand mean
     * @param discounted
     *            Whether a discounted model (rather than average) is to be used.
     */

    @SuppressWarnings("unchecked")
    public InfStochasticDemand(int maxInv, int truckSize, double truckCost, double price, double cost,
            double holdingCost, double intRate, double expDemand, boolean discounted) {
        super(new StatesSet<InvLevel>(new InvLevel(0)));
        this.maxInv = maxInv;
        this.truckSize = truckSize;
        this.truckCost = truckCost;
        this.price = price;
        this.cost = cost;
        this.holdingCost = holdingCost;
        this.expDemand = expDemand;
        initializeProbabilities();
        this.isdisc = discounted;
        this.discountFactor = intRate;
        if (discounted)
            setSolver(new ValueIterationSolver(this, DiscountFactor));
        else
            setSolver(new RelativeValueIterationSolver(this));
    }
    
    /**
     * Creates vectors with the probability density function (pdf), the 
     * complimentary cumulative distribution function (CCDF) and the
     * first order Loss function for a Poisson distribution with parameter 
     * expDemand. A Poisson has infinite support, but the model only requires 
     * this information fro k=0,...,maxInv, so that is all we calculate.
     */
    private void initializeProbabilities() {
        demPMF = new double[maxInv + 1];
        demCCDF = new double[maxInv + 1];
        demandLoss1 = new double[maxInv + 1];
        double cdf, p = Math.exp(-expDemand);
        demPMF[0] = p;
        cdf = 1-p; 
        demCCDF[0]=0;
        demandLoss1[0] = expDemand;
        int maxlevel = maxInv + maxBO;
        for (int i = 1; i <= maxlevel; i++) {
            demCCDF[i] = 1-cdf; // P{demand >= i}
        	demPMF[i] = (p *= expDemand / i); // P{demand = i}
            cdf += p; // P{demand <= i}
            demandLoss1[i] = (expDemand - i) * (1 - cdf) + expDemand * p;  // = E[(D-i)^+]
        }
    }

    /**
     * Determines the feasible actions for state i. Namely, the actions
     * are all order amounts between 0 and maxInv-i.
     * @param i the inventory level at the end of the week
     */
    @Override
    public Actions<Order> feasibleActions(InvLevel i) {
        int max = maxInv - i.getLevel();
        Order[] vec = new Order[max + 1];
        for (int n = 0; n <= max; n++) {
            vec[n] = new Order(n);
        }
        return new ActionsSet<Order>(vec);
    }
    
    /**
     * Determines the active events for state i AFTER action a is taken. 
     * Namely, we have 2 types of events. Those where demand realizations  
     * less than the available, so the demand is between 0 and i+a-1, 
     * and those where demand is  i+a or greater, we define whether an event 
     * is of each type by a boolen variable, and the sales amount by an integer,
     * where sales are (i+a-demand)^+. 
     * @param i the inventory level before ordering
     * @param a the order amount
     */
    @Override
    public Events<DemandEvent> activeEvents ( InvLevel i , Order a) {
    	EventsSet<DemandEvent> eventSet = new EventsSet<DemandEvent>();
    	eventSet .add(new DemandEvent( i . getLevel () + a. getSize () , true ) ) ;
    	for ( int n = 0; n < i . getLevel () + a. getSize ( ) ; n++) {
    		eventSet .add(new DemandEvent(n, false ) ) ;
    	}
    	return eventSet ;
    }
    
    /**
     * Determines the active atates that are reacheble from state i, when
     * action a is taken and event e occurs. Namely, the only reachable state is
     * j= i + a - e.
     * @param i the inventory level before ordering
     * @param a the order amount
     * @param e the demand event
     */
    @Override
    public States<InvLevel> reachable ( InvLevel i , Order a, DemandEvent e) {
    	StatesSet<InvLevel> stSet = new StatesSet<InvLevel >();
    	if (e . getGreaterThan ( ) )
    		stSet .add(new InvLevel (0));
    	else
    		stSet .add(new InvLevel ( i . getLevel () + a. getSize () - e .getDemand ( ) ) ) ;
    	return stSet ;    
    }
    
    /**
     * Determines the active atates that are reacheble from state i, when
     * action a is taken and event e occurs. Namely, the only reachable state is
     * j= i + a - e.
     * @param i the inventory level before ordering
     * @param a the order amount
     * @param e the demand event
     */
    @Override
    public double prob( InvLevel i , DemandEvent e) {
    	if (e . getGreaterThan ( ) )
    	return demCCDF[e .getDemand ( ) ] ;
    	return demPMF[e .getDemand ( ) ] ; 
    }
    
    // These 
    public double immediateCost( InvLevel i , Order a) {
    	int maxSale = i . getLevel () + a. getSize ( ) ;
    	double expectedSales = expDemand − demandLoss1[maxSale ];
    	double netProfit = price*expectedSales − orderCost (a. getSize())− holdCost * i . getLevel ( ) ;
    	return −netProfit ;
    	}
    	double orderCost ( int x) {
    	return truckCost*Math. ceil ((double) x /truckSize ) + x * cost ;
    }


    /**
     * Simple test Program.
     * 
     * @param a
     * @throws SolverException
     */
    public static void main(String a[]) throws SolverException {
        int maxInventory = 25;
        int maxBackorders = 0;
        int truckSize = 4;
        int truckCost = 1000;
        double holdCost = 50;
        double intRate = Math.pow(1.3, 1 / 52);
        double theta = 20;
        double price = 1100;// 22000;
        double cost = 500;// 20000;

        InfStochasticDemand prob = new InfStochasticDemand(maxInventory,
                maxBackorders, truckSize, truckCost, b, price, cost, holdCost, intRate, theta,
                false);

        RelativeValueIterationSolver<InvLevel, Order> solv = new RelativeValueIterationSolver<InvLevel, Order>(
                prob);

        prob.setSolver(solv);
        prob.getSolver().setPrintValueFunction(true);
        prob.solve();
        prob.printSolution();


    }

}