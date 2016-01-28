package examples.jmarkov;

import jmarkov.GeomProcess;
import jmarkov.GeomRelState;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.StatesSet;
import jphase.DenseContPhaseVar;
import jphase.PhaseVar;

import static examples.jmarkov.PriorityQueueMPHPHEvent.Type.ARRIVAL_HI;
import static examples.jmarkov.PriorityQueueMPHPHEvent.Type.SERVICE_END_HI;
import static examples.jmarkov.PriorityQueueMPHPHEvent.Type.SERVICE_PHASECHG_HI;
import static examples.jmarkov.PriorityQueueMPHPHEvent.Type.ARRIVAL_LOW;
import static examples.jmarkov.PriorityQueueMPHPHEvent.Type.SERVICE_END_LOW; 
import static examples.jmarkov.PriorityQueueMPHPHEvent.Type.SERVICE_PHASECHG_LOW; 


public class PriorityQueueMPHPH extends GeomProcess<PriorityQueueMPHPHState, PriorityQueueMPHPHEvent>{

	/**
	 * High-priority arrival rate
	 */
	double lambda_hi;
	
	/**
	 * Low-priority arrival rate
	 */
    double lambda_low;
    
    /**
     * High-priority service time PH variable
     */
    PhaseVar servTime_hi;
    
    /**
     * Low-priority service time PH variable
     */
    PhaseVar servTime_low;
    
    /**
     * 
     */
    int bufferCapacity;  
    
    
    /**
     * 
     * @param lambda_hi
     * @param lambda_low
     * @param servTime_hi
     * @param servTime_low
     */
    public  PriorityQueueMPHPH(double lambda_hi, double lambda_low, PhaseVar servTime_hi, PhaseVar servTime_low, int bufferCapacity) {
        super(new PriorityQueueMPHPHState(0,0), 
        		PriorityQueueMPHPHEvent.getAllEvents(servTime_hi, servTime_low));
        this.lambda_hi = lambda_hi;
        this.lambda_low = lambda_low;
        this.servTime_hi = servTime_hi;
        this.servTime_low = servTime_low;
        this.bufferCapacity = bufferCapacity; 
    }

    
    /**
     * Used by GUI.
     */
    public  PriorityQueueMPHPH() {
        this(1.0, 0.5, 
        		DenseContPhaseVar.HyperExpo(new double[] { 5.0, 8.0 }, new double[] { 0.5, 0.5 }),
        		DenseContPhaseVar.HyperExpo(new double[] { 3.0, 5.0 }, new double[] { 0.5, 0.5 }), 
        		10 );
    }
	
    
    @Override
    public boolean active(PriorityQueueMPHPHState state, int absLevel, PriorityQueueMPHPHEvent event) {

        boolean result = false;
        switch (event.eventType) {
	        case ARRIVAL_HI:
	        	if ( state.getNumberHiJobs() < bufferCapacity )
	        		result = true;
            	break;
	        case SERVICE_END_HI:
	            result =  (state.getNumberHiJobs() > 0 &&  state.getServicePhase() == event.eventPhase); 
	            break;
	        case SERVICE_PHASECHG_HI:
	            result =  (state.getNumberHiJobs() > 0 &&  state.getServicePhase() == event.eventPhase); 
	            break;
			case ARRIVAL_LOW:
				result = true; 
				break;
			case SERVICE_END_LOW:
				result =  (state.getNumberHiJobs() == 0 &&  state.getServicePhase() == event.eventPhase);
				break;
			case SERVICE_PHASECHG_LOW:
				result =  (state.getNumberHiJobs() == 0 &&  state.getServicePhase() == event.eventPhase);
				break;
        }
        return result;
    }

    
    @Override
    public GeomRelState<PriorityQueueMPHPHState>[] dests(PriorityQueueMPHPHState state,
            											int absLevel, PriorityQueueMPHPHEvent event) {
    	
        StatesSet<GeomRelState<PriorityQueueMPHPHState>> destStates = new StatesSet<GeomRelState<PriorityQueueMPHPHState>>();

        int curPhase = state.getServicePhase(); 
        int curNumHiJobs = state.getNumberHiJobs(); 
        int rLevel = 0;

        switch (event.eventType) {
	        case ARRIVAL_LOW:
	            rLevel = +1;
	            if (absLevel == 0 && curNumHiJobs == 0) 
            		// low-priority arrival in an empty system: level increases and starts service 
	                addDestsServiceEnd(rLevel, destStates, servTime_low, curNumHiJobs);
	            else // low-priority arrival in non-empty system: only modifies the level by 1
	                destStates.add(
	                			new GeomRelState<PriorityQueueMPHPHState>(
	                				new PriorityQueueMPHPHState(curNumHiJobs, curPhase), rLevel));
	            break;
	            
	        case SERVICE_END_LOW:
	            rLevel = -1;
	            if (absLevel == 1){
	            	if( curNumHiJobs == 0){ // low-priority service completion that leaves the system empty
	            		destStates.add(new GeomRelState<PriorityQueueMPHPHState>(
                						new PriorityQueueMPHPHState(0, 0))
	            						);
	            	}else{	// low-priority service completion that allows a high-priority job to start service 
	            		addDestsServiceEnd(rLevel, destStates, servTime_hi, curNumHiJobs);
	            	}
	            }
	            else{
	            	if( curNumHiJobs == 0){ // low-priority service completion, a new low-priority job starts service 
	            		addDestsServiceEnd(rLevel, destStates, servTime_low, curNumHiJobs);
	            	}else{ // low-priority service completion, a new high-priority job starts service
	            		addDestsServiceEnd(rLevel, destStates, servTime_hi, curNumHiJobs);
	            	}
	            }
	            break;
	            
	        case SERVICE_PHASECHG_LOW:
	            addDestsServiceChg(rLevel, destStates, servTime_low, curNumHiJobs, curPhase); 
	            break;
	            
			case ARRIVAL_HI:
				if (absLevel == 0 && curNumHiJobs == 0) 
            		// hi-priority arrival in an empty system: number of high jobs increases and starts service 
	                addDestsServiceEnd(rLevel, destStates, servTime_hi, curNumHiJobs+1);
	            else // hi-priority arrival in non-empty system: only increases the number of high-priority jobs
	                destStates.add(
	                			new GeomRelState<PriorityQueueMPHPHState>(
	                				new PriorityQueueMPHPHState(curNumHiJobs+1, curPhase), rLevel));
				break;
			case SERVICE_END_HI:
				if (absLevel == 0){
	            	if( curNumHiJobs == 1){ // hi-priority service completion that leaves the system empty
	            		destStates.add(new GeomRelState<PriorityQueueMPHPHState>(
                						new PriorityQueueMPHPHState(0, 0)) );
	            	}else{	// hi-priority service completion that allows a high-priority job to start service 
	            		addDestsServiceEnd(rLevel, destStates, servTime_hi, curNumHiJobs-1);
	            	}
	            }
	            else{
	            	if( curNumHiJobs == 1){ // hi-priority service completion, a new low-priority job starts service 
	            		addDestsServiceEnd(rLevel, destStates, servTime_low, curNumHiJobs-1);
	            	}else{ // hi-priority service completion, a new high-priority job starts service
	            		addDestsServiceEnd(rLevel, destStates, servTime_hi, curNumHiJobs-1);
	            	}
	            }
				
				break;
			case SERVICE_PHASECHG_HI:
	            addDestsServiceChg(rLevel, destStates, servTime_hi, curNumHiJobs, curPhase); 
	            break;
				
        }
        return destStates.toStateArray();
    }

    
    
    /**
     * Adds to the state set all destinations generated by a service initiation. 
     * @param relLevel: change in the level
     * @param destStates
     * @param serVar: service-time PH variable of the new job 
     */
    void addDestsServiceEnd(int relLevel, StatesSet<GeomRelState<PriorityQueueMPHPHState>> destStates, PhaseVar serVar, int numHiJobs) {
        double alpha[] = serVar.getVectorArray(); // initial probability distribution 
        int m = serVar.getNumPhases();
        for (int n = 0; n < m; n++) {
            if (alpha[n] > 0) {
                destStates.add( 
                		new GeomRelState<PriorityQueueMPHPHState>(new PriorityQueueMPHPHState(numHiJobs, n+1), relLevel)
                		);
            }
        }
    }
    
    /**
     * Adds to the state set all destinations generated by a service phase transition. 
     * @param relLevel: change in the level 
     * @param destStates
     * @param serVar: service-time PH variable of the new job 
     */
    void addDestsServiceChg(int relLevel, StatesSet<GeomRelState<PriorityQueueMPHPHState>> destStates, PhaseVar serVar, int numHiJobs, int currPhase) {
        double transMatrix[][] = serVar.getMatrixArray(); // sub-generator matrix 
        currPhase = currPhase - 1; 
        int m = serVar.getNumPhases();
        for (int n = 0; n < m; n++) {
            if (transMatrix[currPhase][n] > 0) {
                destStates.add( 
                		new GeomRelState<PriorityQueueMPHPHState>(new PriorityQueueMPHPHState(numHiJobs, n+1), relLevel)
                		);
            }
        }
    }
    
    
    
    @Override
    public double rate(PriorityQueueMPHPHState currState, int currLevel, 
    						PriorityQueueMPHPHState destState, int destLevel, PriorityQueueMPHPHEvent event) { {
    	
        StatesSet<GeomRelState<PriorityQueueMPHPHState>> destStates = new StatesSet<GeomRelState<PriorityQueueMPHPHState>>();

        int curPhase = state.getServicePhase(); 
        int curNumHiJobs = state.getNumberHiJobs(); 
        int rLevel = 0;

        switch (event.eventType) {
	        case ARRIVAL_LOW:
	            rLevel = +1;
	            if (absLevel == 0 && curNumHiJobs == 0) 
            		// low-priority arrival in an empty system: level increases and starts service 
	                addDestsServiceEnd(rLevel, destStates, servTime_low, curNumHiJobs);
	            else // low-priority arrival in non-empty system: only modifies the level by 1
	                destStates.add(
	                			new GeomRelState<PriorityQueueMPHPHState>(
	                				new PriorityQueueMPHPHState(curNumHiJobs, curPhase), rLevel));
	            break;
	            
	        case SERVICE_END_LOW:
	            rLevel = -1;
	            if (absLevel == 1){
	            	if( curNumHiJobs == 0){ // low-priority service completion that leaves the system empty
	            		destStates.add(new GeomRelState<PriorityQueueMPHPHState>(
                						new PriorityQueueMPHPHState(0, 0))
	            						);
	            	}else{	// low-priority service completion that allows a high-priority job to start service 
	            		addDestsServiceEnd(rLevel, destStates, servTime_hi, curNumHiJobs);
	            	}
	            }
	            else{
	            	if( curNumHiJobs == 0){ // low-priority service completion, a new low-priority job starts service 
	            		addDestsServiceEnd(rLevel, destStates, servTime_low, curNumHiJobs);
	            	}else{ // low-priority service completion, a new high-priority job starts service
	            		addDestsServiceEnd(rLevel, destStates, servTime_hi, curNumHiJobs);
	            	}
	            }
	            break;
	            
	        case SERVICE_PHASECHG_LOW:
	            addDestsServiceChg(rLevel, destStates, servTime_low, curNumHiJobs, curPhase); 
	            break;
	            
			case ARRIVAL_HI:
				if (absLevel == 0 && curNumHiJobs == 0) 
            		// hi-priority arrival in an empty system: number of high jobs increases and starts service 
	                addDestsServiceEnd(rLevel, destStates, servTime_hi, curNumHiJobs+1);
	            else // hi-priority arrival in non-empty system: only increases the number of high-priority jobs
	                destStates.add(
	                			new GeomRelState<PriorityQueueMPHPHState>(
	                				new PriorityQueueMPHPHState(curNumHiJobs+1, curPhase), rLevel));
				break;
			case SERVICE_END_HI:
				if (absLevel == 0){
	            	if( curNumHiJobs == 1){ // hi-priority service completion that leaves the system empty
	            		destStates.add(new GeomRelState<PriorityQueueMPHPHState>(
                						new PriorityQueueMPHPHState(0, 0)) );
	            	}else{	// hi-priority service completion that allows a high-priority job to start service 
	            		addDestsServiceEnd(rLevel, destStates, servTime_hi, curNumHiJobs-1);
	            	}
	            }
	            else{
	            	if( curNumHiJobs == 1){ // hi-priority service completion, a new low-priority job starts service 
	            		addDestsServiceEnd(rLevel, destStates, servTime_low, curNumHiJobs-1);
	            	}else{ // hi-priority service completion, a new high-priority job starts service
	            		addDestsServiceEnd(rLevel, destStates, servTime_hi, curNumHiJobs-1);
	            	}
	            }
				
				break;
			case SERVICE_PHASECHG_HI:
	            addDestsServiceChg(rLevel, destStates, servTime_hi, curNumHiJobs, curPhase); 
	            break;
				
        }
        return destStates.toStateArray();
    }

    
    
    
    
	public static void main(String[] a) {
        double lambda_hi = 1.0;
        double lambda_low = 0.5;
        PhaseVar servTime_hi = DenseContPhaseVar.HyperExpo(new double[] { 5.0, 8.0 },
                new double[] { 0.5, 0.5 });
        PhaseVar servTime_low = DenseContPhaseVar.HyperExpo(new double[] { 3.0, 5.0 },
                new double[] { 0.5, 0.5 });
        int bufferCapacity = 10; 
        
        PriorityQueueMPHPH model = new PriorityQueueMPHPH(lambda_hi, lambda_low, servTime_hi, servTime_low, bufferCapacity);
        
        
        model.showGUI();
        model.generate();
        model.setDebugLevel(0);
        model.printAll();
    }
}






class PriorityQueueMPHPHState extends PropertiesState {
	
	public PriorityQueueMPHPHState(int numberHiJobs, int servicePhase) {
		super(2);
        setProperty(0, numberHiJobs);
        setProperty(0, servicePhase);
    }
	
	
	public int getNumberHiJobs() {
        return this.prop[0];
    }
	
	public int getServicePhase() {
        return this.prop[1];
    }
	
}

class PriorityQueueMPHPHEvent extends Event {
	
	
	public enum Type {
        /** High-priority arrivals. */
        ARRIVAL_HI,
        /** High-priority service completion. */
        SERVICE_END_HI,
        /** High-priority service phase change. */
        SERVICE_PHASECHG_HI,
        /** Low-priority arrivals. */
        ARRIVAL_LOW,
        /** Low-priority service completion. */
        SERVICE_END_LOW,
        /** Low-priority service phase change. */
        SERVICE_PHASECHG_LOW
        
    }
	
	/**
	 * Type of event
	 */
	Type eventType; 
	
	/**
	 * Phase where event occurs
	 */
	int eventPhase;
	
	
	/** Arrival event */
	PriorityQueueMPHPHEvent(Type eventType) {
		if( eventType == ARRIVAL_HI || eventType == ARRIVAL_LOW){
			this.eventType = eventType;
		}
    }

    /** Service event */
    PriorityQueueMPHPHEvent(Type type, int phase) {
        this.eventType = type;
        this.eventPhase = phase;
    }

	
	static EventsSet<PriorityQueueMPHPHEvent> getAllEvents(PhaseVar servTime_hi, PhaseVar servTime_low) {
		
        EventsSet<PriorityQueueMPHPHEvent> E = new EventsSet<PriorityQueueMPHPHEvent>();
        
        // high-priority arrival event 
        E.add(new PriorityQueueMPHPHEvent(ARRIVAL_HI));
        
        // high-priority service completion and service change in each phase 
        int numPhases_hi = servTime_hi.getNumPhases();
        for (int n = 1; n <= numPhases_hi; n++) {
            // high-priority service completion in phase n
            E.add(new PriorityQueueMPHPHEvent(SERVICE_END_HI, n));
            E.add(new PriorityQueueMPHPHEvent(SERVICE_PHASECHG_HI, n));
        }
        
        // low-priority arrival event 
        E.add(new PriorityQueueMPHPHEvent(ARRIVAL_LOW));
        
        // high-priority service completion and service change in each phase 
        int numPhases_low = servTime_low.getNumPhases();
        for (int n = 1; n <= numPhases_low; n++) {
            // high-priority service completion in phase n
            E.add(new PriorityQueueMPHPHEvent(SERVICE_END_LOW, n));
            E.add(new PriorityQueueMPHPHEvent(SERVICE_PHASECHG_LOW, n));
        }
        
        return E;
    }
	
}