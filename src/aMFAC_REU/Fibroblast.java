package aMFAC_REU;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.parameter.Parameters;
import repast.simphony.query.space.grid.GridCell;
import repast.simphony.query.space.grid.GridCellNgh;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.valueLayer.GridValueLayer;


/**
 * Repast Fibroblast class for Coupled Model
 * @author Tommy Athey
 * Aug 2017
 */
public class Fibroblast {
	private AMFACSpace space;
	Grid grid;

	
	private double[] networkState;
	
	private Parameters p = RunEnvironment.getInstance().getParameters();
	private boolean movement = (Boolean) p.getValue("Fibroblast_move"); //whether fibroblasts should move or not
	
	//These 11 layers are the 11 inputs to the Saucerman network model
	//you will find identical variables at the top of AMFACSpace
	private String[] inputLayerNames = {"TGFB", "Interleukin6", "Interleukin1",
			"TNFalpha"};
	private int[] inputNetworkIndices = {19,38,41,43};
	private ArrayList<GridValueLayer> inputLayers = new ArrayList<GridValueLayer>();
	private int tnfIdx = 4; //used for migration up tnf gradient
	
	
	//network inputs that stay constant
	//note that they only start at constant, if there is feedback, these inputs may grow
	private int[] constantIdxs = {0,31,11,25,6,27,13};
	private double[] constantVals = {0.25,0.25,0.25,0.25,0.25,0.25,0.25};


	private GridValueLayer collagen; //collagen layer
	private GridValueLayer TGFB;
	private GridValueLayer LatentTGFB;
	private GridValueLayer IL6;
	private GridValueLayer IL1;
	private GridValueLayer TNFalpha;
	
	private double kgen = 1; //area fraction/hr
	private double kdeg = 1; //area fraction/hr		
	
	private GridPoint pt;
	
	private int cellsPerGrid;
	
	
	/**
	 * Constructor reads the parameters from the GUI
	 *
	 * @param sp the space in which the cell is placed
	 */
	public Fibroblast(AMFACSpace sp, int cells) {
		
		space = sp;
		grid = (Grid) space.getProjection("grid");
		

		//set the network state to 0 just like the y0 in the model
		networkState = new double[91];
		
		
		collagen = (GridValueLayer) space.getValueLayer("collagen");

		
		cellsPerGrid = cells;
	}
	
	/**
	 * Gets initial location
	 * (can't be in the constructor because the cell
	 * has not yet been added to the world at that point)
	 */
	
	@ScheduledMethod(start = 1, interval = 1, priority = 1)
	public void goThird(){
		//live();
		updateCollagen();
		
		if (movement == true){
			move();
		}
		//System.out.println("Third");

	}
	
	@ScheduledMethod(start = 1, interval = 1, priority = 3)
	public void goFirst(){
		getCellNetwork();
		//System.out.println("First");
	}
	
	@ScheduledMethod(start = 0, priority = 0)
	public void initialize() {
		pt = grid.getLocation(this);
		
		for (int i = 0; i<inputLayerNames.length; i++) {
			inputLayers.add((GridValueLayer) space.getValueLayer(inputLayerNames[i]));
		}
		
		
		for (int i = 0; i < constantIdxs.length; i++) {
			networkState[constantIdxs[i]] = constantVals[i];
		}
		
		int x = pt.getX();
		int y = pt.getY();
			
		
/*		GridValueLayer layer; //for each layer
			for (int i = 0; i<inputLayerNames.length; i++) {
				layer = inputLayers.get(i);
				networkState[inputNetworkIndices[i]] = layer.get(pt.getX(),pt.getY());
			}*/
		
		//System.out.println("Initialize");
		
	}
	
	/**
	 * getter method for the current network state
	 * Used by the Space while processing network simulations
	 * @return network
	 */
	public double[] getNetworkState() {
		return networkState;
	}
	
	/**
	 * Setter method for the network state
	 * Used by the Space while processing network simulations
	 * @return true if the network was the same size as the original network (network was successfully updated)
	 * 			false if the old and new networks are different lengths (network not set successfully)
	 */
	public boolean setNetworkState(double[] n) {
		if (n.length != networkState.length) {
			return false;
		} else {
			networkState = n;
			return true;
		}
	}
	
	/**
	 * Updates the cell's network based on local cytokine values
	 * Does not reset the "constant" inputs back to 0.25
	 */
	
	public void getCellNetwork() {
		pt = grid.getLocation(this);
		
		int x = pt.getX();
		int y = pt.getY();
			
		
			GridValueLayer layer; //for each layer
			for (int i = 0; i<inputLayerNames.length; i++) {
				layer = inputLayers.get(i);
				networkState[inputNetworkIndices[i]] = layer.get(pt.getX(),pt.getY());
			}
		
		//System.out.println("Get Cell Network");
		
	}
	
	
	/**
	 * Returns the fibroblasts location
	 * @return the GridPoint of the cell's location. Used when the space is calculating metrics across the grid
	 */
	public GridPoint getPoint() {
		return pt;
	}
	

	/**
	 * Deals with mitosis and apoptosis of the cell
	 * 
	 * Apoptosis occurs at a constant probability
	 */
	
	/*public void live() {
		pt = grid.getLocation(this);
		
		//System.out.println("Live");

		double mitosisTime = (minMitRate - maxMitRate) * (1 - networkState[69]) + maxMitRate;

		double choice = RandomHelper.nextDoubleFromTo(0, 1);

		if (choice < mitosisTime) {
			
			List<GridPoint> empty = this.findEmptySites();
			if (!empty.isEmpty()) {
				pt = empty.get(0);
				space.addFibroblast(pt, this.networkState);
			}
			
		} else if (choice < mitosisTime + apopProb) { // apoptose
			space.removeFibroblast(this);
		}

	}*/


	public void updateCollagen() {
		pt = grid.getLocation(this);
		
		//System.out.println("Deposit");
		
		try {
		int x = pt.getX();
		int y = pt.getY();
			
		
			collagen = (GridValueLayer) space.getValueLayer("collagen");
		
			double currentCollagen = collagen.get(pt.getX(), pt.getY());
			
			double depLevel = (networkState[87] + networkState[88])/2; //average of CmRNAs
			double degLevel = (networkState[81] + networkState[82] + networkState[83])/3; //average of MMP1, 2, and 9

			double dcdt = depLevel - (degLevel*currentCollagen);
			
			//System.out.println(depLevel);
			
			collagen.set(dcdt, x, y);
		}
		catch  (Exception e) {
			System.out.println(e);
			e.printStackTrace();
		}
		
	}

	
	/**
	 * Direction is chosen by the gradient of tnfalpha
	 * If there is a 10% gradient-the direction will be biased toward the highest gradient
	 */
	
	public void move() {

		pt = grid.getLocation(this);
		
		//System.out.println("Move");
		
		try {
		int x = pt.getX();
		int y = pt.getY();
		
		List<GridPoint> emptySites = findEmptySites();
		GridPoint destination;
		destination = emptySites.get(0);
		grid.moveTo(this, destination.getX(), destination.getY());
		
/*		TNFalpha = (GridValueLayer) space.getValueLayer("TNFalpha");
		double maxTNFalpha = 0;
		GridPoint destination = null;
		
		List<GridPoint> emptySites = findEmptySites();
		for (GridPoint emptySite : emptySites) {
			double TNFvalue = TNFalpha.get(emptySite.getX(), emptySite.getY());
			if (TNFvalue > maxTNFalpha){
				destination = emptySite;
				maxTNFalpha = TNFvalue;
			}
		}
		
		grid.moveTo(this, destination.getX(), destination.getY());*/
		
		}
		
		catch  (Exception e) {
			System.out.println(e);
		}
	}
		/*
		double speed = (maxSpeed - minSpeed) * (networkState[66] - 1) + maxSpeed;
			//cells only move one grid at a time so this warning says that this restriction might not be desirable
			if (speed > 1) { 
				System.out.println("Warning: Calculated Fibroblast speed > 1");
			}
			double mv = RandomHelper.nextDoubleFromTo(0, 1);
			if (mv < speed) {
				double pctThreshold = 0.1; //gradient percentage threshold needed to indicate a strong gradient
				
				TNFalpha = (GridValueLayer) space.getValueLayer("TNFalpha");
				double currVal = TNFalpha.get(x, y);
			
				List<GridPoint> emptySites = findEmptySites();
				List<GridPoint> possibleDestinations = new ArrayList<GridPoint>();
				List<Double> grads = new ArrayList<Double>();
				double gradient;
				for (GridPoint emptySite : emptySites) {
					double otherVal = TNFalpha.get(emptySite.getX(), emptySite.getY());
					gradient = otherVal - currVal;
					//gradient must be substantial, but there also must not be extreme
					//(negligible or saturating) values of chemokines
				if (gradient > pctThreshold*currVal && currVal < 0.9 && otherVal > 0.1) {
					possibleDestinations.add(emptySite);
					grads.add(new Double(gradient));
					}
				}
			
			GridPoint destination;
			if (!possibleDestinations.isEmpty()) {
				destination = chooseDestination(possibleDestinations, grads);
				grid.moveTo(this, destination.getX(), destination.getY());
			} else if (!emptySites.isEmpty()) { //random movement
				destination = emptySites.get(0);
				grid.moveTo(this, destination.getX(), destination.getY());
			}
			pt = grid.getLocation(this);
			}
		}
		catch  (Exception e) {
			System.out.println(e);
			e.printStackTrace();
		}
	
	}

	

	/**
	 * Finds the cells in the adjacent Moore neighborhood that have less than 36 cells
	 * @return list of grid points that have less than 36 cells
	 */
	private List<GridPoint> findEmptySites() {
		List<GridPoint> emptySites = new ArrayList<GridPoint>();
		try {

			GridCellNgh<Fibroblast> nghCreator = new GridCellNgh<Fibroblast>(grid, pt, Fibroblast.class, 1, 1);
			List<GridCell<Fibroblast>> gridCells = nghCreator.getNeighborhood(true);
			for (GridCell<Fibroblast> cell : gridCells) {
				if (cell.size() < cellsPerGrid)
					emptySites.add(cell.getPoint());
			}

			Collections.shuffle(emptySites);
		} catch (NullPointerException e) {
		}
		return emptySites;
	}

	
	/**
	 * Chooses a direction by giving each direction a probability proportional to its gradient value
	 * @param possibleDestinations list of gridpoints that have a 10% gradient
	 * @param grads chemokine gradients (difference of tnf values)
	 * @return a grid point chosen with probability biased toward grad values
	 */
	private GridPoint chooseDestination(List<GridPoint> possibleDestinations, List<Double> grads) {
		double[] probIntervals = new double[grads.size()];
		for (int i = 0; i < grads.size(); i++) {
			if (i == 0) {
				probIntervals[i] = Math.abs(grads.get(i));
			} else {
				probIntervals[i] = Math.abs(grads.get(i)) + grads.get(i-1); //cumulative probability
			}
		}
		
		GridPoint pt = null;
		double c = RandomHelper.nextDoubleFromTo(0, probIntervals[probIntervals.length - 1]);
		for (int i = 0; i < probIntervals.length; i++) {
			if (c <= probIntervals[i]) { //the first cumulative probability that the random number is less than
				pt =  possibleDestinations.get(i);
			}
		}
		return pt;
	}

}
