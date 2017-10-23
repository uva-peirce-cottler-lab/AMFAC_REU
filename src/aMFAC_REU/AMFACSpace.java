/**
 * Repast Space class for Coupled Model
 * @author Tommy Athey
 * Aug 2017
 */
package aMFAC_REU;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import com.mathworks.engine.MatlabEngine;

import repast.simphony.scenario.data.Classpath;

import repast.simphony.context.DefaultContext;
import repast.simphony.context.space.grid.GridFactoryFinder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.parameter.Parameters;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridBuilderParameters;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.space.grid.RandomGridAdder;
import repast.simphony.space.grid.WrapAroundBorders;
import repast.simphony.util.ClassPathEntry;
import repast.simphony.valueLayer.GridValueLayer;
import repast.simphony.valueLayer.ValueLayerDiffuser;

public class AMFACSpace extends DefaultContext<Object> {
	
	Grid<Object> grid;
	
	private Parameters p = RunEnvironment.getInstance().getParameters();
	private int gridWidth = (Integer) p.getValue("gridWidth");
	private int gridHeight = (Integer) p.getValue("gridHeight");
	int initialFibroblastCount = (Integer) p.getValue("initialFibroblastCount");
	
	//arraylist that will hold all fibroblasts in the world (useful for iteration)
	private ArrayList<Fibroblast> fibroblasts = new ArrayList<Fibroblast>();

	//to keep track of all the layers.
	//Entries are in the same order as the input entries in the Saucerman model file
	private ArrayList<GridValueLayer> inputLayers = new ArrayList<GridValueLayer>();
	private ArrayList<ValueLayerDiffuser> inputDiffuseLayers = new ArrayList<ValueLayerDiffuser>();
	
	private ArrayList<GridValueLayer> otherExtracellularLayers = new ArrayList<GridValueLayer>();;
	private ArrayList<ValueLayerDiffuser> otherExtracellularDiffuseLayers = new ArrayList<ValueLayerDiffuser>();
	
	private GridValueLayer collagen;
	
	//These layers will be added to inputLayers
	private String[] inputLayerNames = {"TGFB","Interleukin6", "Interleukin1",
				"TNFalpha"};
	private int[] inputDiffuseIdxs = {}; //numerical indices of inputLayerNames that should diffuse. MUST BE IN INCREASING ORDER
	
	//ints alternate from network indices to indices in the inputLayerNames arraylist
	//for example, the 38th species in the network model is has idx 1 in inputLayerNames (TGFb)
	//remember that java indexes at 0
	private int[] networkLayerInputIdxs = {38,1, 19,0}; 
	
	//indices of inputLayerNames that are inflammatory or anti-inflammatory/fibrotic
	//used for gradient orientation (e.g. the inflammatory cytokines are a gradient from left to right
	private int[] inflam = {1,2,3};
	private int[] antiInflam = {0};
	
	//other cytokines that you might want to keep track of such as latentTGFb etc.
	private String[] otherExtracellularNames = {}; //analagous to inputLayerNames above
	private int[] otherExtracellularDiffuseIdxs = {}; //analagous to inputDiffuseIdxs. MUST BE IN INCREASING ORDER
	//ints alternate from network indices to indices in the inputLayer arraylist
	private int[] networkLayerOtherIdxs = {}; //analagous to networkLayerInputIdxs
	
	
	private MatlabEngine eng;
	
	//these variable are used for cytokines that change over time
	private ArrayList<Double> relChemVals; //the array of relative levels
	private int relChemFactorIdx;
	

	private int cellsPerGrid = 1;
	private int chemokineFeedbackTimeConstant = 36; //how much a cell's network influences local concentrations of chemicals in 
														//networkLayerInputIdxs and networkLayerOtherIdxs

	double[] initialNet;
	
	/**
	 * Create all fields necessary for this class
	 */
	public AMFACSpace() {
		super("AMFACSpace");

		// Define the Grid Space
		grid = GridFactoryFinder.createGridFactory(null).createGrid("grid", this,
				new GridBuilderParameters<Object>(new WrapAroundBorders(), new RandomGridAdder<Object>(), true,
						gridWidth, gridHeight));

		// Build all GridValueLayers
		collagen = new GridValueLayer("collagen", 1.0, true, gridWidth, gridHeight);
		this.addValueLayer(collagen);
		
		//create grid value layers
		GridValueLayer layer;
		for (int i=0; i < inputLayerNames.length; i++) {
			layer = new GridValueLayer(inputLayerNames[i], 1.0, true, gridWidth, gridHeight);
			inputLayers.add(layer);
			this.addValueLayer(layer);
		}
		
		for (int i=0; i < otherExtracellularNames.length; i++) {
			layer = new GridValueLayer(otherExtracellularNames[i], 1.0, true, gridWidth, gridHeight);
			otherExtracellularLayers.add(layer);
			this.addValueLayer(layer);
		}
		
		//read the time signal for the chemokines
		String csvFile = "RelativeChemValues.csv";
		BufferedReader br = null;
		String line = "";
		String csvSplitBy = ",";
		relChemVals = new ArrayList<Double>();
		try {
			br = new BufferedReader(new FileReader(csvFile));
			while ((line = br.readLine()) != null) {
				String[] in = line.split(csvSplitBy);
				try {
					relChemVals.add(new Double(Double.parseDouble(in[0])));
				} catch (NumberFormatException e) {
					System.out.println("Non numerical Value found in RelativeChemValues: " + in[0] + "...");
				}
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		relChemFactorIdx = 0;

		
		
		//create diffusers
		Parameters p = RunEnvironment.getInstance().getParameters();
		double diffEvap = (Double) p.getValue("diffEvap");
		double diffCoeff = (Double) p.getValue("diffCoeff");
		ValueLayerDiffuser diffuse;
		int difIdx;
		difIdx = 0;
		
		for (int i=0; i < inputLayers.size(); i++) {
			if (difIdx < inputDiffuseIdxs.length && i==inputDiffuseIdxs[difIdx]) { //if this layer needs a diffuser
				diffuse = new ValueLayerDiffuser(inputLayers.get(inputDiffuseIdxs[i]), diffEvap, diffCoeff, true);
				inputDiffuseLayers.add(diffuse);
				difIdx++;
			} else {
				inputDiffuseLayers.add(null); //add null elements to fill in the gaps
			}
		}
		
		
		for (int i=0; i < otherExtracellularDiffuseLayers.size(); i++) {
			if (difIdx < otherExtracellularDiffuseIdxs.length && i==otherExtracellularDiffuseIdxs[difIdx]) {
				diffuse = new ValueLayerDiffuser(inputLayers.get(otherExtracellularDiffuseIdxs[i]), diffEvap, diffCoeff, true);
				otherExtracellularDiffuseLayers.add(diffuse);
				difIdx++;
			} else {
				otherExtracellularDiffuseLayers.add(null);
			}
		}
		
		//start the matlab connection
		try {
			eng = MatlabEngine.startMatlab();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		//add fibroblasts
		Fibroblast fibroblast;
		for (int i = 0; i < initialFibroblastCount; i++) {
			fibroblast = new Fibroblast(this, cellsPerGrid);
			this.add(fibroblast);
			fibroblasts.add(fibroblast);
		}
	}

	/**
	 * Sets the collagen layer values as uniform
	 */
	@ScheduledMethod(start = 0, priority = 1)
	public void initializeCollagenLayer(){
		for (int i = 1; i <= gridWidth; i++) {
			for (int j = 1; j <= gridHeight; j++) {
				collagen.set(0.25, i - 1, j - 1);
			}
		}
	}
	

	@ScheduledMethod(start = 0, priority = 2)
	public void initializeChemokineLayer() {
		double wdub = (double) gridWidth;
		double hdub = (double) gridHeight;
		for (int x = 0; x < gridWidth; x++) {
			for (int y = 0; y < gridHeight; y++) {
				//inflammatory chemokines
				for (int i = 0; i < inflam.length; i++) {
					inputLayers.get(inflam[i]).set(((double)x)/wdub, x, y);
				}
				//fibrotic chemokines
				for (int i = 0; i < antiInflam.length; i++) {
					inputLayers.get(antiInflam[i]).set(((double)y)/hdub, x, y);
				}
			}
		}

		// load the mat file only once
		try {
			eng.eval("load network.mat");
			eng.eval("load initialNet.mat");
			initialNet = eng.getVariable("initialNet");
			for (Fibroblast f: fibroblasts) {
				f.setNetworkState(initialNet);
			}
			initialNet = eng.getVariable("initialNet");
			//uncomment if you want to fool areound with parallel computing
			/*eng.eval(" myCluster = parcluster('local')");
			eng.eval(" myCluster.NumWorkers = 3");
			eng.eval("parpool('local',3);");*/
		} catch (Exception e) {
			System.out.println(e);
		}
		
	}
	
	/**
	 * Initialize the layers of the other extracellular species
	 */
	@ScheduledMethod(start = 0, priority = 2)
	public void initializeOtherExLayer() {
		
		GridValueLayer layer;
		for (int l=0; l < otherExtracellularNames.length; l++) {
			layer = otherExtracellularLayers.get(l);
			for (int i = 0; i < gridWidth; i++) {
				for (int j = 0; j < gridHeight; j++) {
					layer.set(0, i,j);
				}
			}
		}
	}
	
	/**
	 * Diffuses each diffusable layer
	 */
	@ScheduledMethod(start = 0, interval = 1, priority = 1)
	public void diffuse() {
		for (ValueLayerDiffuser layer : inputDiffuseLayers) {
			if (layer != null) {
				layer.diffuse();
			}
			
		}
		for (ValueLayerDiffuser layer : otherExtracellularDiffuseLayers) {
			if (layer != null) {
				layer.diffuse();
			}
		}
	}
	
	/**
	 * Updates chemokine layers to show long term time courses
	 * uses data from RelativeChemLevels.csv
	 * WARNING: This only serves as example code that could be used to make long term time courses
	 * I do not actually use this method in my simulations and have not debugged all possible errors
	 * such as a proper initial cytokine profile
	 */
	//@ScheduledMethod(start = 0, interval = 24, priority = 1) //commented out because I don't actually want it to run
	public void relChemVal() {
		if (relChemFactorIdx < relChemVals.size()-1) {
			relChemFactorIdx++;
		}
		double coeff = relChemVals.get(relChemFactorIdx)/relChemVals.get(relChemFactorIdx-1);
		
		for (int i=0; i < inputLayers.size(); i++) {
			GridValueLayer layer = inputLayers.get(i);
			for (int x =0; x < gridWidth; x++) {
				for (int y = 0; y < gridHeight; y++) {
					layer.set(coeff*layer.get(x,y), x,y);
				}
			}
		}
	}
	
	/**
	 * Writes data every 24 ticks
	 */
	@ScheduledMethod(start = 1, interval=24, priority = 1)
	public void writeCollagenData() {
		//collagen
		String fileName = "C:\\Users\\Michaela\\Documents\\col.csv";
		String delim = ",";
		String newline = "\n";
		FileWriter fileWriter;
		try {
			fileWriter = new FileWriter(fileName, true);
			for (int x = 0; x < gridWidth; x++) {
				for (int y = 0; y < gridHeight; y++) {
					fileWriter.append(Double.toString(collagen.get(x,y)));
					//if not at the end
					if (!(x+1 == gridWidth && y+1 == gridHeight)) {
						fileWriter.append(delim);
					}
				}
			}
			fileWriter.append(newline);
			try {
				fileWriter.flush();
				fileWriter.close();
			} catch (IOException e) {
				System.out.println("Error while flushing/closing fileWriter !!!");
				e.printStackTrace();
			}
			
		} catch (Exception e) {
			System.out.println("Error in CsvFileWriter !!!");
			e.printStackTrace();
		}
		
		
		
		
		//tgfb
		fileName = "C:\\Users\\Michaela\\Documents\\tgfb.csv";
		try {
			fileWriter = new FileWriter(fileName, true);
			for (int x = 0; x < gridWidth; x++) {
				for (int y = 0; y < gridHeight; y++) {
					fileWriter.append(Double.toString(inputLayers.get(0).get(x,y))); //index 0 is hardcoded in here
					fileWriter.append(delim);
				}
			}
			fileWriter.append(newline);
			try {
				fileWriter.flush();
				fileWriter.close();
			} catch (IOException e) {
				System.out.println("Error while flushing/closing fileWriter !!!");
				e.printStackTrace();
			}
			
		} catch (Exception e) {
			System.out.println("Error in CsvFileWriter !!!");
			e.printStackTrace();
		}
		
		
		
		//il6
		fileName = "C:\\Users\\Michaela\\Documents\\il6.csv";
		try {
			fileWriter = new FileWriter(fileName, true);
			for (int x = 0; x < gridWidth; x++) {
				for (int y = 0; y < gridHeight; y++) {
					fileWriter.append(Double.toString(inputLayers.get(1).get(x,y))); //index 1 is hardcoded in here
					fileWriter.append(delim);
				}
			}
			fileWriter.append(newline);
			try {
				fileWriter.flush();
				fileWriter.close();
			} catch (IOException e) {
				System.out.println("Error while flushing/closing fileWriter !!!");
				e.printStackTrace();
			}
			
		} catch (Exception e) {
			System.out.println("Error in CsvFileWriter !!!");
			e.printStackTrace();
		}
	}
	
	/**
	 * Writes data every 24 ticks
	 */
	@ScheduledMethod(start = 1, interval=24, priority = 1)
	public void writeCellData() {
		double[] cellNetwork;
		//averages for every quadrant
		double[][] cellCount = new double[2][2];
		double[][] prolifQuad = new double[2][2];
		double[][] migQuad = new double[2][2];
		double[][] depQuad = new double[2][2];
		double[][] degQuad = new double[2][2];
		//network for every grid
		double[][][] avgNet = new double[gridWidth][gridHeight][91];
		
		Fibroblast f;
		//totals for every grid square
		double counter, mig, prolif, dep, deg;
		double[] gridNetwork;

		int xQuad, yQuad;
		//iterate through every grid element
		for (int x = 0; x < gridWidth; x++) {
			//make sure you know which quadrant to be in
			if ( x < gridWidth/2) {
				xQuad = 0;
			} else {
				xQuad = 1;
			}
			
			for (int y = 0; y < gridHeight; y++) {
				
				if ( y < gridHeight/2) {
					yQuad = 0;
				} else {
					yQuad = 1;
				}
				
				//these variables represent a particular grid coordinate
				counter = 0;
				mig = 0;
				prolif = 0;
				dep = 0;
				deg = 0;
				gridNetwork = new double[91];
				//for every cell in this grid coordinate
				for (Object cell : grid.getObjectsAt(x, y)) {
					counter++;
					f = (Fibroblast) cell; //this will need to change if you have different agent types in the world
											//perhaps a solution is to have a try catch when casting the object
					
					//net is for a particular cell
					cellNetwork = f.getNetworkState();
					//add to the grid coordinate variables
					for (int i = 0; i < cellNetwork.length; i++) { //for every node in the network
						gridNetwork[i] += cellNetwork[i];
					}
					mig += cellNetwork[66];
					prolif += cellNetwork[69];
					dep += (cellNetwork[87] + cellNetwork[88]) / 2;
					deg += (cellNetwork[81] + cellNetwork[82] + cellNetwork[83] + cellNetwork[84]) / 4;
				}
				//normalize the whole network array
				for (int i = 0; i < gridNetwork.length; i++) {
					if (counter == 0) {
						gridNetwork[i] = -1; //no cells here
					} else {
						gridNetwork[i] = gridNetwork[i] / counter; //average of networks at that grid coordinate
					}
					
				}
				//set the average network for that coordinate
				avgNet[x][y] = gridNetwork;
				//add to the appropriate quadrant position
				cellCount[xQuad][yQuad] += counter;
				migQuad[xQuad][yQuad] += mig;
				prolifQuad[xQuad][yQuad] += prolif;
				depQuad[xQuad][yQuad] += dep;
				degQuad[xQuad][yQuad] += deg;
			}
		}

		//turn the quad data into averages (or -1 if there were no cells there)
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				if (cellCount[i][j] == 0) {
					migQuad[i][j] = -1;
					prolifQuad[i][j] = -1;
					depQuad[i][j] = -1;
					degQuad[i][j] = -1;
				} else {
					migQuad[i][j] = migQuad[i][j] / cellCount[i][j];
					prolifQuad[i][j] = prolifQuad[i][j] / cellCount[i][j];
					depQuad[i][j] = depQuad[i][j] / cellCount[i][j];
					degQuad[i][j] = degQuad[i][j] / cellCount[i][j];
				}
			}
		}
		
		//write the data
		String fileNameCompleteNet = "C:\\Users\\cards\\workspace\\AMFAC_REU\\completenetwork.csv";
		String fileNameMig = "C:\\Users\\cards\\workspace\\AMFAC_REU\\mig.csv";
		String fileNameProlif = "C:\\Users\\cards\\workspace\\AMFAC_REU\\prolif.csv";
		String fileNameDep = "C:\\Users\\cards\\workspace\\AMFAC_REU\\dep.csv";
		String fileNameDeg = "C:\\Users\\cards\\workspace\\AMFAC_REU\\deg.csv";
		String fileNameCount = "C:\\Users\\cards\\workspace\\AMFAC_REU\\cellcount.csv";
		
		String delim = ",";
		String newline = "\n";
		FileWriter fileNet, fileMig, fileProlif, fileDep, fileDeg, fileCount;
		try {
			fileNet = new FileWriter(fileNameCompleteNet, true);
			fileMig = new FileWriter(fileNameMig, true);
			fileProlif = new FileWriter(fileNameProlif, true);
			fileDep = new FileWriter(fileNameDep, true);
			fileDeg = new FileWriter(fileNameDeg, true);
			fileCount = new FileWriter(fileNameCount, true);
			
			//order of quadrants is 3, 2, 4, 1
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					fileMig.append(Double.toString(migQuad[i][j]));
					fileProlif.append(Double.toString(prolifQuad[i][j]));
					fileDep.append(Double.toString(depQuad[i][j]));
					fileDeg.append(Double.toString(degQuad[i][j]));
					fileCount.append(Double.toString(cellCount[i][j]));
					
					if (!(i+1==2 && j+1==2)) {
						fileMig.append(delim);
						fileProlif.append(delim);
						fileDep.append(delim);
						fileDeg.append(delim);
						fileCount.append(delim);
					}
				}
			}
			fileMig.append(newline);
			fileProlif.append(newline);
			fileDep.append(newline);
			fileDeg.append(newline);
			fileCount.append(newline);
			
			for (int x = 0; x < gridWidth; x++) {
				for (int y = 0; y < gridHeight; y++) {
					for (int i = 0; i < avgNet[1][1].length; i++) {
						fileNet.append(Double.toString(avgNet[x][y][i]));
						if (!(i+1==avgNet[1][1].length)) {
							fileNet.append(delim);
						}
					}
					fileNet.append(newline);
				}
			}
			fileNet.append(newline);
			try {
				fileNet.flush();
				fileMig.flush();
				fileProlif.flush();
				fileDep.flush();
				fileDeg.flush();
				fileCount.flush();
				fileNet.close();
				fileMig.close();
				fileProlif.close();
				fileDep.close();
				fileDeg.close();
				fileCount.close();
			} catch (IOException e) {
				System.out.println("Error while flushing/closing fileWriter !!!");
				e.printStackTrace();
			}
			
		} catch (Exception e) {
			System.out.println("Error in CsvFileWriter !!!");
			e.printStackTrace();
		}
		
	}
	
	/**
	 * Adding new fibroblasts via mitosis is done here so the fibroblast can be cataloged in the arraylist
	 * network state of parent cell is given to daughter cell
	 * Note that it does not check if the point is full/inhabited-this should be done by objects that call addFibroblast
	 * @param pt the grid point to which the fibroblast should be added
	 */
	public void addFibroblast(GridPoint pt, double[] network) {
		Fibroblast f = new Fibroblast(this, cellsPerGrid);
		this.add(f);
		f.setNetworkState(initialNet);
		//f.setNetworkState(network); //if you want daughter cells to inherit network state of the parent cell
		grid.moveTo(f, pt.getX(), pt.getY());
		f.initialize();
		fibroblasts.add(f);
	}
	
	/**
	 * remove a fibroblast from the array list and the world itself
	 * @param f
	 */
	public void removeFibroblast(Fibroblast f) {
		fibroblasts.remove(f); //will remove the first instance of this cell (assuming it does it correctly)
		this.remove(f); //remove from the context
	}
	
	
	/**
	 * Iterates through all the cells in the array list then creates a large array to pass to matlab to process all at once
	 */
	@ScheduledMethod(start = 1, interval=1, priority = 4)
	public void processCellBehavior() {
		
		if (fibroblasts.size() == 0) {
			return;
		}
		
		//get every fibroblast's network state
		double[][] states = new double[fibroblasts.size()][fibroblasts.get(0).getNetworkState().length];
		for (int i=0; i < fibroblasts.size(); i++) {
			states[i] = fibroblasts.get(i).getNetworkState();
		}
		try {
			eng.putVariable("states", states);
			eng.eval("processCellBehavior");
			states =  eng.getVariable("states"); //since states is a 2d array, there must be at least 2 fibroblasts
			GridPoint pt;
			Fibroblast f;
			
			
			for (int i=0; i < fibroblasts.size(); i++) {
				f = fibroblasts.get(i);
				f.setNetworkState(states[i]);
				pt = f.getPoint();
				GridValueLayer layer;
				
				double orig, dvdt;	
				int x,y; 
				for (int j=0; j < networkLayerInputIdxs.length - 1; j+=2) {
					x = pt.getX();
					y = pt.getY();
					layer = inputLayers.get(networkLayerInputIdxs[j+1]);
					orig = layer.get(x,y);
					dvdt = (states[i][networkLayerInputIdxs[j]]-orig)/chemokineFeedbackTimeConstant;
					layer.set(orig+dvdt, x,y);
				}
				
				for (int j=0; j < networkLayerOtherIdxs.length - 1; j+=2) {
					x = pt.getX();
					y = pt.getY();
					layer = otherExtracellularLayers.get(networkLayerOtherIdxs[j+1]);
					orig = layer.get(x,y);
					dvdt = (states[i][networkLayerOtherIdxs[j]]-orig)/chemokineFeedbackTimeConstant;
					layer.set(orig+dvdt, x,y);
				}
				
			}
		} catch (Exception e) {
			System.out.println("Error in Space-processCellBehavior");
			System.out.println(e);
		}
		
	}
	
	public double avColQ1() {
		double total = 0;
		double counter = 0;
		for (int x = gridWidth/2+1; x < gridWidth; x++) {
			for (int y = gridHeight/2+1; y < gridHeight; y++) {
				total += collagen.get(x,y);
				counter++;
			}
		}
		return total/counter;
	}

	public double avColQ2() {
		double total = 0;
		double counter = 0;
		for (int x = 0; x < gridWidth/2; x++) {
			for (int y = gridHeight/2+1; y < gridHeight; y++) {
				total += collagen.get(x,y);
				counter++;
			}
		}
		return total/counter;
	}

	public double avColQ3() {
		double total = 0;
		double counter = 0;
		for (int x = 0; x < gridWidth/2; x++) {
			for (int y = 0; y < gridHeight/2; y++) {
				total += collagen.get(x,y);
				counter++;
			}
		}
		return total/counter;
	}

	public double avColQ4() {
		double total = 0;
		double counter = 0;
		for (int x = gridWidth/2+1; x < gridWidth; x++) {
			for (int y = 0; y < gridHeight/2; y++) {
				total += collagen.get(x,y);
				counter++;
			}
		}
		return total/counter;
	}
}
