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
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;

import com.mathworks.engine.MatlabEngine;

import repast.simphony.scenario.data.Classpath;

import repast.simphony.context.DefaultContext;
import repast.simphony.context.space.grid.GridFactoryFinder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ISchedule;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.parameter.Parameters;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridBuilderParameters;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.space.grid.RandomGridAdder;
import repast.simphony.space.grid.StrictBorders;
import repast.simphony.space.grid.WrapAroundBorders;
import repast.simphony.space.grid.BouncyBorders;
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
	public ArrayList<GridValueLayer> inputLayers = new ArrayList<GridValueLayer>();
	private ArrayList<ValueLayerDiffuser> inputDiffuseLayers = new ArrayList<ValueLayerDiffuser>();
	
	private GridValueLayer collagen;
	
	//ints alternate from network indices to indices in the inputLayerNames arraylist
	//for example, the 38th species in the network model is has idx 1 in inputLayerNames (TGFb)
	//remember that java indexes at 0
	
	private int[] networkLayerInputIdxs = {}; //no feedback
	private int[] networkLayerOutputIdxs = {}; //no feedback
	
	//These layers will be added to inputLayers
	public String[] inputLayerNames = {"TGFB", "LatentTGFB", "Interleukin6", "Interleukin1",
				"TNFalpha"};
	private int[] inputDiffuseIdxs = {}; //numerical indices of inputLayerNames that should diffuse. MUST BE IN INCREASING ORDER
	
	private double latentdegradationRate = 0.0; //constant degradation rate for latent TGFB
	private double activedegradationRate = (Double) p.getValue("activedegRate"); //constant degradation rate for active TGFB
	
	//indices of inputLayerNames that are inflammatory or anti-inflammatory/fibrotic
	//used for gradient orientation (e.g. the inflammatory cytokines are a gradient from left to right
	private int[] inflam = {2,3,4};
	private int[] antiInflam = {0,1};
	
	private MatlabEngine eng;
	
	//define saturating concentrations for each of the chemokines in order to calculate weights for the network model
	private double TGFBsat = 1;
	private double IL1sat = 1;
	private double IL6sat = 1;
	private double TNFasat = 1;
	
	//define max values that would saturate the receptor. This value will result in a weight of 1. Should be the same as the saturating values, but is used to initialize the value layer.
	private double[] inflamMax = {1, 1, 1}; //IL-6 max, IL-1 max, TNFa max - correspond to inflam indexes
	private double[] antiInflamMax = {1, 1}; //TGFB max, latent TGF-B max - correspond to antiInflam indexes

	
	//these variable are used for cytokines that change over time
	private ArrayList<Double> relChemVals; //the array of relative levels
	private int relChemFactorIdx;
	

	private int cellsPerGrid = 1;

	double[] initialNet = new double[91];


	
	/**
	 * Create all fields necessary for this class
	 */
	public AMFACSpace() {
		super("AMFACSpace");

		// Define the Grid Space
		grid = GridFactoryFinder.createGridFactory(null).createGrid("grid", this,
				new GridBuilderParameters<Object>(new BouncyBorders(), new RandomGridAdder<Object>(), false,
						gridWidth, gridHeight));

		// Build all GridValueLayers
		//initialize collagen layer at 3% area fraction
		collagen = new GridValueLayer("collagen", 3, false, gridWidth, gridHeight);
		this.addValueLayer(collagen);
		
		//create grid value layers
		GridValueLayer layer;
		for (int i=0; i < inputLayerNames.length; i++) {
			layer = new GridValueLayer(inputLayerNames[i], 1.0, false, gridWidth, gridHeight);
			inputLayers.add(layer);
			this.addValueLayer(layer);
		}
		

		Parameters p = RunEnvironment.getInstance().getParameters();
		boolean feedback = (Boolean) p.getValue("TGFB_feedback");
		
			if (feedback == true) {
				networkLayerInputIdxs = ArrayUtils.add(networkLayerInputIdxs, 19); //active TGF-B
				networkLayerInputIdxs = ArrayUtils.add(networkLayerInputIdxs, 0); //active TGF-B
				networkLayerOutputIdxs = ArrayUtils.add(networkLayerOutputIdxs, 23); //latent TGF-B
				networkLayerOutputIdxs = ArrayUtils.add(networkLayerOutputIdxs, 1); //latent TGF-B
			}
			
		System.out.println(Arrays.toString(networkLayerInputIdxs));
		System.out.println(Arrays.toString(networkLayerOutputIdxs));
		
	}
	
	public void createDiffusers() {
		
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
			} //else {
				//inputDiffuseLayers.add(null); //add null elements to fill in the gaps
			//}
		}
		
		//System.out.println("Create Diffusers");
	}
	
	@ScheduledMethod(start = 0, priority = 1)
	public void initialize() {
		loadMatlab();
		initializeFibroblasts();
		initializeChemokineLayer();
		initializeNetworkState();
		//createDiffusers();
		writeOutputData();
	}
	
	
	@ScheduledMethod(start = 1, interval = 1, priority = 2)
	public void goSecond() {
		processCellBehavior();
		//System.out.println("Second");
	}
	
	@ScheduledMethod(start = 1, interval = 1, priority = 0)
	public void goLast(){
		Parameters p = RunEnvironment.getInstance().getParameters();
		boolean feedback = (Boolean) p.getValue("TGFB_feedback");
		
			if (feedback == true){
				TGFBactivation();
			}
			
		writeOutputData();
		//System.out.println("Last");
	}
		
	public void loadMatlab() {
		//start the matlab connection
		try {
			eng = MatlabEngine.startMatlab();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		//System.out.println("Load Matlab");
	}
	
	public void initializeFibroblasts() {	
		//add fibroblasts
		Fibroblast fibroblast;
		for (int i = 0; i < initialFibroblastCount; i++) {
			fibroblast = new Fibroblast(this, cellsPerGrid);
			this.add(fibroblast);
			fibroblasts.add(fibroblast);
		}
		
		//System.out.println("Initialize Fibroblasts");
	}

	public void initializeChemokineLayer() {
		double wdub = (double) gridWidth;
		double hdub = (double) gridHeight;
		for (int x = 0; x < gridWidth; x++) {
			for (int y = 0; y < gridHeight; y++) {
				//inflammatory chemokines
				for (int i = 0; i < inflam.length; i++) {
					inputLayers.get(inflam[i]).set(((double)x)/wdub*inflamMax[i], x, y);
				}
				//fibrotic chemokines
				for (int i = 0; i < antiInflam.length; i++) {
					inputLayers.get(antiInflam[i]).set(((double)y)/hdub*antiInflamMax[i], x, y);
				}
			}
		}
		
		//System.out.println("Initialize Chemokines");

	}
	
	public void initializeNetworkState() {
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
		
		//System.out.println("Initialize Network State");
		
	}
	

	public void writeOutputData() {
		//collagen
		
		String fileName = "C:\\Users\\smr2we\\Documents\\collagen.csv";
		String delim = ",";
		String newline = "\n";
		FileWriter fileWriter;
		try {
			fileWriter = new FileWriter(fileName, true);
			for (int y = 0; y < gridWidth; y++) {
				for (int x = 0; x < gridHeight; x++) {
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
		fileName = "C:\\Users\\smr2we\\Documents\\tgfb.csv";
		try {
			fileWriter = new FileWriter(fileName, true);
			for (int y = 0; y < gridWidth; y++) {
				for (int x = 0; x < gridHeight; x++) {
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
		
		//latent tgfb
		fileName = "C:\\Users\\smr2we\\Documents\\LatentTgfb.csv";
		try {
			fileWriter = new FileWriter(fileName, true);
			for (int y = 0; y < gridWidth; y++) {
				for (int x = 0; x < gridHeight; x++) {
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
		
		//deposition (average of Col I and III mRNA)
		fileName = "C:\\Users\\smr2we\\Documents\\deposition.csv";
		try {
			fileWriter = new FileWriter(fileName, true);
			for (int y = 0; y < gridWidth; y++) {
				for (int x = 0; x < gridHeight; x++) {
					for( Object object: grid.getObjectsAt(x,y)) {
						
						Fibroblast f = (Fibroblast)object;
						double[] states = f.getNetworkState();
						double deposition = (states[87]+states[88])/2;
					

					fileWriter.append(Double.toString(deposition)); //index 1 is hardcoded in here
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
		
		//degradation (average of MMP 1, 2, and 9)
		fileName = "C:\\Users\\smr2we\\Documents\\degradation.csv";
		try {
			fileWriter = new FileWriter(fileName, true);
			for (int y = 0; y < gridWidth; y++) {
				for (int x = 0; x < gridHeight; x++) {
					for( Object object: grid.getObjectsAt(x,y)) {
						
						Fibroblast f = (Fibroblast)object;
						double[] states = f.getNetworkState();
						double degradation = (states[81]+states[82]+states[83])/3;
					

					fileWriter.append(Double.toString(degradation)); //index 1 is hardcoded in here
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
		
		//MMP1
		fileName = "C:\\Users\\smr2we\\Documents\\MMP1.csv";

		try {
			fileWriter = new FileWriter(fileName, true);
			for (int y = 0; y < gridWidth; y++) {
				for (int x = 0; x < gridHeight; x++) {
					for( Object object: grid.getObjectsAt(x,y)) {
						
						Fibroblast f = (Fibroblast)object;
						double[] states = f.getNetworkState();
						double degradation = (states[81]);
					

					fileWriter.append(Double.toString(degradation)); //index 1 is hardcoded in here
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
		
		//MMP2
		fileName = "C:\\Users\\smr2we\\Documents\\MMP2.csv";

		try {
			fileWriter = new FileWriter(fileName, true);
			for (int y = 0; y < gridWidth; y++) {
				for (int x = 0; x < gridHeight; x++) {
					for( Object object: grid.getObjectsAt(x,y)) {
						
						Fibroblast f = (Fibroblast)object;
						double[] states = f.getNetworkState();
						double degradation = (states[82]);
					

					fileWriter.append(Double.toString(degradation)); //index 1 is hardcoded in here
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
		
		//MMP9
		fileName = "C:\\Users\\smr2we\\Documents\\MMP9.csv";

		try {
			fileWriter = new FileWriter(fileName, true);
			for (int y = 0; y < gridWidth; y++) {
				for (int x = 0; x < gridHeight; x++) {
					for( Object object: grid.getObjectsAt(x,y)) {
						
						Fibroblast f = (Fibroblast)object;
						double[] states = f.getNetworkState();
						double degradation = (states[83]);
					

					fileWriter.append(Double.toString(degradation)); //index 1 is hardcoded in here
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
		
		//MMP14
		fileName = "C:\\Users\\smr2we\\Documents\\MMP14.csv";

		try {
			fileWriter = new FileWriter(fileName, true);
			for (int y = 0; y < gridWidth; y++) {
				for (int x = 0; x < gridHeight; x++) {
					for( Object object: grid.getObjectsAt(x,y)) {
						
						Fibroblast f = (Fibroblast)object;
						double[] states = f.getNetworkState();
						double degradation = (states[84]);
					

					fileWriter.append(Double.toString(degradation)); //index 1 is hardcoded in here
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
		
		//ColI
				fileName = "C:\\Users\\smr2we\\Documents\\ColI.csv";
				try {
					fileWriter = new FileWriter(fileName, true);
					for (int y = 0; y < gridWidth; y++) {
						for (int x = 0; x < gridHeight; x++) {
							for( Object object: grid.getObjectsAt(x,y)) {
								
								Fibroblast f = (Fibroblast)object;
								double[] states = f.getNetworkState();
								double deposition = (states[89]);
							

							fileWriter.append(Double.toString(deposition)); //index 1 is hardcoded in here
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
				
			//ColIII
				fileName = "C:\\Users\\smr2we\\Documents\\ColIII.csv";
				try {
					fileWriter = new FileWriter(fileName, true);
					for (int y = 0; y < gridWidth; y++) {
						for (int x = 0; x < gridHeight; x++) {
							for( Object object: grid.getObjectsAt(x,y)) {
								
								Fibroblast f = (Fibroblast)object;
								double[] states = f.getNetworkState();
								double deposition = (states[90]);
							

							fileWriter.append(Double.toString(deposition)); //index 1 is hardcoded in here
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
				
				//final network state
	/*			ISchedule schedule = RunEnvironment.getInstance().getCurrentSchedule();
				double tick = schedule.getTickCount();
				
				if (tick == 168){
					fileName = "C:\\Users\\smr2we\\Documents\\SteadyStateNetwork.csv";
					
					try {
						fileWriter = new FileWriter(fileName, true);
						for (int y = 0; y < gridWidth; y++) {
							for (int x = 0; x < gridHeight; x++) {
								for( Object object: grid.getObjectsAt(x,y)) {
									
									Fibroblast f = (Fibroblast)object;
									double[] states = f.getNetworkState();
								
									for (int i = 0; i < states.length; i++){
										fileWriter.append(Double.toString(states[i]));
										fileWriter.append(delim);
									}

								}
								fileWriter.append(newline);
						}
						}

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
				}*/
				
				//network state of all fibroblasts
				/*if (fibroblasts.size() != 0) {
				
				for (int i=0; i < fibroblasts.size(); i++) {

					fileName = "C:\\Users\\smr2we\\Documents\\Fibroblast" + Integer.toString(i) + ".csv";
					//System.out.println(fileName);
				
					
					try {
						fileWriter = new FileWriter(fileName, true);
						
						double[] networkstate = fibroblasts.get(i).getNetworkState();
						GridPoint pt = fibroblasts.get(i).getPoint();
						double y = pt.getY();
						double x = pt.getX();
						
								
						for (int j = 0; j < networkstate.length; j++){
								fileWriter.append(Double.toString(networkstate[j]));
								fileWriter.append(delim);
						}
						fileWriter.append(Double.toString(x));
						fileWriter.append(delim);
						fileWriter.append(Double.toString(y));
						fileWriter.append(delim);

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
				}*/
					
					
				
		//System.out.println("Write Output Data");
		
		
		
/*		//il6
		fileName = "C:\\Users\\smr2we\\Documents\\il6.csv";
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
		}*/
		
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
		double[] initial = new double[91];
		f.setNetworkState(initial);
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

	public void processCellBehavior() {
		if (fibroblasts.size() == 0) {
			return;
		}
		
		//get every fibroblast's network state
		
		double[][] states = new double[fibroblasts.size()][fibroblasts.get(0).getNetworkState().length];
		double[] TGFBweights = new double[fibroblasts.size()];
		double[] IL1weights = new double[fibroblasts.size()];
		double[] IL6weights = new double[fibroblasts.size()];
		double[] TNFaweights = new double[fibroblasts.size()];
		
		GridValueLayer TGFB = inputLayers.get(0); //hardcoded indexes, may need to be changed in the future
		GridValueLayer IL1 = inputLayers.get(3);
		GridValueLayer IL6 = inputLayers.get(2);
		GridValueLayer TNFa = inputLayers.get(4);
		
		for (int i=0; i < fibroblasts.size(); i++) {
			states[i] = fibroblasts.get(i).getNetworkState();
			Fibroblast f = fibroblasts.get(i);
			GridPoint pt = f.getPoint();
			double y = pt.getY();
			double x = pt.getY();
			
			TGFBweights[i] = TGFB.get(x,y)/TGFBsat;
			IL1weights[i] = IL1.get(x,y)/IL1sat;
			IL6weights[i] = IL6.get(x,y)/IL6sat;
			TNFaweights[i] = TNFa.get(x,y)/TNFasat;
			//TGFBweights[i] = y/gridHeight;
		}
		
			//System.out.println(Arrays.toString(IL6weights));
			
		try {
			eng.putVariable("states", states);
			eng.putVariable("TGFBweights", TGFBweights);
			eng.putVariable("IL1weights", IL1weights);
			eng.putVariable("IL6weights", IL6weights);
			eng.putVariable("TNFaweights", TNFaweights);
			
			//System.out.println(System.currentTimeMillis());
			eng.eval("processCellBehavior");
			//System.out.println(System.currentTimeMillis());
			
			states =  eng.getVariable("states"); //since states is a 2d array, there must be at least 2 fibroblasts
			
			GridPoint pt;
			Fibroblast f;
					
			for (int i=0; i < fibroblasts.size(); i++) {
				f = fibroblasts.get(i);
				f.setNetworkState(states[i]);
				pt = f.getPoint();
				GridValueLayer layer;
				
				double orig, dvdt, updated;	
				int x,y; 
				
				for (int j=0; j < networkLayerOutputIdxs.length - 1; j+=2) {
					x = pt.getX();
					y = pt.getY();
					layer = inputLayers.get(networkLayerOutputIdxs[j+1]);
					orig = layer.get(x,y);
					dvdt = (states[i][networkLayerOutputIdxs[j]]-orig);
					layer.set(orig+dvdt, x,y);
				}
				
				
			}
		} catch (Exception e) {
			System.out.println("Error in Space-processCellBehavior");
			System.out.println(e);
		}
		
		//System.out.println("Process Cell Behavior");
	}
	
	public void TGFBactivation() {
		for (int x=0; x < gridWidth; x++){
			for (int y=0; y < gridHeight; y++ ){
				for( Object object: grid.getObjectsAt(x,y)) {
				
			Fibroblast f = (Fibroblast)object;
					
			double[] states = f.getNetworkState();
			double activation = (states[82] + states[83])/2; //average of MMP2 and MMP9
			
			GridValueLayer TGFB;
			GridValueLayer latentTGFB;
			double origTGFB, origLatentTGFB;
			
			//get TGFB value layer at every grid point
			TGFB = inputLayers.get(networkLayerInputIdxs[1]);
			origTGFB = TGFB.get(x, y);
			
			//get latent TGFB value layer at every grid point
			latentTGFB = inputLayers.get(networkLayerOutputIdxs[1]);
			origLatentTGFB = latentTGFB.get(x, y);
			
			//activate TGFB based on activation rate and latent TGFB present
			TGFB.set(origTGFB + activation*origLatentTGFB, x, y);
			
			//remove the same amount of latent TGFB from that grid point
			latentTGFB.set(origLatentTGFB - activation*origLatentTGFB - latentdegradationRate*origLatentTGFB, x, y);
			
			//degrade TGFB at a constant rate
			TGFB.set(TGFB.get(x, y)*(1-activedegradationRate), x, y);
			
			}
		}
		}
		
		//System.out.println("TGFB Activation");
	}
}
