/**
 * Used for chemokine colors (gradient of green to red)
 * @author Tommy Athey
 * Aug 2017
 */

package aMFAC_REU;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;

import repast.simphony.valueLayer.ValueLayer;
import repast.simphony.visualizationOGL2D.ValueLayerStyleOGL;

public class ChemokineStyle implements ValueLayerStyleOGL {

	private ValueLayer cglayer;
	Map<Integer, Color> colorMap;
	
	public ChemokineStyle() {
		colorMap = new HashMap<Integer,Color>();
		for (int i = 0; i < 100; i++) {
			double v = (double) i;
			colorMap.put(i, new Color((int) (255 * v /99.0), (int) (255*(1 - v /99.0)) , 0));
			
		}
	}
	
	@Override
	public Color getColor(double... coordinates) {
		double val = (double) cglayer.get(coordinates)*100;
		int roundVal = (int) val;
		Color color = colorMap.get(roundVal);
		
		if (roundVal >= 100) {
			color = Color.red;
		} else if (roundVal < 0) {
			color = Color.green;
		} else if (color==null) {
			color = Color.gray;
		}
		
		return color;
	}

	@Override
	public float getCellSize() {
		return 15.0f;
	}

	@Override
	public void init(ValueLayer layer) {
		this.cglayer = layer;
	}

}
