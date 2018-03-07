/**
 * Used for collagen color (gradient of black to red)
 * @author Tommy Athey
 * Aug 2017
 */
package aMFAC_REU;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;

import repast.simphony.valueLayer.ValueLayer;
import repast.simphony.visualizationOGL2D.ValueLayerStyleOGL;

public class CollagenStyle implements ValueLayerStyleOGL {

	private ValueLayer collagen;
	private Map<Integer, Color> colorMap = new HashMap<Integer, Color>();

	public CollagenStyle() {
		for (int i = 0; i < 100; i++) {
			double d = (double) i;
			int v = (int) (d*255.0/99.0);
			colorMap.put(i, new Color(v, 0, 0));
		}
	}

	@Override
	public Color getColor(double... coordinates) {
		int roundVal = (int) (collagen.get(coordinates)*100);
		Color color = colorMap.get(roundVal);
		//Color color = new Color(roundVal/100*255+100, 0, 0);
		
		if ( roundVal >= 100) {
			color = Color.red;
		} else if (roundVal < 0) {
			color = Color.black;
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
	public void init(ValueLayer collagen) {
		this.collagen = collagen;

	}
	
	

}
