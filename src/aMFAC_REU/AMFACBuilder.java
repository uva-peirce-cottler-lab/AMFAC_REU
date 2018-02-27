/**
 * Builder used for Coupled Model
 * @author Tommy Athey
 * Aug 2017
 */
package aMFAC_REU;

import repast.simphony.context.Context;
import repast.simphony.dataLoader.ContextBuilder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.parameter.Parameters;

//Build context "AMFACBuilder"
public class AMFACBuilder implements ContextBuilder<Object> {
	public Context<Object> build(Context<Object> context) {
		Parameters p = RunEnvironment.getInstance().getParameters();
		int endTick = (Integer) p.getValue("endTick");

		// Build subcontext "AMFACSpace"
		AMFACSpace amfacSpace = new AMFACSpace();
		context.addSubContext(amfacSpace);
		context.add(amfacSpace);


		RunEnvironment.getInstance().endAt(endTick);
		return context;
	}
}
