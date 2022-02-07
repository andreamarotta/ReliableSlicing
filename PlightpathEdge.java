package com.net2plan.examples.general.offline.nfv;

import org.jgrapht.graph.DefaultWeightedEdge;

public class PlightpathEdge extends DefaultWeightedEdge {

	private static Integer current_id=0;
	public Integer id;
	
	

	public PlightpathEdge() {
		// TODO Auto-generated constructor stub
		id=current_id;
		current_id++;
		
	}

	public Integer getId() {
		return id;
	}
	
	
	
	

}
