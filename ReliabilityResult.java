package com.net2plan.examples.general.offline.nfv;

import java.util.ArrayList;
import java.util.Collections;

//% 1 isolation, 2 nslices, 3 time, 4 tot_wl, 5 w_wl, 6 p_wl,
//% 7 total_bandwidth, 8 generated_traffic, 9 act_nodes, 10 load, 11 virtual links,
//% 12 avg_wl_per_link, 13 avg_wl_load,14 kmax

public class ReliabilityResult {
	public int isolation;
	public int nslices;
	public double time;
    public int tot_wl;
    public int w_wl;
    public int p_wl; 
    
    public double total_bandwidth;
    public double w_bandwidth;
    public double p_bandwidth;
    
    public double tot_generated_traffic;
    public double w_generated_traffic;
    public double p_generated_traffic;
    
    public int tot_connections;
    public int w_connections;
    public int p_connections;
    
    public int tot_act_nodes;
    public int tot_act_nodes_without_RUs;
    public int w_act_nodes;
    public int w_act_nodes_without_RUs;
    public int p_act_nodes;
    public int p_act_nodes_without_RUs;
    
    public int tot_functions;
    public int tot_functions_without_RUs;
    public int w_functions;
    public int w_functions_without_RUs;
    public int p_functions;
    public int p_functions_without_RUs;
    
    public double tot_load;
    public double tot_load_without_RUs;
    public double w_load;
    public double w_load_without_RUs;
    public double p_load;
    public double p_load_without_RUs;
    
    public int tot_grooming_ports;
    public int w_grooming_ports;
    public int p_grooming_ports;

    public ArrayList<Double> tot_load_nodes= new ArrayList<Double>();
    public ArrayList<Double> tot_load_nodes_noRUs= new ArrayList<Double>();
    public ArrayList<Double> w_load_nodes= new ArrayList<Double>();
    public ArrayList<Double> w_load_nodes_noRUs= new ArrayList<Double>();
    public ArrayList<Double> p_load_nodes= new ArrayList<Double>();
    public ArrayList<Double> p_load_nodes_noRUs= new ArrayList<Double>();
    
    
    public ReliabilityResult() {
    }
    
//    public ReliabilityResult(int isolation, int nslices, double time, int wlchann, int wwlchann, int pwlchann, double total_bandwidth, double generated_traffic,  int act_nodes, double load, int conn, int functions, double load_without_RUs, int act_nodes_without_RUs,int functions_without_RUs, ArrayList<Double> load_nodes,ArrayList<Double> load_nodes_noRUs ) {
//        super();
//        this.time = time;
//        this.tot_wl = wlchann;
//        this.w_wl= wwlchann;
//        this.p_wl= pwlchann;
//        this.total_bandwidth =total_bandwidth;
//        this.generated_traffic= generated_traffic;
//        this.isolation=isolation;
//        this.nslices=nslices;
//        this.act_nodes=act_nodes;
//        this.load=load;
//        this.connections=conn;
//        this.load_without_RUs=load_without_RUs;
//        this.functions=functions;
//        this.functions_without_RUs=functions_without_RUs;
//        this.act_nodes_without_RUs=act_nodes_without_RUs;
//        this.load_nodes=load_nodes;
//        this.load_nodes_noRUs=load_nodes_noRUs;
//        
//    }
    
    public String toString()
    {
    	
    	String s = String.format("%1$d,%2$d,%3$f,%4$d,%5$d,%6$d,%7$f,%8$f,%9$f,%10$f,%11$f, %12$f, %13$d, %14$d, %15$d, %16$d, %17$d, %18$d, %19$d, %20$d, %21$d, %22$d, %23$d, %24$d, %25$d, %26$d, %27$d, %28$f, %29$f, %30$f, %31$f, %32$f, %33$f, %34$d, %35$d, %36$d", isolation,nslices, time, tot_wl,w_wl, p_wl,total_bandwidth, w_bandwidth, p_bandwidth, tot_generated_traffic, w_generated_traffic, p_generated_traffic, tot_connections, w_connections,p_connections, tot_act_nodes, tot_act_nodes_without_RUs, w_act_nodes, w_act_nodes_without_RUs, p_act_nodes , p_act_nodes_without_RUs, tot_functions, tot_functions_without_RUs, w_functions, w_functions_without_RUs, p_functions, p_functions_without_RUs, tot_load, tot_load_without_RUs, w_load,w_load_without_RUs, p_load, p_load_without_RUs,tot_grooming_ports,w_grooming_ports,p_grooming_ports);
    	System.out.println(s);
    	return s;
    }
    
    public String tot_load_ValuesAsString()
    {
    	Collections.sort(this.tot_load_nodes);
    	return this.tot_load_nodes.toString().replace("[", "").replace("]", "").trim(); 
    }
    
    public String tot_load_noRUs_ValuesAsString()
    {
    	Collections.sort(this.tot_load_nodes_noRUs);
    	return this.tot_load_nodes_noRUs.toString().replace("[", "").replace("]", "").trim(); 
    }
    
    public String w_load_ValuesAsString()
    {
    	Collections.sort(this.w_load_nodes);
    	return this.w_load_nodes.toString().replace("[", "").replace("]", "").trim(); 
    }
    
    public String w_load_noRUs_ValuesAsString()
    {
    	Collections.sort(this.w_load_nodes_noRUs);
    	return this.w_load_nodes_noRUs.toString().replace("[", "").replace("]", "").trim(); 
    }
    
    public String p_load_ValuesAsString()
    {
    	Collections.sort(this.p_load_nodes);
    	return this.p_load_nodes.toString().replace("[", "").replace("]", "").trim(); 
    }
    
    public String p_load_noRUs_ValuesAsString()
    {
    	Collections.sort(this.p_load_nodes_noRUs);
    	return this.p_load_nodes_noRUs.toString().replace("[", "").replace("]", "").trim(); 
    }
    
    
    
}