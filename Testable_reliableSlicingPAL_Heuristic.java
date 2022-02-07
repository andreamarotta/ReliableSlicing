package com.net2plan.examples.general.offline.nfv;

import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import cern.colt.matrix.DoubleFactory3D;
//import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tint.IntMatrix1D;
import cern.colt.matrix.tint.IntMatrix2D;
import cern.colt.matrix.tint.IntFactory1D;
import cern.colt.matrix.tint.IntFactory2D;

import com.jom.DoubleMatrixND;
import com.jom.Expression;
import com.jom.OptimizationProblem;
import com.mxgraph.layout.mxCircleLayout;
import com.mxgraph.layout.mxIGraphLayout;
import com.mxgraph.layout.mxOrganicLayout;
import com.mxgraph.swing.mxGraphComponent;
import com.net2plan.examples.general.offline.Offline_ipOverWdm_routingSpectrumAndModulationAssignmentILPNotGrooming;
import com.net2plan.interfaces.networkDesign.*;
import com.net2plan.libraries.GraphUtils;
import com.net2plan.utils.Constants.RoutingType;
import com.net2plan.utils.InputParameter;
import com.net2plan.utils.Pair;
import com.net2plan.utils.Triple;
import com.sun.javafx.geom.Edge;

import org.apache.commons.collections15.map.HashedMap;
import org.apache.commons.lang3.StringUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import javax.swing.JFrame;

import java.sql.Timestamp;
import java.text.SimpleDateFormat;

import org.jgrapht.*;
import org.jgrapht.alg.DijkstraShortestPath;
import org.jgrapht.alg.*;
import org.jgrapht.alg.KShortestPaths;
import org.jgrapht.alg.interfaces.KShortestPathAlgorithm;
import org.jgrapht.alg.shortestpath.JohnsonShortestPaths;
import org.jgrapht.alg.shortestpath.KShortestSimplePaths;
import org.jgrapht.alg.shortestpath.EppsteinKShortestPath;
import org.jgrapht.alg.shortestpath.EppsteinShortestPathIterator;
import org.jgrapht.graph.*;
import org.jgrapht.alg.shortestpath.SuurballeKDisjointShortestPaths;
import org.jgrapht.alg.shortestpath.YenKShortestPath;
import org.jgrapht.alg.shortestpath.YenShortestPathIterator;
import org.jgrapht.ext.JGraphXAdapter;
import org.jgrapht.util.SupplierUtil;

//import org.jgrapht.io.*;
//import org.jgrapht.traverse.*

// ASSUMES YOU CAN NOT SPLIT ONE TRAFFIC REQUEST OVER DIFFERENT WAVELENGHTS
// DOES IT MAKES SENSE?

/**
 *
 * 
 * @net2plan.keywords JOM, NFV
 * @net2plan.inputParameters
 * @author Andrea Marotta
 */




public class Testable_reliableSlicingPAL_Heuristic {

	
	
	//graphedge-> list of <Link,wl,capacity>
	public class Lightpath extends  ArrayList<Triple<Link,Integer,Double>> {};
	public class PLightpath extends  ArrayList<ArrayList<Triple<Link,Integer,Double>>> {};
	public class PLightpathMap extends HashMap<DefaultWeightedEdge,PLightpath> {};
	
	private InputParameter C = new InputParameter("C", (double) 40, "Wavelength linerate");

	private InputParameter W = new InputParameter("W", (int) 40, "Wavelength per link");

	private InputParameter isolation = new InputParameter("isolation", (int) 0, "Slice isolation level");

	//private InputParameter SI = new InputParameter("SI", (int) 2, "Number of common segments");

	private InputParameter S = new InputParameter("S", (int) 1, "Number of slice requests to be generated");

	private InputParameter Cn = new InputParameter("Cn", (int) 4, "Node capacity");

	private InputParameter K = new InputParameter("K", (int) 1, "Cost of function placement");

	private InputParameter solverName = new InputParameter("solverName", "#select# glpk ipopt xpress cplex",
			"The solver name to be used by JOM. GLPK and IPOPT are free, XPRESS and CPLEX commercial. GLPK, XPRESS and CPLEX solve linear problems w/w.o integer contraints. IPOPT is can solve nonlinear problems (if convex, returns global optimum), but cannot handle integer constraints");
	private InputParameter solverLibraryName = new InputParameter("solverLibraryName", "",
			"The solver library full or relative path, to be used by JOM. Leave blank to use JOM default.");
	private InputParameter maxSolverTimeInSeconds = new InputParameter("maxSolverTimeInSeconds", (double) -1,
			"Maximum time granted to the solver to solve the problem. If this time expires, the solver returns the best solution found so far (if a feasible solution is found)");
	private NetworkLayer wdmLayer, ipLayer;

	private InputParameter virt= new InputParameter("Virtualization", "SF", "#SF #VN");
	
	private InputParameter OP= new InputParameter("Optical path", "VWP", "#VWP #WP");
	private InputParameter bandwidth= new InputParameter("Bandwidth", (int) 100, "System Bandwidth");
	private InputParameter mean= new InputParameter("Mean", (double) 2, "Average slice capacity");
	private InputParameter std= new InputParameter("Std", (double) 2, "Std slice capacity");
	private InputParameter seed= new InputParameter("seed", (int) 1, "Seed of simulation");
	private InputParameter grooming= new InputParameter("grooming", (boolean) true, "Grooming capable network");
	private InputParameter onlyplacement= new InputParameter("onlyplacement", (boolean) true, "Only placement simulation");
	
	
	private InputParameter networkName= new InputParameter("network", "Default", "Name of the network");
	private InputParameter caplim= new InputParameter("caplim", (int)0, "Capacity limits 0: no limit; 1: 500000 GOPS; 2: 250000 GOPS");
	
	private DoubleMatrixND M_sun;// = new DoubleMatrixND(new int[] { SI.getInt(), Nf, N });
	private DoubleMatrixND Omega_SUsu;
	private DoubleMatrixND b_suv;// = new DoubleMatrixND(new int[] { SI.getInt(), Nf, Nf });
	private DoubleMatrixND Y_sun;//
	private ArrayList<Double> Capn= new ArrayList<Double>();
	private ArrayList<Double> Capn_wRUs= new ArrayList<Double>();
	private ArrayList<Double> Capn_native= new ArrayList<Double>();
	private ArrayList<Integer> D;//= new ArrayList<Integer>();
	private ArrayList<Integer> R;//= new ArrayList<Integer>();
	private ArrayList <Pair<Integer,Integer>> slice_cseg;
	
	private Map<Integer, Pair<Integer,Double>> slice_Requests; 
	private DoubleMatrix2D K_su;//
	private int N;
	private int Nf;
	private int w_wl;
	private int p_wl;
	
	private int tot_grooming_ports=0;
	private int w_grooming_ports=0;
	private int p_grooming_ports=0;
	
	private double gentraff;
	private DoubleMatrix2D linkcap,linkcap_temp;
	private ArrayList<Integer> antennas= new ArrayList<Integer>();

	private int number_of_nodes_sfc;
	private ArrayList<Integer> core_capable= new ArrayList<Integer>();
	private Set<Integer> activeNodes= new HashSet<Integer>();
	private Set<Integer> activeNodes_withoutRUs= new HashSet<Integer>();
	private double load_withoutRUs=0;
	ArrayList<Triple<Pair<Integer,Integer>,PLightpath,Double>> existing_plightpaths = new ArrayList<Triple<Pair<Integer,Integer>,PLightpath,Double>>();

	private int kmax=0;
	private int K_sur_placement=100;
	
	
	private int plightpathcounter=0;
	

	private int Ns,Nc;

	private boolean DEBUG = false;

	
	
	
		
	public static void main(String[] args) {
		Locale.setDefault(Locale.US);
		String netname="Milan02";
		Map<String, String> parameters = new HashMap<>();
		parameters.put("network",netname);
		
		NetPlan netPlan = NetPlan.loadFromFile(new File(netname+".n2p"));
		//NetPlan netPlan = NetPlan.loadFromFile(new File("N6_2.n2p"));
		
		//NetPlan netPlan = NetPlan.loadFromFile(new File("N60.n2p"));
		//NetPlan netPlan = NetPlan.loadFromFile(new File("example7nodes_withTraffic.n2p"));
		//NetPlan netPlan = NetPlan.loadFromFile(new File("eon_N18_E66_withTraffic.n2p"));

		
		parameters.put("solverName", "cplex");
		parameters.put("solverLibraryName","C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio129\\opl\\bin\\x64_win64\\cplex1290.dll");
		parameters.put("maxSolverTimeInSeconds", String.valueOf((double) -1));
		//parameters.put("solverName", "glpk" );
		//parameters.put("solverLibraryName","C:\\WINDOWS\\system32\\glpk_4_48.dll");
		
		parameters.put("Virtualization","VN");
		
		//parameters.put("Virtualization","VN4");
		//parameters.put("Bandwidth", "100");;
		parameters.put("C", String.valueOf((double) 40));
		parameters.put("W", String.valueOf((int) 80));
		parameters.put("isolation", String.valueOf((int) 0));
		
		parameters.put("S", String.valueOf((int) 200));
		
		parameters.put("Cn", String.valueOf((int) 255));
		parameters.put("K", String.valueOf((int) 1));
		parameters.put("Bandwidth", "100");
		parameters.put("Mean", "2");
		parameters.put("Std", "1");
//		parameters.put("Bandwidth", "50");
//		parameters.put("Mean", "1");
//		parameters.put("Std", "0.5");
		parameters.put("Optical path", "VWP");
		
		parameters.put("grooming", String.valueOf((boolean) false));
		parameters.put("caplim", String.valueOf((int) 0));
		
		parameters.put("seed", "1");
		
//		double a=5.3234234/Double.POSITIVE_INFINITY;
//		double b=6.3453453453/Double.POSITIVE_INFINITY;
//
//		System.out.println("Infinity");
//		System.out.println(a);
//		System.out.println(b);
		
		Testable_reliableSlicingPAL_Heuristic Algorithm = new Testable_reliableSlicingPAL_Heuristic();
		Algorithm.getParameters();
		ReliabilityResult res=Algorithm.executeAlgorithm(netPlan, parameters, null);
		
		System.out.println(res.tot_load_ValuesAsString());
		System.out.println(res.toString());
		
		netPlan.saveToFile(new File("RESULT.n2p"));


	}


	public ReliabilityResult executeAlgorithm(NetPlan netPlan, Map<String, String> algorithmParameters,
			Map<String, String> net2planParameters) {

		/*
		 * Initialize all InputParameter objects defined in this object (this uses Java
		 * reflection)
		 */
		InputParameter.initializeAllInputParameterFieldsOfObject(this, algorithmParameters);

		final double PRECISION_FACTOR = 0.00001;

		/* Initialize variables */
		final int E = netPlan.getNumberOfLinks();
		final int N = netPlan.getNumberOfNodes();
		final int Nf=4;
		this.N=N;
		this.Nf=Nf;

		DoubleMatrix2D Fmn= DoubleFactory2D.sparse.make(N,N);
		for (Link e : netPlan.getLinks())
		{
			Fmn.set(e.getOriginNode().getIndex(), e.getDestinationNode().getIndex(), 1);
			Fmn.set( e.getDestinationNode().getIndex(),e.getOriginNode().getIndex(), 1);
		}


		SimpleDateFormat sdf = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss");

//
//		for (int i=0;i<10;i++)
//			generateAndSaveInput(netPlan, mean.getDouble(), std.getDouble(), mean.getDouble()*2, 1000,i);
		
		loadInput(netPlan, mean.getDouble(), std.getDouble(), S.getInt(),seed.getInt());
		createInput(netPlan, S.getInt(), bandwidth.getInt(), isolation.getInt(),grooming.getBoolean());

		



		Timestamp tic = new Timestamp(System.currentTimeMillis());

		ArrayList<Double> CapnOriginal= new ArrayList<Double>( this.Capn);




		//================================================================================
		// PRINT INPUT
		//================================================================================

		IntArrayList ind= new IntArrayList();
		IntArrayList ind2= new IntArrayList();

		DoubleArrayList val= new DoubleArrayList();

		int[] pos= new int[]{Ns,Nf,Nf};

		//b
		//System.out.println("");
		//System.out.println("b_suv: "+ String.valueOf(b_suv.zSum()));
		b_suv.getNonZeros(ind, val);
		for (int i=0; i< ind.size();i++)
		{	 
			//System.out.println("b_suv("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+")="+String.valueOf(val.get(i))); 
		}

		//M
		//System.out.println("");
		//System.out.println("M_sun: "+ String.valueOf(M_sun.zSum()));
		M_sun.getNonZeros(ind, val);
		pos= new int[]{Ns,Nf,N};
		for (int i=0; i< ind.size();i++)
		{	 
			//System.out.println("M_sun("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+")="+String.valueOf(val.get(i))); 
		}

		//K
		//System.out.println("");
		//System.out.println("K_su: "+ String.valueOf(K_su.zSum()));
		K_su.getNonZeros(ind,ind2, val);
		//pos= new int[]{Ns,Nf};
		for (int i=0; i< ind.size();i++)
		{	 
			//System.out.println("K_su("+String.valueOf(ind.get(i))+","+String.valueOf(ind2.get(i))+")="+String.valueOf(val.get(i))); 
		}


		
		System.out.println("//----------------NODE MAPPING----------------//");
//		if (virt.getString().contentEquals("VN"))
//		{ 
//			if (isolation.getInt()==2 || isolation.getInt()==4 || isolation.getInt()==3)
//			NodeMappingVN4( netPlan);
//			if (isolation.getInt()==0 || isolation.getInt()==1|| isolation.getInt()==5)
//			NodeMappingVN2( netPlan);
//				
//		}else {
//			NodeMappingSFC( netPlan);
//		}

		if (virt.getString().contentEquals("VN"))
			NodeMappingVNM_sur(netPlan);
		else
			NodeMappingSFCM_sur(netPlan);
		
		//else if (virt.getString().contentEquals("VN2"))
//			NodeMappingVN2(netPlan);
//		else if (virt.getString().contentEquals("VN3"))
//			NodeMappingVN3(netPlan);
//		else if (virt.getString().contentEquals("VN4"))
//			NodeMappingVN4(netPlan);
//		else if (virt.getString().contentEquals("SFC"))
//			NodeMappingSFC(netPlan);
//		else if (virt.getString().contentEquals("SFC2"))
//			NodeMappingSFC2(netPlan);
//		else if (virt.getString().contentEquals("SFC3"))
//			NodeMappingSFC3(netPlan);

		System.out.println("//----------------LINK MAPPING----------------//");
		
		//if (!onlyplacement.getBoolean())
		LinkMapping3(netPlan);
		//LinkMapping(netPlan);
				

		//Compute wavelenght channels
		int wl_channels=0;
		for(int i=0;i<E;i++)
			for (int j=0; j<W.getInt();j++)
				if (linkcap.get(i, j) < C.getDouble())
					{
						int n=netPlan.getLink(i).getOriginNode().getIndex();
						int m=netPlan.getLink(i).getDestinationNode().getIndex();
						//System.out.println("("+String.valueOf(n)+","+String.valueOf(m)+"):"+String.valueOf(j));
						//System.out.println("Used wavelength "+j+" on link ("+n+","+m+")");
						wl_channels++;
					}



		if (w_wl+p_wl != wl_channels)
		{
			System.out.println("Link Mapping Error");
			throw new IllegalArgumentException("Something wrong in link mapping");
			}
	
			
		


		//================================================================================
		// PRINT OUTPUT
		//================================================================================

		//y
		//System.out.println("");
		//System.out.println("y_sun: "+ String.valueOf(Y_sun.zSum()));
		Y_sun.getNonZeros(ind, val);
		pos= new int[]{Ns,Nf,N};
		int act[];    //declaring array
		act = new int[N];  // allocating memory to array
		int functions_without_RUs=0;
		int functions_RUs=0;
		Set<Integer> activeNodes_RUs= new HashSet<Integer>();
		Set<Integer> activeNodes_DUs= new HashSet<Integer>();
		Set<Integer> activeNodes_CUs= new HashSet<Integer>();
		Set<Integer> activeNodes_NGCs= new HashSet<Integer>();

		int n,u;

		
		for (int i=0; i< ind.size();i++)
		{	 
			//System.out.println("y_sun("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+")="+String.valueOf(val.get(i))); 
			n=(int)DoubleMatrixND.ind2sub(ind.get(i), pos).get(2);
			act[(int)DoubleMatrixND.ind2sub(ind.get(i), pos).get(2)]=1;
			
			u=DoubleMatrixND.ind2sub(ind.get(i), pos).get(1);
			if (DoubleMatrixND.ind2sub(ind.get(i), pos).get(1) == 0) {
				functions_RUs++;
				activeNodes_RUs.add(n);
			} else if (u == 1)
				activeNodes_DUs.add(n);
			else if (u == 2)
				activeNodes_CUs.add(n);
			else if (u == 3)
				activeNodes_NGCs.add(n);
		}
		if (isolation.getInt()==1)
			functions_RUs=functions_RUs-(Ns-Nc)
			;


		System.out.println("RUs");
		System.out.println(activeNodes_RUs.toString());
		System.out.println("DUs");
		System.out.println(activeNodes_DUs.toString());
		System.out.println("CUs");
		System.out.println(activeNodes_CUs.toString());
		System.out.println("NGCs");
		System.out.println(activeNodes_NGCs.toString());
		
		
		//z
		//System.out.println("");
		//System.out.println("z: "+ String.valueOf(wl_channels));



		//calculate active nodes
		int act_n=0;
		for (int i =0; i<N;i++)
		{
			if (act[i]==1)
				act_n++;
		}
		
		//if (act_n!=activeNodes.size()) //in active_nodes you have also
			//throw new IllegalArgumentException("Something wrong in active nodes calculation");
		
		int functions=0;
				if (isolation.getInt()!=0)
					functions= Y_sun.getNumberOfNonZeros()-(Ns-Nc); //Ns-Nc is number of dedicated slices
				else
					functions= Y_sun.getNumberOfNonZeros();
				
		functions_without_RUs=functions-functions_RUs;		
				
		//System.out.println("functions: "+String.valueOf(functions));
		//System.out.println("functions without RUs: "+String.valueOf(functions_without_RUs));
		//System.out.println("active nodes: "+String.valueOf(act_n));

		//load_withoutRUs+=10000*functions_without_RUs;
		load_withoutRUs += 10000.0 * activeNodes_withoutRUs.size();
		
		//calculate load
		double gops_with_RUs=0;
		for (int i =0; i<N;i++)
		{
			gops_with_RUs=gops_with_RUs+Capn_native.get(i)-Capn_wRUs.get(i);
			//System.out.println(String.valueOf(pu));
		}
		//pu+= functions*10000;
		gops_with_RUs+= act_n*10000;

		
		Timestamp tac = new Timestamp(System.currentTimeMillis());

		double diff= (tac.getTime()-tic.getTime())/1000;
		
		double tot_b= b_suv.zSum();
		
		int conn = b_suv.getNumberOfNonZeros();	

		ArrayList<Double> load= new ArrayList<Double>();
		ArrayList<Double> load_noRUs= new ArrayList<Double>();
		
		
		for (int i=0; i<N; i++)
			{
			load.add(Capn_native.get(i)-Capn_wRUs.get(i));
			load_noRUs.add(Capn_native.get(i)-Capn.get(i));
			
			}
		
		ReliabilityResult res= new ReliabilityResult();
		res.isolation=isolation.getInt();
		res.nslices=S.getInt();
		res.time=diff;
		res.tot_wl=wl_channels;
		res.w_wl=w_wl;
		res.p_wl=p_wl;
		res.total_bandwidth=tot_b;
		res.w_bandwidth=tot_b;
		res.p_bandwidth=0;
		
		res.tot_generated_traffic=gentraff;
		res.w_generated_traffic=gentraff;
		res.p_generated_traffic=0;
		
		res.tot_connections=conn;
		res.w_connections=conn;
		res.p_connections=0;
		
		res.tot_act_nodes=act_n;
		res.tot_act_nodes_without_RUs=activeNodes_withoutRUs.size();
		res.w_act_nodes=act_n;
		res.w_act_nodes=activeNodes_withoutRUs.size();
		res.p_act_nodes=0;
		res.p_act_nodes_without_RUs=0;
		
		res.tot_functions=functions;
		res.tot_functions_without_RUs=functions_without_RUs;
		res.w_functions=functions;
		res.w_functions_without_RUs=functions_without_RUs;
		res.p_functions=0;
		res.p_functions_without_RUs=0;
		
		res.tot_load=gops_with_RUs;
		res.tot_load_without_RUs=load_withoutRUs;
		res.w_load=gops_with_RUs;
		res.w_load_without_RUs=load_withoutRUs;
		res.p_load=0;
		res.p_load_without_RUs=0;
		
		res.tot_load_nodes=load;
		res.tot_load_nodes_noRUs=load_noRUs;
		res.w_load_nodes=load;
		res.w_load_nodes_noRUs=load_noRUs;
		//res.p_load_nodes=null;
		//res.p_load_nodes_noRUs=null;
		
		res.tot_grooming_ports=existing_plightpaths.size()*2;
		res.w_grooming_ports=existing_plightpaths.size();
		res.p_grooming_ports=existing_plightpaths.size();
		
		
		return res;

	}

	public ArrayList<Integer> finduv(int s)
	{
		ArrayList<Integer> uv= new ArrayList<Integer>();
		if(b_suv.get(new int[] {s,0,1})>0)
		{
			uv.add(0);
			uv.add(1);
		}
		if (b_suv.get(new int[] {s,1,2})>0)
		{
			uv.add(1);
			uv.add(2);
		}
		if (b_suv.get(new int[] {s,2,3})>0)
		{
			uv.add(2);
			uv.add(3);
		}

		Object[] aux=Arrays.stream(uv.toArray()).distinct().toArray();

		Integer[] in = Arrays.copyOf(aux, aux.length, Integer[].class);

		uv.clear();
		uv.addAll(Arrays.asList(in));



		return uv;
		//return uv;

	}

	public ArrayList<Integer> findCandidateLocation(int s, int u) {
		ArrayList<Integer> candi= new ArrayList<Integer>();
		boolean constrained =false;
		int is= isolation.getInt();
		if ((is==2 && u==1) ||(is==3 && u==2) ||(is==4 && u==2) )
			constrained = true;//only for first function in dedicated slices
		if (s>=Nc && isolation.getInt()<5 && isolation.getInt()>1 && constrained)
		{//dedicated slices
			//check if there are constraints
			int cindex=0;/// THIS NEEDS TO BE UPDATED
			
			for(int i=0;i<slice_cseg.size();i++)
			{
				//if (slice_cseg.get(i).getFirst()+Nc==s) //now we have it already with right indexes
					if (slice_cseg.get(i).getFirst()==s)
					cindex=slice_cseg.get(i).getSecond();
			}
			
			if (Omega_SUsu.get(new int[] {cindex,u,s,u})==1)
			{
				for (int n=0;n<N;n++) {
					if (Y_sun.get(new int[] {cindex,u,n})==1)
					{
						candi.add(n);
					}
				}
			}
			return candi;	
		}

		//common slices or non constrained functions
		for (int n=0; n<this.N;n++)
			if (M_sun.get(new int[] {s,u,n})==1)
				candi.add(n);

		//look for compatibles

		return candi;

	}

	public boolean isCapable(int slice,int function,int node)
	{
		boolean net= M_sun.get(new int[] {slice,function,node})==1.0;
		
		//System.out.println(M_sun.get(new int[] {slice,function,node}));
		boolean cap=true;//if it is RU capacity is not a proolem
		if (function!=0)
		cap= (this.Capn_wRUs.get(node)> K_su.get(slice, function) );
		//cap= (this.Capn.get(node)> K_su.get(slice, function) );
		
//		if (!cap)// && slice==123)
//			System.out.println("WARNING, node "+node+" is full");
		
//		if (!net)// && slice==123)
//			System.out.println("WARNING, node "+node+" not compatible with function "+ function);
		
		return (net && cap);
	}

	public ArrayList<Pair<Integer,Integer>> orderCandidatesByDCdistance(int slice, int function, NetPlan netplan, ArrayList<Integer> used)
	{
		ArrayList<Pair<Integer,Integer>> node_dist_dc =new ArrayList<Pair<Integer,Integer>>();
		//find DCs
		ArrayList<Integer> DCs=new ArrayList<Integer>();
		for (int i=0;i<N;i++)
		{
			for (int s=0;s<Ns;s++)
			if (isCapable(s, 3, i))
			{DCs.add(i);}
		}//make distinct
		
		SortedMap<Pair<Node,Node>,List<List<Link>>> cpl = netplan.computeUnicastCandidatePathList(null , 1, -1, -1, -1, -1, -1, -1 , null );
		
		for (int i=0;i<N;i++)
		{
			int distmin=100000;
			if ((isCapable(slice, function, i)) && !DCs.contains(i))
			{
				for (int d=0;d<DCs.size();d++)
				{
					Pair<Node,Node> nodePair;
					nodePair = Pair.of(netplan.getNode(i) , netplan.getNode(DCs.get(d)));
					int dist=cpl.get(nodePair).get(0).size();
					if (dist<distmin)
					{
						distmin=dist;
					}
				}
				
				node_dist_dc.add(Pair.of(i, distmin));
			}

		}

		//now we have to sort //check
		node_dist_dc.sort(new Comparator<Pair<Integer, Integer>>() {
			@Override
			public int compare(Pair<Integer, Integer> o1, Pair<Integer, Integer> o2) {
				if (o1.getSecond() < o2.getSecond()) {
					return -1;
				} else if (o1.getSecond().equals(o2.getSecond())) {
					return 0; // You can change this to make it then look at the
					//words alphabetical order
				} else {
					return 1;
				}
			}
		});
		
		return node_dist_dc;
			
	}

	public ArrayList<Pair<Integer,Integer>> computeBetweenness (int source, int destination, NetPlan netplan)
	{
		SortedMap<Pair<Node,Node>,List<List<Link>>> cpl = netplan.computeUnicastCandidatePathList(null , 1, -1, -1, -1, -1, -1, -1 , null );

		int d1,d2;
		Pair<Node,Node> nodePair;

		ArrayList<Pair<Integer,Integer>> bw = new ArrayList<Pair<Integer,Integer>>();

		for (int i=0; i<N;i++)
		{
			if (i!=source && i!=destination)
			{
				nodePair = Pair.of(netplan.getNode(i) , netplan.getNode(source));
				d1=cpl.get(nodePair).get(0).size();
	
				nodePair = Pair.of(netplan.getNode(i) , netplan.getNode(destination));
				d2=cpl.get(nodePair).get(0).size();
	
				bw.add(Pair.of(i,d1+d2));
			}

		}

		bw.sort(new Comparator<Pair<Integer, Integer>>() {
			@Override
			public int compare(Pair<Integer, Integer> o1, Pair<Integer, Integer> o2) {
				if (o1.getSecond() < o2.getSecond()) {
					return -1;
				} else if (o1.getSecond().equals(o2.getSecond())) {
					return 0; // You can change this to make it then look at the
					//words alphabetical order
				} else {
					return 1;
				}
			}
		});


		return bw;

	}


	public void loadInput(NetPlan netPlan, double mean, double std, int size, int seed)
	{

		this.slice_Requests= new HashMap<Integer, Pair<Integer,Double>>();
		//Load input for comparisons
		try {
			//String filename=  ("./InputN"+String.valueOf(N)+"B"+String.valueOf(mean)+"S"+String.valueOf(std)+"seed"+String.valueOf(seed));
			
			String filename=  (networkName.getString()+"B"+String.valueOf(mean)+"S"+String.valueOf(std)+"seed"+String.valueOf(seed));
			
			//System.out.println(filename);
			
			FileInputStream fi = new FileInputStream(new File(filename));
			ObjectInputStream oi = new ObjectInputStream(fi);
			
			int id;
			int antenna;
			double b;
			
			for (int i=0;i<size;i++)
			{
				id=(int) oi.readObject();
				antenna=(int) oi.readObject();
				b= (double) oi.readObject();
				slice_Requests.put(id, Pair.of(antenna, b));
			}
			//System.out.println(slice_Requests.toString());
			oi.close();
			fi.close();

		} catch (FileNotFoundException e) {
			System.out.println("File not found");
		} catch (IOException e) {
			System.out.println("Error initializing stream");
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	

	public void generateAndSaveInput(NetPlan netPlan, double mean, double std, double max, int size,int seed)
	{
		//================================================================================
		// INPUT GENERATION
		//================================================================================
		
		int N= netPlan.getNodes().size();
		
		//Store indexes of RUs,DUs,CUs,NGCs
				//we will use for slice requests generation
				ArrayList<Integer> RUs = new ArrayList<Integer>(); 


		
		
		for (int i = 0; i < N; i++) 
		{
			if (netPlan.getNode(i).getAttribute("RU").toString().equals("1")) 
			{
				RUs.add(i); 
			}
		}
		
		try {
		FileOutputStream fileOut = new FileOutputStream(networkName.getString()+"B"+String.valueOf(mean)+"S"+String.valueOf(std)+"seed"+String.valueOf(seed));
		ObjectOutputStream objectOut = new ObjectOutputStream(fileOut);		
		Random r = new Random();
		double b= 0;
		int antenna=0;
		
		for (int i=0;i<size;i++)
		{
			boolean ok =false;
			
			while (!ok)
			{
				b= r.nextGaussian()*std+mean;
				if (b>0 && b< max)
					ok=true;
			}
			
			//b=mean;
			
			antenna=RUs.get(r.nextInt(RUs.size()));
			objectOut.writeObject(i);
			objectOut.writeObject(antenna);
			objectOut.writeObject(b);

		}
		}catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	
	public void createInput(NetPlan netPlan, int nslices, int system_bandwidth,int isolation, boolean grooming)
	{
		
		//Store indexes of RUs,DUs,CUs,NGCs
				//we will use for slice requests generation
				ArrayList<Integer> RUs = new ArrayList<Integer>(); 
				ArrayList<Integer> DUs = new ArrayList<Integer>(); 
				ArrayList<Integer> CUs = new ArrayList<Integer>(); 
				ArrayList<Integer> NGCs = new ArrayList<Integer>(); 

				//for (int i=0;i<N;i++)
				//this.Capn.add((double) Cn.getInt());
				double capacity_limit;
				switch (caplim.getInt()){
				case 0:
					capacity_limit=1000000000;
					break;
				case 1:
					capacity_limit=1500000;
					break;
				case 2:
					capacity_limit=500000;
					break;
				default:
					capacity_limit=1000000000;
						
				}
				
				
				for (int i = 0; i < N; i++) {
					this.Capn.add(0.0);

					try {

						if (netPlan.getNode(i).getAttribute("RU").toString().equals("1")) 
						{
							RUs.add(i);
//							if (caplim.getBoolean()==false)
//							this.Capn.set(i, 1000000000.0); 
//							else
								this.Capn.set(i, capacity_limit); 
							antennas.add(i);
							//activeNodes.add(i);
						}

						if (netPlan.getNode(i).getAttribute("DU").toString().equals("1")) 
						{DUs.add(i);
//						if (caplim.getBoolean()==false)
//							this.Capn.set(i, 1000000000.0); 
//							else
								this.Capn.set(i, capacity_limit); 
						
						}

						if (netPlan.getNode(i).getAttribute("CU").toString().equals("1")) 
						{CUs.add(i);
						
//						if (caplim.getBoolean()==false)
//							this.Capn.set(i, 1000000000.0); 
//							else
								this.Capn.set(i, capacity_limit); 
						}

						if (netPlan.getNode(i).getAttribute("NGC").toString().equals("1")) 
						{NGCs.add(i); 
						
//						if (caplim.getBoolean()==false)
//							this.Capn.set(i, 1000000000.0); 
//							else
								this.Capn.set(i, capacity_limit); 
						
						
						core_capable.add(i);
						}
					}
					catch (NullPointerException e) {
						if (Math.random() < 1)
						{RUs.add(i);this.Capn.set(i, 1000000000.0); 
						antennas.add(i);}
						if (Math.random() < 0.5)
						{DUs.add(i);this.Capn.set(i, 1000000000.0); }
						if (Math.random() < 0.5)
						{CUs.add(i); this.Capn.set(i, 1000000000.0); }
						if (Math.random() < 0.5)
						{NGCs.add(i);this.Capn.set(i, 1000000000.0); 
						core_capable.add(i);
						}

					}


				}
				
		Capn_native= new ArrayList<Double>(Capn);
		Capn_wRUs= new ArrayList<Double>(Capn);
		
				
		double RU_cost=0;
		//if (caplim.getBoolean()==false)
			RU_cost=57600.0;
		double DU_cost=89200.0;
		double CU_cost=3200.0;
		double NGC_cost=3200.0;
		
		double max_b=0;
		double fronthaul=0;
		double midhaul=0;
		double backhaul=0;
		switch (system_bandwidth)
		{
		case 200:
			max_b=8;
			fronthaul=44.4;
			midhaul=8.032;
			backhaul=8;
			break;
		case 100:
			max_b=4;
			fronthaul=22.2;
			midhaul=4.016;
			backhaul=4;
			break;
		case 50:
			max_b=2;
			fronthaul=11.1;
			midhaul=2.008;
			backhaul=2;
			break;
		
		}

		
		DoubleMatrixND M_sun = new DoubleMatrixND(new int[] { nslices, Nf, N });
		DoubleMatrixND b_suv = new DoubleMatrixND(new int[] { nslices, Nf, Nf });
		ArrayList<Integer> D= new ArrayList<Integer>();
		ArrayList<Integer> R= new ArrayList<Integer>();
		DoubleMatrix2D K_su=  DoubleFactory2D.sparse.make(nslices, Nf);
		ArrayList <Triple< Integer, Integer,Double>> cseg=new ArrayList<Triple<Integer,Integer,Double>>(); //List<Id,Antenna,capacity>>
		slice_cseg = new ArrayList<Pair<Integer,Integer>>();
		
		for (int i=0;i<nslices;i++)
		{
//			if (i==123)
//				System.out.println("stop");
			boolean foundcseg;
			double slice_weight= slice_Requests.get(i).getSecond()/backhaul;
			switch (isolation)
			{
				case 0:
					//set deployability
					M_sun.set(new int[] { i, 0, slice_Requests.get(i).getFirst()}, 1 );
					for (int j =0 ; j< DUs.size();j++)
					{
						M_sun.set(new int[] { i, 1, DUs.get(j)}, 1 );
					}
					for (int j =0 ; j< CUs.size();j++)
					{
						M_sun.set(new int[] { i, 2, CUs.get(j)}, 1 );
						//System.out.println("M_sun("+i+",2,"+CUs.get(j)+")");
					}
					for (int j =0 ; j< NGCs.size();j++)
					{
						M_sun.set(new int[] { i, 3, NGCs.get(j)}, 1 );
					}
					K_su.set(i, 0, RU_cost*slice_weight);
					K_su.set(i, 1, DU_cost*slice_weight);
					K_su.set(i, 2, CU_cost*slice_weight);
					K_su.set(i, 3, NGC_cost*slice_weight);
					
					
					
					//set 3 links bandwidth
					b_suv.set(new int[] { i, 0, 1 }, fronthaul*slice_weight);
					//b_suv.set(new int[] { i, 1, 0 }, 20);//biderctionality
					b_suv.set(new int[] { i, 1, 2 }, midhaul*slice_weight);
					//b_suv.set(new int[] { i, 2, 1 }, 20);//biderctionality
					b_suv.set(new int[] { i, 2, 3 }, backhaul*slice_weight);
					
					//dedicated transport is only on the non common part
					if(grooming) D.add(0); else D.add(1);

					//reliability?
					R.add(1);
				break;
				case 1:
					M_sun.set(new int[] { i, 0, slice_Requests.get(i).getFirst()}, 1 );
					//set deployability
					for (int j =0 ; j< DUs.size();j++)
					{
						M_sun.set(new int[] { i, 1, DUs.get(j)}, 1 );
					}
					for (int j =0 ; j< CUs.size();j++)
					{
						M_sun.set(new int[] { i, 2, CUs.get(j)}, 1 );
					}
					for (int j =0 ; j< NGCs.size();j++)
					{
						M_sun.set(new int[] { i, 3, NGCs.get(j)}, 1 );
					}
					
					//look for cseg
					foundcseg= false;
									
					for (int c=0; c<cseg.size() && !foundcseg;c++)
					{
						if (cseg.get(c).getSecond()==slice_Requests.get(i).getFirst() && cseg.get(c).getThird()>= slice_Requests.get(i).getSecond())
						{
							foundcseg=true;
							cseg.get(c).setThird(cseg.get(c).getThird()-slice_Requests.get(i).getSecond());
							//add bsuv
							//set 3 links bandwidth
							b_suv.set(new int[] { i, 0, 1 }, fronthaul*slice_weight);
							//b_suv.set(new int[] { i, 1, 0 }, 20);//biderctionality
							b_suv.set(new int[] { i, 1, 2 }, midhaul*slice_weight);
							//b_suv.set(new int[] { i, 2, 1 }, 20);//biderctionality
							b_suv.set(new int[] { i, 2, 3 }, backhaul*slice_weight);
							
							//ad ksu
							K_su.set(i, 0, RU_cost*slice_weight); //the cost is added only first time
							K_su.set(i, 1, DU_cost*slice_weight);
							K_su.set(i, 2, CU_cost*slice_weight);
							K_su.set(i, 3, NGC_cost*slice_weight);
						}
					
					}
					
					if (!foundcseg)
					{
						//create new cseg
						cseg.add(Triple.of(cseg.size(),slice_Requests.get(i).getFirst(),max_b-slice_Requests.get(i).getSecond() ));
						
						b_suv.set(new int[] { i, 0, 1 }, fronthaul*slice_weight);
						//b_suv.set(new int[] { i, 1, 0 }, 20);//biderctionality
						b_suv.set(new int[] { i, 1, 2 }, midhaul*slice_weight);
						//b_suv.set(new int[] { i, 2, 1 }, 20);//biderctionality
						b_suv.set(new int[] { i, 2, 3 }, backhaul*slice_weight);
						
						//ad ksu
						K_su.set(i, 0, RU_cost*slice_weight); //the cost is added only first time
						K_su.set(i, 1, DU_cost*slice_weight);
						K_su.set(i, 2, CU_cost*slice_weight);
						K_su.set(i, 3, NGC_cost*slice_weight);
					}
					
					//dedicated transport is only on the non common part
					if(grooming) D.add(0); else D.add(1);

					//reliability?
					R.add(1);
					
				break;
				case 2:
					
					//set deployability basing on network
					for (int j =0 ; j< DUs.size();j++)
					{
						M_sun.set(new int[] { i, 1, DUs.get(j)}, 1 );
					}

					for (int j =0 ; j< CUs.size();j++)
					{
						M_sun.set(new int[] { i, 2, CUs.get(j)}, 1 );
					}
					for (int j =0 ; j< NGCs.size();j++)
					{
						M_sun.set(new int[] { i, 3, NGCs.get(j)}, 1 );
					}
					
					//create or use common segment
					//look for cseg
					foundcseg= false;
									
					for (int c=0; c<cseg.size() && !foundcseg;c++)
					{
						if (cseg.get(c).getSecond()==slice_Requests.get(i).getFirst() && cseg.get(c).getThird()>= slice_Requests.get(i).getSecond())
						{
							foundcseg=true;
							cseg.get(c).setThird(cseg.get(c).getThird()-slice_Requests.get(i).getSecond());
							//add bsuv
							//set 3 links bandwidth
							//b_suv.set(new int[] { i, 0, 1 }, fronthaul*slice_weight);
							//b_suv.set(new int[] { i, 1, 0 }, 20);//biderctionality
							b_suv.set(new int[] { i, 1, 2 }, midhaul*slice_weight);
							//b_suv.set(new int[] { i, 2, 1 }, 20);//biderctionality
							b_suv.set(new int[] { i, 2, 3 }, backhaul*slice_weight);
							
							//ad ksu
							//K_su.set(i, 0, 0); 
							K_su.set(i, 1, 0);//the cost is at the common segment
							K_su.set(i, 2, CU_cost*slice_weight);
							K_su.set(i, 3, NGC_cost*slice_weight);
							
							slice_cseg.add(Pair.of(i, cseg.get(c).getFirst()));
						}
					
					}
					
					if (!foundcseg)
					{
						//create new cseg
						int new_seg= cseg.size();
						cseg.add(Triple.of(new_seg,slice_Requests.get(i).getFirst(),max_b-slice_Requests.get(i).getSecond() ));
						
						slice_cseg.add(Pair.of(i, new_seg));
						
						//b_suv.set(new int[] { i, 0, 1 }, fronthaul*slice_weight);
						//b_suv.set(new int[] { i, 1, 0 }, 20);//biderctionality
						b_suv.set(new int[] { i, 1, 2 }, midhaul*slice_weight);
						//b_suv.set(new int[] { i, 2, 1 }, 20);//biderctionality
						b_suv.set(new int[] { i, 2, 3 }, backhaul*slice_weight);
						
						//ad ksu
						//K_su.set(i, 0, 1); 
						K_su.set(i, 1, 0); //the cost is at the common segment
						K_su.set(i, 2, CU_cost*slice_weight);
						K_su.set(i, 3, NGC_cost*slice_weight);
					}
					
					
				
					
					
					
					//dedicated transport is only on the non common part
					if(grooming) D.add(0); else D.add(1);

					//reliability?
					R.add(1);
					
				break;
				case 3:
					
					//set deployability basing on network

					for (int j =0 ; j< CUs.size();j++)
					{
						M_sun.set(new int[] { i, 2, CUs.get(j)}, 1 );
					}

					for (int j =0 ; j< NGCs.size();j++)
					{
						M_sun.set(new int[] { i, 3, NGCs.get(j)}, 1 );
					}
					
					//create or use common segment
					//look for cseg
					foundcseg= false;
									
					for (int c=0; c<cseg.size() && !foundcseg;c++)
					{
						if (cseg.get(c).getSecond()==slice_Requests.get(i).getFirst() && cseg.get(c).getThird()>= slice_Requests.get(i).getSecond())
						{
							foundcseg=true;
							if (cseg.get(c).getThird()-slice_Requests.get(i).getSecond()<0)
								System.out.println("OUT!");
							cseg.get(c).setThird(cseg.get(c).getThird()-slice_Requests.get(i).getSecond());
							//add bsuv
							//set 3 links bandwidth
							//b_suv.set(new int[] { i, 0, 1 }, fronthaul*slice_weight);
							//b_suv.set(new int[] { i, 1, 0 }, 20);//biderctionality
							//b_suv.set(new int[] { i, 1, 2 }, midhaul*slice_weight);
							//b_suv.set(new int[] { i, 2, 1 }, 20);//biderctionality
							b_suv.set(new int[] { i, 2, 3 }, backhaul*slice_weight);
							
							//ad ksu
							//K_su.set(i, 0, 0); 
							//K_su.set(i, 1, 0);//the cost is added only first time
							K_su.set(i, 2, 0); //the cost is at the common segment
							K_su.set(i, 3, NGC_cost*slice_weight);
							
							slice_cseg.add(Pair.of(i, cseg.get(c).getFirst()));
						}
					
					}
					
					if (!foundcseg)
					{
						//create new cseg
						int new_seg= cseg.size();
						if (max_b-slice_Requests.get(i).getSecond()<0)
							System.out.println("OUT!");
						cseg.add(Triple.of(new_seg,slice_Requests.get(i).getFirst(),max_b-slice_Requests.get(i).getSecond() ));
						
						slice_cseg.add(Pair.of(i, new_seg));
						
						//b_suv.set(new int[] { i, 0, 1 }, fronthaul*slice_weight);
						//b_suv.set(new int[] { i, 1, 0 }, 20);//biderctionality
						//b_suv.set(new int[] { i, 1, 2 }, midhaul*slice_weight);
						//b_suv.set(new int[] { i, 2, 1 }, 20);//biderctionality
						b_suv.set(new int[] { i, 2, 3 }, backhaul*slice_weight);
						
						//ad ksu
						//K_su.set(i, 0, 1); 
						//K_su.set(i, 1, DU_cost);
						K_su.set(i, 2, 0); //the cost is at the common segment
						K_su.set(i, 3, NGC_cost*slice_weight);
					}
					
					
				
					
					
					
					//dedicated transport is only on the non common part
					D.add(1);

					//reliability?
					R.add(1);
					
					
					break;
					
					
					case 4:
					
					//set deployability basing on network

					for (int j =0 ; j< CUs.size();j++)
					{
						M_sun.set(new int[] { i, 2, CUs.get(j)}, 1 );
					}

					for (int j =0 ; j< NGCs.size();j++)
					{
						M_sun.set(new int[] { i, 3, NGCs.get(j)}, 1 );
					}
					
					//create or use common segment
					//look for cseg
					foundcseg= false;
									
					for (int c=0; c<cseg.size() && !foundcseg;c++)
					{
						if (cseg.get(c).getSecond()==slice_Requests.get(i).getFirst() && cseg.get(c).getThird()>= slice_Requests.get(i).getSecond())
						{
							foundcseg=true;
							cseg.get(c).setThird(cseg.get(c).getThird()-slice_Requests.get(i).getSecond());
							//add bsuv
							//set 3 links bandwidth
							//b_suv.set(new int[] { i, 0, 1 }, fronthaul*slice_weight);
							//b_suv.set(new int[] { i, 1, 0 }, 20);//biderctionality
							//b_suv.set(new int[] { i, 1, 2 }, midhaul*slice_weight);
							//b_suv.set(new int[] { i, 2, 1 }, 20);//biderctionality
							b_suv.set(new int[] { i, 2, 3 }, backhaul*slice_weight);
							
							//ad ksu
							//K_su.set(i, 0, 0); 
							//K_su.set(i, 1, 0);//the cost is added only first time
							K_su.set(i, 2, 0); //the cost is at the common segment
							K_su.set(i, 3, NGC_cost*slice_weight);
							
							slice_cseg.add(Pair.of(i, cseg.get(c).getFirst()));
						}
					
					}
					
					if (!foundcseg)
					{
						//create new cseg
						int new_seg= cseg.size();
						cseg.add(Triple.of(new_seg,slice_Requests.get(i).getFirst(),max_b-slice_Requests.get(i).getSecond() ));
						
						slice_cseg.add(Pair.of(i, new_seg));
						
						//b_suv.set(new int[] { i, 0, 1 }, fronthaul*slice_weight);
						//b_suv.set(new int[] { i, 1, 0 }, 20);//biderctionality
						//b_suv.set(new int[] { i, 1, 2 }, midhaul*slice_weight);
						//b_suv.set(new int[] { i, 2, 1 }, 20);//biderctionality
						b_suv.set(new int[] { i, 2, 3 }, backhaul*slice_weight);
						
						//ad ksu
						//K_su.set(i, 0, 1); //the cost is added only first time
						//K_su.set(i, 1, DU_cost);
						K_su.set(i, 2, 0); //the cost is at the common segment
						K_su.set(i, 3, NGC_cost*slice_weight);
					}
					
					
					//dedicated transport is only on the non common part
					if(grooming) D.add(0); else D.add(1);

					//reliability?
					R.add(1);
					
					
					break;
					case 5:
						//NO TRAFFIC IN THE ISOLATED SEGMENT
						
						//create or use common segment
						//look for cseg
						foundcseg= false;
										
						for (int c=0; c<cseg.size() && !foundcseg;c++)
						{
							if (cseg.get(c).getSecond()==slice_Requests.get(i).getFirst() && cseg.get(c).getThird()>= slice_Requests.get(i).getSecond())
							{
								foundcseg=true;
								cseg.get(c).setThird(cseg.get(c).getThird()-slice_Requests.get(i).getSecond());
								slice_cseg.add(Pair.of(i, cseg.get(c).getFirst()));
							}
						
						}
						
						if (!foundcseg)
						{
							//create new cseg
							int new_seg= cseg.size();
							cseg.add(Triple.of(new_seg,slice_Requests.get(i).getFirst(),max_b-slice_Requests.get(i).getSecond() ));
							
							slice_cseg.add(Pair.of(i, new_seg));
						}
						
						
						//dedicated transport is only on the non common part
						D.add(0);

						//reliability?
						R.add(1);
						
					break;
			}
	
		}
		
		//now that we processed slices we should process common segments and create input datastructures
		this.Nc= cseg.size();
		switch (isolation) {
		case 0:
		case 1:
			Ns=nslices;
			break;
		case 5:
			Ns=Nc;
			break;
		default:
			Ns =nslices+Nc;
			break;
		}
		
		
		//COMMON SEGMENT
		if (isolation !=0)
		{
		DoubleMatrixND Mc_sun = new DoubleMatrixND(new int[] { Nc, Nf, N });
		DoubleMatrixND bc_suv = new DoubleMatrixND(new int[] { Nc, Nf, Nf });
		ArrayList<Integer> cD= new ArrayList<Integer>();
		ArrayList<Integer> cR= new ArrayList<Integer>();
		DoubleMatrix2D Kc_su=  DoubleFactory2D.sparse.make(Nc, Nf);

		
		switch(isolation)
		{
		case 0:
		case 1:
			break;
		case 2:
			
			for (int cs =0; cs<cseg.size();cs++)
			{
				double slice_weight= (max_b-cseg.get(cs).getThird())/max_b;
				
				bc_suv.set(new int[] {cseg.get(cs).getFirst(),0,1}, (max_b-cseg.get(cs).getThird())/max_b *fronthaul);
				
				Mc_sun.set(new int[] {cseg.get(cs).getFirst(),0,cseg.get(cs).getSecond() }, 1);
				
				for (int j =0 ; j< DUs.size();j++)
				{
					Mc_sun.set(new int[] { cseg.get(cs).getFirst(), 1, DUs.get(j)}, 1 );
				}
				
				Kc_su.set(cseg.get(cs).getFirst(), 0, RU_cost*slice_weight);
				Kc_su.set(cseg.get(cs).getFirst(), 1, DU_cost*slice_weight);
				//dedicated transport is only on the non common part
				if(grooming) cD.add(0); else cD.add(1);
				

				//reliability?
				cR.add(1);

			}
			
			break;
		case 3:
		case 4:
			for (int cs =0; cs<cseg.size();cs++)
			{
				double slice_weight= (max_b-cseg.get(cs).getThird())/max_b;
				bc_suv.set(new int[] {cseg.get(cs).getFirst(),0,1}, (max_b-cseg.get(cs).getThird())/max_b *fronthaul);
				bc_suv.set(new int[] {cseg.get(cs).getFirst(),1,2}, (max_b-cseg.get(cs).getThird())/max_b *midhaul);
				
				Mc_sun.set(new int[] {cseg.get(cs).getFirst(),0,cseg.get(cs).getSecond() }, 1);
				
				for (int j =0 ; j< DUs.size();j++)
				{
					Mc_sun.set(new int[] { cseg.get(cs).getFirst(), 1, DUs.get(j)}, 1 );
				}
				
				for (int j =0 ; j< CUs.size();j++)
				{
					Mc_sun.set(new int[] { cseg.get(cs).getFirst(), 2, CUs.get(j)}, 1 );
				}
				
				
				Kc_su.set(cseg.get(cs).getFirst(), 0, RU_cost*slice_weight);
				Kc_su.set(cseg.get(cs).getFirst(), 1, DU_cost*slice_weight);
				Kc_su.set(cseg.get(cs).getFirst(), 2, CU_cost*slice_weight);
				//dedicated transport is only on the non common part
				if(grooming) cD.add(0); else cD.add(1);

				//reliability?
				cR.add(1);

			}
			
			break;
		case 5:
			for (int cs =0; cs<cseg.size();cs++)
			{
				double slice_weight= (max_b-cseg.get(cs).getThird())/max_b;
				bc_suv.set(new int[] {cseg.get(cs).getFirst(),0,1}, (max_b-cseg.get(cs).getThird())/max_b *fronthaul);
				bc_suv.set(new int[] {cseg.get(cs).getFirst(),1,2}, (max_b-cseg.get(cs).getThird())/max_b *midhaul);
				bc_suv.set(new int[] {cseg.get(cs).getFirst(),2,3}, (max_b-cseg.get(cs).getThird())/max_b *backhaul);
				
				Mc_sun.set(new int[] {cseg.get(cs).getFirst(),0,cseg.get(cs).getSecond() }, 1);
				
				for (int j =0 ; j< DUs.size();j++)
				{
					Mc_sun.set(new int[] { cseg.get(cs).getFirst(), 1, DUs.get(j)}, 1 );
				}
				
				for (int j =0 ; j< CUs.size();j++)
				{
					Mc_sun.set(new int[] { cseg.get(cs).getFirst(), 2, CUs.get(j)}, 1 );
				}
				
				for (int j =0 ; j< NGCs.size();j++)
				{
					Mc_sun.set(new int[] { cseg.get(cs).getFirst(), 3, NGCs.get(j)}, 1 );
				}
				
				Kc_su.set(cseg.get(cs).getFirst(), 0, RU_cost*slice_weight);
				Kc_su.set(cseg.get(cs).getFirst(), 1, DU_cost*slice_weight);
				Kc_su.set(cseg.get(cs).getFirst(), 2, CU_cost*slice_weight);
				Kc_su.set(cseg.get(cs).getFirst(), 3, NGC_cost*slice_weight);
				//dedicated transport is only on the non common part
				if(grooming) cD.add(0); else cD.add(1);

				//reliability?
				cR.add(1);

			}
			
			break;
		
		}
		
		//MERGE 
		Omega_SUsu = new DoubleMatrixND(new int[] {Ns,Nf,Ns,Nf});
		// OMEGA Cslice,U,Dslice,U
		
		if (isolation>1 && isolation<5)//merge only for isolation 2,3,4
		{
			cD.addAll(D);
			D=cD;
			cR.addAll(R);
			R=cR;
			
			DoubleMatrixND M_sun_new = new DoubleMatrixND(new int[] { Ns, Nf, N });
			DoubleMatrixND b_suv_new = new DoubleMatrixND(new int[] { Ns, Nf, Nf });
			DoubleMatrix2D K_su_new=  DoubleFactory2D.sparse.make(Ns, Nf);
			
			for(int i=0; i<nslices;i++)
				for(int j=0; j<Nf; j++)
					for(int k=0;k<N;k++)
					{
						M_sun_new.set(new int[] { i+Nc, j, k }, M_sun.get(new int[] { i, j, k }) );
					}
			for(int i=0; i<Nc;i++)
				for(int j=0; j<Nf; j++)
					for(int k=0;k<N;k++)
					{
						M_sun_new.set(new int[] { i, j, k }, Mc_sun.get(new int[] { i, j, k }) );
					}
			for(int i=0; i<nslices;i++)
				for(int j=0; j<Nf; j++)
					for(int k=0;k<Nf;k++)
					{
						b_suv_new.set(new int[] { i+Nc, j, k }, b_suv.get(new int[] { i, j, k }) );
					}
			for(int i=0; i<Nc;i++)
				for(int j=0; j<Nf; j++)
					for(int k=0;k<Nf;k++)
					{
						b_suv_new.set(new int[] { i, j, k }, bc_suv.get(new int[] { i, j, k }) );
					}
			for(int i=0; i<nslices;i++)
				for(int j=0; j<Nf; j++)
					K_su_new.set(i+Nc, j, K_su.get(i, j));

			for(int i=0; i<Nc;i++)
				for(int j=0; j<Nf; j++)
					K_su_new.set(i, j, Kc_su.get(i, j));

			b_suv=b_suv_new;
			M_sun=M_sun_new;
			K_su=K_su_new;
			
			// WHAT ABOUT OMEGA? OMEGA EXISTS ONLY FOR ISOLATION 234
			
			for (int i=0; i< slice_cseg.size();i++)
			{
				switch (isolation) {
				case 2:
					Omega_SUsu.set(new int[] {slice_cseg.get(i).getSecond(),1,slice_cseg.get(i).getFirst()+Nc,1},1);
					break;
				case 3:
				case 4:
					Omega_SUsu.set(new int[] {slice_cseg.get(i).getSecond(),2,slice_cseg.get(i).getFirst()+Nc,2},1);
					break;
				
				}
			}
			
			//UPDATE CSEG WITH NEW IDs
			for (int i =0; i< slice_cseg.size();i++ )
			{
				slice_cseg.get(i).setFirst(slice_cseg.get(i).getFirst()+Nc);
			}
		
		}
		
		if(isolation==5)
		{
			b_suv=bc_suv;
			M_sun=Mc_sun;
			K_su=Kc_su;
			D=cD;
			R=cR;
		}
		
		}
		
		this.b_suv=b_suv;
		this.M_sun=M_sun;
		this.K_su=K_su;
		this.R=R;
		this.D=D;
		this.Ns=Ns;
		this.Omega_SUsu=Omega_SUsu; //should it be bidirectional?;
		
	}
	
	
	public void NodeMappingVNM(NetPlan netPlan)
	{
		Y_sun= new DoubleMatrixND(new int[] { Ns, Nf, N });
		
		if (isolation.getInt()==1 || isolation.getInt()==0 || isolation.getInt()==5)
			NodeMappingVN015 (netPlan);
			else
				NodeMappingVN234 (netPlan);
		
		
	}
	
	public void NodeMappingVNM_sur(NetPlan netPlan)
	{
		Y_sun= new DoubleMatrixND(new int[] { Ns, Nf, N });
		
		if (isolation.getInt()==1 || isolation.getInt()==0 || isolation.getInt()==5)
			NodeMappingVN015_sur (netPlan);
			else
				NodeMappingVN234_sur (netPlan);
		
		
	}
	
	public void NodeMappingVN234(NetPlan netPlan)
	{
		//================================================================================
		// NODE MAPPING VN
		// based on shortest path + already active nodes are considered in the cost 
		//================================================================================
		
		
		//create the structure with the mapping
		// map <cseg,List<dedicated>
		// for each entry
		// find h, the ru
		// find k the core
		// do the allocation form k to h
		// for each element remaining in the list do the allocation from k' to h'
		// 
		// 
		
		//create the map
		
		Map<Integer,ArrayList<Integer>> csegmapping= new HashedMap<Integer, ArrayList<Integer>>();
		//cseg, List of ded_slices
		for (int i =0; i< Nc; i++)
		{
			csegmapping.put(i, new ArrayList<Integer>() );
		}
		
		for (int i =0; i< slice_cseg.size() ; i++)
		{
			csegmapping.get(slice_cseg.get(i).getSecond()).add(slice_cseg.get(i).getFirst()) ;
		}
		
		ArrayList<Double> Capn_original= new ArrayList<>(Capn); //CHECK IF IT IS CORRECT HERE OR LATER
		
		for (Map.Entry<Integer,ArrayList<Integer>> entry: csegmapping.entrySet())
		{ 
			int sc= entry.getKey();
			int sd = entry.getValue().get(0);
			int antenna = slice_Requests.get(sd-Nc).getFirst();
			
			
			//do the mapping
			ArrayList<Integer> toServe= new ArrayList<Integer>();  
			toServe.add(3);toServe.add(2);toServe.add(1);toServe.add(0);
			
			int k= antenna;
			
			ArrayList<Integer> h =  findCandidateLocation (sd,3);//find cores
			SortedSet<Pair<Node,Node>> nodepairs= new TreeSet<Pair<Node,Node>>();

			
			for (int psp=0;psp<h.size();psp++)
			{
				if (k!=h.get(psp))
					nodepairs.add(Pair.of(netPlan.getNode(h.get(psp)),netPlan.getNode(k)));
			}

			SortedMap<Pair<Node,Node>,List<List<Link>>> cpl = netPlan.computeUnicastCandidatePathList(null , 10, -1, -1, -1, -1, -1, -1 , nodepairs );
			//sort ALL by cost
			//Arrailist<costo,<Coppia,Path>>
			ArrayList<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>> sortedcpl=new ArrayList<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>>();
			//compute cost of kshortestpaths
			double cost=0;

			for (Map.Entry<Pair<Node,Node>,List<List<Link>>> entry_cand: cpl.entrySet())
			{
				for (int i=0;i<entry_cand.getValue().size();i++)//for each path
				{
					//cost=entry.getValue().get(i).size();// lets try to do something better
					cost=computeCostVNMapping(entry_cand.getValue().get(i), Capn_original, Capn);//// -> TODO HERE IS THE WRONG CAPN_ORIGINAL!!!!!!
					sortedcpl.add(Pair.of(cost,Pair.of(entry_cand.getKey(), entry_cand.getValue().get(i))));
				}
			}
			//now we have candidate paths and cost
			//sort

			sortedcpl.sort(new Comparator<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>>() {
				@Override
				public int compare(Pair<Double,Pair<Pair<Node,Node>,List<Link>>> o1, Pair<Double,Pair<Pair<Node,Node>,List<Link>>> o2) {
					if (o1.getFirst() > o2.getFirst()) {
						return 1;
					} else if (o1.getFirst().equals(o2.getFirst())) {
						return 0; // You can change this to make it then look at the
						//words alphabetical order
					} else {
						return -1;
					}
				}
			});
			//now we have our set of candidate paths sorted
			
			//now we have our set of candidate paths sorted
			ArrayList<Integer> toServeOriginal=new ArrayList<Integer>(toServe);
			boolean foundpath =false;
			for (int sp=0;sp<sortedcpl.size() && !foundpath;sp++) //for each candidate path
			{
				toServe=new ArrayList<Integer>(toServeOriginal);
				ArrayList<Pair<Integer,Integer>> mapping = new ArrayList<Pair<Integer,Integer>>();
				//analyzes each shortest path
				ArrayList<Integer> path= new ArrayList<Integer>();
				for (Link l: sortedcpl.get(sp).getSecond().getSecond() )// NB HOW ARE YOU SURE THAT YOUR MAPPING ENDS IN K?
				{
					path.add(l.getOriginNode().getIndex());
					path.add(l.getDestinationNode().getIndex());
				}

				Object[] aux=Arrays.stream(path.toArray()).distinct().toArray();
				Integer[] in = Arrays.copyOf(aux, aux.length, Integer[].class);
				path.clear();
				path.addAll(Arrays.asList(in));
				//remove u from path
				//path.remove(0);//NOT FOR SFC //NOW WE DO THE MAPPING IN ONE SHOT AT THE END
				//path.remove(path.size()-1);
				
				//reverse path and t
		        //Collections.reverse(path); 
		        //Collections.reverse(toServe); 

				ArrayList<Integer> s_t = new ArrayList<Integer>();
				s_t.add(sc);s_t.add(sc);s_t.add(sd);s_t.add(sd);//check if it is valid always

				
				
				for (int t=0; t< toServe.size();t++)//maybe we can emprove with sequential mapping
				{
					
					
					boolean found=false;
					for (int i=0;i<path.size() && !found;i++)
					{
						int slice= s_t.get(toServe.get(t));
						if (isCapable(slice,toServe.get(0) , path.get(i))) //HERE IS THE MISTAKE S_T(t?, toserve(t)?, toserve(0)? ) 
						{

							mapping.add(Pair.of(toServe.get(t), path.get(i))); //CHECK: t or 0?
							//just temporary store toServefunction and physical node
							found=true;
							path.remove(i);//DIFFERENT FOR SFC
							i--;
							toServe.remove(t);
							t--;
						}
					}
				}
				//if still to serve, try new kpath
				if (toServe.size()==0)
				{
					//deplo
					//distribute among slices
							for (int m=0;m<mapping.size();m++)
							{
								int u= mapping.get(m).getFirst();
								int n= mapping.get(m).getSecond();
								
								switch (u)
								{
								case 0:
									Y_sun.set(new int[] {sc,u ,n},1);
									Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
									break;
								case 1:
									if (isolation.getInt()==2)
									{
										Y_sun.set(new int[] {sc,u ,n},1);
										//Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
										Y_sun.set(new int[] {sd,u ,n},1);
										Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
										
									}
									else
									{
										Y_sun.set(new int[] {sc,u ,n},1);
										Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 	
									}
									
									break;
								case 2:
									if (isolation.getInt()!=2)
									{
										Y_sun.set(new int[] {sc,u ,n},1);
										//Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
										Y_sun.set(new int[] {sd,u ,n},1);
										Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
										
									}
									else
									{
										Y_sun.set(new int[] {sd,u ,n},1);
										Capn.set(n,Capn.get(n)-K_su.get(sd,u)); 	
									}
									
									break;
								case 3:
									Y_sun.set(new int[] {sd,u ,n},1);
									Capn.set(n,Capn.get(n)-K_su.get(sd,u)); 
									break;
								
								}
								
							}
					foundpath=true;
				}
			
		}
			if (toServe.size()>0)
			{System.out.println("Mapping Error");throw new IllegalArgumentException("Unfeasible mapping");} 
		
			
			
			
			
			
			///////////////////////////////////////////////////////////////////////////////////////// DEDICATED SLICES
			
			//HERE WE have to map remaining dedicated slices
			
			for (int j=1; j< entry.getValue().size();j++)
			{
				int s = entry.getValue().get(j);
				 toServe= finduv(s);
				k=findCandidateLocation(s, toServe.get(0)).get(0);//source (epc) is many locations
				//Capn.set(h,Capn.get(h)-K_su.get(s, toServe.get(0)));
				//Y_sun.set(new int[] {s,toServe.get(0),h },1);
				//toServe.remove(0);
				//FOR SFC HERE I SHOULD ADD THE CASE H IS CAPABLE TO BE ENDPOINT

				nodepairs= new TreeSet<Pair<Node,Node>>();
				h= findCandidateLocation(s, toServe.get(toServe.size()-1));//look for possible destinations

				for (int psp=0;psp<h.size();psp++)
				{
					if (k!=h.get(psp))
						nodepairs.add(Pair.of(netPlan.getNode(k),netPlan.getNode(h.get(psp))));
				}

				cpl = netPlan.computeUnicastCandidatePathList(null , 10, -1, -1, -1, -1, -1, -1 , nodepairs );

				//sort ALL by cost
				//Arrailist<costo,<Coppia,Path>>
				sortedcpl=new ArrayList<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>>();
				//compute cost of kshortestpaths
				cost=0;

				for (Map.Entry<Pair<Node,Node>,List<List<Link>>> entry_cand: cpl.entrySet())
				{
					for (int i=0;i<entry.getValue().size();i++)//for each path
					{
						//cost=entry.getValue().get(i).size();// lets try to do something better
						cost=computeCostVNMapping(entry_cand.getValue().get(i), Capn_original, Capn);
						sortedcpl.add(Pair.of(cost,Pair.of(entry_cand.getKey(), entry_cand.getValue().get(i))));
					}
				}
				//now we have candidate paths and cost
				//sort

				sortedcpl.sort(new Comparator<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>>() {
					@Override
					public int compare(Pair<Double,Pair<Pair<Node,Node>,List<Link>>> o1, Pair<Double,Pair<Pair<Node,Node>,List<Link>>> o2) {
						if (o1.getFirst() > o2.getFirst()) {
							return 1;
						} else if (o1.getFirst().equals(o2.getFirst())) {
							return 0; // You can change this to make it then look at the
							//words alphabetical order
						} else {
							return -1;
						}
					}
				});
				//now we have our set of candidate paths sorted
				toServeOriginal=new ArrayList<Integer>(toServe);
				foundpath =false;
				for (int sp=0;sp<sortedcpl.size() && !foundpath;sp++) //for each candidate path
				{
					toServe=new ArrayList<Integer>(toServeOriginal);
					ArrayList<Pair<Integer,Integer>> mapping = new ArrayList<Pair<Integer,Integer>>();
					//analyzes each shortest path
					ArrayList<Integer> path= new ArrayList<Integer>();
					for (Link l: sortedcpl.get(sp).getSecond().getSecond() )// NB HOW ARE YOU SURE THAT YOUR MAPPING ENDS IN K?
					{
						path.add(l.getOriginNode().getIndex());
						path.add(l.getDestinationNode().getIndex());
					}

					Object[] aux=Arrays.stream(path.toArray()).distinct().toArray();
					Integer[] in = Arrays.copyOf(aux, aux.length, Integer[].class);
					path.clear();
					path.addAll(Arrays.asList(in));
					//remove u from path
					//path.remove(0);//NOT FOR SFC
					//path.remove(path.size()-1);
					
					//reverse path and t
			        //Collections.reverse(path); 
			        //Collections.reverse(toServe); 

					

					for (int t=0; t< toServe.size();t++)//maybe we can emprove we sequential mapping
					{
						boolean found=false;
						for (int i=0;i<path.size() && !found;i++)
						{
							if (isCapable(s,toServe.get(0) , path.get(i)))
							{

								mapping.add(Pair.of(toServe.get(t), path.get(i)));
								//just temporary store toServefunction and physical node
								found=true;
								path.remove(i);//DIFFERENT FOR SFC
								i--;
								toServe.remove(t);
								t--;
							}
						}
					}
					//if still to serve, try new kpath
					if (toServe.size()==0)
					{
						//deploy
						for (int m=0;m<mapping.size();m++)
						{
							int n= mapping.get(m).getSecond();
							int u= mapping.get(m).getFirst();
							Y_sun.set(new int[] {s,u ,n},1);
							Capn.set(n,Capn.get(n)-K_su.get(s,u)); 
						}
						foundpath=true;
					}
				}
				
				if (toServe.size()>0)
				{System.out.println("Mapping Error");throw new IllegalArgumentException("Unfeasible mapping");} 		
			}///END OF DEDICATED SLICES
		}/// END OF THE MAP common<-> dedicated
		
	}

	public void NodeMappingVN234_sur(NetPlan netPlan)
	{
		//================================================================================
		// NODE MAPPING VN
		// based on shortest path + already active nodes are considered in the cost 
		//================================================================================
		
		
		//create the structure with the mapping
		// map <cseg,List<dedicated>
		// for each entry
		// find h, the ru
		// find k the core
		// do the allocation form k to h
		// for each element remaining in the list do the allocation from k' to h'
		// 
		// 
		
		//create the map
		
		Map<Integer,ArrayList<Integer>> csegmapping= new HashedMap<Integer, ArrayList<Integer>>();
		//cseg, List of ded_slices
		for (int i =0; i< Nc; i++)
		{
			csegmapping.put(i, new ArrayList<Integer>() );
		}
		
		for (int i =0; i< slice_cseg.size() ; i++)
		{
			csegmapping.get(slice_cseg.get(i).getSecond()).add(slice_cseg.get(i).getFirst()) ;
		}
		
		ArrayList<Double> Capn_original= new ArrayList<>(Capn); //CHECK IF IT IS CORRECT HERE OR LATER
		SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k_without_placement=  makeConnectivityGraph(netPlan);
		SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k_with_placement= updateConnectivityGraphWithPlacement(netPlan, g_k_without_placement);
		
		for (Map.Entry<Integer,ArrayList<Integer>> entry: csegmapping.entrySet())
		{ 
			int sc= entry.getKey();
			int sd = entry.getValue().get(0);
			int antenna = slice_Requests.get(sd-Nc).getFirst();
			
			//System.out.println("Antenna "+antenna);
			
			
			
			//do the mapping
			ArrayList<Integer> toServe= new ArrayList<Integer>();  
			toServe.add(3);toServe.add(2);toServe.add(1);toServe.add(0);
			
			int k= antenna;
			
			//ArrayList<Integer> h =  findCandidateLocation (sd,3);//find cores
			int h = findSurNearestCore(k, netPlan);
			
//			SortedSet<Pair<Node,Node>> nodepairs= new TreeSet<Pair<Node,Node>>();
//
//			
//			for (int psp=0;psp<h.size();psp++)
//			{
//				if (k!=h.get(psp))
//					nodepairs.add(Pair.of(netPlan.getNode(h.get(psp)),netPlan.getNode(k)));
//			}
//
//			SortedMap<Pair<Node,Node>,List<List<Link>>> cpl = netPlan.computeUnicastCandidatePathList(null , 10, -1, -1, -1, -1, -1, -1 , nodepairs );
//			
//			
			//*
//			//sort ALL by cost
//			//Arrailist<costo,<Coppia,Path>>
//			//ArrayList<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>> sortedcpl=new ArrayList<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>>();
//			int K=K_sur_placement;
//
//			
//			KShortestSimplePaths<Integer, DefaultWeightedEdge> ksp = new KShortestSimplePaths<Integer, DefaultWeightedEdge>(g_k_sur, toServe.size()-1);
//			List<GraphPath<Integer, DefaultWeightedEdge>> sps= ksp.getPaths(netPlan.getNode(h.get(0)).getIndex(), netPlan.getNode(k).getIndex(), K);
//			System.out.println("Found"+ sps.size()+" paths over k="+K);
//			
//			ArrayList<Pair<Double,GraphPath<Integer,DefaultWeightedEdge>>> sortedcpl=new ArrayList<Pair<Double,GraphPath<Integer,DefaultWeightedEdge>>>();
//			
//			
//			//compute cost of kshortestpaths
//			double cost=0;
//
//			for (int i=0; i < sps.size();i++)
//			{
//					//cost=entry.getValue().get(i).size();// lets try to do something better
//					cost=computeCostVNMapping_sur(sps.get(i), g_k_sur, Capn_original, Capn);
//					sortedcpl.add(Pair.of(cost,sps.get(i)));
//			}
//			//now we have candidate paths and cost
//			//sort
//	
//			sortedcpl.sort(new Comparator<Pair<Double,GraphPath<Integer,DefaultWeightedEdge>>>() {
//				@Override
//				public int compare(Pair<Double,GraphPath<Integer,DefaultWeightedEdge>> o1, Pair<Double,GraphPath<Integer,DefaultWeightedEdge>> o2) {
//					if (o1.getFirst() > o2.getFirst()) {
//						return 1;
//					} else if (o1.getFirst().equals(o2.getFirst())) {
//						return 0; // You can change this to make it then look at the
//						//words alphabetical order
//					} else {
//						return -1;
//					}
//				}
//			});
//			//now we have our set of candidate paths sorted
//						
//			//now we have our set of candidate paths sorted
//			
//			//*
			
			
			g_k_with_placement=updateConnectivityGraphWithPlacement(netPlan, g_k_without_placement);
			
			//System.out.println(g_k_with_placement.getEdgeWeight(g_k_with_placement.getEdge(4, 1)));
			
			//EppsteinShortestPathIterator<Integer,DefaultWeightedEdge> it = new EppsteinShortestPathIterator<Integer,DefaultWeightedEdge>(g_k_with_placement, netPlan.getNode(h.get(0)).getIndex(), netPlan.getNode(k).getIndex());
			EppsteinShortestPathIterator<Integer,DefaultWeightedEdge> it = new EppsteinShortestPathIterator<Integer,DefaultWeightedEdge>(g_k_with_placement, netPlan.getNode(h).getIndex(), netPlan.getNode(k).getIndex());
			//YenShortestPathIterator<Integer,DefaultWeightedEdge> it = new YenShortestPathIterator<Integer, DefaultWeightedEdge>(g_k_with_placement, netPlan.getNode(h).getIndex(), netPlan.getNode(k).getIndex());
			
			
			ArrayList<Integer> toServeOriginal=new ArrayList<Integer>(toServe);
			boolean foundpath =false;
			
			int sp=0;
			//for (sp=0;sp>=0  && !foundpath;sp++) //for each candidate path
			while(it.hasNext() && !foundpath)
			{
				toServe=new ArrayList<Integer>(toServeOriginal);
				ArrayList<Pair<Integer,Integer>> mapping = new ArrayList<Pair<Integer,Integer>>();
				//analyzes each shortest path
				List<Integer> path= new ArrayList<Integer>();
//				for (Link l: sortedcpl.get(sp).getSecond().getSecond() )// NB HOW ARE YOU SURE THAT YOUR MAPPING ENDS IN K?
//				{
//					path.add(l.getOriginNode().getIndex());
//					path.add(l.getDestinationNode().getIndex());
//				}
	
				GraphPath<Integer, DefaultWeightedEdge> a = it.next();
				path= a.getVertexList();
				
//				System.out.println(a);
//				System.out.println(path);
//				System.out.println(a.getWeight());
//				System.out.println(a.getLength());
				
				
				Object[] aux=Arrays.stream(path.toArray()).distinct().toArray();
				Integer[] in = Arrays.copyOf(aux, aux.length, Integer[].class);
				path.clear();
				path.addAll(Arrays.asList(in));
				
				List<Integer> path_aux= new ArrayList<Integer>(path);
	
				//remove u from path
				//path.remove(0);//NOT FOR SFC //NOW WE DO THE MAPPING IN ONE SHOT AT THE END
				//path.remove(path.size()-1);
				
				//reverse path and t
		        //Collections.reverse(path); 
		        //Collections.reverse(toServe); 

				ArrayList<Integer> s_t = new ArrayList<Integer>();
				s_t.add(sc);s_t.add(sc);s_t.add(sd);s_t.add(sd);//check if it is valid always

				
				
				for (int t=0; t< toServe.size();t++)//maybe we can emprove with sequential mapping
				{
					
					
					boolean found=false;
					for (int i=0;i<path.size() && !found;i++)
					{
						int slice= s_t.get(toServe.get(t));
						if (isCapable(slice,toServe.get(0) , path.get(i))) //HERE IS THE MISTAKE S_T(t?, toserve(t)?, toserve(0)? ) 
						{

							mapping.add(Pair.of(toServe.get(t), path.get(i))); //CHECK: t or 0?
							//just temporary store toServefunction and physical node
							found=true;
							path.remove(i);//DIFFERENT FOR SFC
							i--;
							toServe.remove(t);
							t--;
						}
					}
				}
				//if still to serve, try new kpath
				
				
				
				if (toServe.size()==0)
				{
//					int c = countMappingNodes(mapping);
//					System.out.println("Counting mapping nodes: "+c);
					//deplo
					//distribute among slices
							for (int m=0;m<mapping.size();m++)
							{
								int u= mapping.get(m).getFirst();
								int n= mapping.get(m).getSecond();
								activeNodes.add(n);
								
								switch (u)
								{
								case 0:
									Y_sun.set(new int[] {sc,u ,n},1);
									//Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
									//only capnwithrus
									Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(sc,u)); 
									
									break;
								case 1:
									if (isolation.getInt()==2)
									{
										Y_sun.set(new int[] {sc,u ,n},1);
										//Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
										Y_sun.set(new int[] {sd,u ,n},1);
										Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
										Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(sc,u)); 
										load_withoutRUs+=K_su.get(sc,u);
										activeNodes_withoutRUs.add(n);
									}
									else
									{
										Y_sun.set(new int[] {sc,u ,n},1);
										Capn.set(n,Capn.get(n)-K_su.get(sc,u));
										Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(sc,u)); 
										load_withoutRUs+=K_su.get(sc,u);
										activeNodes_withoutRUs.add(n);
									}
									
									break;
								case 2:
									if (isolation.getInt()!=2)
									{
										Y_sun.set(new int[] {sc,u ,n},1);
										//Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
										Y_sun.set(new int[] {sd,u ,n},1);
										Capn.set(n,Capn.get(n)-K_su.get(sc,u));
										Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(sc,u)); 
										load_withoutRUs+=K_su.get(sc,u);
										activeNodes_withoutRUs.add(n);
										
									}
									else
									{
										Y_sun.set(new int[] {sd,u ,n},1);
										Capn.set(n,Capn.get(n)-K_su.get(sd,u)); 
										Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(sd,u)); 
										load_withoutRUs+=K_su.get(sd,u);
										activeNodes_withoutRUs.add(n);
									}
									
									break;
								case 3:
									Y_sun.set(new int[] {sd,u ,n},1);
									Capn.set(n,Capn.get(n)-K_su.get(sd,u));
									Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(sd,u)); 
									load_withoutRUs+=K_su.get(sd,u);
									activeNodes_withoutRUs.add(n);
									break;
								
								}
								
							}
					foundpath=true;
					
					//System.out.println("Selected path");
					//System.out.println(path_aux.toString());
				}
				
				sp++;
				if (sp>10000)//EMERGENCY EXIT since eppstein has cycles
				break;
		}
			
		
			
			
			if (toServe.size()>0)
			{
				System.out.println("Mapping Error");throw new IllegalArgumentException("Unfeasible mapping");
			} 
		
			
			
			
			
			
			///////////////////////////////////////////////////////////////////////////////////////// DEDICATED SLICES
			
			//HERE WE have to map remaining dedicated slices
			
			for (int j=1; j< entry.getValue().size();j++)
			{
				int s = entry.getValue().get(j);
		
				 toServe= finduv(s);
				 Collections.reverse(toServe); 
				
				 
				 //Capn.set(h,Capn.get(h)-K_su.get(s, toServe.get(0)));
				//Y_sun.set(new int[] {s,toServe.get(0),h },1);
				//toServe.remove(0);
				//FOR SFC HERE I SHOULD ADD THE CASE H IS CAPABLE TO BE ENDPOINT

				//nodepairs= new TreeSet<Pair<Node,Node>>();
				k= findCandidateLocation(s, toServe.get(toServe.size()-1)).get(0);//destination of a dedicated slice is one location

				//k=findCandidateLocation(s, toServe.get(0)).get(0);//source (epc) is many locations //source of a dedicated slice is always a core
				h=findSurNearestCore(k, netPlan);//source (epc) is many locations //source of a dedicated slice is always a core
				
				
//				for (int psp=0;psp<h.size();psp++)
//				{
//					if (k!=h.get(psp))
//						nodepairs.add(Pair.of(netPlan.getNode(k),netPlan.getNode(h.get(psp))));
//				}
//
//				cpl = netPlan.computeUnicastCandidatePathList(null , 10, -1, -1, -1, -1, -1, -1 , nodepairs );
//
//				//sort ALL by cost
//				//Arrailist<costo,<Coppia,Path>>
//				sortedcpl=new ArrayList<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>>();
				//compute cost of kshortestpaths
				//*
//				ksp = new KShortestSimplePaths<Integer, DefaultWeightedEdge>(g_k_sur, toServe.size()-1);
//				sps= ksp.getPaths(netPlan.getNode(h.get(0)).getIndex(), netPlan.getNode(k).getIndex(), K);
//				System.out.println("Found"+ sps.size()+" paths over k="+K);
//				
//				cost=0;
//
//				for (int i=0; i < sps.size();i++)
//				{
//						//cost=entry.getValue().get(i).size();// lets try to do something better
//						cost=computeCostVNMapping_sur(sps.get(i), g_k_sur, Capn_original, Capn);
//						sortedcpl.add(Pair.of(cost,sps.get(i)));
//				}
//				//now we have candidate paths and cost
//				//sort
//		
//				sortedcpl.sort(new Comparator<Pair<Double,GraphPath<Integer,DefaultWeightedEdge>>>() {
//					@Override
//					public int compare(Pair<Double,GraphPath<Integer,DefaultWeightedEdge>> o1, Pair<Double,GraphPath<Integer,DefaultWeightedEdge>> o2) {
//						if (o1.getFirst() > o2.getFirst()) {
//							return 1;
//						} else if (o1.getFirst().equals(o2.getFirst())) {
//							return 0; // You can change this to make it then look at the
//							//words alphabetical order
//						} else {
//							return -1;
//						}
//					}
//				});
		//*
				//now we have our set of candidate paths sorted
				
				g_k_with_placement=updateConnectivityGraphWithPlacement(netPlan, g_k_without_placement);
				
				//it = new EppsteinShortestPathIterator<Integer,DefaultWeightedEdge>(g_k_with_placement,  netPlan.getNode(k).getIndex(),netPlan.getNode(h.get(0)).getIndex());
				it = new EppsteinShortestPathIterator<Integer,DefaultWeightedEdge>(g_k_with_placement,  netPlan.getNode(h).getIndex(),netPlan.getNode(k).getIndex());
				//it = new YenShortestPathIterator<Integer, DefaultWeightedEdge>(g_k_with_placement, netPlan.getNode(h).getIndex(), netPlan.getNode(k).getIndex());

				
				toServeOriginal=new ArrayList<Integer>(toServe);
				foundpath =false;
				
				sp=0;
				//for (sp=0;sp>=0  && !foundpath;sp++) //for each candidate path
				while(it.hasNext() && !foundpath)
				{
					toServe=new ArrayList<Integer>(toServeOriginal);
					ArrayList<Pair<Integer,Integer>> mapping = new ArrayList<Pair<Integer,Integer>>();
					//analyzes each shortest path
					List<Integer> path= new ArrayList<Integer>();
//					for (Link l: sortedcpl.get(sp).getSecond().getSecond() )// NB HOW ARE YOU SURE THAT YOUR MAPPING ENDS IN K?
//					{
//						path.add(l.getOriginNode().getIndex());
//						path.add(l.getDestinationNode().getIndex());
//					}
		
					path= it.next().getVertexList();
					Object[] aux=Arrays.stream(path.toArray()).distinct().toArray();
					Integer[] in = Arrays.copyOf(aux, aux.length, Integer[].class);
					path.clear();
					path.addAll(Arrays.asList(in));
					//remove u from path
					//path.remove(0);//NOT FOR SFC
					//path.remove(path.size()-1);
					
					//reverse path and t
			        //Collections.reverse(path); 
			        Collections.reverse(toServe); 
			        
			        List<Integer> path_aux= new ArrayList<Integer>(path);

					

					for (int t=0; t< toServe.size();t++)//maybe we can emprove we sequential mapping
					{
						boolean found=false;
						for (int i=0;i<path.size() && !found;i++)
						{
							if (isCapable(s,toServe.get(0) , path.get(i)))
							{

//								if (path.get(i)==5)
//								{
//									System.out.println("Path");
//									System.out.println(path_aux.toString());
//									
//								}
								mapping.add(Pair.of(toServe.get(t), path.get(i)));
								//just temporary store toServefunction and physical node
								found=true;
								path.remove(i);//DIFFERENT FOR SFC
								i--;
								toServe.remove(t);
								t--;
							}
						}
					}
					//if still to serve, try new kpath
					if (toServe.size()==0)
					{
						//deploy
						for (int m=0;m<mapping.size();m++)
						{
							int n= mapping.get(m).getSecond();
							int u= mapping.get(m).getFirst();
							activeNodes.add(n);
							Y_sun.set(new int[] {s,u ,n},1);
							Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(s,u)); 
							if (u!=0)
							{
								Capn.set(n,Capn.get(n)-K_su.get(s,u)); 
								load_withoutRUs+=K_su.get(s,u);
								activeNodes_withoutRUs.add(n);
							}
							
						}
						foundpath=true;
					}
					sp++;
					if (sp>10000)//EMERGENCY EXIT since eppstein has cycles
						break;
				}
				
				if (toServe.size()>0)
				{System.out.println("Mapping Error");throw new IllegalArgumentException("Unfeasible mapping");} 		
			}///END OF DEDICATED SLICES
		}/// END OF THE MAP common<-> dedicated
		
	}

	
	public void NodeMappingVN015(NetPlan netPlan)
	{
		//================================================================================
		// NODE MAPPING VN
		// based on shortest path + already active nodes are considered in the cost 
		//================================================================================
		IntArrayList ind= new IntArrayList();
		IntArrayList ind2= new IntArrayList();

		int[] pos= new int[]{Ns,Nf,Nf};

		DoubleArrayList val= new DoubleArrayList();
		
		//System.out.println(Cn.getInt());
		ArrayList<Double> Capn_original= new ArrayList<>(Capn);

		for (int s =0; s<Ns; s++)//
		{
			ArrayList<Integer> toServe= finduv(s);
			int h=findCandidateLocation(s, toServe.get(0)).get(0);//source is always one location
			Capn.set(h,Capn.get(h)-K_su.get(s, toServe.get(0)));
			Y_sun.set(new int[] {s,toServe.get(0),h },1);
			
			toServe.remove(0);
			//FOR SFC HERE I SHOULD ADD THE CASE H IS CAPABLE TO BE ENDPOINT

			SortedSet<Pair<Node,Node>> nodepairs= new TreeSet<Pair<Node,Node>>();
			ArrayList<Integer> k= findCandidateLocation(s, toServe.get(toServe.size()-1));//look for possible destinations

			for (int psp=0;psp<k.size();psp++)
			{
				if (h!=k.get(psp))
					nodepairs.add(Pair.of(netPlan.getNode(h),netPlan.getNode(k.get(psp))));
			}

			SortedMap<Pair<Node,Node>,List<List<Link>>> cpl = netPlan.computeUnicastCandidatePathList(null , 10, -1, -1, -1, -1, -1, -1 , nodepairs );

			//sort ALL by cost
			//Arrailist<costo,<Coppia,Path>>
			ArrayList<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>> sortedcpl=new ArrayList<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>>();
			//compute cost of kshortestpaths
			double cost=0;

			for (Map.Entry<Pair<Node,Node>,List<List<Link>>> entry: cpl.entrySet())
			{
				for (int i=0;i<entry.getValue().size();i++)//for each path
				{
					//cost=entry.getValue().get(i).size();// lets try to do something better
					cost=computeCostVNMapping(entry.getValue().get(i), Capn_original, Capn);
					sortedcpl.add(Pair.of(cost,Pair.of(entry.getKey(), entry.getValue().get(i))));
				}
			}
			//now we have candidate paths and cost
			//sort

			sortedcpl.sort(new Comparator<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>>() {
				@Override
				public int compare(Pair<Double,Pair<Pair<Node,Node>,List<Link>>> o1, Pair<Double,Pair<Pair<Node,Node>,List<Link>>> o2) {
					if (o1.getFirst() > o2.getFirst()) {
						return 1;
					} else if (o1.getFirst().equals(o2.getFirst())) {
						return 0; // You can change this to make it then look at the
						//words alphabetical order
					} else {
						return -1;
					}
				}
			});
			//now we have our set of candidate paths sorted
			ArrayList<Integer> toServeOriginal=new ArrayList<Integer>(toServe);
			boolean foundpath =false;
			for (int sp=0;sp<sortedcpl.size() && !foundpath;sp++) //for each candidate path
			{
				toServe=new ArrayList<Integer>(toServeOriginal);
				ArrayList<Pair<Integer,Integer>> mapping = new ArrayList<Pair<Integer,Integer>>();
				//analyzes each shortest path
				ArrayList<Integer> path= new ArrayList<Integer>();
				for (Link l: sortedcpl.get(sp).getSecond().getSecond() )// NB HOW ARE YOU SURE THAT YOUR MAPPING ENDS IN K?
				{
					path.add(l.getOriginNode().getIndex());
					path.add(l.getDestinationNode().getIndex());
				}

				Object[] aux=Arrays.stream(path.toArray()).distinct().toArray();
				Integer[] in = Arrays.copyOf(aux, aux.length, Integer[].class);
				path.clear();
				path.addAll(Arrays.asList(in));
				//remove u from path
				path.remove(0);//NOT FOR SFC
				//path.remove(path.size()-1);
				
				//reverse path and t
		        Collections.reverse(path); 
		        Collections.reverse(toServe); 

				

				for (int t=0; t< toServe.size();t++)//maybe we can emprove we sequential mapping
				{
					boolean found=false;
					for (int i=0;i<path.size() && !found;i++)
					{
						if (isCapable(s,toServe.get(0) , path.get(i)))
						{

							mapping.add(Pair.of(toServe.get(t), path.get(i)));
							//just temporary store toServefunction and physical node
							found=true;
							path.remove(i);//DIFFERENT FOR SFC
							i--;
							toServe.remove(t);
							t--;
						}
					}
				}
				//if still to serve, try new kpath
				if (toServe.size()==0)
				{
					//deploy
					for (int m=0;m<mapping.size();m++)
					{
						int n= mapping.get(m).getSecond();
						int u= mapping.get(m).getFirst();
						Y_sun.set(new int[] {s,u ,n},1);
						Capn.set(n,Capn.get(n)-K_su.get(s,u)); 
					}


					foundpath=true;
				}
			}
			if (toServe.size()>0)
			{System.out.println("Mapping Error");throw new IllegalArgumentException("Unfeasible mapping");} 
			//else
				//System.out.println("Mapping Done");
		}

		

		//y
		//System.out.println("");
		//System.out.println("y_sun: "+ String.valueOf(Y_sun.zSum()));
		Y_sun.getNonZeros(ind, val);
		pos= new int[]{Ns,Nf,N};
		for (int i=0; i< ind.size();i++)
		{	 
			//System.out.println("y_sun("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+")="+String.valueOf(val.get(i))); 
		}
	}
	
	
	public void NodeMappingVN015_sur(NetPlan netPlan)
	{
		//================================================================================
		// NODE MAPPING VN
		// based on shortest path + already active nodes are considered in the cost 
		//================================================================================
		IntArrayList ind= new IntArrayList();
		IntArrayList ind2= new IntArrayList();

		int[] pos= new int[]{Ns,Nf,Nf};

		DoubleArrayList val= new DoubleArrayList();
		
		//System.out.println(Cn.getInt());
		ArrayList<Double> Capn_original= new ArrayList<>(Capn);
		SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k_without_placement=  makeConnectivityGraph(netPlan);
		SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k_with_placement= new SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		for (int s =0; s<Ns; s++)//
		{
			//System.out.println("New slice "+s);
//			if(s==123)
//				System.out.println("Stop");
			ArrayList<Integer> toServe= finduv(s);
			Collections.reverse(toServe); 
			//ArrayList<Integer> k= findCandidateLocation(s, toServe.get(toServe.size()-1));//antenna
			int k= findCandidateLocation(s, toServe.get(toServe.size()-1)).get(0);//antenna is always unique
			
			//System.out.println("Antenna "+k);
			
			//int h=findCandidateLocation(s, toServe.get(0)).get(0);// possible cores
			int h=findSurNearestCore(k, netPlan);
			
//			Capn.set(h,Capn.get(h)-K_su.get(s, toServe.get(0)));
//			Y_sun.set(new int[] {s,toServe.get(0),h },1);
//			if (toServe.get(0)!=0)
//			{
//				load_withoutRUs+=K_su.get(s,toServe.get(0));
//				activeNodes_withoutRUs.add(h);
//			}
//			
//			
//			toServe.remove(0);
			//FOR SFC HERE I SHOULD ADD THE CASE H IS CAPABLE TO BE ENDPOINT

			//SortedSet<Pair<Node,Node>> nodepairs= new TreeSet<Pair<Node,Node>>();
			
			
			g_k_with_placement=updateConnectivityGraphWithPlacement(netPlan, g_k_without_placement);
			
			//EppsteinShortestPathIterator<Integer,DefaultWeightedEdge> it = new EppsteinShortestPathIterator<Integer,DefaultWeightedEdge>(g_k_with_placement, netPlan.getNode(h).getIndex(), netPlan.getNode(k.get(0)).getIndex());
			EppsteinShortestPathIterator<Integer,DefaultWeightedEdge> it = new EppsteinShortestPathIterator<Integer,DefaultWeightedEdge>(g_k_with_placement, netPlan.getNode(h).getIndex(), netPlan.getNode(k).getIndex());
			//YenShortestPathIterator<Integer,DefaultWeightedEdge> it = new YenShortestPathIterator<Integer, DefaultWeightedEdge>(g_k_with_placement, netPlan.getNode(h).getIndex(), netPlan.getNode(k).getIndex());
			
			ArrayList<Integer> toServeOriginal=new ArrayList<Integer>(toServe);
			
			boolean foundpath =false;
			int sp=0;
			//for (sp=0;sp>=0  && !foundpath;sp++) //for each candidate path
			while(it.hasNext() && !foundpath)
			{
				//System.out.println("New path "+sp);
				toServe=new ArrayList<Integer>(toServeOriginal);
				ArrayList<Pair<Integer,Integer>> mapping = new ArrayList<Pair<Integer,Integer>>();
				//analyzes each shortest path
				List<Integer> path= new ArrayList<Integer>();
//				for (Link l: sortedcpl.get(sp).getSecond().getSecond() )// NB HOW ARE YOU SURE THAT YOUR MAPPING ENDS IN K?
//				{
//					path.add(l.getOriginNode().getIndex());
//					path.add(l.getDestinationNode().getIndex());
//				}
	
				//path=sortedcpl.get(sp).getSecond().getVertexList();
				GraphPath<Integer, DefaultWeightedEdge> a = it.next();
				path= a.getVertexList();
				
//				System.out.println(a);
//				System.out.println(path);
//				System.out.println(a.getWeight());
//				System.out.println(a.getLength());
				
				//System.out.println(path.toString());
				Object[] aux=Arrays.stream(path.toArray()).distinct().toArray();
				Integer[] in = Arrays.copyOf(aux, aux.length, Integer[].class);
				path.clear();
				path.addAll(Arrays.asList(in));
				
				List<Integer> path_aux= new ArrayList<Integer>(path);
				//remove u from path
				//path.remove(0);//NOT FOR SFC
				//path.remove(path.size()-1);
				
				//reverse path and t
		        //Collections.reverse(path); 
		        //Collections.reverse(toServe); 

				for (int t=0; t< toServe.size();t++)//maybe we can emprove we sequential mapping
				{
					boolean found=false;
					for (int i=0;i<path.size() && !found;i++)
					{
						if (isCapable(s,toServe.get(0) , path.get(i)))
						{
							mapping.add(Pair.of(toServe.get(t), path.get(i)));
							//just temporary store toServefunction and physical node
							found=true;
							path.remove(i);//DIFFERENT FOR SFC
							i--;
							toServe.remove(t);
							t--;
						}
					}
				}
				//if still to serve, try new kpath
				if (toServe.size()==0)
				{
					
					//System.out.println(mapping.toString());
					//deploy
					for (int m=0;m<mapping.size();m++)
					{
						int n= mapping.get(m).getSecond();
						int u= mapping.get(m).getFirst();
						//activeNodes.add(n);
						Y_sun.set(new int[] {s,u ,n},1);
						Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(s,u)); 
						activeNodes.add(n);
						if (u!=0)
						{
							Capn.set(n,Capn.get(n)-K_su.get(s,u)); 
							load_withoutRUs+=K_su.get(s,u);
							activeNodes_withoutRUs.add(n);
						}
						
					}

					//here we have to update the graph
					
					foundpath=true;
					//System.out.println("Selected path");
					//System.out.println(path_aux.toString());
					
					
				}
				sp++;
				
			if (sp>10000)//EMERGENCY EXIT since eppstein has cycles
				break;
			}
			
			if (toServe.size()>0)
			{
				System.out.println("Mapping Error");
				throw new IllegalArgumentException("Unfeasible mapping");
			} 
			//else
				//System.out.println("Mapping Done");
		}

		

		//y
//		System.out.println("");
//		System.out.println("y_sun: "+ String.valueOf(Y_sun.zSum()));
//		Y_sun.getNonZeros(ind, val);
//		pos= new int[]{Ns,Nf,N};
//		for (int i=0; i< ind.size();i++)
//		{	 
//			System.out.println("y_sun("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+")="+String.valueOf(val.get(i))); 
//		}
		
	}
	
	
	public void NodeMappingSFCM(NetPlan netPlan)
	{
		Y_sun= new DoubleMatrixND(new int[] { Ns, Nf, N });
		
		if (isolation.getInt()==1 || isolation.getInt()==0 || isolation.getInt()==5)
			NodeMappingSFC015 (netPlan);
			else
				NodeMappingSFC234 (netPlan);
		
		
	}

	
	public void NodeMappingSFC234(NetPlan netPlan)
		{
			//================================================================================
			// NODE MAPPING VN
			// based on shortest path + already active nodes are considered in the cost 
			//================================================================================
			
			
			//create the structure with the mapping
			// map <cseg,List<dedicated>
			// for each entry
			// find h, the ru
			// find k the core
			// do the allocation form k to h
			// for each element remaining in the list do the allocation from k' to h'
			// 
			// 
			
			//create the map
			
			Map<Integer,ArrayList<Integer>> csegmapping= new HashedMap<Integer, ArrayList<Integer>>();
			//cseg, List of ded_slices
			for (int i =0; i< Nc; i++)
			{
				csegmapping.put(i, new ArrayList<Integer>() );
			}
			
			for (int i =0; i< slice_cseg.size() ; i++)
			{
				csegmapping.get(slice_cseg.get(i).getSecond()).add(slice_cseg.get(i).getFirst()) ;
			}
			ArrayList<Double> Capn_original= new ArrayList<>(Capn); //CHECK IF IT IS CORRECT HERE OR LATE
			for (Map.Entry<Integer,ArrayList<Integer>> entry: csegmapping.entrySet())
			{ 
				int sc= entry.getKey();
				int sd = entry.getValue().get(0);
				int antenna = slice_Requests.get(sd-Nc).getFirst();

				
				//do the mapping
				ArrayList<Integer> toServe= new ArrayList<Integer>();  
				toServe.add(3);toServe.add(2);toServe.add(1);toServe.add(0);
				
				int k= antenna;
				
				ArrayList<Integer> h =  findCandidateLocation (sd,3);//find cores
				SortedSet<Pair<Node,Node>> nodepairs= new TreeSet<Pair<Node,Node>>();
		
				
				for (int psp=0;psp<h.size();psp++)
				{
					if (k!=h.get(psp))
						nodepairs.add(Pair.of(netPlan.getNode(h.get(psp)),netPlan.getNode(k)));
				}
		
				SortedMap<Pair<Node,Node>,List<List<Link>>> cpl = netPlan.computeUnicastCandidatePathList(null , 10, -1, -1, -1, -1, -1, -1 , nodepairs );
				//sort ALL by cost
				//Arrailist<costo,<Coppia,Path>>
				ArrayList<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>> sortedcpl=new ArrayList<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>>();
				//compute cost of kshortestpaths
				double cost=0;
		
				for (Map.Entry<Pair<Node,Node>,List<List<Link>>> entry_cand: cpl.entrySet())
				{
					for (int i=0;i<entry_cand.getValue().size();i++)//for each path
					{
						//cost=entry.getValue().get(i).size();// lets try to do something better
						cost=computeCostVNMapping(entry_cand.getValue().get(i), Capn_original, Capn);
						sortedcpl.add(Pair.of(cost,Pair.of(entry_cand.getKey(), entry_cand.getValue().get(i))));
					}
				}
				//now we have candidate paths and cost
				//sort
		
				sortedcpl.sort(new Comparator<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>>() {
					@Override
					public int compare(Pair<Double,Pair<Pair<Node,Node>,List<Link>>> o1, Pair<Double,Pair<Pair<Node,Node>,List<Link>>> o2) {
						if (o1.getFirst() > o2.getFirst()) {
							return 1;
						} else if (o1.getFirst().equals(o2.getFirst())) {
							return 0; // You can change this to make it then look at the
							//words alphabetical order
						} else {
							return -1;
						}
					}
				});
				//now we have our set of candidate paths sorted
				
				//now we have our set of candidate paths sorted
				ArrayList<Integer> toServeOriginal=new ArrayList<Integer>(toServe);
				boolean foundpath =false;
				for (int sp=0;sp<sortedcpl.size() && !foundpath;sp++) //for each candidate path
				{
					toServe=new ArrayList<Integer>(toServeOriginal);
					ArrayList<Pair<Integer,Integer>> mapping = new ArrayList<Pair<Integer,Integer>>();
					//analyzes each shortest path
					ArrayList<Integer> path= new ArrayList<Integer>();
					for (Link l: sortedcpl.get(sp).getSecond().getSecond() )// NB HOW ARE YOU SURE THAT YOUR MAPPING ENDS IN K?
					{
						path.add(l.getOriginNode().getIndex());
						path.add(l.getDestinationNode().getIndex());
					}
		
					Object[] aux=Arrays.stream(path.toArray()).distinct().toArray();
					Integer[] in = Arrays.copyOf(aux, aux.length, Integer[].class);
					path.clear();
					path.addAll(Arrays.asList(in));
					//remove u from path
					//path.remove(0);//NOT FOR SFC //NOW WE DO THE MAPPING IN ONE SHOT AT THE END
					//path.remove(path.size()-1);
					
					//reverse path and t
			        //Collections.reverse(path); 
			        //Collections.reverse(toServe); 
		
					ArrayList<Integer> s_t = new ArrayList<Integer>();
					s_t.add(sc);s_t.add(sc);s_t.add(sd);s_t.add(sd);//check if it is valid always
		
					
					
	//				for (int t=0; t< toServe.size();t++)//maybe we can emprove with sequential mapping
	//				{
	//					
	//					
	//					boolean found=false;
	//					for (int i=0;i<path.size() && !found;i++)
	//					{
	//						int slice= s_t.get(toServe.get(t));
	//						if (isCapable(slice,toServe.get(0) , path.get(i))) //HERE IS THE MISTAKE S_T(t?, toserve(t)?, toserve(0)? ) 
	//						{
	//	
	//							mapping.add(Pair.of(toServe.get(t), path.get(i))); //CHECK: t or 0?
	//							//just temporary store toServefunction and physical node
	//							found=true;
	//							path.remove(i);//DIFFERENT FOR SFC
	//							i--;
	//							toServe.remove(t);
	//							t--;
	//						}
	//					}
	//				}
					
					
			    	int p=0;
					while (toServe.size()>0 && p<path.size())
					{
						int slice= s_t.get(toServe.get(0));
						if (isCapable(slice, toServe.get(0), path.get(p)))
						{
							mapping.add(Pair.of(toServe.get(0), path.get(p)));
							toServe.remove(0);
						}
						else
						{
							p++;
						}
	
					} 
					
					
					//if still to serve, try new kpath
					if (toServe.size()==0) /// CHECK HERE SOMETHING IS STRANGE
					{
						//deplo
						//distribute among slices
								for (int m=0;m<mapping.size();m++)
								{
									int u= mapping.get(m).getFirst();
									int n= mapping.get(m).getSecond();
									
									switch (u)
									{
									case 0:
										Y_sun.set(new int[] {sc,u ,n},1);
										Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
										break;
									case 1:
										if (isolation.getInt()==2)
										{
											Y_sun.set(new int[] {sc,u ,n},1);
											//Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
											Y_sun.set(new int[] {sd,u ,n},1);
											Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
											
										}
										else
										{
											Y_sun.set(new int[] {sc,u ,n},1);
											Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 	
										}
										
										break;
									case 2:
										if (isolation.getInt()==2)
										{
											Y_sun.set(new int[] {sd,u ,n},1);
											Capn.set(n,Capn.get(n)-K_su.get(sd,u)); 	
										}
										else
										{
											Y_sun.set(new int[] {sc,u ,n},1);
											//Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
											Y_sun.set(new int[] {sd,u ,n},1);
											Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
										}
										
										break;
									case 3:
										Y_sun.set(new int[] {sd,u ,n},1);
										Capn.set(n,Capn.get(n)-K_su.get(sd,u)); 
										break;
									
									}
									
								}
						foundpath=true;
					}
				
			}
				if (toServe.size()>0)
				{System.out.println("Mapping Error");throw new IllegalArgumentException("Unfeasible mapping");} 
			
				
				
				
				
				
				///////////////////////////////////////////////////////////////////////////////////////// DEDICATED SLICES
				
				//HERE WE have to map remaining dedicated slices
				
				for (int j=1; j< entry.getValue().size();j++)
				{
					int s = entry.getValue().get(j);
					 toServe= finduv(s);
					k=findCandidateLocation(s, toServe.get(0)).get(0);//source (epc) is many locations
					//Capn.set(h,Capn.get(h)-K_su.get(s, toServe.get(0)));
					//Y_sun.set(new int[] {s,toServe.get(0),h },1);
					//toServe.remove(0);
					//FOR SFC HERE I SHOULD ADD THE CASE H IS CAPABLE TO BE ENDPOINT
		
					nodepairs= new TreeSet<Pair<Node,Node>>();
					h= findCandidateLocation(s, toServe.get(toServe.size()-1));//look for possible destinations
		
					for (int psp=0;psp<h.size();psp++)
					{
						if (k!=h.get(psp))
							nodepairs.add(Pair.of(netPlan.getNode(k),netPlan.getNode(h.get(psp))));
					}
		
					cpl = netPlan.computeUnicastCandidatePathList(null , 10, -1, -1, -1, -1, -1, -1 , nodepairs );
		
					//sort ALL by cost
					//Arrailist<costo,<Coppia,Path>>
					sortedcpl=new ArrayList<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>>();
					//compute cost of kshortestpaths
					cost=0;
					
					//if (sortedcpl.size()==0 && k==h.get(0)) //means everything could be in the same node
						//sortedcpl.add(Pair.of(0,Pair.of(Pair.of(netPlan.getNode(k),netPlan.getNode(h.get(0))), entry_cand.getValue().get(i))));
		
					for (Map.Entry<Pair<Node,Node>,List<List<Link>>> entry_cand: cpl.entrySet())
					{
						for (int i=0;i<entry.getValue().size();i++)//for each path
						{
							//cost=entry.getValue().get(i).size();// lets try to do something better
							cost=computeCostVNMapping(entry_cand.getValue().get(i), Capn_original, Capn);
							sortedcpl.add(Pair.of(cost,Pair.of(entry_cand.getKey(), entry_cand.getValue().get(i))));
						}
					}
					//now we have candidate paths and cost
					//sort
		
					sortedcpl.sort(new Comparator<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>>() {
						@Override
						public int compare(Pair<Double,Pair<Pair<Node,Node>,List<Link>>> o1, Pair<Double,Pair<Pair<Node,Node>,List<Link>>> o2) {
							if (o1.getFirst() > o2.getFirst()) {
								return 1;
							} else if (o1.getFirst().equals(o2.getFirst())) {
								return 0; // You can change this to make it then look at the
								//words alphabetical order
							} else {
								return -1;
							}
						}
					});
					//now we have our set of candidate paths sorted
					toServeOriginal=new ArrayList<Integer>(toServe);
					foundpath =false;
					for (int sp=0;sp<sortedcpl.size() && !foundpath;sp++) //for each candidate path
					{
						toServe=new ArrayList<Integer>(toServeOriginal);
						ArrayList<Pair<Integer,Integer>> mapping = new ArrayList<Pair<Integer,Integer>>();
						//analyzes each shortest path
						ArrayList<Integer> path= new ArrayList<Integer>();
						for (Link l: sortedcpl.get(sp).getSecond().getSecond() )// NB HOW ARE YOU SURE THAT YOUR MAPPING ENDS IN K?
						{
							path.add(l.getOriginNode().getIndex());
							path.add(l.getDestinationNode().getIndex());
						}
		
						Object[] aux=Arrays.stream(path.toArray()).distinct().toArray();
						Integer[] in = Arrays.copyOf(aux, aux.length, Integer[].class);
						path.clear();
						path.addAll(Arrays.asList(in));
						//remove u from path
						//path.remove(0);//NOT FOR SFC
						//path.remove(path.size()-1);
						
						//reverse path and t
				        //Collections.reverse(path); 
				        //Collections.reverse(toServe); 
		
						
		
	//					for (int t=0; t< toServe.size();t++)//maybe we can emprove we sequential mapping
	//					{
	//						boolean found=false;
	//						for (int i=0;i<path.size() && !found;i++)
	//						{
	//							if (isCapable(s,toServe.get(0) , path.get(i)))
	//							{
	//	
	//								mapping.add(Pair.of(toServe.get(t), path.get(i)));
	//								//just temporary store toServefunction and physical node
	//								found=true;
	//								path.remove(i);//DIFFERENT FOR SFC
	//								i--;
	//								toServe.remove(t);
	//								t--;
	//							}
	//						}
	//					}
	//					
		
						
						
				    	int p=0;
						while (toServe.size()>0 && p<path.size())
						{
							if (isCapable(s, toServe.get(0), path.get(p)))
							{
								mapping.add(Pair.of(toServe.get(0), path.get(p)));
								toServe.remove(0);
							}
							else
							{
								p++;
							}
	
						} 
						
						
						
						//if still to serve, try new kpath
						if (toServe.size()==0)
						{
							//deploy
							for (int m=0;m<mapping.size();m++)
							{
								int n= mapping.get(m).getSecond();
								int u= mapping.get(m).getFirst();
								Y_sun.set(new int[] {s,u ,n},1);
								Capn.set(n,Capn.get(n)-K_su.get(s,u)); 
							}
							foundpath=true;
						}
					}
					
					if (toServe.size()>0 && sortedcpl.size()==0 && k==h.get(0)) //means everything could be in the same node
					{
						toServe=new ArrayList<Integer>(toServeOriginal);
						ArrayList<Pair<Integer,Integer>> mapping = new ArrayList<Pair<Integer,Integer>>();
						//analyzes each shortest path
						ArrayList<Integer> path= new ArrayList<Integer>();
						path.add(k);
						int p=0;
						while (toServe.size()>0 && p<path.size())
						{
							if (isCapable(s, toServe.get(0), path.get(p)))
							{
								mapping.add(Pair.of(toServe.get(0), path.get(p)));
								toServe.remove(0);
							}
							else
							{
								p++;
							}
	
						} 
							
						//if still to serve, try new kpath
						if (toServe.size()==0)
						{
							//deploy
							for (int m=0;m<mapping.size();m++)
							{
								int n= mapping.get(m).getSecond();
								int u= mapping.get(m).getFirst();
								Y_sun.set(new int[] {s,u ,n},1);
								Capn.set(n,Capn.get(n)-K_su.get(s,u)); 
							}
							foundpath=true;
						}
						
						
					}
						
					
					if (toServe.size()>0)
					{System.out.println("Mapping Error");throw new IllegalArgumentException("Unfeasible mapping");} 		
				}///END OF DEDICATED SLICES
			}/// END OF THE MAP common<-> dedicated
			
		}


	public void NodeMappingSFC015(NetPlan netPlan)
	{
		//================================================================================
		// NODE MAPPING SFC
		// based on shortest path + already active nodes are considered in the cost 
		//================================================================================
		IntArrayList ind= new IntArrayList();
		IntArrayList ind2= new IntArrayList();
	
		int[] pos= new int[]{Ns,Nf,Nf};
	
		DoubleArrayList val= new DoubleArrayList();
		
		//System.out.println(Cn.getInt());
		ArrayList<Double> Capn_original= new ArrayList<>(Capn);
	
		for (int s =0; s<Ns; s++)//
		{
			ArrayList<Integer> toServe= finduv(s);
			int h=findCandidateLocation(s, toServe.get(0)).get(0);//source is always one location
			Capn.set(h,Capn.get(h)-K_su.get(s, toServe.get(0)));
			Y_sun.set(new int[] {s,toServe.get(0),h },1);
			toServe.remove(0);
			//FOR SFC HERE I SHOULD ADD THE CASE H IS CAPABLE TO BE ENDPOINT
	
			SortedSet<Pair<Node,Node>> nodepairs= new TreeSet<Pair<Node,Node>>();
			ArrayList<Integer> k= findCandidateLocation(s, toServe.get(toServe.size()-1));//look for possible destinations
	
			for (int psp=0;psp<k.size();psp++)
			{
				if (h!=k.get(psp))
					nodepairs.add(Pair.of(netPlan.getNode(h),netPlan.getNode(k.get(psp))));
			}
	
			SortedMap<Pair<Node,Node>,List<List<Link>>> cpl = netPlan.computeUnicastCandidatePathList(null , 10, -1, -1, -1, -1, -1, -1 , nodepairs );
	
			//sort ALL by cost
			//Arrailist<costo,<Coppia,Path>>
			ArrayList<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>> sortedcpl=new ArrayList<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>>();
			//compute cost of kshortestpaths
			double cost=0;
	
			for (Map.Entry<Pair<Node,Node>,List<List<Link>>> entry: cpl.entrySet())
			{
				for (int i=0;i<entry.getValue().size();i++)//for each path
				{
					//cost=entry.getValue().get(i).size();// lets try to do something better
					cost=computeCostVNMapping(entry.getValue().get(i), Capn_original, Capn);
					sortedcpl.add(Pair.of(cost,Pair.of(entry.getKey(), entry.getValue().get(i))));
				}
			}
			//now we have candidate paths and cost
			//sort
	
			sortedcpl.sort(new Comparator<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>>() {
				@Override
				public int compare(Pair<Double,Pair<Pair<Node,Node>,List<Link>>> o1, Pair<Double,Pair<Pair<Node,Node>,List<Link>>> o2) {
					if (o1.getFirst() > o2.getFirst()) {
						return 1;
					} else if (o1.getFirst().equals(o2.getFirst())) {
						return 0; // You can change this to make it then look at the
						//words alphabetical order
					} else {
						return -1;
					}
				}
			});
			//now we have our set of candidate paths sorted
			ArrayList<Integer> toServeOriginal=new ArrayList<Integer>(toServe);
			boolean foundpath =false;
			for (int sp=0;sp<sortedcpl.size() && !foundpath;sp++) //for each candidate path
			{
				toServe=new ArrayList<Integer>(toServeOriginal);
				ArrayList<Pair<Integer,Integer>> mapping = new ArrayList<Pair<Integer,Integer>>();
				//analyzes each shortest path
				ArrayList<Integer> path= new ArrayList<Integer>();
				for (Link l: sortedcpl.get(sp).getSecond().getSecond() )// NB HOW ARE YOU SURE THAT YOUR MAPPING ENDS IN K?
				{
					path.add(l.getOriginNode().getIndex());
					path.add(l.getDestinationNode().getIndex());
				}
	
				Object[] aux=Arrays.stream(path.toArray()).distinct().toArray();
				Integer[] in = Arrays.copyOf(aux, aux.length, Integer[].class);
				path.clear();
				path.addAll(Arrays.asList(in));
				//remove u from path
				path.remove(0);//NOT FOR SFC
				//path.remove(path.size()-1);
				
				//reverse path and t
		        Collections.reverse(path); 
		        Collections.reverse(toServe); 
	
				
	
		    	int p=0;
				while (toServe.size()>0 && p<path.size())
				{
					if (isCapable(s, toServe.get(0), path.get(p)))
					{
						mapping.add(Pair.of(toServe.get(0), path.get(p)));
						toServe.remove(0);
					}
					else
					{
						p++;
					}
	
				} 
				
				
				//if still to serve, try new kpath
				if (toServe.size()==0)
				{
					//deploy
					for (int m=0;m<mapping.size();m++)
					{
						int n= mapping.get(m).getSecond();
						int u= mapping.get(m).getFirst();
						Y_sun.set(new int[] {s,u ,n},1);
						Capn.set(n,Capn.get(n)-K_su.get(s,u)); 
					}
	
	
					foundpath=true;
				}
			}
			
			if (toServe.size()>0 && sortedcpl.size()==0 && h==k.get(0)) //means everything could be in the same node
			{
				toServe=new ArrayList<Integer>(toServeOriginal);
				ArrayList<Pair<Integer,Integer>> mapping = new ArrayList<Pair<Integer,Integer>>();
				//analyzes each shortest path
				ArrayList<Integer> path= new ArrayList<Integer>();
				path.add(h);
				int p=0;
				while (toServe.size()>0 && p<path.size())
				{
					if (isCapable(s, toServe.get(0), path.get(p)))
					{
						mapping.add(Pair.of(toServe.get(0), path.get(p)));
						toServe.remove(0);
					}
					else
					{
						p++;
					}

				} 
					
				//if still to serve, try new kpath
				if (toServe.size()==0)
				{
					//deploy
					for (int m=0;m<mapping.size();m++)
					{
						int n= mapping.get(m).getSecond();
						int u= mapping.get(m).getFirst();
						Y_sun.set(new int[] {s,u ,n},1);
						Capn.set(n,Capn.get(n)-K_su.get(s,u)); 
					}
					foundpath=true;
				}
				
				
			}
			
			if (toServe.size()>0)
			{System.out.println("Mapping Error");throw new IllegalArgumentException("Unfeasible mapping");} 
			//else
				//System.out.println("Mapping Done");
		}
	
		
	
		//y
		//System.out.println("");
		//System.out.println("y_sun: "+ String.valueOf(Y_sun.zSum()));
		Y_sun.getNonZeros(ind, val);
		pos= new int[]{Ns,Nf,N};
		for (int i=0; i< ind.size();i++)
		{	 
			//System.out.println("y_sun("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+")="+String.valueOf(val.get(i))); 
		}
	}

	public void NodeMappingSFCM_sur(NetPlan netPlan)
	{
		Y_sun= new DoubleMatrixND(new int[] { Ns, Nf, N });
		number_of_nodes_sfc= Character.getNumericValue(virt.getString().charAt(3));
		if (isolation.getInt()==1 || isolation.getInt()==0 || isolation.getInt()==5)
			NodeMappingSFC015_sur (netPlan);
			else
				NodeMappingSFC234_sur (netPlan);
		
		
	}
	
	public void NodeMappingSFC015_sur(NetPlan netPlan)
	{
		//================================================================================
		// NODE MAPPING SFC
		// based on shortest path on Surballes connectivity graph + already active nodes are considered in the cost 
		//================================================================================
		IntArrayList ind= new IntArrayList();
		IntArrayList ind2= new IntArrayList();
	
		int[] pos= new int[]{Ns,Nf,Nf};
	
		DoubleArrayList val= new DoubleArrayList();
		
		//System.out.println(Cn.getInt());
		ArrayList<Double> Capn_original= new ArrayList<>(Capn);
	
		SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k_without_placement=  makeConnectivityGraph(netPlan);
		SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k_with_placement= new SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		
		for (int s =0; s<Ns; s++)//
		{
			//System.out.println("Slice: "+s);
			//if (s==46)
			//	System.out.println("stop");
			ArrayList<Integer> toServe= finduv(s);
			Collections.reverse(toServe); 
			//ArrayList<Integer> k= findCandidateLocation(s, toServe.get(toServe.size()-1));//antenna is always one
			int k= findCandidateLocation(s, toServe.get(toServe.size()-1)).get(0);//antenna is always unique
			
			
			//int h=findCandidateLocation(s, toServe.get(0)).get(0);//source are possible cores
			int h= findSurNearestCore(k, netPlan);//source are possible cores

			//			Capn.set(h,Capn.get(h)-K_su.get(s, toServe.get(0)));
//			Y_sun.set(new int[] {s,toServe.get(0),h },1);
//			if (toServe.get(0)!=0)
//			{
//				load_withoutRUs+=K_su.get(s,toServe.get(0));
//				activeNodes_withoutRUs.add(h);
//			}
			
			//toServe.remove(0);
			//FOR SFC HERE I SHOULD ADD THE CASE H IS CAPABLE TO BE ENDPOINT
	 
			//*
			//SortedSet<Pair<Node,Node>> nodepairs= new TreeSet<Pair<Node,Node>>();
//			ArrayList<Integer> k= findCandidateLocation(s, toServe.get(toServe.size()-1));//look for possible destinations
//	
//			int K=K_sur_placement;
//
//			KShortestSimplePaths<Integer, DefaultWeightedEdge> ksp = new KShortestSimplePaths<Integer, DefaultWeightedEdge>(g_k_sur, toServe.size()-1+1);
//			List<GraphPath<Integer, DefaultWeightedEdge>> sps= ksp.getPaths(netPlan.getNode(h).getIndex(), netPlan.getNode(k.get(0)).getIndex(), K);
//			System.out.println("Found"+ sps.size()+" paths over k="+K);
//			
//			//List<DefaultWeightedEdge> edges;//=sps.get(0).getEdgeList();
//
//
//			
//			//SuurballeKDisjointShortestPaths<Integer,DefaultWeightedEdge> sur= new SuurballeKDisjointShortestPaths<Integer,DefaultWeightedEdge>(g_k);
//			//we take only the nearest core
//			//List<GraphPath<Integer,DefaultWeightedEdge>> ksp=sur.getPaths(netPlan.getNode(h).getIndex() ,netPlan.getNode(k.get(0)).getIndex() , 999); 
//
//			
//			//SortedMap<Pair<Node,Node>,List<List<Link>>> cpl = netPlan.computeUnicastCandidatePathList(null , 10, -1, -1, -1, -1, -1, -1 , nodepairs );
//	
//			//sort ALL by cost
//			//Arrailist<costo,<Coppia,Path>>
//			ArrayList<Pair<Double,GraphPath<Integer,DefaultWeightedEdge>>> sortedcpl=new ArrayList<Pair<Double,GraphPath<Integer,DefaultWeightedEdge>>>();
//			//compute cost of kshortestpaths
//			double cost=0;
//			
//			
//			for (int i=0; i < sps.size();i++)
//			{
//					//cost=entry.getValue().get(i).size();// lets try to do something better
//					cost=computeCostVNMapping_sur(sps.get(i), g_k_sur, Capn_original, Capn);
//					sortedcpl.add(Pair.of(cost,sps.get(i)));
//			}
//			//now we have candidate paths and cost
//			//sort
//	
//			sortedcpl.sort(new Comparator<Pair<Double,GraphPath<Integer,DefaultWeightedEdge>>>() {
//				@Override
//				public int compare(Pair<Double,GraphPath<Integer,DefaultWeightedEdge>> o1, Pair<Double,GraphPath<Integer,DefaultWeightedEdge>> o2) {
//					if (o1.getFirst() > o2.getFirst()) {
//						return 1;
//					} else if (o1.getFirst().equals(o2.getFirst())) {
//						return 0; // You can change this to make it then look at the
//						//words alphabetical order
//					} else {
//						return -1;
//					}
//				}
//			});
//			//now we have our set of candidate paths sorted
//			ArrayList<Integer> toServeOriginal=new ArrayList<Integer>(toServe);
			
			//*
			
			
			g_k_with_placement=updateConnectivityGraphWithPlacement(netPlan, g_k_without_placement);
			
			//EppsteinShortestPathIterator<Integer,DefaultWeightedEdge> it = new EppsteinShortestPathIterator<Integer,DefaultWeightedEdge>(g_k_with_placement, netPlan.getNode(h).getIndex(),netPlan.getNode(k.get(0)).getIndex() );
			EppsteinShortestPathIterator<Integer,DefaultWeightedEdge> it = new EppsteinShortestPathIterator<Integer,DefaultWeightedEdge>(g_k_with_placement, netPlan.getNode(h).getIndex(),netPlan.getNode(k).getIndex() );
			//YenShortestPathIterator<Integer,DefaultWeightedEdge> it = new YenShortestPathIterator<Integer, DefaultWeightedEdge>(g_k_with_placement, netPlan.getNode(h).getIndex(), netPlan.getNode(k).getIndex());

			
			ArrayList<Integer> toServeOriginal=new ArrayList<Integer>(toServe);
			//it.next();
			
			boolean foundpath =false;
			int sp=0;
			//for (sp=0;sp>=0  && !foundpath;sp++) //for each candidate path
			while(it.hasNext() && !foundpath)
			{
				toServe=new ArrayList<Integer>(toServeOriginal);
				ArrayList<Pair<Integer,Integer>> mapping = new ArrayList<Pair<Integer,Integer>>();
				//analyzes each shortest path
				List<Integer> path= new ArrayList<Integer>();
//				for (Link l: sortedcpl.get(sp).getSecond().getSecond() )// NB HOW ARE YOU SURE THAT YOUR MAPPING ENDS IN K?
//				{
//					path.add(l.getOriginNode().getIndex());
//					path.add(l.getDestinationNode().getIndex());
//				}
	
				GraphPath<Integer, DefaultWeightedEdge> gp= it.next();
				
				path= gp.getVertexList();
//				if (s==123) {
//				System.out.println(path.toString());
//					System.out.println("stop");
//				//System.out.println(gp.getWeight());
//				}
				if (path.size()==0)
					break;
				
				Object[] aux=Arrays.stream(path.toArray()).distinct().toArray();
				Integer[] in = Arrays.copyOf(aux, aux.length, Integer[].class);
				if (path.size()>1)
				{
				path.clear();
				path.addAll(Arrays.asList(in));
				}
				//remove u from path
				//path.remove(0);//NOT FOR SFC
				//path.remove(path.size()-1);
				
				//reverse path and t
		        //Collections.reverse(path); 
		        //Collections.reverse(toServe); 
	
				
	
		    	int p=0;
				while (toServe.size()>0 && p<path.size())
				{
					if (isCapable(s, toServe.get(0), path.get(p)))
					{
						mapping.add(Pair.of(toServe.get(0), path.get(p)));
						toServe.remove(0);
						
						if (number_of_nodes_sfc==4)
						{
							p++;
						}
						
					}
					else
					{
						p++;
					}
	
				} 
				
				int c = countMappingNodes(mapping);
				
				//if still to serve, try new kpath
				if ((toServe.size()==0 && c==number_of_nodes_sfc) || (toServe.size()==0 && number_of_nodes_sfc==1) )
				{
					//System.out.println(mapping.toString());
					//deploy
					for (int m=0;m<mapping.size();m++)
					{
						int n= mapping.get(m).getSecond();
						int u= mapping.get(m).getFirst();
						
						activeNodes.add(n);
						Y_sun.set(new int[] {s,u,n},1);
						Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(s,u)); 
						
						if (u!=0)
						{
							Capn.set(n,Capn.get(n)-K_su.get(s,u)); 
							load_withoutRUs+=K_su.get(s,u);
							activeNodes_withoutRUs.add(n);
						}
					}
	
	
					foundpath=true;
				}
				
				sp++;
				if (sp>10000)//EMERGENCY EXIT since eppstein has cycles
					break;
			}
			
			//if (toServe.size()>0 && sp==0 && h==k.get(0)) //means everything could be in the same node (sp==1??)
			if (toServe.size()>0 && sp==0 && h==k) //means everything could be in the same node (sp==1??)
			{
				toServe=new ArrayList<Integer>(toServeOriginal);
				ArrayList<Pair<Integer,Integer>> mapping = new ArrayList<Pair<Integer,Integer>>();
				//analyzes each shortest path
				ArrayList<Integer> path= new ArrayList<Integer>();
				path.add(h);
				int p=0;
				while (toServe.size()>0 && p<path.size())
				{
					if (isCapable(s, toServe.get(0), path.get(p)))
					{
						mapping.add(Pair.of(toServe.get(0), path.get(p)));
						toServe.remove(0);
					}
					else
					{
						p++;
					}

				} 
					
				//if still to serve, try new kpath
				if (toServe.size()==0)
				{
					//deploy
					for (int m=0;m<mapping.size();m++)
					{
						int n= mapping.get(m).getSecond();
						int u= mapping.get(m).getFirst();
						activeNodes.add(n);
						Y_sun.set(new int[] {s,u ,n},1);
						Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(s,u)); 
						if (u!=0)
						{
							Capn.set(n,Capn.get(n)-K_su.get(s,u)); 
							load_withoutRUs+=K_su.get(s,u);
							activeNodes_withoutRUs.add(n);
						}
					}
					foundpath=true;
				}
				
				
			}
			
			if (toServe.size()>0)
			{
				System.out.println("Mapping Error");
				throw new IllegalArgumentException("Unfeasible mapping");
			} 
			//else
				//System.out.println("Mapping Done");
		}
	
		
	
		//y
		//System.out.println("");
		//System.out.println("y_sun: "+ String.valueOf(Y_sun.zSum()));
//		Y_sun.getNonZeros(ind, val);
//		pos= new int[]{Ns,Nf,N};
//		for (int i=0; i< ind.size();i++)
//		{	 
//			//System.out.println("y_sun("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+")="+String.valueOf(val.get(i))); 
//		}
	}

	private int countMappingNodes(ArrayList<Pair<Integer, Integer>> mapping) {

		Set<Integer> nodes=  new HashSet<Integer>();
		
		for (int i=0; i< mapping.size();i++)
		{
			nodes.add(mapping.get(i).getSecond());
		}
		
		return nodes.size();
	}


	public void NodeMappingSFC234_sur(NetPlan netPlan)
	{
		//================================================================================
		// NODE MAPPING VN
		// based on shortest path + already active nodes are considered in the cost 
		//================================================================================
		
		
		//create the structure with the mapping
		// map <cseg,List<dedicated>
		// for each entry
		// find h, the ru
		// find k the core
		// do the allocation form k to h
		// for each element remaining in the list do the allocation from k' to h'
		// 
		// 
		
		//create the map
		
		Map<Integer,ArrayList<Integer>> csegmapping= new HashedMap<Integer, ArrayList<Integer>>();
		//cseg, List of ded_slices
		for (int i =0; i< Nc; i++)
		{
			csegmapping.put(i, new ArrayList<Integer>() );
		}
		
		for (int i =0; i< slice_cseg.size() ; i++)
		{
			csegmapping.get(slice_cseg.get(i).getSecond()).add(slice_cseg.get(i).getFirst()) ;
		}
		ArrayList<Double> Capn_original= new ArrayList<>(Capn); //CHECK IF IT IS CORRECT HERE OR LATE
		SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k_without_placement=  makeConnectivityGraph(netPlan);
		SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k_with_placement= updateConnectivityGraphWithPlacement(netPlan, g_k_without_placement);
		
		
		//try save and load
		FileOutputStream fileOut;
		try {
			fileOut = new FileOutputStream("Graph_Epp");
			ObjectOutputStream objectOut = new ObjectOutputStream(fileOut);	
			objectOut.writeObject(g_k_with_placement);
			objectOut.close();
			fileOut.close();
			
			FileInputStream fi = new FileInputStream(new File("Graph_Epp"));
			ObjectInputStream oi = new ObjectInputStream(fi);
			SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k_with_placement_loaded= (SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>) oi.readObject();
			g_k_with_placement=g_k_with_placement_loaded;
			oi.close();
			fi.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
		
		for (Map.Entry<Integer,ArrayList<Integer>> entry: csegmapping.entrySet())
		{ 
			int sc= entry.getKey();
			int sd = entry.getValue().get(0);
			int antenna = slice_Requests.get(sd-Nc).getFirst();

			
			//do the mapping
			ArrayList<Integer> toServe= new ArrayList<Integer>();  
			toServe.add(3);toServe.add(2);toServe.add(1);toServe.add(0);
			
			int k= antenna;
			
			//ArrayList<Integer> h =  findCandidateLocation (sd,3);//find cores
			int h = findSurNearestCore(k, netPlan);
			//SortedSet<Pair<Node,Node>> nodepairs= new TreeSet<Pair<Node,Node>>();
	
			
//			for (int psp=0;psp<h.size();psp++)
//			{
//				if (k!=h.get(psp))
//					nodepairs.add(Pair.of(netPlan.getNode(h.get(psp)),netPlan.getNode(k)));
//			}
	
			//SortedMap<Pair<Node,Node>,List<List<Link>>> cpl = netPlan.computeUnicastCandidatePathList(null , 10, -1, -1, -1, -1, -1, -1 , nodepairs );
			//sort ALL by cost
			//Arrailist<costo,<Coppia,Path>>
			//ArrayList<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>> sortedcpl=new ArrayList<Pair<Double,Pair<Pair<Node,Node>,List<Link>>>>();
			//compute cost of kshortestpaths
			
//			int K=K_sur_placement;
//
//			KShortestSimplePaths<Integer, DefaultWeightedEdge> ksp = new KShortestSimplePaths<Integer, DefaultWeightedEdge>(g_k_sur, toServe.size()-1);
//			List<GraphPath<Integer, DefaultWeightedEdge>> sps= ksp.getPaths(netPlan.getNode(h.get(0)).getIndex(), netPlan.getNode(k).getIndex(), K);
//
//			System.out.println("Found"+ sps.size()+" paths over k="+K);
//			ArrayList<Pair<Double,GraphPath<Integer,DefaultWeightedEdge>>> sortedcpl=new ArrayList<Pair<Double,GraphPath<Integer,DefaultWeightedEdge>>>();
//			//compute cost of kshortestpaths
//			double cost=0;
//	
//			for (int i=0; i < sps.size();i++)
//			{
//					//cost=entry.getValue().get(i).size();// lets try to do something better
//					cost=computeCostVNMapping_sur(sps.get(i), g_k_sur, Capn_original, Capn);
//					sortedcpl.add(Pair.of(cost,sps.get(i)));
//			}
//			//now we have candidate paths and cost
//			//sort
//	
//			sortedcpl.sort(new Comparator<Pair<Double,GraphPath<Integer,DefaultWeightedEdge>>>() {
//				@Override
//				public int compare(Pair<Double,GraphPath<Integer,DefaultWeightedEdge>> o1, Pair<Double,GraphPath<Integer,DefaultWeightedEdge>> o2) {
//					if (o1.getFirst() > o2.getFirst()) {
//						return 1;
//					} else if (o1.getFirst().equals(o2.getFirst())) {
//						return 0; // You can change this to make it then look at the
//						//words alphabetical order
//					} else {
//						return -1;
//					}
//				}
//			});
			//now we have our set of candidate paths sorted
			
			
			g_k_with_placement=updateConnectivityGraphWithPlacement(netPlan, g_k_without_placement);
			
			//System.out.println(g_k_with_placement.getEdgeWeight(g_k_with_placement.getEdge(4, 1)));
			
			//EppsteinShortestPathIterator<Integer,DefaultWeightedEdge> it = new EppsteinShortestPathIterator<Integer,DefaultWeightedEdge>(g_k_with_placement, netPlan.getNode(h.get(0)).getIndex(), netPlan.getNode(k).getIndex());
			EppsteinShortestPathIterator<Integer,DefaultWeightedEdge> it = new EppsteinShortestPathIterator<Integer,DefaultWeightedEdge>(g_k_with_placement, netPlan.getNode(h).getIndex(), netPlan.getNode(k).getIndex());
			//YenShortestPathIterator<Integer,DefaultWeightedEdge> it = new YenShortestPathIterator<Integer, DefaultWeightedEdge>(g_k_with_placement, netPlan.getNode(h).getIndex(), netPlan.getNode(k).getIndex());

			
//			for (int rovap=0;rovap<=5  ;rovap++)
//			{
//				GraphPath<Integer, DefaultWeightedEdge> a = it.next();
//				List<Integer> path = a.getVertexList();
//				
//				System.out.println(a);
//				System.out.println(path);
//				System.out.println(a.getWeight());
//				System.out.println(a.getLength());
//			}
//			
			//it.next();
			
			//now we have our set of candidate paths sorted
			ArrayList<Integer> toServeOriginal=new ArrayList<Integer>(toServe);
			boolean foundpath =false;
			int sp=0;
			//for (sp=0;sp>=0  && !foundpath;sp++) //for each candidate path
			while(it.hasNext() && !foundpath)
			{
				toServe=new ArrayList<Integer>(toServeOriginal);
				ArrayList<Pair<Integer,Integer>> mapping = new ArrayList<Pair<Integer,Integer>>();
				//analyzes each shortest path
				List<Integer> path= new ArrayList<Integer>();
//				for (Link l: sortedcpl.get(sp).getSecond().getSecond() )// NB HOW ARE YOU SURE THAT YOUR MAPPING ENDS IN K?
//				{
//					path.add(l.getOriginNode().getIndex());
//					path.add(l.getDestinationNode().getIndex());
//				}
	
				
				GraphPath<Integer, DefaultWeightedEdge> a = it.next();
				path= a.getVertexList();
				
//				System.out.println(a);
//				System.out.println(path);
//				System.out.println(a.getWeight());
//				System.out.println(a.getLength());
				
				Object[] aux=Arrays.stream(path.toArray()).distinct().toArray();
				Integer[] in = Arrays.copyOf(aux, aux.length, Integer[].class);
				path.clear();
				path.addAll(Arrays.asList(in));
				//remove u from path
				//path.remove(0);//NOT FOR SFC //NOW WE DO THE MAPPING IN ONE SHOT AT THE END
				//path.remove(path.size()-1);
				
				//reverse path and t
		        //Collections.reverse(path); 
		        //Collections.reverse(toServe); 
	
				ArrayList<Integer> s_t = new ArrayList<Integer>();
				s_t.add(sc);s_t.add(sc);s_t.add(sd);s_t.add(sd);//check if it is valid always
				
				
				
				//HERE I HAVE TO BE CAREFUL OF SFC1(Unconstrained), SFC2, SFC3, SFC4
		    	int p=0;
				while (toServe.size()>0 && p<path.size())
				{
					int slice= s_t.get(toServe.get(0));
					if (isCapable(slice, toServe.get(0), path.get(p)))
					{
						mapping.add(Pair.of(toServe.get(0), path.get(p)));
						toServe.remove(0);
						//here p++ or not basing on FCx and function
						
						if (number_of_nodes_sfc==4)
						{
							p++;
						}
					}
					else
					{
						p++;
					}

				} 
				

				int c = countMappingNodes(mapping);
				//System.out.println("Counting mapping nodes: "+c);
				//if still to serve, try new kpath
				if ((toServe.size()==0 && c==number_of_nodes_sfc) || (toServe.size()==0 && number_of_nodes_sfc==1) ) /// CHECK HERE SOMETHING IS STRANGE
				{
					
					//deplo
					//distribute among slices
							for (int m=0;m<mapping.size();m++)
							{
								int u= mapping.get(m).getFirst();
								int n= mapping.get(m).getSecond();
								activeNodes.add(n);
								switch (u)
								{
								case 0:
									Y_sun.set(new int[] {sc,u ,n},1);
									Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(sc,u)); 
									break;
								case 1:
									if (isolation.getInt()==2)
									{
										Y_sun.set(new int[] {sc,u ,n},1);
										//Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
										Y_sun.set(new int[] {sd,u ,n},1);
										Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
										Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(sc,u)); 
										load_withoutRUs+=K_su.get(sc,u);
										activeNodes_withoutRUs.add(n);
										
									}
									else
									{
										Y_sun.set(new int[] {sc,u ,n},1);
										Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
										Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(sc,u)); 
										load_withoutRUs+=K_su.get(sc,u);
										activeNodes_withoutRUs.add(n);
									}
									
									break;
								case 2:
									if (isolation.getInt()==2)
									{
										Y_sun.set(new int[] {sd,u ,n},1);
										Capn.set(n,Capn.get(n)-K_su.get(sd,u)); 
										Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(sd,u)); 
										load_withoutRUs+=K_su.get(sd,u);
										activeNodes_withoutRUs.add(n);
									}
									else
									{
										Y_sun.set(new int[] {sc,u ,n},1);
										//Capn.set(n,Capn.get(n)-K_su.get(sc,u)); 
										Y_sun.set(new int[] {sd,u ,n},1);
										Capn.set(n,Capn.get(n)-K_su.get(sc,u));
										Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(sc,u)); 
										load_withoutRUs+=K_su.get(sc,u);
										activeNodes_withoutRUs.add(n);
									}
									
									break;
								case 3:
									Y_sun.set(new int[] {sd,u ,n},1);
									Capn.set(n,Capn.get(n)-K_su.get(sd,u)); 
									Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(sd,u)); 
									load_withoutRUs+=K_su.get(sd,u);
									activeNodes_withoutRUs.add(n);
									break;
								
								}
								
							}
					foundpath=true;
				}
				
				sp++;
				if (sp>10000)//EMERGENCY EXIT since eppstein has cycles
					break;
			
		}
			if (toServe.size()>0)
			{System.out.println("Mapping Error");throw new IllegalArgumentException("Unfeasible mapping");} 
		
			
			
			
			
			
			///////////////////////////////////////////////////////////////////////////////////////// DEDICATED SLICES
			
			//HERE WE have to map remaining dedicated slices
			
			for (int j=1; j< entry.getValue().size();j++)
			{
				int s = entry.getValue().get(j);
				toServe= finduv(s);
				Collections.reverse(toServe);
				
				//h= findCandidateLocation(s, toServe.get(toServe.size()-1));//destination of a dedicated slice is one location
				k= findCandidateLocation(s, toServe.get(toServe.size()-1)).get(0);//destination of a dedicated slice is unique location
				
				
				//h=findCandidateLocation(s, toServe.get(0)).get(0);//source (epc) is many locations// source of a dedicated slice is a core
				//h=findCandidateLocation(s, toServe.get(0)).get(0);//source (epc) is many locations// source of a dedicated slice is a core
				h = findSurNearestCore(k, netPlan);
				//Capn.set(h,Capn.get(h)-K_su.get(s, toServe.get(0)));
				//Y_sun.set(new int[] {s,toServe.get(0),h },1);
				//toServe.remove(0);
				//FOR SFC HERE I SHOULD ADD THE CASE H IS CAPABLE TO BE ENDPOINT
	
				//nodepairs= new TreeSet<Pair<Node,Node>>();
				
//				for (int psp=0;psp<h.size();psp++)
//				{
//					if (k!=h.get(psp))
//						nodepairs.add(Pair.of(netPlan.getNode(k),netPlan.getNode(h.get(psp))));
//				}
//	
//				cpl = netPlan.computeUnicastCandidatePathList(null , 10, -1, -1, -1, -1, -1, -1 , nodepairs );
//	
//				int K=5;

//				ksp = new KShortestSimplePaths<Integer, DefaultWeightedEdge>(g_k_sur, toServe.size()-1);
//				sps= ksp.getPaths(netPlan.getNode(h.get(0)).getIndex(), netPlan.getNode(k).getIndex(), K);
//
//				System.out.println("Found"+ sps.size()+" paths over k="+K);
//				
//				cost=0;
//				
//				for (int i=0; i < sps.size();i++)
//				{
//						//cost=entry.getValue().get(i).size();// lets try to do something better
//						cost=computeCostVNMapping_sur(sps.get(i), g_k_sur, Capn_original, Capn);
//						sortedcpl.add(Pair.of(cost,sps.get(i)));
//				}
//				//now we have candidate paths and cost
//				//sort
//		
//				sortedcpl.sort(new Comparator<Pair<Double,GraphPath<Integer,DefaultWeightedEdge>>>() {
//					@Override
//					public int compare(Pair<Double,GraphPath<Integer,DefaultWeightedEdge>> o1, Pair<Double,GraphPath<Integer,DefaultWeightedEdge>> o2) {
//						if (o1.getFirst() > o2.getFirst()) {
//							return 1;
//						} else if (o1.getFirst().equals(o2.getFirst())) {
//							return 0; // You can change this to make it then look at the
//							//words alphabetical order
//						} else {
//							return -1;
//						}
//					}
//				});
					//now we have our set of candidate paths sorted
				

				g_k_with_placement=updateConnectivityGraphWithPlacement(netPlan, g_k_without_placement);
				
				//it = new EppsteinShortestPathIterator<Integer,DefaultWeightedEdge>(g_k_with_placement,  netPlan.getNode(k).getIndex(),netPlan.getNode(h.get(0)).getIndex());
				it = new EppsteinShortestPathIterator<Integer,DefaultWeightedEdge>(g_k_with_placement,  netPlan.getNode(h).getIndex(),netPlan.getNode(k).getIndex());
				//it = new YenShortestPathIterator<Integer, DefaultWeightedEdge>(g_k_with_placement, netPlan.getNode(h).getIndex(), netPlan.getNode(k).getIndex());

				
				toServeOriginal=new ArrayList<Integer>(toServe);
				foundpath =false;
				sp=0;
				//for (sp=0;sp>=0  && !foundpath;sp++) //for each candidate path
				while(it.hasNext() && !foundpath)
				{
					toServe=new ArrayList<Integer>(toServeOriginal);
					ArrayList<Pair<Integer,Integer>> mapping = new ArrayList<Pair<Integer,Integer>>();
					//analyzes each shortest path
					List<Integer> path= new ArrayList<Integer>();
//					for (Link l: sortedcpl.get(sp).getSecond().getSecond() )// NB HOW ARE YOU SURE THAT YOUR MAPPING ENDS IN K?
//					{
//						path.add(l.getOriginNode().getIndex());
//						path.add(l.getDestinationNode().getIndex());
//					}
		
					path= it.next().getVertexList();
					Object[] aux=Arrays.stream(path.toArray()).distinct().toArray();
					Integer[] in = Arrays.copyOf(aux, aux.length, Integer[].class);
					
					if (path.size()>1)
					{
					path.clear();
					path.addAll(Arrays.asList(in));
					}
					
						//remove u from path
					//path.remove(0);//NOT FOR SFC
					//path.remove(path.size()-1);
					
					//reverse path and t
			        //Collections.reverse(path); 
			        //Collections.reverse(toServe); 
	
					
	
//					for (int t=0; t< toServe.size();t++)//maybe we can emprove we sequential mapping
//					{
//						boolean found=false;
//						for (int i=0;i<path.size() && !found;i++)
//						{
//							if (isCapable(s,toServe.get(0) , path.get(i)))
//							{
//	
//								mapping.add(Pair.of(toServe.get(t), path.get(i)));
//								//just temporary store toServefunction and physical node
//								found=true;
//								path.remove(i);//DIFFERENT FOR SFC
//								i--;
//								toServe.remove(t);
//								t--;
//							}
//						}
//					}
//					
	
					
					
			    	int p=0;
					while (toServe.size()>0 && p<path.size())
					{
						if (isCapable(s, toServe.get(0), path.get(p)))
						{
							mapping.add(Pair.of(toServe.get(0), path.get(p)));
							toServe.remove(0);
							
							if (number_of_nodes_sfc==4)
							{
								p++;
							}
							
						}
						else
						{
							p++;
						}

					} 
					
					
					
					//if still to serve, try new kpath
					if (toServe.size()==0)
					{
						//deploy
						for (int m=0;m<mapping.size();m++)
						{
							int n= mapping.get(m).getSecond();
							int u= mapping.get(m).getFirst();
							activeNodes.add(n);
							Y_sun.set(new int[] {s,u ,n},1);
							Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(s,u)); 
							
							if (u!=0)
							{
								Capn.set(n,Capn.get(n)-K_su.get(s,u)); 
								load_withoutRUs+=K_su.get(s,u);
								activeNodes_withoutRUs.add(n);
							}
							
						}
						foundpath=true;
					}
					
					sp++;
					if (sp>10000)//EMERGENCY EXIT since eppstein has cycles
						break;
				}
				
				//if (toServe.size()>0 && sp==0 && k==h.get(0)) //means everything could be in the same node
				if (toServe.size()>0 && sp==0 && k==h) //means everything could be in the same node
					
				{
					toServe=new ArrayList<Integer>(toServeOriginal);
					ArrayList<Pair<Integer,Integer>> mapping = new ArrayList<Pair<Integer,Integer>>();
					//analyzes each shortest path
					ArrayList<Integer> path= new ArrayList<Integer>();
					path.add(k);
					int p=0;
					while (toServe.size()>0 && p<path.size())
					{
						if (isCapable(s, toServe.get(0), path.get(p)))
						{
							mapping.add(Pair.of(toServe.get(0), path.get(p)));
							toServe.remove(0);
						}
						else
						{
							p++;
						}

					} 
						
					//if still to serve, try new kpath
					if (toServe.size()==0)
					{
						//deploy
						for (int m=0;m<mapping.size();m++)
						{
							int n= mapping.get(m).getSecond();
							int u= mapping.get(m).getFirst();
							activeNodes.add(n);
							Y_sun.set(new int[] {s,u ,n},1);
							Capn_wRUs.set(n,Capn_wRUs.get(n)-K_su.get(s,u)); 
							
							if (u!=0)
							{
								Capn.set(n,Capn.get(n)-K_su.get(s,u)); 
								load_withoutRUs+=K_su.get(s,u);
								activeNodes_withoutRUs.add(n);
							}
							
						}
						foundpath=true;
					}
					
					
				}
					
				
				if (toServe.size()>0)
				{System.out.println("Mapping Error");throw new IllegalArgumentException("Unfeasible mapping");} 		
			}///END OF DEDICATED SLICES
		}/// END OF THE MAP common<-> dedicated
		
	}

	
	
	public void LinkMapping(NetPlan netPlan)
	{
		//prerequisite: have b_suv ready
		
		//IN V6 USES K-SURBALLES ALGORITHM
		IntArrayList ind= new IntArrayList();
		IntArrayList ind2= new IntArrayList();

		int[] pos= new int[]{Ns,Nf,N};
		DoubleArrayList val= new DoubleArrayList();

		w_wl=0;
		p_wl=0;
		gentraff=0;

		//Convert Y_sun in Y(s,u)
		DoubleMatrix2D Y_su=DoubleFactory2D.sparse.make(Ns,Nf);
		Y_sun.getNonZeros(ind, val);

		for (int i=0; i< ind.size();i++)
		{	 
			Y_su.set(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0), DoubleMatrixND.ind2sub(ind.get(i), pos).get(1), DoubleMatrixND.ind2sub(ind.get(i), pos).get(2));
		}

		//netPlan.getVectorLinkLengthInKm(netPlan.getNetworkLayerDefault());
		
		
		//COMMENTED FOR V5
		//System.out.println("CALCULATING DISJOINT PATHS");
		//SortedMap<Pair<Node,Node>,List<List<Link>>> cpl = netPlan.computeUnicastCandidatePathList( null, 200, -1, -1, -1, -1, -1, -1 ,null,netPlan.getNetworkLayerDefault());
		//SortedMap<Pair<Node,Node>,List<Pair<List<Link>,List<Link>>>> cpl11 = NetPlan.computeUnicastCandidate11PathList(cpl, 2);
		
		
		//GraphUtils.getTwoLinkDisjointPaths(nodes, links, originNode, destinationNode, linkCostMap)
		
		
		//System.out.println("DISJOINT PATHS CALCULATED");
		Pair<Node,Node> nodePair;

		//create JGraph for plightpath and grooming
		//SimpleWeightedGraph<Integer, DefaultWeightedEdge> g = new SimpleWeightedGraph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g= new SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);

		//index difference between electrical and plightpathnodes
		// electrical 0,1,2,3,...
		// plighpaths 1000,1001,1002,1003,...
		int offset=1000;
		
		//System.out.println("CREATING TRANSPODER EDGES");
		
		

		for (Node n:netPlan.getNodes())		 
		{
			//pl.addVertex(n.getIndex());
			g.addVertex(n.getIndex());
			g.addVertex(n.getIndex()+offset);
			g.addEdge(n.getIndex(), n.getIndex()+offset);
			g.addEdge(n.getIndex()+offset,n.getIndex());
			g.setEdgeWeight(g.getEdge(n.getIndex(), n.getIndex()+offset),999999);
			g.setEdgeWeight(g.getEdge(n.getIndex()+offset,n.getIndex()),999999);
		}

		linkcap= DoubleFactory2D.dense.make(netPlan.getNumberOfLinks(), W.getInt(), C.getDouble());
		
		//System.out.println("CREATING PLIGHTPATHS");

		for (int i =0; i<N;i++)
			for (int j =0;j<N;j++)
				if (i!=j)
				{
					if (i==0 && j==18)
					{
						//System.out.println("Ehy");
					}

					nodePair = Pair.of(netPlan.getNode(i) , netPlan.getNode(j));
					//List<Pair<List<Link>, List<Link>>> prova= cpl11.get(nodePair);
					//Pair<List<Link>, List<Link>> prova2=prova.get(0);

					//System.out.println("Disjoint paths between ("+String.valueOf(i)+","+String.valueOf(j)+")");
					
					Map<Link,Double> linkCostMap = new HashedMap<Link, Double>();
					for (int l=0; l<netPlan.getLinks(netPlan.getNetworkLayerDefault()).size();l++)
					{
						linkCostMap.put(netPlan.getLink(l, netPlan.getNetworkLayerDefault()), 1.0);
					}
					
					List<List<Link>> disj = GraphUtils.getTwoLinkDisjointPaths(netPlan.getNodes(), netPlan.getLinks(netPlan.getNetworkLayerDefault()), nodePair.getFirst(), nodePair.getSecond(), linkCostMap);

					//int a= cpl11.get(nodePair).get(0).getFirst().size();
					//int b= cpl11.get(nodePair).get(0).getSecond().size();
					
					int a=disj.get(0).size();
					int b=disj.get(1).size();
					
					
					if (a>0 && b>0)
					{ //there is a plightpath, add it in grooming graph
						g.addEdge(i+offset,j+offset);
						//g.addEdge(j+offset,i+offset);

						double weight1= a+b;
						//double weight2= cpl11.get(nodePair).get(0).getFirst().size()+cpl11.get(nodePair).get(0).getSecond().size();
						g.setEdgeWeight(g.getEdge(i+offset, j+offset), weight1);
						//g.setEdgeWeight(g.getEdge(j+offset,i+offset), cpl11.get(nodePair).get(0).getFirst().size()+cpl11.get(nodePair).get(0).getSecond().size());
						//store link in plightpaths, no need, they are in cpl and cpl11
						
					}

				}
		//now we have the grooming graph ready

		//System.out.println("PLIGHTPATHS CREATED");
		
		//Lets prepare a structure to store wavelenghts of plightpaths

		//<we,List<<Link,<wl,cap>>>
		Map<DefaultWeightedEdge,ArrayList<Pair<Link,ArrayList<Pair<Integer,Double>>>>> plinkcapw=new HashMap<DefaultWeightedEdge,ArrayList<Pair<Link,ArrayList<Pair<Integer,Double>>>>>();
		Map<DefaultWeightedEdge,ArrayList<Pair<Link,ArrayList<Pair<Integer,Double>>>>> plinkcapb=new HashMap<DefaultWeightedEdge,ArrayList<Pair<Link,ArrayList<Pair<Integer,Double>>>>>();
		//Generate traffic requests


		
		//use b_suv and Y_sun to understand if there is a traffic request
		b_suv.getNonZeros(ind, val);
		pos= new int[]{Ns,Nf,Nf};
		for (int i=0; i< ind.size();i++) //traffic requests
		{	
			int s=DoubleMatrixND.ind2sub(ind.get(i), pos).get(0);
			int u=DoubleMatrixND.ind2sub(ind.get(i), pos).get(1);
			int v=DoubleMatrixND.ind2sub(ind.get(i), pos).get(2);
			int un=(int)Y_su.get(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0),DoubleMatrixND.ind2sub(ind.get(i), pos).get(1));
			int vn=(int)Y_su.get(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0),DoubleMatrixND.ind2sub(ind.get(i), pos).get(2));
			if (un!=vn)//traffic is generated
			{   

				//if (s==19 && u==2 && v==3)
					//System.out.println("ehy");
				nodePair = Pair.of(netPlan.getNode(un) , netPlan.getNode(vn));
				//increase genttraff
				gentraff=gentraff+ b_suv.get(new int[] {s,u,v});
				//we know startin and end points un e vn, we know bandwidth b_suv(s,u,v)
				//we can groom it
				//find shortest path on grooming graph

				System.out.println( String.format("Serving slice b(%1$d,%2$d,%3$d)",s,u,v));

				KShortestPaths<Integer, DefaultWeightedEdge> ksp = new KShortestPaths<Integer, DefaultWeightedEdge>(g, un,100);
				List<GraphPath<Integer, DefaultWeightedEdge>> sps= ksp.getPaths(vn);

				//DijkstraShortestPath<Integer, DefaultWeightedEdge> alg = new DijkstraShortestPath<Integer, DefaultWeightedEdge>(g,un,vn);
				//List<DefaultWeightedEdge> list= alg.getPathEdgeList();


				boolean found=false;
				if (sps.size()==0)
					throw new IllegalArgumentException("Something went wrong in finding disjoint paths");

				for (int p=0; p< sps.size() && !found;p++) //paths in grooming graph
				{
					List<DefaultWeightedEdge> edges=sps.get(p).getEdgeList();

					ArrayList<Pair<DefaultWeightedEdge,ArrayList<Pair<Integer,Integer>>>> wlw=new ArrayList<Pair<DefaultWeightedEdge,ArrayList<Pair<Integer,Integer>>>>(); //<plightpat,<link,wavel>
					ArrayList<Pair<DefaultWeightedEdge,ArrayList<Pair<Integer,Integer>>>> wlb=new ArrayList<Pair<DefaultWeightedEdge,ArrayList<Pair<Integer,Integer>>>>(); //<plightpat,<link,wavel>
					Map<DefaultWeightedEdge,ArrayList<Pair<Link,ArrayList<Pair<Integer,Double>>>>> plinkcapw_original= new HashMap<DefaultWeightedEdge,ArrayList<Pair<Link,ArrayList<Pair<Integer,Double>>>>>(plinkcapw);
					Map<DefaultWeightedEdge,ArrayList<Pair<Link,ArrayList<Pair<Integer,Double>>>>> plinkcapb_original= new HashMap<DefaultWeightedEdge,ArrayList<Pair<Link,ArrayList<Pair<Integer,Double>>>>>(plinkcapb);
					DoubleMatrix2D linkcap_original= DoubleFactory2D.dense.make(linkcap.toArray());
					int w_wl_original=w_wl;
					int p_wl_original=p_wl;

					//DirectedGraph<Integer, DefaultWeightedEdge> revGraph = new EdgeReversedGraph<Integer, DefaultWeightedEdge>(g);
					//DirectedGraph<Integer, DefaultWeightedEdge> g_original = new EdgeReversedGraph<Integer, DefaultWeightedEdge>(revGraph);

					SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_original = (SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>) g.clone();


					boolean enoughres=true;
					for (int pl=0; pl<edges.size() && enoughres;pl++ )// link (plightpaths) in grooming graph
					{	

						//check if pl is a plightpath or a tpedge
						int source= g.getEdgeSource(edges.get(pl));
						int dest= g.getEdgeTarget(edges.get(pl));
						ArrayList<Pair<Integer,Integer>> wlsw = null; //linkindex,wavelength
						ArrayList<Pair<Integer,Integer>> wlsb= new ArrayList<Pair<Integer,Integer>>();//wavelengths in each link
						if ((source >=offset && dest >=offset) ||(source <offset && dest <offset)) // it is not a transponder edge
						{
							// check available resources in the plightpath
							//take lightpath
							boolean existing=false;
							if (source >=offset) //new plightpath
							{
								nodePair = Pair.of(netPlan.getNode(source-offset) , netPlan.getNode(dest-offset));
								existing =false;
							}
							else //existing plightpath
							{
								nodePair = Pair.of(netPlan.getNode(source) , netPlan.getNode(dest));
								existing =true;
							}


							if (existing) //existing plightpath // EXISTING PLIGHTPATH //Existing plightpath
							{
								boolean availablew = true;
								boolean availableb = true;
								wlsw= new ArrayList<Pair<Integer,Integer>>();//List<link,wl>
								wlsb= new ArrayList<Pair<Integer,Integer>>();

								DefaultWeightedEdge newedge=g.getEdge(source, dest);

								//scan working path of plightpath

								for (int h=0;h<plinkcapw.get(newedge).size() && availablew;h++) //sono link fisici
								{
									boolean foundwl=false;
									for (int w=0;w<plinkcapw.get(newedge).get(h).getSecond().size() &&!foundwl;w++) //sono wavelenghts
									{
										if (plinkcapw.get(newedge).get(h).getSecond().get(w).getSecond()>=b_suv.get(new int[] {s,u,v}) && D.get(s)==0)
										{
											foundwl= true;
											//i have to store link and wavelength of the plink
											wlsw.add(Pair.of((Integer)plinkcapw.get(newedge).get(h).getFirst().getIndex(),plinkcapw.get(newedge).get(h).getSecond().get(w).getFirst()));										
										}
										if (plinkcapw.get(newedge).get(h).getSecond().get(w).getSecond()>=C.getDouble() && D.get(s)==1)
										{
											foundwl= true;
											//i have to store link and wavelength of the plink
											wlsw.add(Pair.of((Integer)plinkcapw.get(newedge).get(h).getFirst().getIndex(),plinkcapw.get(newedge).get(h).getSecond().get(w).getFirst()));										
										}

									}
									if (!foundwl)//we exit the for l cycle
										availablew=false;
									else
									{ //we can scan next link l

									}
								}



								//scan backup path of plightpath
								if (availablew)//otherwise no need to scan backup
									for (int h=0;h<plinkcapb.get(newedge).size() && availableb;h++) //sono link fisici
									{
										boolean foundwl=false;
										for (int w=0;w<plinkcapb.get(newedge).get(h).getSecond().size() &&!foundwl;w++) //sono wavelenghts
										{
											if (plinkcapb.get(newedge).get(h).getSecond().get(w).getSecond()>=b_suv.get(new int[] {s,u,v}) && D.get(s)==0)
											{
												foundwl= true;
												//i have to store link and wavelength of the plink
												wlsb.add(Pair.of((Integer)plinkcapb.get(newedge).get(h).getFirst().getIndex(),plinkcapb.get(newedge).get(h).getSecond().get(w).getFirst()));										
											}
											if (plinkcapb.get(newedge).get(h).getSecond().get(w).getSecond()>=C.getDouble() && D.get(s)==1)
											{
												foundwl= true;
												//i have to store link and wavelength of the plink
												wlsb.add(Pair.of((Integer)plinkcapb.get(newedge).get(h).getFirst().getIndex(),plinkcapb.get(newedge).get(h).getSecond().get(w).getFirst()));										
											}

										}
										if (!foundwl)//we exit the for l cycle
											availableb=false;
										else
										{ //we can scan next link l

										}
									}

								if (availablew && availableb) //substract resources from plightpath
								{
									//workingpath
									// no need to create the plinkcap, just subtract capacity

									for (int lin=0;lin<wlsw.size();lin++)
									{
										Link linknet=netPlan.getLink(wlsw.get(lin).getFirst());
										for (int h=0;h<plinkcapw.get(newedge).size();h++) //sono link fisici
										{
											if (plinkcapw.get(newedge).get(h).getFirst()==linknet)
											{
												for (int w=0;w<plinkcapw.get(newedge).get(h).getSecond().size();w++)
												{
													if (wlsw.get(lin).getSecond()== plinkcapw.get(newedge).get(h).getSecond().get(w).getFirst())
													{
														if (D.get(s)==0)
															plinkcapw.get(newedge).get(h).getSecond().get(w).setSecond(plinkcapw.get(newedge).get(h).getSecond().get(w).getSecond()-b_suv.get(new int[] {s,u,v}));
														else if (D.get(s)==1)
															plinkcapw.get(newedge).get(h).getSecond().get(w).setSecond(plinkcapw.get(newedge).get(h).getSecond().get(w).getSecond()-C.getDouble());

													}
												}
											}
										}

									}

									//backuppath
									// no need to create the plinkcap, just subtract capacity

									for (int lin=0;lin<wlsb.size();lin++)
									{
										Link linknet=netPlan.getLink(wlsb.get(lin).getFirst());
										for (int h=0;h<plinkcapb.get(newedge).size();h++) //sono link fisici
										{
											if (plinkcapb.get(newedge).get(h).getFirst()==linknet)
											{
												for (int w=0;w<plinkcapb.get(newedge).get(h).getSecond().size();w++)
												{
													if (wlsb.get(lin).getSecond()== plinkcapb.get(newedge).get(h).getSecond().get(w).getFirst())
													{
														if (D.get(s)==0)
															plinkcapb.get(newedge).get(h).getSecond().get(w).setSecond(plinkcapb.get(newedge).get(h).getSecond().get(w).getSecond()-b_suv.get(new int[] {s,u,v}));
														else if (D.get(s)==1)
															plinkcapb.get(newedge).get(h).getSecond().get(w).setSecond(plinkcapb.get(newedge).get(h).getSecond().get(w).getSecond()-C.getDouble());


													}
												}
											}
										}

									}

								}
								else //we have to consider another shortest path in grooming graph
								{
									enoughres=false;
									//restore stato precedente della rete
									linkcap=linkcap_original;
									plinkcapb=plinkcapb_original;
									plinkcapw=plinkcapw_original;
									g= g_original;
									p_wl=p_wl_original;
									w_wl=w_wl_original;


								}	

							}
							else //new plightpath 
							{
								// IF NOT EXISTING IS OK
								
								//List<Link> working= cpl11.get(nodePair).get(0).getFirst();
								//List<Link> backup= cpl11.get(nodePair).get(0).getSecond();
								
								
								//calc disj
								Map<Link,Double> linkCostMap = new HashedMap<Link, Double>();
								
								//netPlan.getlin
								for (int l=0; l<netPlan.getLinks(netPlan.getNetworkLayerDefault()).size();l++)
								{
									linkCostMap.put(netPlan.getLink(l, netPlan.getNetworkLayerDefault()), 1.0);
								}
								
								List<List<Link>> disj = GraphUtils.getTwoLinkDisjointPaths(netPlan.getNodes(), netPlan.getLinks(netPlan.getNetworkLayerDefault()), nodePair.getFirst(), nodePair.getSecond(), linkCostMap);
								List<Link> working = disj.get(0);
								List<Link> backup= disj.get(1);
								
								
								
								
								boolean availablew = true;
								boolean availableb = true;
								wlsw= new ArrayList<Pair<Integer,Integer>>();//List<link,wl>
								wlsb= new ArrayList<Pair<Integer,Integer>>();

								//scan working path
								for (int l=0; l<working.size() && availablew ; l++) // for each link in working path
								{
									
									System.out.println("Scanning for working link: ("+String.valueOf(working.get(l).getIndex())+")");

									boolean foundwl=false;
									for (int w=0; w< W.getInt() &&!foundwl ;w++)
									{

										if(working.get(l).getIndex()==17 && w==1)
										{System.out.println("Ehy b_suv="+String.valueOf(b_suv.get(new int[] {s,u,v})));
										System.out.println(String.format("Ehy actual capacity of l %1$d w %2$d is %3$f",working.get(l).getIndex(),w,linkcap.get(working.get(l).getIndex(), w)));
										}

										if(linkcap.get(working.get(l).getIndex(),w)>= b_suv.get(new int[] {s,u,v}) && D.get(s)==0)//add case D(s)=1
										{
											foundwl= true;
											wlsw.add(Pair.of(working.get(l).getIndex(), w));
										}
										if(linkcap.get(working.get(l).getIndex(),w)>= C.getDouble() && D.get(s)==1)//add case D(s)=1
										{
											foundwl= true;
											wlsw.add(Pair.of(working.get(l).getIndex(), w));
										}
									}

									if (!foundwl)//we exit the for l cycle
										availablew=false;
									else
									{ //we can scan next link l

									}
								}

								//scan backup path
								if (availablew)//otherwise no need to scan backup path
									for (int l=0; l<backup.size() && availableb ; l++) // for each link in working path
									{

										boolean foundwl=false;
										
										System.out.println("Scanning for backup link: ("+String.valueOf(backup.get(l).getIndex())+")");
										
										for (int w=0; w< W.getInt() &&!foundwl ;w++)
										{
											if(backup.get(l).getIndex()==17 && w==1)
											{System.out.println("Ehy b_suv="+String.valueOf(b_suv.get(new int[] {s,u,v})));
											System.out.println(String.format("Ehy actual capacity of l %1$d w %2$d is %3$f",backup.get(l).getIndex(),w,linkcap.get(backup.get(l).getIndex(), w)));
											}

											if(linkcap.get(backup.get(l).getIndex(),w)>= b_suv.get(new int[] {s,u,v})&& D.get(s)==0)
											{
												foundwl= true;
												wlsb.add(Pair.of(backup.get(l).getIndex(), w));
											}
											if(linkcap.get(backup.get(l).getIndex(),w)>= C.getDouble() && D.get(s)==1)//add case D(s)=1
											{

												foundwl= true;
												wlsb.add(Pair.of(backup.get(l).getIndex(), w));
											}
										}

										if (!foundwl)//we exit the for l cycle 
											availableb=false;
										else
										{ //we can scan next link l

										}
									}
								//WLW AND WLB ARE SUFFICIENT? MAYBE FOR NEW PLIGHTPATH YES
								if (availablew && availableb)
								{
									wlw.add(Pair.of(edges.get(pl), wlsw));
									wlb.add(Pair.of(edges.get(pl), wlsb));
									//do the magic please!

									//create the plightpath
									g.addEdge(source-offset, dest-offset);
									//g.setEdgeWeight(g.getEdge(source-offset, dest-offset), g.getEdgeWeight(g.getEdge(source-offset, dest-offset))/1000);
									double newedgeweight = g.getEdgeWeight(g.getEdge(source, dest))- (g.getEdgeWeight(g.getEdge(source, dest))/1000);
									System.out.println(String.format("Added plightpath (%1$d,%2$d):%3$f ",source-offset,dest-offset,newedgeweight ));

									//wrong
									//g.setEdgeWeight(g.getEdge(source-offset, dest-offset),g.getEdgeWeight(g.getEdge(source-offset, dest-offset))- (g.getEdgeWeight(g.getEdge(source-offset, dest-offset))/1000));

									//correct
									g.setEdgeWeight(g.getEdge(source-offset, dest-offset),g.getEdgeWeight(g.getEdge(source, dest))- (g.getEdgeWeight(g.getEdge(source, dest))/1000));

									DefaultWeightedEdge newedge=g.getEdge(source-offset, dest-offset);

									//workingpath
									//create the plinkcap
									ArrayList<Pair<Link,ArrayList<Pair<Integer,Double>>>> toadd =new ArrayList<Pair<Link,ArrayList<Pair<Integer,Double>>>>();
									for (int lin=0;lin<wlsw.size();lin++)
									{
										Link linknet=netPlan.getLink(wlsw.get(lin).getFirst());
										ArrayList<Pair<Integer,Double>> toaddwl=new ArrayList<Pair<Integer,Double>>();
										if (D.get(s)==0)
										{toaddwl.add(Pair.of(wlsw.get(lin).getSecond(), C.getDouble()-b_suv.get(new int[] {s,u,v})));
										w_wl++;
										System.out.println("Working uses: ("+String.valueOf(linknet.getIndex())+","+String.valueOf(wlsw.get(lin).getSecond())+")");

										}
										else if (D.get(s)==1)
										{	toaddwl.add(Pair.of(wlsw.get(lin).getSecond(), C.getDouble()-C.getDouble()));
										w_wl++;
										System.out.println("Working uses: ("+String.valueOf(linknet.getIndex())+","+String.valueOf(wlsw.get(lin).getSecond())+")");

										}

										toadd.add(Pair.of(linknet,toaddwl));
										//remove capacity from substratenetwork						
										linkcap.set(wlsw.get(lin).getFirst(), wlsw.get(lin).getSecond(), 0);

										if (wlsw.get(lin).getFirst()==17 && wlsw.get(lin).getSecond()==1)
											System.out.println(String.format("Ehy actual capacity of l %1$d w %2$d is %3$f",wlsw.get(lin).getFirst(),wlsw.get(lin).getSecond(), linkcap.get(wlsw.get(lin).getFirst(), wlsw.get(lin).getSecond()) ));

									}
									plinkcapw.put(newedge, toadd);

									//backuppath
									//create the plinkcap
									toadd =new ArrayList<Pair<Link,ArrayList<Pair<Integer,Double>>>>();
									for (int lin=0;lin<wlsb.size();lin++)
									{
										Link linknet=netPlan.getLink(wlsb.get(lin).getFirst());
										ArrayList<Pair<Integer,Double>> toaddwl=new ArrayList<Pair<Integer,Double>>();
										if (D.get(s)==0)
										{toaddwl.add(Pair.of(wlsb.get(lin).getSecond(), C.getDouble()-b_suv.get(new int[] {s,u,v})));
										p_wl++;
										System.out.println("Backup uses: ("+String.valueOf(linknet.getIndex())+","+String.valueOf(wlsb.get(lin).getSecond())+")");

										}
										else if (D.get(s)==1)
										{	toaddwl.add(Pair.of(wlsb.get(lin).getSecond(), C.getDouble()-C.getDouble()));
										p_wl++;
										System.out.println("Backup uses: ("+String.valueOf(linknet.getIndex())+","+String.valueOf(wlsb.get(lin).getSecond())+")");

										}
										toadd.add(Pair.of(linknet,toaddwl));
										//remove capacity from substratenetwork
										if (wlsb.get(lin).getFirst()==17 && wlsb.get(lin).getSecond()==1)
											System.out.println("Ehy");
										linkcap.set(wlsb.get(lin).getFirst(), wlsb.get(lin).getSecond(), 0);
									}
									plinkcapb.put(newedge, toadd);

								}
								else //we have to consider another shortest path in grooming graph
								{
									enoughres=false;
									//restore stato precedente della rete
									linkcap=linkcap_original;
									plinkcapb=plinkcapb_original;
									plinkcapw=plinkcapw_original;
									g=(SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>) g_original;
									p_wl=p_wl_original;
									w_wl=w_wl_original;
								}	
							}

						}

						//end of plightpath analysis


					}

					if (enoughres)
					{
						found=true;
					}

				}
				if (!found)
				{System.out.println("Link Mapping error"); throw new IllegalArgumentException("Something wrong in link mapping");}


			}
		}

	}
	
	
	
	public void LinkMapping2(NetPlan netPlan)
	{
		//prerequisite: have b_suv ready

		//finds k solutions per traffic request... do we need them?

		w_wl=0;
		p_wl=0;
		gentraff=0;

		//IN V6 USES K-SURBALLES ALGORITHM
		IntArrayList ind= new IntArrayList();
		IntArrayList ind2= new IntArrayList();

		int[] pos= new int[]{Ns,Nf,N};
		DoubleArrayList val= new DoubleArrayList();

		//Convert Y_sun in Y(s,u)
		DoubleMatrix2D Y_su=DoubleFactory2D.sparse.make(Ns,Nf);
		Y_sun.getNonZeros(ind, val);

		for (int i=0; i< ind.size();i++)
		{	 
			Y_su.set(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0), DoubleMatrixND.ind2sub(ind.get(i), pos).get(1), DoubleMatrixND.ind2sub(ind.get(i), pos).get(2));
		}

		//System.out.println("DISJOINT PATHS CALCULATED");
		Pair<Node,Node> nodePair;

		linkcap= DoubleFactory2D.dense.make(netPlan.getNumberOfLinks(), W.getInt(), C.getDouble());

		//use b_suv and Y_sun to understand if there is a traffic request
		b_suv.getNonZeros(ind, val);
		pos= new int[]{Ns,Nf,Nf};
		for (int i=0; i< ind.size();i++) //traffic requests
		{				
			//			
			//if (i==1)
				//System.out.print("Error?");
			
			System.out.print(i+"-th request");
			PLightpathMap pligthpathsmap = new PLightpathMap();

			int s=DoubleMatrixND.ind2sub(ind.get(i), pos).get(0);
			int u=DoubleMatrixND.ind2sub(ind.get(i), pos).get(1);
			int v=DoubleMatrixND.ind2sub(ind.get(i), pos).get(2);
			//System.out.println( String.format("Serving slice b(%1$d,%2$d,%3$d)",s,u,v));
			int un=(int)Y_su.get(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0),DoubleMatrixND.ind2sub(ind.get(i), pos).get(1));
			int vn=(int)Y_su.get(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0),DoubleMatrixND.ind2sub(ind.get(i), pos).get(2));
			if (un!=vn)//traffic is generated
			{   
				//increase genttraff
				gentraff=gentraff+ b_suv.get(new int[] {s,u,v});
				double cap_req= b_suv.get(new int[] {s,u,v});

				//boolean request_served=false;

				//now we have to create the plightpathgrooming graph
				//for each nodepair we run ksurballe and look for the minimum cost plightpath with enough capacity

				//create the graph suitable for the algorithm
				SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k= new SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);

				//SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k_debug= new SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);

				//SimpleDirectedGraph<Integer, DefaultEdge> g_k2=new SimpleDirectedGraph<Integer,DefaultEdge>()
				//SimpleDirectedGraph<Integer, DefaultEdge> g_k= new SimpleDirectedGraph<Integer, DefaultEdge>(DefaultEdge.class);
				for (int l=0; l<netPlan.getLinks(netPlan.getNetworkLayerDefault()).size();l++)
				{

					//System.out.println("Link: "+ netPlan.getLink(l).getIndex());

					int sourcenode=(int) netPlan.getLink(l).getOriginNode().getIndex();
					//Node sourcenode= netPlan.getLink(l).getOriginNode();
					g_k.addVertex(sourcenode);

					int destnode=(int) netPlan.getLink(l).getDestinationNode().getIndex();
					//Node destnode= netPlan.getLink(l).getDestinationNode();
					g_k.addVertex(destnode);

					//System.out.println("Source: "+sourcenode+ " Destination: "+ destnode);
					g_k.addEdge(sourcenode,destnode);
					//g_k.addEdge(destnode,sourcenode);
					g_k.setEdgeWeight(g_k.getEdge(sourcenode, destnode), 1.0);
				}

				//DEBUG: facciamo un grafo completo
				//				for (int h=0; h < N; h++)
				//					for (int k=0; k < N; k++)
				//					{
				//						if (h!=k)
				//						//run ksuurballe
				//						{
				//							g_k_debug.addVertex(h);
				//							g_k_debug.addVertex(k);
				//							g_k_debug.addEdge(h,k);
				//						}
				//					}
				//				


				//TO PRINT THE GRAPH
				//				JFrame frame = new JFrame("NewGraph");
				//				frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
				//				JGraphXAdapter<Integer, DefaultWeightedEdge> graphAdapter = new JGraphXAdapter<Integer, DefaultWeightedEdge>(g_k);
				//				//JGraphXAdapter<Integer, DefaultWeightedEdge> graphAdapter = new JGraphXAdapter<Integer, DefaultWeightedEdge>(g_k_debug);
				//				mxIGraphLayout layout = new mxCircleLayout(graphAdapter);
				//				layout.execute(graphAdapter.getDefaultParent());
				//				frame.add(new mxGraphComponent(graphAdapter));
				//				frame.pack();
				//				frame.setLocationByPlatform(true);
				//				frame.setVisible(true);

				//create the plightpathgraph
				SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> plightpath_graph= new SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);
				//add transponder edges
				int offset=1000;
				//we need a temporary link capacity occupation to avoid conflicts between lightpaths
				linkcap_temp= DoubleFactory2D.dense.make(linkcap.toArray());
				for (int h=0; h < N; h++)
					for (int k=0; k < N; k++)
					{
						if (h!=k)
							//run ksuurballe
						{
							plightpath_graph.addVertex(h);
							plightpath_graph.addVertex(k);
							plightpath_graph.addVertex(h+offset);
							plightpath_graph.addVertex(k+offset);
							plightpath_graph.addEdge(h, h+offset);
							plightpath_graph.addEdge(h+offset,h);
							plightpath_graph.addEdge(k, k+offset);
							plightpath_graph.addEdge(k+offset,k);
							plightpath_graph.setEdgeWeight(plightpath_graph.getEdge(h, h+offset), 9999);
							plightpath_graph.setEdgeWeight(plightpath_graph.getEdge(k, k+offset), 9999);
							plightpath_graph.setEdgeWeight(plightpath_graph.getEdge( h+offset,h), 9999);
							plightpath_graph.setEdgeWeight(plightpath_graph.getEdge(k+offset,k), 9999);

//							if (h==34 && k==27)
//							{
//								System.out.println("34,27");
//							}
							PLightpath pl;
							//if (OP.getString().contentEquals("VWP"))
							pl= feasiblePlightpath(h, k, cap_req, g_k, netPlan);
							
							

							if (pl!=null) //at least working and backup
							{
								plightpath_graph.addEdge(h, k);
								double weight= pl.get(0).size() +pl.get(1).size();
//								if ((h==34 && k==28) || (h==28 && k==27 ) || (h==34 && k==27) )
//									System.out.println("Check weight of ("+h+","+k+"):"+weight);
								plightpath_graph.setEdgeWeight(plightpath_graph.getEdge(h, k), weight);								
								//we have to store somewhere the plightpath
								pligthpathsmap.put(plightpath_graph.getEdge(h, k), pl);
							}

							//scan for existing plightpaths
							PLightpath existing_pl = feasibleExistingPligthpath (h+offset,k+offset, cap_req, netPlan);
							if (existing_pl!=null)
							{
								
								double weight = (existing_pl.get(0).size()+existing_pl.get(1).size())-(existing_pl.get(0).size()+existing_pl.get(1).size())/1000;
//								if ((h==34 && k==28) || (h==28 && k==27 ) || (h==34 && k==27) )
//									System.out.println("Check weight of ("+h+","+k+"):"+weight);
								plightpath_graph.addEdge(h+offset, k+offset);
								plightpath_graph.setEdgeWeight(plightpath_graph.getEdge(h+offset, k+offset), weight);
								pligthpathsmap.put(plightpath_graph.getEdge(h+offset, k+offset), existing_pl);
							}
						}
					}
				//now we have the plightpath graph

				//calculate shortest path in plightpathgraph
				//KShortestPaths<Integer, DefaultWeightedEdge> ksp = new KShortestPaths<Integer, DefaultWeightedEdge>(plightpath_graph, un+1000,2);
				//List<GraphPath<Integer, DefaultWeightedEdge>> sps= ksp.getPaths(vn+1000);
				//in theory we do not need k shortest paths, we know that there is capacity, sure?

				int K=3;

				KShortestSimplePaths<Integer, DefaultWeightedEdge> ksp = new KShortestSimplePaths<Integer, DefaultWeightedEdge>(plightpath_graph);
				List<GraphPath<Integer, DefaultWeightedEdge>> sps= ksp.getPaths(un+1000, vn+1000, K);
				List<DefaultWeightedEdge> edges;//=sps.get(0).getEdgeList();


				//List<DefaultWeightedEdge> edges= DijkstraShortestPath.findPathBetween(plightpath_graph, un+1000, vn+1000);

				boolean served=false;
				DoubleMatrix2D linkcap_original= DoubleFactory2D.dense.make(linkcap.toArray());
				ArrayList<Triple<Pair<Integer, Integer>, PLightpath, Double>>  existing_plightpaths_original= (ArrayList<Triple<Pair<Integer, Integer>, PLightpath, Double>> ) existing_plightpaths.clone();
				
				
				
				for (int k=0; k<K && !served; k++)
				{
					//stat
					if(k>kmax)
						kmax=k;

					edges=sps.get(k).getEdgeList();
					
					boolean feasible=true;
					for(int e=0; e<edges.size() && feasible;e++) 
					//(DefaultWeightedEdge de:edges)//DONE put the index and the condition on served here othewise next edge will put served=true after served=false
					{
						DefaultWeightedEdge de=edges.get(e);
						//check if pl is a plightpath or a tpedge
						int source= plightpath_graph.getEdgeSource(de);
						int dest= plightpath_graph.getEdgeTarget(de);

						if ((source >=offset && dest >=offset) ||(source <offset && dest <offset)) // it is not a transponder edge
						{
							// check available resources in the plightpath
							//take lightpath
							boolean existing=false;
							if (source <offset) //new plightpath 
							{
								nodePair = Pair.of(netPlan.getNode(source) , netPlan.getNode(dest));
								existing =false;
							}
							else //existing plightpath
							{
								nodePair = Pair.of(netPlan.getNode(source-offset) , netPlan.getNode(dest-offset));
								existing =true;
							}

							IntArrayList ind3= new IntArrayList();
							IntArrayList ind4= new IntArrayList();

							DoubleArrayList val2= new DoubleArrayList();


							if (!existing) //New plightpath
							{
								//we have to create a new existing plightpaht
								PLightpath npl=pligthpathsmap.get(de);
								if (!checkConflictOfLightpath(npl))
									{
									reservePLightpathResources(npl);
									}
								else {
									PLightpath new_npl=feasiblePlightpath(source , dest, cap_req, g_k, netPlan);
									pligthpathsmap.put(de, new_npl);
									reservePLightpathResources(new_npl);
									npl=new_npl;
								}

								//create existing lightpath
								if (D.get(s)==0)//Grooming
								{//existing_plightpaths.add(Triple.of(Pair.of(source-offset, dest-offset), npl,C.getDouble()-cap_req));
									existing_plightpaths.add(Triple.of(Pair.of(source+offset, dest+offset), npl,C.getDouble()-cap_req));
									w_wl=w_wl+ npl.get(0).size();
									p_wl=p_wl+ npl.get(1).size();
									linkcap.getNonZeros(ind3, ind4, val2);
									int totwl=linkcap.columns()*linkcap.rows();
									int wl_used= totwl-ind3.size();
									if (w_wl+p_wl != wl_used )
										feasible=false; //avoid conflicts among plightpaths in the same shortest path of the pligthpathgraph
									else
										feasible=true;
								}
								else
								{
									//existing_plightpaths.add(Triple.of(Pair.of(source-offset, dest-offset), npl,0.0));
									existing_plightpaths.add(Triple.of(Pair.of(source+offset, dest+offset), npl,0.0));
									w_wl=w_wl+ npl.get(0).size();
									p_wl=p_wl+ npl.get(1).size();
									linkcap.getNonZeros(ind3, ind4, val2);
									int totwl=linkcap.columns()*linkcap.rows();
									int wl_used= totwl-ind3.size();
									if (w_wl+p_wl != wl_used )
										feasible=false;//avoid conflicts among plightpaths in the same shortest path of the pligthpathgraph
									else
										feasible=true;
								}


							}
							else //Existing plightpath
							{
								//subtract capacity from plightpath
								PLightpath npl=pligthpathsmap.get(de);
								if (D.get(s)==0)//Grooming
								{
									for (int p=0; p< existing_plightpaths.size();p++)
									{
										if (existing_plightpaths.get(p).getSecond()==npl)
										{
											//add check on capacity for feasibility to avoid conflicts among plightpaths in the same shortest path of the pligthpathgraph
											if (existing_plightpaths.get(p).getThird()-cap_req>=0.0 )
											{
											existing_plightpaths.set(p, Triple.of(Pair.of(source, dest), npl, existing_plightpaths.get(p).getThird()-cap_req));
											feasible=true;
											}
											else
											{
												feasible=false;
											}
										}

									}
								}
								//existing_plightpaths.set(Triple.of(Pair.of(source-offset, dest-offset), npl,C.getDouble()-cap_req));
								else
								{
									for (int p=0; p< existing_plightpaths.size();p++)
									{
										if (existing_plightpaths.get(p).getSecond()==npl)
										{
											//add check on capacity for feasibility to avoid conflicts among plightpaths in the same shortest path of the pligthpathgraph
											if (existing_plightpaths.get(p).getThird()!=0.0 )
											{
											existing_plightpaths.set(p, Triple.of(Pair.of(source, dest), npl, 0.0));
											feasible=true;
											}
											else
											{
												feasible=false;
											}
										}

									}
									//existing_plightpaths.add(Triple.of(Pair.of(source-offset, dest-offset), npl,0.0));
								}


							}
						}

					}

					if (!feasible)
					{
						linkcap= linkcap_original;
						existing_plightpaths=existing_plightpaths_original;
						//throw new IllegalArgumentException("Request unserved");
					}
					else
					{
						served=true;
						//System.out.println("Served using the"+k+"-th shortest path in plightpath grooming graph");
						
					}


					//put the end of the for here
				}
				
				if (!served)
				{
					throw new IllegalArgumentException("Request unserved");
				}
			}
		}

	}

	
	public void LinkMapping3(NetPlan netPlan)
	{
		//prerequisite: have b_suv ready

		//finds only 1 solution (NOT K) per traffic request... 

		w_wl=0;
		p_wl=0;
		gentraff=0;

		//IN V6 USES K-SURBALLES ALGORITHM
		IntArrayList ind= new IntArrayList();
		IntArrayList ind2= new IntArrayList();

		int[] pos= new int[]{Ns,Nf,N};
		DoubleArrayList val= new DoubleArrayList();

		//Convert Y_sun in Y(s,u)
		DoubleMatrix2D Y_su=DoubleFactory2D.sparse.make(Ns,Nf);
		Y_sun.getNonZeros(ind, val);

		for (int i=0; i< ind.size();i++)
		{	 
			Y_su.set(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0), DoubleMatrixND.ind2sub(ind.get(i), pos).get(1), DoubleMatrixND.ind2sub(ind.get(i), pos).get(2));
		}

		//System.out.println("DISJOINT PATHS CALCULATED");
		Pair<Node,Node> nodePair;

		linkcap= DoubleFactory2D.dense.make(netPlan.getNumberOfLinks(), W.getInt(), C.getDouble());

		//use b_suv and Y_sun to understand if there is a traffic request
		b_suv.getNonZeros(ind, val);
		pos= new int[]{Ns,Nf,Nf};
		for (int i=0; i< ind.size();i++) //traffic requests
		{				
			//			
			//if (i==78)
			//	System.out.println(i+"-th request");
			//if (i>=75)
				//System.out.println("Stop");
				
			PLightpathMap pligthpathsmap = new PLightpathMap();

			int s=DoubleMatrixND.ind2sub(ind.get(i), pos).get(0);
			int u=DoubleMatrixND.ind2sub(ind.get(i), pos).get(1);
			int v=DoubleMatrixND.ind2sub(ind.get(i), pos).get(2);
			//System.out.println( String.format("Serving slice b(%1$d,%2$d,%3$d)",s,u,v));
			
			//if (s==107 && u== 2 && v==3)
				//System.out.println("Debug");
			
			int un=(int)Y_su.get(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0),DoubleMatrixND.ind2sub(ind.get(i), pos).get(1));
			int vn=(int)Y_su.get(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0),DoubleMatrixND.ind2sub(ind.get(i), pos).get(2));
			if (un!=vn)//traffic is generated
			{   
				//increase genttraff
				gentraff=gentraff+ b_suv.get(new int[] {s,u,v});
				double cap_req= b_suv.get(new int[] {s,u,v});

				//boolean request_served=false;

				//now we have to create the plightpathgrooming graph
				//for each nodepair we run ksurballe and look for the minimum cost plightpath with enough capacity

				//create the graph suitable for the algorithm
				SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k=  makeAvailableGraph(netPlan,cap_req);
				

				//DEBUG: facciamo un grafo completo
				//				for (int h=0; h < N; h++)
				//					for (int k=0; k < N; k++)
				//					{
				//						if (h!=k)
				//						//run ksuurballe
				//						{
				//							g_k_debug.addVertex(h);
				//							g_k_debug.addVertex(k);
				//							g_k_debug.addEdge(h,k);
				//						}
				//					}
				//				



				//create the plightpathgraph
				SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> plightpath_graph= new SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);
				//add transponder edges
				int offset=1000;
				//we need a temporary link capacity occupation to avoid conflicts between lightpaths
				linkcap_temp= DoubleFactory2D.dense.make(linkcap.toArray());
				for (int h=0; h < N; h++)
					for (int k=0; k < N; k++)
					{
						if (h!=k)
							//run ksuurballe
						{
							plightpath_graph.addVertex(h);
							plightpath_graph.addVertex(k);
							plightpath_graph.addVertex(h+offset);
							plightpath_graph.addVertex(k+offset);
							plightpath_graph.addEdge(h, h+offset);
							plightpath_graph.addEdge(h+offset,h);
							plightpath_graph.addEdge(k, k+offset);
							plightpath_graph.addEdge(k+offset,k);
							plightpath_graph.setEdgeWeight(plightpath_graph.getEdge(h, h+offset), 9999);
							plightpath_graph.setEdgeWeight(plightpath_graph.getEdge(k, k+offset), 9999);
							plightpath_graph.setEdgeWeight(plightpath_graph.getEdge( h+offset,h), 9999);
							plightpath_graph.setEdgeWeight(plightpath_graph.getEdge(k+offset,k), 9999);

//							if (h==44 && k==1 && i>=75)
//							{
//								System.out.println("44,1");
//							}
							PLightpath pl;
							//if (OP.getString().contentEquals("VWP"))
//							if ( (h==28 && k==27 && i== 327))
//								{
//								System.out.println("Stop!");
//								JFrame frame = new JFrame("NewGraph");
//								frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//								JGraphXAdapter<Integer, DefaultWeightedEdge> graphAdapter = new JGraphXAdapter<Integer, DefaultWeightedEdge>(g_k);
//								//JGraphXAdapter<Integer, DefaultWeightedEdge> graphAdapter = new JGraphXAdapter<Integer, DefaultWeightedEdge>(g_k_debug);
//								//mxIGraphLayout layout = new mxCircleLayout(graphAdapter);
//								mxIGraphLayout layout = new mxOrganicLayout(graphAdapter); 
//								layout.execute(graphAdapter.getDefaultParent());
//								frame.add(new mxGraphComponent(graphAdapter));
//								frame.pack();
//								frame.setLocationByPlatform(true);
//								frame.setVisible(true);
//								
//								
//								}
							pl= feasiblePlightpath(h, k, cap_req, g_k, netPlan);
							if (pl!=null) //at least working and backup
							{
								plightpath_graph.addEdge(h, k);
								
								double weight= pl.get(0).size() +pl.get(1).size();
//								if ((h==34 && k==28) || (h==28 && k==27 ) || (h==34 && k==27) )
//									System.out.println("Check weight of ("+h+","+k+"):"+weight);
								plightpath_graph.setEdgeWeight(plightpath_graph.getEdge(h, k), weight);								
								//we have to store somewhere the plightpath
								pligthpathsmap.put(plightpath_graph.getEdge(h, k), pl);
							}

							//scan for existing plightpaths
							PLightpath existing_pl = feasibleExistingPligthpath (h+offset,k+offset, cap_req, netPlan);
							if (existing_pl!=null)
							{
								
								double weight = (existing_pl.get(0).size()+existing_pl.get(1).size())-(existing_pl.get(0).size()+existing_pl.get(1).size())/1000;
//								if ((h==34 && k==28) || (h==28 && k==27 ) || (h==34 && k==27) )
//									System.out.println("Check weight of ("+h+","+k+"):"+weight);
								plightpath_graph.addEdge(h+offset, k+offset);
								plightpath_graph.setEdgeWeight(plightpath_graph.getEdge(h+offset, k+offset), weight);
								pligthpathsmap.put(plightpath_graph.getEdge(h+offset, k+offset), existing_pl);
							}
						}
					}
				//now we have the plightpath graph
				//calculate shortest path in plightpathgraph
				//KShortestPaths<Integer, DefaultWeightedEdge> ksp = new KShortestPaths<Integer, DefaultWeightedEdge>(plightpath_graph, un+1000,2);
				//List<GraphPath<Integer, DefaultWeightedEdge>> sps= ksp.getPaths(vn+1000);
				//in theory we do not need k shortest paths, we know that there is capacity, sure?
				//int K=3;
				//KShortestSimplePaths<Integer, DefaultWeightedEdge> ksp = new KShortestSimplePaths<Integer, DefaultWeightedEdge>(plightpath_graph);
				//List<GraphPath<Integer, DefaultWeightedEdge>> sps= ksp.getPaths(un+1000, vn+1000, K);
				//List<DefaultWeightedEdge> edges;//=sps.get(0).getEdgeList();
//				//TO PRINT THE GRAPH
//				if (un==28 && vn==27 && i==327)
//				{
//								JFrame frame = new JFrame("NewGraph");
//								frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//								JGraphXAdapter<Integer, DefaultWeightedEdge> graphAdapter = new JGraphXAdapter<Integer, DefaultWeightedEdge>(plightpath_graph);
//								//JGraphXAdapter<Integer, DefaultWeightedEdge> graphAdapter = new JGraphXAdapter<Integer, DefaultWeightedEdge>(g_k_debug);
//								mxIGraphLayout layout = new mxCircleLayout(graphAdapter);
//								//mxIGraphLayout layout = new mxOrganicLayout(graphAdapter); 
//								layout.execute(graphAdapter.getDefaultParent());
//								frame.add(new mxGraphComponent(graphAdapter));
//								frame.pack();
//								frame.setLocationByPlatform(true);
//								frame.setVisible(true);
//								
//				}
				List<DefaultWeightedEdge> edges= DijkstraShortestPath.findPathBetween(plightpath_graph, un+1000, vn+1000);
				
				//KShortestSimplePaths<Integer, DefaultWeightedEdge> ksp = new KShortestSimplePaths<Integer, DefaultWeightedEdge>(plightpath_graph);
				//List<GraphPath<Integer, DefaultWeightedEdge>> sps= ksp.getPaths(un+1000, vn+1000, 1);
				//List<DefaultWeightedEdge> edges=sps.get(0).getEdgeList();
				
				//boolean served=false;
				//DoubleMatrix2D linkcap_original= DoubleFactory2D.dense.make(linkcap.toArray());
				//ArrayList<Triple<Pair<Integer, Integer>, PLightpath, Double>>  existing_plightpaths_original= (ArrayList<Triple<Pair<Integer, Integer>, PLightpath, Double>> ) existing_plightpaths.clone();
				
				
				
//				for (int k=0; k<K && !served; k++)
//				{
//					//stat
//					if(k>kmax)
//						kmax=k;
//
//					edges=sps.get(k).getEdgeList();
//					
//				if (i>=78)
//					System.out.println("SP between "+un+" and "+vn);
				
				boolean feasible=true;
					for(int e=0; e<edges.size() && feasible;e++) 
					//(DefaultWeightedEdge de:edges)//DONE put the index and the condition on served here othewise next edge will put served=true after served=false
					{
						DefaultWeightedEdge de=edges.get(e);
						//check if pl is a plightpath or a tpedge
						int source= plightpath_graph.getEdgeSource(de);
						int dest= plightpath_graph.getEdgeTarget(de);
//						if(i==79)// && source==2 && dest==1)
//							System.out.println("STOP?");
						if ((source >=offset && dest >=offset) ||(source <offset && dest <offset)) // it is not a transponder edge
						{
							// check available resources in the plightpath
							//take lightpath
							boolean existing=false;
							if (source <offset) //new plightpath 
							{
								nodePair = Pair.of(netPlan.getNode(source) , netPlan.getNode(dest));
								existing =false;
							}
							else //existing plightpath
							{
								nodePair = Pair.of(netPlan.getNode(source-offset) , netPlan.getNode(dest-offset));
								existing =true;
							}
							IntArrayList ind3= new IntArrayList();
							IntArrayList ind4= new IntArrayList();
							DoubleArrayList val2= new DoubleArrayList();
							if (!existing) //New plightpath
							{
								//we have to create a new existing plightpaht
								PLightpath npl=pligthpathsmap.get(de);
								if (!checkConflictOfLightpath(npl))
									{
									reservePLightpathResources(npl);
									}
								else {
									g_k=  makeAvailableGraph(netPlan,cap_req);

//									//TO PRINT THE GRAPH
//									if (i>=78)
//									{
//													JFrame frame = new JFrame("NewGraph");
//													frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//													JGraphXAdapter<Integer, DefaultWeightedEdge> graphAdapter = new JGraphXAdapter<Integer, DefaultWeightedEdge>(g_k);
//													//JGraphXAdapter<Integer, DefaultWeightedEdge> graphAdapter = new JGraphXAdapter<Integer, DefaultWeightedEdge>(g_k_debug);
//													//mxIGraphLayout layout = new mxCircleLayout(graphAdapter);
//													mxIGraphLayout layout = new mxOrganicLayout(graphAdapter); 
//													layout.execute(graphAdapter.getDefaultParent());
//													frame.add(new mxGraphComponent(graphAdapter));
//													frame.pack();
//													frame.setLocationByPlatform(true);
//													frame.setVisible(true);
//									}
									PLightpath new_npl=feasiblePlightpath(source , dest, cap_req, g_k, netPlan);
									pligthpathsmap.put(de, new_npl);
									reservePLightpathResources(new_npl);
									npl=new_npl;
								}
								//create existing lightpath
								if (D.get(s)==0)//Grooming
								{//existing_plightpaths.add(Triple.of(Pair.of(source-offset, dest-offset), npl,C.getDouble()-cap_req));
									existing_plightpaths.add(Triple.of(Pair.of(source+offset, dest+offset), npl,C.getDouble()-cap_req));
									w_wl=w_wl+ npl.get(0).size();
									p_wl=p_wl+ npl.get(1).size();
									linkcap.getNonZeros(ind3, ind4, val2);
									int totwl=linkcap.columns()*linkcap.rows();
									int wl_used= totwl-ind3.size();
									if (w_wl+p_wl != wl_used )
										feasible=false; //avoid conflicts among plightpaths in the same shortest path of the pligthpathgraph
									else
										feasible=true;
								}
								else
								{
									//existing_plightpaths.add(Triple.of(Pair.of(source-offset, dest-offset), npl,0.0));
									existing_plightpaths.add(Triple.of(Pair.of(source+offset, dest+offset), npl,0.0));
									w_wl=w_wl+ npl.get(0).size();
									p_wl=p_wl+ npl.get(1).size();
									linkcap.getNonZeros(ind3, ind4, val2);
									int totwl=linkcap.columns()*linkcap.rows();
									int wl_used= totwl-ind3.size();
									if (w_wl+p_wl != wl_used )
										feasible=false;//avoid conflicts among plightpaths in the same shortest path of the pligthpathgraph
									else
										feasible=true;
								}
							}
							else //Existing plightpath
							{
								//subtract capacity from plightpath
								PLightpath npl=pligthpathsmap.get(de);
								if (D.get(s)==0)//Grooming
								{
									for (int p=0; p< existing_plightpaths.size();p++)
									{
										if (existing_plightpaths.get(p).getSecond()==npl)
										{
											//add check on capacity for feasibility to avoid conflicts among plightpaths in the same shortest path of the pligthpathgraph
											if (existing_plightpaths.get(p).getThird()-cap_req>=0.0 )
											{
											existing_plightpaths.set(p, Triple.of(Pair.of(source, dest), npl, existing_plightpaths.get(p).getThird()-cap_req));
											feasible=true;
											}
											else
											{
												feasible=false;
											}
										}

									}
								}
								//existing_plightpaths.set(Triple.of(Pair.of(source-offset, dest-offset), npl,C.getDouble()-cap_req));
								else
								{
									for (int p=0; p< existing_plightpaths.size();p++)
									{
										if (existing_plightpaths.get(p).getSecond()==npl)
										{
											//add check on capacity for feasibility to avoid conflicts among plightpaths in the same shortest path of the pligthpathgraph
											if (existing_plightpaths.get(p).getThird()!=0.0 )
											{
											existing_plightpaths.set(p, Triple.of(Pair.of(source, dest), npl, 0.0));
											feasible=true;
											}
											else
											{
												feasible=false;
											}
										}

									}
									//existing_plightpaths.add(Triple.of(Pair.of(source-offset, dest-offset), npl,0.0));
								}


							}
						}

					}

					if (!feasible)
					{
						//linkcap= linkcap_original;
						//existing_plightpaths=existing_plightpaths_original;
						throw new IllegalArgumentException("Request unserved");
					}
//					else
//					{
//						served=true;
//						//System.out.println("Served using the"+k+"-th shortest path in plightpath grooming graph");
//						
//					}
//
//
//					//put the end of the for here
//				}
//				
//				if (!served)
//				{
//					throw new IllegalArgumentException("Request unserved");
//				}
			}
		}

	}

	public SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> updateConnectivityGraphWithPlacement (NetPlan netPlan, SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k_without_placement)
	{
		SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k_with_placement = new SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		
		//g_k_with_placement= (SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>) g_k_without_placement.clone();
		g_k_with_placement= cloneGraph(netPlan, g_k_without_placement);
		
		for (int i=0; i<N ;i++)
		{
			//System.out.println("Updating incoming links of node "+i);
			if (activeNodes.contains(i))
			{
				//System.out.println("I active node "+i);
			Set<DefaultWeightedEdge> incoming= g_k_with_placement.incomingEdgesOf(i);
			Iterator it =incoming.iterator();
			while (it.hasNext())
			{
				//incoming.iterator().next();
				
				DefaultWeightedEdge edge=(DefaultWeightedEdge) it.next();
			
				//System.out.println("Incoming edge "+ edge.toString());
				
				int source=g_k_with_placement.getEdgeSource(edge);
				int target=g_k_with_placement.getEdgeTarget(edge);
				//System.out.println(g_k_with_placement.getEdgeWeight(edge));
				DefaultWeightedEdge edgetoset= g_k_with_placement.getEdge(source, target);
				double hasfunctions= (Capn_native.get(i)-Capn.get(i))>0?1.0:0.0;
				double weight= g_k_with_placement.getEdgeWeight(edge)-(hasfunctions/2.0);
				//weight=weight*1000;
				Double truncatedDouble = BigDecimal.valueOf(weight).setScale(3, RoundingMode.HALF_UP).doubleValue();
				//System.out.println(edgetoset);
				//System.out.println(String.valueOf(truncatedDouble));
				
				g_k_with_placement.setEdgeWeight(edgetoset, weight);
				//g_k_with_placement.setEdgeWeight(edgetoset, g_k_with_placement.getEdgeWeight(edge)-(1.0/10.0)-(hasfunctions/10.0));
				//g_k_with_placement.setEdgeWeight(edgetoset, 999);
				
				//g_k_with_placement.setEdgeWeight(g_k_with_placement.getEdge(source, target), g_k_with_placement.getEdgeWeight(edge)-(1/10));
				
				//System.out.println(g_k_with_placement.getEdgeWeight(edgetoset));

				
			}	
			}
		}
		//g_k_with_placement.removeEdge(4, 1);
		
		
		return g_k_with_placement;
	}
	

	
	
	public SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> cloneGraph (NetPlan netPlan, SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> old)
	{
		SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> newgraph = new SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		
		Iterator vert_it = old.vertexSet().iterator();
		while (vert_it.hasNext())
		{
			Integer v= (Integer) vert_it.next();
			newgraph.addVertex(v);
		}
		
		Iterator edge_it = old.edgeSet().iterator();
		while (edge_it.hasNext())
		{
			DefaultWeightedEdge edge = (DefaultWeightedEdge) edge_it.next();
			int source= old.getEdgeSource(edge);
			int target= old.getEdgeTarget(edge);
			double weight= old.getEdgeWeight(edge);
			newgraph.addEdge(source,target);
			newgraph.setEdgeWeight(newgraph.getEdge(source, target), weight);
			
		}
		
		
		return newgraph;
	} 

	
	public SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> makeConnectivityGraph (NetPlan netPlan)
	{
		SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k= new SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);

		
		for (int l=0; l<netPlan.getLinks(netPlan.getNetworkLayerDefault()).size();l++)
		{
			//System.out.println("Link: "+ netPlan.getLink(l).getIndex());
			int sourcenode=(int) netPlan.getLink(l).getOriginNode().getIndex();
			//Node sourcenode= netPlan.getLink(l).getOriginNode();
			g_k.addVertex(sourcenode);
			int destnode=(int) netPlan.getLink(l).getDestinationNode().getIndex();
			//Node destnode= netPlan.getLink(l).getDestinationNode();
			g_k.addVertex(destnode);	
//			if (checkLinkCapacity(netPlan.getLink(l),cap_req))
//			{
//				//System.out.println("Source: "+sourcenode+ " Destination: "+ destnode);
				g_k.addEdge(sourcenode,destnode);
//				//g_k.addEdge(destnode,sourcenode);
//				g_k.setEdgeWeight(g_k.getEdge(sourcenode, destnode), 1.0);
//			}
		}
		
		SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k_sur= new SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);

//		FileOutputStream fileOut;
//		try {
//			fileOut = new FileOutputStream("Graph_Sur");
//			ObjectOutputStream objectOut = new ObjectOutputStream(fileOut);	
//			objectOut.writeObject(g_k);
//			objectOut.close();
//			fileOut.close();
//			
//			FileInputStream fi = new FileInputStream(new File("Graph_Sur"));
//			ObjectInputStream oi = new ObjectInputStream(fi);
//			SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k_loaded= (SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>) oi.readObject();
//			g_k=g_k_loaded;
//			oi.close();
//			fi.close();
//
//		} catch (FileNotFoundException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (ClassNotFoundException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		
		
		
		for (int i=0; i<N;i++)
			for (int j=0;j<N;j++)
			{
				if (i!=j)
				{
					g_k_sur.addVertex(netPlan.getNode(i).getIndex());
					g_k_sur.addVertex(netPlan.getNode(j).getIndex());
					
					
					SuurballeKDisjointShortestPaths<Integer,DefaultWeightedEdge> sur= new SuurballeKDisjointShortestPaths<Integer,DefaultWeightedEdge>(g_k);
					List<GraphPath<Integer,DefaultWeightedEdge>> ksp=sur.getPaths(netPlan.getNode(i).getIndex() ,netPlan.getNode(j).getIndex() , 999); 
					
					assert i==netPlan.getNode(i).getIndex() : "Node indexing error!";
					
					
					if (ksp.size()>1)
					{
						DefaultWeightedEdge toadd=new DefaultWeightedEdge();
						g_k_sur.addEdge(i, j, toadd );
						Double weight= Double.valueOf(ksp.get(0).getEdgeList().size()+ksp.get(1).getEdgeList().size());
						
						
						//weight=weight*1000;
						Double truncatedDouble = BigDecimal.valueOf(weight).setScale(3, RoundingMode.HALF_UP).doubleValue();
//						if ((i==27 && j==34) || (i==34 && j==38) )
//						{
//							System.out.println("--------PLACEMENT-----");
//							System.out.println("("+i+","+j+") weights: "+truncatedDouble.toString());
//							System.out.println("Working:");
//							System.out.println(ksp.get(0).getVertexList());
//							System.out.println(ksp.get(0).getWeight());
//							System.out.println("Backup:");
//							System.out.println(ksp.get(1).getVertexList());
//							System.out.println(ksp.get(1).getWeight());
//						}
						g_k_sur.setEdgeWeight(toadd,truncatedDouble );						
						
					}
				}
				
				
			}
		
		
		//introduce distance to nearest core as a metric
		for (int i=0; i<N ;i++)
		{ 
			if(! core_capable.contains(i)) 
			{
			double dist_min=99999999999.0;
			double dist=99999999999.0;
			for (int j=0; j<core_capable.size();j++)
			{
				dist= calculate_distance_sur(g_k_sur,i,core_capable.get(j));
				if (dist<dist_min)
					dist_min=dist;
			}
			Set<DefaultWeightedEdge> incoming= g_k_sur.incomingEdgesOf(i);
			Iterator it =incoming.iterator();
			
			//System.out.println(dist_min);
			while (it.hasNext())
			{
				DefaultWeightedEdge edge=(DefaultWeightedEdge) it.next();
				int source=g_k_sur.getEdgeSource(edge);
				int target=g_k_sur.getEdgeTarget(edge);
				DefaultWeightedEdge edgetoset= g_k_sur.getEdge(source, target);
//				if ( g_k_sur.getEdgeWeight(edge)-(1/dist_min)<0.0)
//					System.out.println("Negative edge");
				double weight=g_k_sur.getEdgeWeight(edge)-(1/dist_min);
				
				//double weight=g_k_sur.getEdgeWeight(edge);//-(0.1/dist_min);
				
				
				//System.out.println("Dj \\ 20 = "+String.valueOf(dist_min/10));
				//double weight=1*g_k_sur.getEdgeWeight(edge)+g_k_sur.getEdgeWeight(edge)*(dist_min);
				//weight=weight*1000;
				Double truncatedDouble = BigDecimal.valueOf(weight).setScale(3, RoundingMode.HALF_UP).doubleValue();
				//System.out.println(edge);
				//System.out.println(String.valueOf(truncatedDouble));
				g_k_sur.setEdgeWeight(edgetoset, weight);
			}
		}
		}
		
		
		return g_k_sur;
		
	}
	
	public double calculate_distance_sur (SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g, int source, int dest)
	{
		
		DijkstraShortestPath<Integer, DefaultWeightedEdge> dsp= new DijkstraShortestPath<Integer, DefaultWeightedEdge>(g,source,dest);
		dsp.getPathLength();
		
		return dsp.getPathLength();
	
	}
	
	public SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> makeAvailableGraph (NetPlan netPlan, double cap_req)
	{
		
		SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k= new SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);

		
		for (int l=0; l<netPlan.getLinks(netPlan.getNetworkLayerDefault()).size();l++)
		{
			//System.out.println("Link: "+ netPlan.getLink(l).getIndex());
			int sourcenode=(int) netPlan.getLink(l).getOriginNode().getIndex();
			//Node sourcenode= netPlan.getLink(l).getOriginNode();
			g_k.addVertex(sourcenode);
			int destnode=(int) netPlan.getLink(l).getDestinationNode().getIndex();
			//Node destnode= netPlan.getLink(l).getDestinationNode();
			g_k.addVertex(destnode);	
			if (checkLinkCapacity(netPlan.getLink(l),cap_req))
			{
				//System.out.println("Source: "+sourcenode+ " Destination: "+ destnode);
				g_k.addEdge(sourcenode,destnode);
				//g_k.addEdge(destnode,sourcenode);
				g_k.setEdgeWeight(g_k.getEdge(sourcenode, destnode), 1.0);
			}
		}

		
		return g_k;
	}
	
	public boolean checkConflictOfLightpath(PLightpath pl)
	{
		
		for (int i=0; i<pl.size();i++)
			for (int l=0; l<pl.get(i).size();l++)
			{
				if (linkcap.get(pl.get(i).get(l).getFirst().getIndex(), pl.get(i).get(l).getSecond())==0.0 )//DONE //it means that there is a conflict in lightpath reservation for the current request
				{	
					return true;	
				}
				
			}
		
		return false;
		
	}
	
	public void reservePLightpathResources(PLightpath pl)
	{
		//DEBUG
		IntArrayList ind3= new IntArrayList();
		IntArrayList ind4= new IntArrayList();
		DoubleArrayList val2= new DoubleArrayList();
		
		for (int i=0; i<pl.size();i++)
			for (int l=0; l<pl.get(i).size();l++)
			{
				//linkcap.getNonZeros(ind3, ind4, val2);
				//int before = ind3.size();
				int n=pl.get(i).get(l).getFirst().getOriginNode().getIndex();
				int m= pl.get(i).get(l).getFirst().getDestinationNode().getIndex();
				//System.out.println("Going to set to 0.0 the value in ( ("+n+","+m+"),"+pl.get(i).get(l).getSecond()+") which actually is "+linkcap.get(pl.get(i).get(l).getFirst().getIndex(), pl.get(i).get(l).getSecond()));
				//linkcap.getNonZeros(ind3, ind4, val2);
				//int after = ind3.size();
				//if (after == before )
				if (linkcap.get(pl.get(i).get(l).getFirst().getIndex(), pl.get(i).get(l).getSecond())==0.0 )//DONE //it means that there is a conflict in lightpath reservation for the current request
				{	
					
					//we have to 
					new IllegalArgumentException("Wavelengths conflict");
					
					
				}
				else
				{
					linkcap.set(pl.get(i).get(l).getFirst().getIndex(), pl.get(i).get(l).getSecond(),0.0);
				}
				//	throw new IllegalArgumentException("Wavelengths mysmatching");
				
			}
		
	}
	
	
	public PLightpath feasibleExistingPligthpath(int source, int dest, double cap_req,  NetPlan netPlan)
	{
		ArrayList<Triple<Pair<Integer, Integer>, PLightpath, Double>>  candidates= new ArrayList<Triple<Pair<Integer,Integer>,PLightpath,Double>>();
		
		for (int i=0; i< existing_plightpaths.size();i++)
		{
			if( existing_plightpaths.get(i).getFirst().getFirst()==source && existing_plightpaths.get(i).getFirst().getSecond()==dest)
			{
				//if (existing_plightpaths.get(i).getSecond().get(0).get(0).getThird()>cap_req)
				if (existing_plightpaths.get(i).getThird()>cap_req)
					candidates.add(existing_plightpaths.get(i));
			}
		}
		
		candidates.sort(new Comparator<Triple<Pair<Integer, Integer>, PLightpath, Double>>() {
			@Override
			public int compare(Triple<Pair<Integer, Integer>, PLightpath, Double> o1, Triple<Pair<Integer, Integer>, PLightpath, Double> o2) {
				if (o1.getThird() > o2.getThird()) {
					return 1;
				} else if (o1.getFirst().equals(o2.getFirst())) {
					return 0; // You can change this to make it then look at the
					//words alphabetical order
				} else {
					return -1;
				}
			}
		});
		
		if (candidates.size()>0)
		return candidates.get(0).getSecond();
		else
			return null;
	}
	
	
	//Plightpath is a list: 0 working, 1 backup
	public PLightpath feasiblePlightpath(int source, int dest, double cap_req , SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k,  NetPlan netPlan)
	{
		PLightpath feas=new Testable_reliableSlicingPAL_Heuristic.PLightpath();
		SuurballeKDisjointShortestPaths<Integer,DefaultWeightedEdge> sur= new SuurballeKDisjointShortestPaths<Integer,DefaultWeightedEdge>(g_k);
		//	SuurballeKDisjointShortestPaths<Integer,DefaultWeightedEdge> sur= new SuurballeKDisjointShortestPaths<Integer,DefaultWeightedEdge>(g_k_debug);

		//System.out.println("Source:"+source+" Dest:"+dest);
		
		List<GraphPath<Integer,DefaultWeightedEdge>> ksp=sur.getPaths(source ,dest , 999);
		//take the first two paths with enough capacity
		//System.out.println("Ksurballe: "+ksp.size());
		//System.out.println("Paths: " +ksp.toString());
		
//		if ((source ==27 && dest ==34) || (source==34 && dest==38) )
//		{
//			
//			System.out.println("--------PLIGHTPATH-----");
//			Double weight= Double.valueOf(ksp.get(0).getEdgeList().size()+ksp.get(1).getEdgeList().size());
//			System.out.println("("+source+","+dest+") weights: "+weight.toString());
//			System.out.println("Working:");
//			System.out.println(ksp.get(0).getVertexList());
//			System.out.println(ksp.get(0).getWeight());
//			System.out.println("Backup:");
//			System.out.println(ksp.get(1).getVertexList());
//			System.out.println(ksp.get(1).getWeight());
//		}
		
		for (int p=0; p< ksp.size() && feas.size()<2;p++ )
		{
			
			if (checkEnoughRes(ksp.get(p),cap_req,g_k,netPlan))//make version wp
			{
				Lightpath l;
				if (OP.getString().contentEquals("VWP"))
				l= lightpathFromPathVWP(ksp.get(p),cap_req,g_k,netPlan);
				else
					l= lightpathFromPathWP(ksp.get(p),cap_req,g_k,netPlan);
				
				if (l!= null)
				{
				feas.add(l);
				}
			}
		}
		
		if (feas.size() <2)
		{
			return null;
		}
		else
		{
		return feas;
		}
	}
	
	public Lightpath lightpathFromPathWP(GraphPath<Integer,DefaultWeightedEdge> path, double cap_req , SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k,  NetPlan netPlan)
	{
		
		Lightpath l= new Lightpath();
		List<DefaultWeightedEdge> links= path.getEdgeList();

		boolean foundwl=false;
		int w;
		for (w=0; w< W.getInt() && !foundwl ;w++)
		{
			boolean wl_available_in_all_links=true;
			for (int i=0;i<links.size()&&wl_available_in_all_links ;i++)	
			{
				int source=g_k.getEdgeSource(links.get(i));
				int dest= g_k.getEdgeTarget(links.get(i));
				Link link= netPlan.getNodePairLinks(netPlan.getNode(source), netPlan.getNode(dest), false).first();
				if(linkcap.get(link.getIndex(), w) <cap_req)
				{
					wl_available_in_all_links=false;
				}
			}
			
			if (wl_available_in_all_links)
			{
				foundwl=true;
			}
			
		}
			
		w=w-1; //because for cycle increaeses before to exit
		if (!foundwl)
		{
			return null;
		}
		
		//now we know the w and the links, we add the to l
		for (int i=0;i<links.size() ;i++)	
		{
			int source=g_k.getEdgeSource(links.get(i));
			int dest= g_k.getEdgeTarget(links.get(i));
			Link link= netPlan.getNodePairLinks(netPlan.getNode(source), netPlan.getNode(dest), false).first();
			l.add(new Triple<Link, Integer, Double>(link, w, C.getDouble(), true));
		}
		
		return l;
		
	}
	
	public Lightpath lightpathFromPathVWP(GraphPath<Integer,DefaultWeightedEdge> path, double cap_req , SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k,  NetPlan netPlan)
	{
		
		Lightpath l= new Lightpath();
		List<DefaultWeightedEdge> links= path.getEdgeList();

		for (int i=0;i<links.size();i++)
		{
			int source=g_k.getEdgeSource(links.get(i));
			int dest= g_k.getEdgeTarget(links.get(i));
			Link link= netPlan.getNodePairLinks(netPlan.getNode(source), netPlan.getNode(dest), false).first();
			boolean wl_in_link=false;
			for (int w=0; w< W.getInt() && !wl_in_link;w++)
			{
				if(linkcap.get(link.getIndex(), w) >=cap_req)
				{
					
					int n=link.getOriginNode().getIndex();
					int m=link.getDestinationNode().getIndex();
					//System.out.println("("+String.valueOf(n)+","+String.valueOf(m)+"):"+String.valueOf(j));
					
					//System.out.println("Going to reserve wavelength "+w+" on link ("+n+","+m+")");
					
					
					l.add(new Triple<Link, Integer, Double>(link, w, C.getDouble(), true));
					wl_in_link=true;
				}
			
			}
			if (!wl_in_link)
			{
				return null;
			}
			
		}
		
		return l;
		
	}
	
	public boolean checkEnoughRes(GraphPath<Integer,DefaultWeightedEdge> path, double cap_req , SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k,  NetPlan netPlan)
	{
		
		List<DefaultWeightedEdge> links= path.getEdgeList();

		for (int i=0;i<links.size();i++)
		{
			int source=g_k.getEdgeSource(links.get(i));
			int dest= g_k.getEdgeTarget(links.get(i));
			Link l= netPlan.getNodePairLinks(netPlan.getNode(source), netPlan.getNode(dest), false).first();
			boolean wl_in_link=false;
			for (int w=0; w< W.getInt() && !wl_in_link;w++)
			{
				if(linkcap.get(l.getIndex(), w) >=cap_req)
				{
					wl_in_link=true;
				}
			
			}
			if (!wl_in_link)
				return false;
		}
		return true;
	}
	
	public boolean checkLinkCapacity(Link l, double cap_req)
	{
		
		
			for (int w=0; w< W.getInt();w++)
			{
				if(linkcap.get(l.getIndex(), w) >=cap_req)
				{
					return true;
				}
			
			}

		return false;
	}
	
	
	public ArrayList<Integer> findAlreadyActive_sametype(int s, int u , ArrayList<Double> CapnOriginal,ArrayList<Double> capn, ArrayList<Integer> used)
	{	//returns active not used by current slice
		//TO DO: FIND ALREADY ACTIVE OF THE SAME TYPE
		ArrayList<Integer> candi=new ArrayList<Integer>();
		IntArrayList ind= new IntArrayList();
		int[] pos= new int[]{Ns,Nf,N};
		DoubleArrayList val= new DoubleArrayList();
		
		Y_sun.getNonZeros(ind, val);
		pos= new int[]{Ns,Nf,N};
		for (int i=0; i< ind.size();i++)
		{
			int s1= DoubleMatrixND.ind2sub(ind.get(i), pos).get(0);
			int u1=DoubleMatrixND.ind2sub(ind.get(i), pos).get(1);
			int n1= DoubleMatrixND.ind2sub(ind.get(i), pos).get(2);
			if (u1==u && !used.contains(n1) && isCapable(s, u, n1))
			{
				candi.add(n1);
			}
		}
			
		return candi;

	}
	
	public int findNextActiveInPath(int s, int function, ArrayList<Integer> path ,int curren_pos)
	{
		for (int i=curren_pos; i<path.size();i++)
		{
			if ( M_sun.get( new int[] {s,function,path.get(i)})==1 )
			{
				if(i>curren_pos+1)
					System.out.println("Jump may occurr, path size is "+path.size());
				
				return i;
			}
		}
		
		return curren_pos+1;
	}
	
	public int findSurNearestCore(int node,NetPlan netPlan)
	{
		
		
		//TODO add check for slice compatibility with selected core node, need slice as a parameter of the method
		SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k= new SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);

		
		for (int l=0; l<netPlan.getLinks(netPlan.getNetworkLayerDefault()).size();l++)
		{
			//System.out.println("Link: "+ netPlan.getLink(l).getIndex());
			int sourcenode=(int) netPlan.getLink(l).getOriginNode().getIndex();
			//Node sourcenode= netPlan.getLink(l).getOriginNode();
			g_k.addVertex(sourcenode);
			int destnode=(int) netPlan.getLink(l).getDestinationNode().getIndex();
			//Node destnode= netPlan.getLink(l).getDestinationNode();
			g_k.addVertex(destnode);	
//			if (checkLinkCapacity(netPlan.getLink(l),cap_req))
//			{
//				//System.out.println("Source: "+sourcenode+ " Destination: "+ destnode);
				g_k.addEdge(sourcenode,destnode);
//				//g_k.addEdge(destnode,sourcenode);
//				g_k.setEdgeWeight(g_k.getEdge(sourcenode, destnode), 1.0);
//			}
		}
		
		
		Iterator it = core_capable.iterator();
		
		SuurballeKDisjointShortestPaths<Integer,DefaultWeightedEdge> sur= new SuurballeKDisjointShortestPaths<Integer,DefaultWeightedEdge>(g_k);
		List<GraphPath<Integer,DefaultWeightedEdge>> ksp;
		
		int pcore=(int) it.next();
		double dist=0;
		double dist_min=0;
		int coremin=0;
		if ((pcore==node)&&virt.getString().equals("SFC"))
			return pcore;
		if ((pcore==node)&&virt.getString().equals("VN")&& it.hasNext())
		{
			pcore=(int) it.next();//skip
		}
		
			ksp=sur.getPaths(pcore ,node, 999); 
			dist= ksp.get(0).getEdgeList().size()+ksp.get(1).getEdgeList().size();
			dist_min=dist;
			//System.out.println("Distance from node "+node+" to core "+pcore+" is "+ dist);
			coremin=pcore;	
		
		
		while (it.hasNext())
		{ 
			
				pcore=(int)it.next();
				if ((pcore==node)&&virt.equals("SFC"))
					return pcore;
				if ((pcore!=node))
				{
					
					ksp=sur.getPaths(pcore ,node, 999); 
					dist= ksp.get(0).getEdgeList().size()+ksp.get(1).getEdgeList().size();
					//System.out.println("Distance from node "+node+" to core "+pcore+" is "+ dist);
					if (dist<dist_min )
					{	
						//System.out.print(" and it is lower, new nearest core is "+ pcore + " with distance "+dist);
						dist_min=dist;
						coremin = pcore;
					}
				}
				
				
				
		}
			
		return coremin;
	}
	
	public ArrayList<Integer> findAlreadyActive_generic(int s, int u , ArrayList<Double> CapnOriginal,ArrayList<Double> capn, ArrayList<Integer> used)
	{	
		//returns active not used by current slice
		//TO DO: FIND ALREADY ACTIVE OF THE SAME TYPE
		ArrayList<Integer> candi=new ArrayList<Integer>();
		for (int i=0;i<N;i++)
		{
			if (CapnOriginal.get(i)-Capn.get(i) >0 && isCapable(s, u, i) && !used.contains(i)) {
				candi.add(i);
			}
		}
		return candi;

	}
	
	
	public double computeCostVNMapping(List<Link> pathlinks, ArrayList<Double> CapnOriginal,ArrayList<Double> capn ) {
		double cost=0;
		
		cost= pathlinks.size();
		double vncost=0;
		ArrayList <Integer> path = new ArrayList<Integer>();
		
		for (Link l: pathlinks )// NB HOW ARE YOU SURE THAT YOUR MAPPING ENDS IN K?
		{
			path.add(l.getOriginNode().getIndex());
			path.add(l.getDestinationNode().getIndex());
		}

		Object[] aux=Arrays.stream(path.toArray()).distinct().toArray();
		Integer[] in = Arrays.copyOf(aux, aux.length, Integer[].class);
		path.clear();
		path.addAll(Arrays.asList(in));
		
		//here we have distinct list of node
		int actn=0;
		
		for (int i=0;i<path.size();i++)
		{
			if (CapnOriginal.get(path.get(i))-Capn.get(path.get(i)) >0) {
				actn++;
			}
		}
		double cost_with_act= cost -(((double) (actn)/1000));
		return  cost_with_act;
		//return -actn*1000+cost;
	}

	public double computeCostVNMapping_sur(GraphPath<Integer, DefaultWeightedEdge> pathlinks,SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> g_k_sur, ArrayList<Double> CapnOriginal,ArrayList<Double> capn ) {
		double cost=0;
		
		cost= 0;
		double vncost=0;
		List <Integer> path = new ArrayList<Integer>();
		
		for (DefaultWeightedEdge l: pathlinks.getEdgeList() )// NB HOW ARE YOU SURE THAT YOUR MAPPING ENDS IN K?
		{
			cost=cost+g_k_sur.getEdgeWeight(l);
//			path.add(l.getOriginNode().getIndex());
//			path.add(l.getDestinationNode().getIndex());
		}

		path=pathlinks.getVertexList();
		Object[] aux=Arrays.stream(path.toArray()).distinct().toArray();
		Integer[] in = Arrays.copyOf(aux, aux.length, Integer[].class);
		path.clear();
		path.addAll(Arrays.asList(in));
		
		//here we have distinct list of node
		int actn=0;
		
		for (int i=0;i<path.size();i++)
		{
			if (CapnOriginal.get(path.get(i))-Capn.get(path.get(i)) >0 || antennas.contains(path.get(i))) {
				actn++;
			}
		}
		double cost_with_act= cost -(((double) (actn)/1000));
		return  cost_with_act;
		//return -actn*1000+cost;
	}
	
	public double computeCostSFCMapping(List<Link> pathlinks, ArrayList<Double> CapnOriginal,ArrayList<Double> capn ) {
		double cost=0;
		
		cost= pathlinks.size();
		double vncost=0;
		ArrayList <Integer> path = new ArrayList<Integer>();
		
		for (Link l: pathlinks )// NB HOW ARE YOU SURE THAT YOUR MAPPING ENDS IN K?
		{
			path.add(l.getOriginNode().getIndex());
			path.add(l.getDestinationNode().getIndex());
		}

		Object[] aux=Arrays.stream(path.toArray()).distinct().toArray();
		Integer[] in = Arrays.copyOf(aux, aux.length, Integer[].class);
		path.clear();
		path.addAll(Arrays.asList(in));
		
		//here we have distinct list of node
		int actn=0;
		
		for (int i=0;i<path.size();i++)
		{
			if (CapnOriginal.get(path.get(i))-Capn.get(path.get(i)) >0) {
				actn++;
			}
		}
		double cost_with_act= cost -(((double) (actn)/1000));
		return  cost_with_act;
		//return -actn*1000+cost;
	}

	public String getDescription() {
		return "Algorithm for reliable slicing, Protection at Lightpath";
	}


	public List<Triple<String, String, String>> getParameters() {
		/*
		 * Returns the parameter information for all the InputParameter objects defined
		 * in this object (uses Java reflection)
		 */
		return InputParameter.getInformationAllInputParameterFieldsOfObject(this);
	}
	
	
	


}