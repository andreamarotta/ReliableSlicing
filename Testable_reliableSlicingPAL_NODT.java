package com.net2plan.examples.general.offline.nfv;

import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
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
import com.net2plan.examples.general.offline.Offline_ipOverWdm_routingSpectrumAndModulationAssignmentILPNotGrooming;
import com.net2plan.interfaces.networkDesign.*;
import com.net2plan.utils.Constants.RoutingType;
import com.net2plan.utils.InputParameter;
import com.net2plan.utils.Pair;
import com.net2plan.utils.Triple;
import org.apache.commons.lang3.StringUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.*;
import java.util.stream.Collectors;
import java.sql.Timestamp;
import java.text.SimpleDateFormat;

// ASSUMES YOU CAN NOT SPLIT ONE TRAFFIC REQUEST OVER DIFFERENT WAVELENGHTS
// DOES IT MAKES SENSE?

/**
 *
 * 
 * @net2plan.keywords JOM, NFV
 * @net2plan.inputParameters
 * @author Andrea Marotta
 */
public class Testable_reliableSlicingPAL_NODT  {

	private InputParameter C = new InputParameter("C", (double) 40, "Wavelength linerate");

	private InputParameter W = new InputParameter("W", (int) 40, "Wavelength per link");

	private InputParameter isolation = new InputParameter("isolation", (int) 0, "Slice isolation level");

	//private InputParameter SI = new InputParameter("SI", (int) 2, "Number of slices per common segment");

	private InputParameter S = new InputParameter("S", (int) 2, "Number of slice requests to be generated");

	//private InputParameter Cn = new InputParameter("Cn", (int) 4, "Node capacity");

	private InputParameter K = new InputParameter("K", (int) 1, "Cost of function placement");

	private InputParameter solverName = new InputParameter("solverName", "#select# glpk ipopt xpress cplex",
			"The solver name to be used by JOM. GLPK and IPOPT are free, XPRESS and CPLEX commercial. GLPK, XPRESS and CPLEX solve linear problems w/w.o integer contraints. IPOPT is can solve nonlinear problems (if convex, returns global optimum), but cannot handle integer constraints");
	private InputParameter solverLibraryName = new InputParameter("solverLibraryName", "",
			"The solver library full or relative path, to be used by JOM. Leave blank to use JOM default.");
	private InputParameter maxSolverTimeInSeconds = new InputParameter("maxSolverTimeInSeconds", (double) -1,
			"Maximum time granted to the solver to solve the problem. If this time expires, the solver returns the best solution found so far (if a feasible solution is found)");
	private NetworkLayer wdmLayer, ipLayer;

	private InputParameter virt= new InputParameter("Virtualization", "SF", "#SF #VN");
	private InputParameter bandwidth= new InputParameter("Bandwidth", (int) 50, "System Bandwidth");
	
	private InputParameter mean= new InputParameter("Mean", (double) 2, "Average slice capacity");
	private InputParameter std= new InputParameter("Std", (double) 2, "Std slice capacity");
	private InputParameter seed= new InputParameter("seed", (int) 1, "Seed of simulation");

	private DoubleMatrixND M_sun;// = new DoubleMatrixND(new int[] { SI.getInt(), Nf, N });
	private DoubleMatrixND Omega_SUsu;
	private DoubleMatrixND b_suv;// = new DoubleMatrixND(new int[] { SI.getInt(), Nf, Nf });
	private DoubleMatrixND Y_sun;//
	//private ArrayList<Double> Capn= new ArrayList<Double>();;
	private ArrayList<Integer> D;//= new ArrayList<Integer>();
	private ArrayList<Integer> R;//= new ArrayList<Integer>();

	private ArrayList <Pair<Integer,Integer>> slice_cseg;

	private Map<Integer, Pair<Integer,Double>> slice_Requests; 
	private DoubleMatrix2D K_su;//
	private int N;
	private int Nf;
	private int Ns,Nc;


	private boolean DEBUG = false;

	private ArrayList<Double> Capn= new ArrayList<Double>();

	private ArrayList<Double> Capn_native= new ArrayList<Double>();

	private ArrayList<Integer> antennas= new ArrayList<Integer>();

	private Set<Integer> activeNodes= new HashSet<Integer>();

	private InputParameter grooming= new InputParameter("grooming", (boolean) true, "Grooming capable network");


	public static void main(String[] args) {
		//NetPlan netPlan = NetPlan.loadFromFile(new File("resil.n2p"));
		//NetPlan netPlan = NetPlan.loadFromFile(new File("eon_N18_E66_withTraffic.n2p"));
		NetPlan netPlan = NetPlan.loadFromFile(new File("N6_3.n2p"));



		Map<String, String> parameters = new HashMap<>();
		parameters.put("solverName", "cplex");
		parameters.put("solverLibraryName","C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio129\\opl\\bin\\x64_win64\\cplex1290.dll");
		parameters.put("maxSolverTimeInSeconds", String.valueOf((double) -1));
		//parameters.put("solverName", "glpk" );
		//parameters.put("solverLibraryName","C:\\WINDOWS\\system32\\glpk_4_48.dll");

		parameters.put("C", String.valueOf((double) 40));
		parameters.put("W", String.valueOf((int) 40));
		parameters.put("isolation", String.valueOf((int) 3));
		//parameters.put("SI", String.valueOf((int) 3));
		parameters.put("S", String.valueOf((int) 3));
		//parameters.put("Cn", String.valueOf((int) 255));
		parameters.put("K", String.valueOf((int) 1));
		parameters.put("Virtualization","VN");
		parameters.put("Bandwidth","100");
		parameters.put("Mean","2");
		parameters.put("Std","1");
		parameters.put("seed", "1");



		Testable_reliableSlicingPAL_NODT Algorithm = new Testable_reliableSlicingPAL_NODT();
		Algorithm.getParameters();
		System.out.println(Algorithm.executeAlgorithm(netPlan, parameters, null));
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
		DoubleMatrix2D Fmn= DoubleFactory2D.sparse.make(N,N);
		for (Link e : netPlan.getLinks())
		{
			Fmn.set(e.getOriginNode().getIndex(), e.getDestinationNode().getIndex(), 1);
			Fmn.set( e.getDestinationNode().getIndex(),e.getOriginNode().getIndex(), 1);
		}
		final int Nf=4;
		this.N=N;
		this.Nf=Nf;
		IntArrayList ind= new IntArrayList();
		IntArrayList ind2= new IntArrayList();
		DoubleArrayList val= new DoubleArrayList();
		int[] pos= new int[]{Ns,Nf,Nf};
		SimpleDateFormat sdf = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss");

		
		

		loadInput(netPlan, mean.getDouble(), std.getDouble(), S.getInt(), seed.getInt());
		createInput(netPlan, S.getInt(), bandwidth.getInt(), isolation.getInt(),grooming.getBoolean());



		//================================================================================
		// OPTIMIZATION PROBLEM
		//================================================================================

		/* Create the optimization problem object (JOM library) */
		OptimizationProblem op = new OptimizationProblem();

		/* Set some input parameters to the problem */
		op.setInputParameter("b_suv", b_suv);
		op.setInputParameter("M_sun", M_sun);
		op.setInputParameter("K_su", K_su);
		op.setInputParameter("N", N);
		op.setInputParameter("Ns", Ns);
		op.setInputParameter("W", W.getInt());
		op.setInputParameter("C", C.getDouble());
		op.setInputParameter("Cn", Capn, "row" );
		op.setInputParameter("M", 9999);
		op.setInputParameter("F_mn", Fmn);
		op.setInputParameter("Ns",Ns);
		op.setInputParameter("N",N);
		op.setInputParameter("Nf",Nf);

		op.addDecisionVariable("Z_sijmn", true, new int[] { Ns , N , N, N, N }, 0, 1); 
		op.addDecisionVariable("Z_ijmn", true, new int[] {  N , N, N, N }); 
		op.addDecisionVariable("Zc_ijmn", true, new int[] {  N , N, N, N },0,1); 
		op.addDecisionVariable("Z_lijmn", true, new int[] { W.getInt() , N , N, N, N }, 0, 1); 
		op.addDecisionVariable("Z_slijmn", true, new int[] { Ns , W.getInt() , N , N, N, N }, 0, 1); 
		op.addDecisionVariable("Zeta_sijmn", true, new int[] { Ns , N , N, N, N }, 0, 1); 
		op.addDecisionVariable("Zeta_ijmn", true, new int[] {  N , N, N, N });
		op.addDecisionVariable("Zetac_ijmn", true, new int[] {  N , N, N, N },0,1);
		op.addDecisionVariable("Zeta_lijmn", true, new int[] { W.getInt() , N , N, N, N }, 0, 1); 
		op.addDecisionVariable("Zeta_slijmn", true, new int[] { Ns , W.getInt() , N , N, N, N }, 0, 1);	
		op.addDecisionVariable("X_lij", true, new int[] { W.getInt() , N, N }, 0, 1); 
		op.addDecisionVariable("X_ij", true, new int[] { N, N }); 
		op.addDecisionVariable("X_slij", true, new int[] {Ns, W.getInt() , N, N }, 0, 1);
		op.addDecisionVariable("X_sij", true, new int[] {Ns , N, N }, 0, 1);
		op.addDecisionVariable("Xi_ij", true, new int[] { N, N }); 
		op.addDecisionVariable("Xi_lij", true, new int[] { W.getInt() , N, N }, 0, 1); 
		op.addDecisionVariable("Xi_slij", true, new int[] {Ns, W.getInt() , N, N }, 0, 1);
		op.addDecisionVariable("Xi_sij", true, new int[] {Ns , N, N }, 0, 1);
		op.addDecisionVariable("p_sbeij", false, new int[] { Ns , N , N, N, N }); 
		op.addDecisionVariable("pi_sbeij", false, new int[] { Ns , N , N, N, N }); 
		op.addDecisionVariable("pc_sbeij", true, new int[] { Ns , N , N, N, N },0,1); 
		op.addDecisionVariable("pic_sbeij", true, new int[] { Ns , N , N, N, N },0,1); 
		op.addDecisionVariable("l_sbe", false, new int[] { Ns , N , N }); 
		op.addDecisionVariable("h_suvbe", true, new int[] { Ns , Nf, Nf, N , N },0,1);
		op.addDecisionVariable("l_suvbe", false, new int[] { Ns , Nf, Nf, N , N });//auxiliary
		op.addDecisionVariable("y_sun", true, new int[] {Ns, Nf,N },0,1); 
		op.addDecisionVariable("K_sun", false, new int[] {Ns, Nf,N}); //auxiliary
		op.addDecisionVariable("t_sbelij", false, new int[] { Ns , N , N, W.getInt() ,N, N }); 
		op.addDecisionVariable("t_sbeij", false, new int[] { Ns , N , N ,N, N });
		op.addDecisionVariable("ti_sbeij", false, new int[] { Ns , N , N ,N, N });
		op.addDecisionVariable("y_n", true, new int[] {1,N},0,1);


		//op.addConstraint("y_n>=y_sun");//WRONG
		op.addConstraint("l_sbe>=0");
		op.addConstraint("pi_sbeij>=0");
		op.addConstraint("p_sbeij>=0");
		op.addConstraint("p_sbeij<=M");
		op.addConstraint("Zeta_ijmn>=0");
		op.addConstraint("Z_ijmn>=0");
		op.addConstraint("t_sbeij>=0");
		op.addConstraint("ti_sbeij>=0");

		//binarization

		//op.addConstraint("Zc_ijmn>=Z_ijmn./(ones([N N N N]).*M)");
		op.addConstraint("Zc_ijmn>=Z_ijmn./M");
		op.addConstraint("Zc_ijmn<=Z_ijmn.*M");
		op.addConstraint("Zetac_ijmn>=Zeta_ijmn./M");
		op.addConstraint("Zetac_ijmn<=Zeta_ijmn.*M");
		//op.addConstraint("pc_sbeij>=p_sbeij/M");
		//op.addConstraint("pc_sbeij<=p_sbeij*M");
		//op.addConstraint("pic_sbeij>=pi_sbeij/M");
		//op.addConstraint("pic_sbeij<=pi_sbeij*M");

		String message;





		// C30 //caso i==j da escludere?
		System.out.println("");
		Timestamp timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C30");
		for (int s=0 ; s< Ns; s++)
			for (int i =0; i<N; i++)
				for (int b=0; b<N; b++)
					for (int e=0; e<N; e++)
					{
						if (i==b)
						{
							message= String.format("sum(p_sbeij(%1$d,%2$d,%3$d,%4$d,all))-sum(p_sbeij(%1$d,%2$d,%3$d,all,%4$d)) == l_sbe(%1$d,%2$d,%3$d)",s,b,e,i );
							if (DEBUG) System.out.println(message);
							op.addConstraint(message);	
						}
						else 
						{
							if (i==e)
							{
								message= String.format("sum(p_sbeij(%1$d,%2$d,%3$d,%4$d,all))-sum(p_sbeij(%1$d,%2$d,%3$d,all,%4$d)) == -l_sbe(%1$d,%2$d,%3$d)",s,b,e,i );
								if (DEBUG) System.out.println(message);
								op.addConstraint(message);
							}
							else
							{
								message= String.format("sum(p_sbeij(%1$d,%2$d,%3$d,%4$d,all))-sum(p_sbeij(%1$d,%2$d,%3$d,all,%4$d)) == 0",s,b,e,i );
								if (DEBUG) System.out.println(message);
								op.addConstraint(message);
							}
						}

					}


		// C31
		System.out.println("");
		timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C31");
		for (int s=0 ; s< Ns; s++)
			for (int b=0; b<N; b++)
				for (int e=0; e<N; e++)
					for (int u=0; u<Nf; u++)
						for (int v=0; v<Nf; v++)
						{
							message= String.format("l_suvbe(%1$d,%2$d,%3$d,%4$d,%5$d)== sum(b_suv(%1$d,%2$d,%3$d)) * sum(h_suvbe(%1$d,%2$d,%3$d,%4$d,%5$d))",s,u,v,b,e);
							if (DEBUG) System.out.println(message);
							op.addConstraint(message);
						}

		//C31
		System.out.println("");
		timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C31bis");
		for (int s=0 ; s< Ns; s++)
			for (int b=0; b<N; b++)
				for (int e=0; e<N; e++)
				{
					message= String.format("l_sbe(%1$d,%2$d,%3$d)== sum(l_suvbe(%1$d,all,all,%2$d,%3$d)) ",s,b,e);
					if (DEBUG) System.out.println(message);
					op.addConstraint(message);
				}

		//C32
		System.out.println("");
		timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C32");
		for (int s=0 ; s< Ns; s++)
			for (int i =0; i<N; i++)
				for (int j =0; j<N; j++)
				{	
					message= String.format("sum(p_sbeij(%1$d,all,all,%2$d,%3$d))<= C*sum(X_sij(%1$d,%2$d,%3$d))",s,i,j);
					if (DEBUG) System.out.println(message);
					op.addConstraint(message);
				}	


		//C33
		System.out.println("");
		timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C33");
		for (int s=0 ; s< Ns; s++)
			for (int b=0; b<N; b++)
				for (int e=0; e<N; e++)		
					for (int i =0; i<N; i++)
						for (int j =0; j<N; j++)
						{
							message= String.format("sum(t_sbeij(%1$d,%2$d,%3$d,%4$d,%5$d))<= M .* sum(X_sij(%1$d,%4$d,%5$d)) ",s,b,e,i,j);
							if (DEBUG) System.out.println(message);
							op.addConstraint(message);
							message= String.format("sum(t_sbeij(%1$d,%2$d,%3$d,%4$d,%5$d))<= sum(p_sbeij(%1$d,%2$d,%3$d,%4$d,%5$d)) ",s,b,e,i,j);
							if (DEBUG) System.out.println(message);
							op.addConstraint(message);
							message= String.format("sum(t_sbeij(%1$d,%2$d,%3$d,%4$d,%5$d))>= sum(p_sbeij(%1$d,%2$d,%3$d,%4$d,%5$d)) - (M*(1-sum(X_sij(%1$d,%4$d,%5$d)))) ",s,b,e,i,j);
							if (DEBUG) System.out.println(message);
							op.addConstraint(message);

							//extra xslij
							//message= String.format("sum(pc_sbeij(%1$d,%2$d,%3$d,%5$d,%6$d))>= sum(X_slij(%1$d,%4$d,%5$d,%6$d)) ",s,b,e,l,i,j);
							//if (DEBUG) System.out.println(message);
							//op.addConstraint(message);


						}

		//C33
		System.out.println("");
		timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C33bis");
		for (int i =0; i<N; i++)
			for (int j =0; j<N; j++)
			{	//1
				message= String.format("sum(t_sbeij(all,all,all,%1$d,%2$d))<= C*sum(X_ij(%1$d,%2$d))",i,j);
				if (DEBUG) System.out.println(message);
				op.addConstraint(message);
			}

		// C37 //caso i==j da escludere?
		System.out.println("");
		timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C37");
		for (int m =0; m<N; m++)
			for (int i=0; i<N; i++)
				for (int j=0; j<N; j++)
				{
					if (i==m)
					{
						message= String.format("sum(Z_ijmn(%1$d,%2$d,%3$d,all))-sum(Z_ijmn(%1$d,%2$d,all,%3$d)) == X_ij(%1$d,%2$d)",i,j,m );
						if (DEBUG) System.out.println(message);
						op.addConstraint(message);	
					}
					else 
					{
						if (m==j)
						{
							message= String.format("sum(Z_ijmn(%1$d,%2$d,%3$d,all))-sum(Z_ijmn(%1$d,%2$d,all,%3$d)) == -X_ij(%1$d,%2$d)",i,j,m );
							if (DEBUG) System.out.println(message);
							op.addConstraint(message);
						}
						else
						{
							message= String.format("sum(Z_ijmn(%1$d,%2$d,%3$d,all))-sum(Z_ijmn(%1$d,%2$d,all,%3$d)) == 0",i,j,m );
							if (DEBUG) System.out.println(message);
							op.addConstraint(message);
						}
					}

				}


		//C39		
		System.out.println("");
		timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C39");
		for (int s=0 ; s< Ns; s++)
		{
			for (int b=0; b<N; b++)
				for (int e=0; e<N; e++)
					for (int u=0; u<Nf; u++)
						for (int v=0; v<Nf; v++)
						{
							if (b!=e) 
							{
								message= String.format("h_suvbe(%1$d,%2$d,%3$d,%4$d,%5$d) >= sum(y_sun(%1$d,%2$d,%4$d)) +sum( y_sun(%1$d,%3$d,%5$d)) -1",s,u,v,b,e);
								if (DEBUG) System.out.println(message);
								op.addConstraint(message);
								message= String.format("h_suvbe(%1$d,%2$d,%3$d,%4$d,%5$d) <= sum(y_sun(%1$d,%2$d,%4$d))",s,u,v,b,e);
								if (DEBUG) System.out.println(message);
								op.addConstraint(message);
								message= String.format("h_suvbe(%1$d,%2$d,%3$d,%4$d,%5$d) <= sum(y_sun(%1$d,%3$d,%5$d))",s,u,v,b,e);
								if (DEBUG) System.out.println(message);
								op.addConstraint(message);
							}
						}
		}

		//C40
		System.out.println("");
		timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C40");
		boolean apply;
		for (int s=0 ; s< Ns; s++)
			for (int u=0; u<Nf; u++)
			{
				apply=false;
				for (int v=0; v<Nf;v++)
				{
					if (b_suv.get(new int[]{s,u,v})>0 ||b_suv.get(new int[]{s,v,u})>0 )
						apply=true;
				}
				if(apply)
				{
					message= String.format("sum(y_sun(%1$d,%2$d,all))==1 ",s,u);
					if (DEBUG) System.out.println(message);
					op.addConstraint(message);
				}
			}


		//CVN
		if (virt.getString().contentEquals("VN"))
		{
			System.out.println("");
			timestamp = new Timestamp(System.currentTimeMillis());
			System.out.println(sdf.format(timestamp));
			System.out.println("CVN");
			for (int s=0 ; s< Ns; s++)
				for (int n=0; n<N; n++)
				{
					message= String.format("sum(y_sun(%1$d,all,%2$d))<=1 ",s,n);
					if (DEBUG) System.out.println(message);
					op.addConstraint(message);
				}
		}

		//C41
		System.out.println("");
		timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C41");
		for (int s=0 ; s< Ns; s++)
			for (int u=0; u<Nf; u++)
				for (int n=0; n<N; n++)
				{
					message= String.format("sum(K_sun(%1$d,%2$d,%3$d))==sum(y_sun(%1$d,%2$d,%3$d))*sum(K_su(%1$d,%2$d))",s,u,n);
					if (DEBUG) System.out.println(message);
					op.addConstraint(message);
				}

		//C41
		System.out.println("");
		timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C41Bis");
		for (int n=0; n<N; n++)
		{
			message= String.format(" sum(K_sun(all,all,%1$d))<=Cn(%1$d)",n);
			if (DEBUG) System.out.println(message);
			op.addConstraint(message);
		}



		//C42
		System.out.println("");
		timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C42");
		message= String.format("y_sun<=M_sun");
		if (DEBUG) System.out.println(message);
		op.addConstraint(message);

		//C42Bis
		int s1,s,u1,u;
		if (isolation.getInt()<5 && isolation.getInt()>1)//only for 2,3,4
		{
		System.out.println("");
		timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C42bis");
		Omega_SUsu.getNonZeros(ind, val);
		pos= new int[]{Ns,Nf,Ns,Nf};
		
		for (int i=0; i< ind.size();i++)
		{	s1= (int) DoubleMatrixND.ind2sub(ind.get(i), pos).get(0);
		u1= (int) DoubleMatrixND.ind2sub(ind.get(i), pos).get(1);
		s= (int) DoubleMatrixND.ind2sub(ind.get(i), pos).get(2);
		u= (int) DoubleMatrixND.ind2sub(ind.get(i), pos).get(3);
		for (int n=0;n<N;n++)
		{
			message= String.format("sum(y_sun(%1$d,%2$d,%3$d))==sum(y_sun(%4$d,%5$d,%3$d))",s1,u1,n,s,u);
			if (DEBUG) System.out.println(message);
			op.addConstraint(message);
		}
		}
		}

		//C42 tris active nodes
		System.out.println("");
		timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C42tris active nodes");
		for ( s=0 ; s< Ns; s++)
			for ( u=0; u<Nf; u++)
				for (int n=0; n<N; n++)
				{
					message= String.format("sum(y_sun(%1$d,%2$d,%3$d))<=sum(y_n(%3$d))",s,u,n);
					if (DEBUG) System.out.println(message);
					op.addConstraint(message);
				}



		//C52 // si puó accorciare nm
		System.out.println("");
		timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C52");
		for (int n=0; n<N;n++)
			for (int m =0; m<N; m++)
			{
				message= String.format("sum(Z_ijmn(all,all,%1$d,%2$d))+sum(Zeta_ijmn(all,all,%1$d,%2$d))<=W*sum(F_mn(%1$d,%2$d))",m,n);
				if (DEBUG) System.out.println(message);
				op.addConstraint(message);
			}



		//================================================================================
		// PROTECTION CONSTRAINTS
		//================================================================================

		//		System.out.println("");
		//		timestamp = new Timestamp(System.currentTimeMillis());
		//		System.out.println(sdf.format(timestamp));
		//		System.out.println("C30");
		//		for (int s=0 ; s< Ns; s++)
		//			if(R.get(s)==1)
		//			for (int i =0; i<N; i++)
		//				for (int b=0; b<N; b++)
		//					for (int e=0; e<N; e++)
		//					{
		//						if (i==b)
		//						{
		//							message= String.format("sum(pi_sbeij(%1$d,%2$d,%3$d,%4$d,all))-sum(pi_sbeij(%1$d,%2$d,%3$d,all,%4$d)) == l_sbe(%1$d,%2$d,%3$d)",s,b,e,i );
		//							if (DEBUG) System.out.println(message);
		//							op.addConstraint(message);	
		//						}
		//						else 
		//						{
		//							if (i==e)
		//							{
		//								message= String.format("sum(pi_sbeij(%1$d,%2$d,%3$d,%4$d,all))-sum(pi_sbeij(%1$d,%2$d,%3$d,all,%4$d)) == -l_sbe(%1$d,%2$d,%3$d)",s,b,e,i );
		//								if (DEBUG) System.out.println(message);
		//								op.addConstraint(message);
		//							}
		//							else
		//							{
		//								message= String.format("sum(pi_sbeij(%1$d,%2$d,%3$d,%4$d,all))-sum(pi_sbeij(%1$d,%2$d,%3$d,all,%4$d)) == 0",s,b,e,i );
		//								if (DEBUG) System.out.println(message);
		//								op.addConstraint(message);
		//							}
		//						}
		//						
		//					}
		//					


		//					
		//
		////C32
		//		System.out.println("");
		//		timestamp = new Timestamp(System.currentTimeMillis());
		//		System.out.println(sdf.format(timestamp));
		//		System.out.println("C32");
		//		for (int s=0 ; s< Ns; s++)
		//			if(R.get(s)==1)
		//			for (int i =0; i<N; i++)
		//				for (int j =0; j<N; j++)
		//				{	
		//					message= String.format("sum(pi_sbeij(%1$d,all,all,%2$d,%3$d))<= C*sum(Xi_sij(%1$d,%2$d,%3$d))",s,i,j);
		//					if (DEBUG) System.out.println(message);
		//					op.addConstraint(message);
		//				}	
		//		
		//		
		//		//C33
		//		System.out.println("");
		//		timestamp = new Timestamp(System.currentTimeMillis());
		//		System.out.println(sdf.format(timestamp));
		//		System.out.println("C33");
		//		for (int s=0 ; s< Ns; s++)
		//			if(R.get(s)==1)
		//			for (int b=0; b<N; b++)
		//				for (int e=0; e<N; e++)		
		//					for (int i =0; i<N; i++)
		//						for (int j =0; j<N; j++)
		//							{
		//								message= String.format("sum(ti_sbeij(%1$d,%2$d,%3$d,%4$d,%5$d))<= M .* sum(Xi_sij(%1$d,%4$d,%5$d)) ",s,b,e,i,j);
		//								if (DEBUG) System.out.println(message);
		//								op.addConstraint(message);
		//								message= String.format("sum(ti_sbeij(%1$d,%2$d,%3$d,%4$d,%5$d))<= sum(pi_sbeij(%1$d,%2$d,%3$d,%4$d,%5$d)) ",s,b,e,i,j);
		//								if (DEBUG) System.out.println(message);
		//								op.addConstraint(message);
		//								message= String.format("sum(ti_sbeij(%1$d,%2$d,%3$d,%4$d,%5$d))>= sum(pi_sbeij(%1$d,%2$d,%3$d,%4$d,%5$d)) - (M*(1-sum(Xi_sij(%1$d,%4$d,%5$d)))) ",s,b,e,i,j);
		//								if (DEBUG) System.out.println(message);
		//								op.addConstraint(message);
		//								
		//								//extra xslij
		//								//message= String.format("sum(pc_sbeij(%1$d,%2$d,%3$d,%5$d,%6$d))>= sum(Xi_slij(%1$d,%4$d,%5$d,%6$d)) ",s,b,e,l,i,j);
		//								//if (DEBUG) System.out.println(message);
		//								//op.addConstraint(message);
		//								
		//								
		//							}
		//		
		//		//C33
		//		System.out.println("");
		//		timestamp = new Timestamp(System.currentTimeMillis());
		//		System.out.println(sdf.format(timestamp));
		//		System.out.println("C33bis");
		//		for (int i =0; i<N; i++)
		//			for (int j =0; j<N; j++)
		//				{	//1
		//					message= String.format("sum(ti_sbeij(all,all,all,%1$d,%2$d))<= C*sum(Xi_ij(%1$d,%2$d))",i,j);
		//					if (DEBUG) System.out.println(message);
		//					op.addConstraint(message);
		//				}

		// C37 //caso i==j da escludere?
		System.out.println("");
		timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C37");
		for (int m =0; m<N; m++)
			for (int i=0; i<N; i++)
				for (int j=0; j<N; j++)
				{
					if (i==m)
					{
						message= String.format("sum(Zeta_ijmn(%1$d,%2$d,%3$d,all))-sum(Zeta_ijmn(%1$d,%2$d,all,%3$d)) == X_ij(%1$d,%2$d)",i,j,m );
						if (DEBUG) System.out.println(message);
						op.addConstraint(message);	
					}
					else 
					{
						if (m==j)
						{
							message= String.format("sum(Zeta_ijmn(%1$d,%2$d,%3$d,all))-sum(Zeta_ijmn(%1$d,%2$d,all,%3$d)) == -X_ij(%1$d,%2$d)",i,j,m );
							if (DEBUG) System.out.println(message);
							op.addConstraint(message);
						}
						else
						{
							message= String.format("sum(Zeta_ijmn(%1$d,%2$d,%3$d,all))-sum(Zeta_ijmn(%1$d,%2$d,all,%3$d)) == 0",i,j,m );
							if (DEBUG) System.out.println(message);
							op.addConstraint(message);
						}
					}

				}



		//C51 si puó accorciare nm
		System.out.println("");
		timestamp = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("C51");
		for (int n=0; n<N;n++)
			for (int m =0; m<N; m++)
				for (int i=0; i<N; i++)
					for (int j=0; j<N; j++)
						//for (int l=0; l<W.getInt();l++)
					{
						message= String.format("sum(Zc_ijmn(%1$d,%2$d,%3$d,%4$d))+sum(Zc_ijmn(%1$d,%2$d,%4$d,%3$d))+sum(Zetac_ijmn(%1$d,%2$d,%3$d,%4$d))+sum(Zetac_ijmn(%1$d,%2$d,%4$d,%3$d))<=1",i,j,m,n);
						if (DEBUG) System.out.println(message);
						op.addConstraint(message);
					}



		//================================================================================
		// OBJECTIVE FUNCTION
		//================================================================================


		//op.setObjectiveFunction("minimize", "sum(y_sun)");
		//op.setObjectiveFunction("minimize", "sum(p_sbeij)");
		//op.setObjectiveFunction("minimize", "sum(X_ij)");
		//op.setObjectiveFunction("minimize", "sum(Z_ijmn)");
		op.setObjectiveFunction("minimize", "sum(Z_ijmn+Zeta_ijmn) + 0.001*sum(y_n)");

		//op.setObjectiveFunction("minimize", "100000*sum(Z_ijmn+Zeta_ijmn) + 10000*sum(Z_ijmn)+ 5000*sum(Zeta_ijmn) + 0.001*sum(y_n)");

		Timestamp tic = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		System.out.println("SOLVING");

		op.solve(solverName.getString(), "solverLibraryName", solverLibraryName.getString());
		//op.solve(solverName.getString(), "solverLibraryName", solverLibraryName.getString(),"msg_lev",3);

		Timestamp tac = new Timestamp(System.currentTimeMillis());
		System.out.println(sdf.format(timestamp));
		double diff= (tac.getTime()-tic.getTime())/1000;
		System.out.println(String.format("SOLVED in %1$f seconds",diff));


		//================================================================================
		// PRINT OUTPUT
		//================================================================================


		if (!op.solutionIsFeasible())
			throw new Net2PlanException("A feasible solution was not found");

		System.out.println("Optimal cost " + String.valueOf(op.getOptimalCost()));

		ind= new IntArrayList();
		val= new DoubleArrayList();
		pos= new int[] {N+Nf, N+Nf, N+Nf, N+Nf};

		//y
		System.out.println("");
		double[][][] ysun= op.getPrimalSolution("y_sun").to3DArray();
		Y_sun= new DoubleMatrixND(ysun);
		System.out.println("y_sun: "+ String.valueOf(Y_sun.zSum()));
		Y_sun.getNonZeros(ind, val);
		pos= new int[]{Ns,Nf,N};
		int act[];    //declaring array
		act = new int[N];  // allocating memory to array
		for (int i=0; i< ind.size();i++)
		{	 
			System.out.println("y_sun("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+")="+String.valueOf(val.get(i))); 
			act[(int)DoubleMatrixND.ind2sub(ind.get(i), pos).get(2)]=1;
		}

		//calculate active nodes
		int act_n=0;
		for (int i =0; i<N;i++)
		{
			if (act[i]==1)
				act_n++;
		}


		//l
		System.out.println("");
		double[][][] lsbe= op.getPrimalSolution("l_sbe").to3DArray();
		DoubleMatrixND LSBE= new DoubleMatrixND(lsbe);
		System.out.println("l_sbe: "+ String.valueOf(LSBE.zSum()));
		LSBE.getNonZeros(ind, val);
		pos= new int[]{Ns,N,N};
		for (int i=0; i< ind.size();i++)
		{	 
			System.out.println("l_sbe("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+")="+String.valueOf(val.get(i))); 
		}

		//b
		System.out.println("");
		System.out.println("b_suv: "+ String.valueOf(b_suv.zSum()));
		b_suv.getNonZeros(ind, val);
		pos= new int[]{Ns,Nf,Nf};
		for (int i=0; i< ind.size();i++)
		{	 
			System.out.println("b_suv("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+")="+String.valueOf(val.get(i))); 
		}

		//h
		System.out.println("");
		System.out.println("h_suvbe: "+ String.valueOf(op.getPrimalSolution("h_suvbe").zSum()));
		op.getPrimalSolution("h_suvbe").getNonZeros(ind, val);
		pos= new int[]{Ns,Nf,Nf,N,N};
		for (int i=0; i< ind.size();i++)
		{	 
			System.out.println("h_suvbe("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(3))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(4))+")="+String.valueOf(val.get(i))); 
		}

		//p
		System.out.println("");
		System.out.println("p_sbeij: "+ String.valueOf(op.getPrimalSolution("p_sbeij").zSum()));
		op.getPrimalSolution("p_sbeij").getNonZeros(ind, val);
		pos= new int[]{Ns,N,N,N,N};
		for (int i=0; i< ind.size();i++)
		{	 
			System.out.println("p_sbeij("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(3))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(4))+")="+String.valueOf(val.get(i))); 
		}

		//pi_sbeij
		System.out.println("");
		System.out.println("pi_sbeij: "+ String.valueOf(op.getPrimalSolution("pi_sbeij").zSum()));
		op.getPrimalSolution("pi_sbeij").getNonZeros(ind, val);
		pos= new int[]{Ns,N,N,N,N};
		for (int i=0; i< ind.size();i++)
		{	 
			System.out.println("pi_sbeij("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(3))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(4))+")="+String.valueOf(val.get(i))); 
		}

		//t
		System.out.println("");
		System.out.println("t_sbeij: "+ String.valueOf(op.getPrimalSolution("t_sbeij").zSum()));
		op.getPrimalSolution("t_sbeij").getNonZeros(ind, val);
		pos= new int[]{Ns,N,N,N,N};
		for (int i=0; i< ind.size();i++)
		{	 
			System.out.println("t_sbeij("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(3))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(4))+")="+String.valueOf(val.get(i))); 
		}

		//ti
		System.out.println("");
		System.out.println("ti_sbeij: "+ String.valueOf(op.getPrimalSolution("ti_sbeij").zSum()));
		op.getPrimalSolution("ti_sbeij").getNonZeros(ind, val);
		pos= new int[]{Ns,N,N,N,N};
		for (int i=0; i< ind.size();i++)
		{	 
			System.out.println("ti_sbeij("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(3))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(4))+")="+String.valueOf(val.get(i))); 
		}





		//X_ij
		System.out.println("");
		System.out.println("X_ij: "+ String.valueOf(op.getPrimalSolution("X_ij").zSum()));
		op.getPrimalSolution("X_ij").getNonZeros(ind, val);
		pos= new int[]{N,N};
		for (int i=0; i< ind.size();i++)
		{	 
			System.out.println("X_ij("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+")="+String.valueOf(val.get(i))); 
		}

		//Xi_ij
		System.out.println("");
		System.out.println("Xi_ij: "+ String.valueOf(op.getPrimalSolution("Xi_ij").zSum()));
		op.getPrimalSolution("Xi_ij").getNonZeros(ind, val);
		pos= new int[]{N,N};
		for (int i=0; i< ind.size();i++)
		{	 
			System.out.println("Xi_ij("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+")="+String.valueOf(val.get(i))); 
		}


		//Z_ijmn
		System.out.println("");
		System.out.println("Z_ijmn: "+ String.valueOf(op.getPrimalSolution("Z_ijmn").zSum()));
		op.getPrimalSolution("Z_ijmn").getNonZeros(ind, val);
		pos= new int[]{N,N,N,N};
		for (int i=0; i< ind.size();i++)
		{	 
			System.out.println("Z_ijmn("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(3))+")="+String.valueOf(val.get(i))); 
		}

		//Zeta_ijmn
		System.out.println("");
		System.out.println("Zeta_ijmn: "+ String.valueOf(op.getPrimalSolution("Zeta_ijmn").zSum()));
		op.getPrimalSolution("Zeta_ijmn").getNonZeros(ind, val);
		pos= new int[]{N,N,N,N};
		for (int i=0; i< ind.size();i++)
		{	 
			System.out.println("Zeta_ijmn("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(3))+")="+String.valueOf(val.get(i))); 
		}

		//K_sun
		System.out.println("");
		System.out.println("K_sun: "+ String.valueOf(op.getPrimalSolution("K_sun").zSum()));
		op.getPrimalSolution("K_sun").getNonZeros(ind, val);
		pos= new int[]{Ns,Nf,N};
		for (int i=0; i< ind.size();i++)
		{	 
			System.out.println("K_sun("+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(0))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(1))+","+String.valueOf(DoubleMatrixND.ind2sub(ind.get(i), pos).get(2))+")="+String.valueOf(val.get(i))); 
		}
		
		
		int wlchann= (int) op.getPrimalSolution("Zeta_ijmn").zSum()+ (int) op.getPrimalSolution("Z_ijmn").zSum();
		int protec_wlchann= (int) op.getPrimalSolution("Zeta_ijmn").zSum();
		int w_wlchann= (int) op.getPrimalSolution("Z_ijmn").zSum();
		
		double tot_b= b_suv.zSum();
		double gen_t= op.getPrimalSolution("l_sbe").zSum();

		//% 1 isolation, 2 nslices, 3 time, 4 tot_wl, 5 w_wl, 6 p_wl,
		//% 7 total_bandwidth, 8 generated_traffic, 9 act_nodes, 10 load, 11 virtual links,
		//% 12 avg_wl_per_link, 13 avg_wl_load,14 kmax
		
		return new ReliabilityResult(isolation.getInt(),S.getInt(),diff,wlchann,w_wlchann, protec_wlchann,tot_b,gen_t, (int)op.getPrimalSolution("y_n").zSum(), op.getPrimalSolution("K_sun").zSum(), b_suv.getNumberOfNonZeros() ,Y_sun.getNumberOfNonZeros(),(double) 0.0,0);

		//  ReliabilityResult(int isolation, int nslices, double time, int wlchann, int wwlchann, int pwlchann, double total_bandwidth, double generated_traffic,  int act_nodes, double load, int conn, double avg_wl_per_link, double avg_wl_load, int kmax) {
        
		// return "Ok!: The solution found is guaranteed to be optimal: " +
		// op.solutionIsOptimal() + ". Number routes = " + netPlan.getNumberOfRoutes();
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


	public void loadInput(NetPlan netPlan, double mean, double std, int size, int seed)
	{
	
		this.slice_Requests= new HashMap<Integer, Pair<Integer,Double>>();
		//Load input for comparisons
		try {
			String filename=  ("./InputN"+String.valueOf(N)+"B"+String.valueOf(mean)+"S"+String.valueOf(std)+"seed"+String.valueOf(seed));
			System.out.println(filename);
			
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
		FileOutputStream fileOut = new FileOutputStream("./InputN"+String.valueOf(N)+"B"+String.valueOf(mean)+"S"+String.valueOf(std)+"seed"+String.valueOf(seed));
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


	public void createInput_old(NetPlan netPlan, int nslices, int system_bandwidth,int isolation)
	{
		
		//Store indexes of RUs,DUs,CUs,NGCs
				//we will use for slice requests generation
				ArrayList<Integer> RUs = new ArrayList<Integer>(); 
				ArrayList<Integer> DUs = new ArrayList<Integer>(); 
				ArrayList<Integer> CUs = new ArrayList<Integer>(); 
				ArrayList<Integer> NGCs = new ArrayList<Integer>(); 
	
				//for (int i=0;i<N;i++)
				//this.Capn.add((double) Cn.getInt());
	
				for (int i = 0; i < N; i++) {
					this.Capn.add(0.0);
	
					try {
	
						if (netPlan.getNode(i).getAttribute("RU").toString().equals("1")) 
						{RUs.add(i);this.Capn.set(i, 10.0); }
	
						if (netPlan.getNode(i).getAttribute("DU").toString().equals("1")) 
						{DUs.add(i);this.Capn.set(i, 10.0); }
	
						if (netPlan.getNode(i).getAttribute("CU").toString().equals("1")) 
						{CUs.add(i);this.Capn.set(i, 10.0); }
	
						if (netPlan.getNode(i).getAttribute("NGC").toString().equals("1")) 
						{NGCs.add(i); this.Capn.set(i, 1000.0); }
					}
					catch (NullPointerException e) {
						if (Math.random() < 1)
						{RUs.add(i);this.Capn.set(i, 10.0); }
						if (Math.random() < 0.5)
						{DUs.add(i);this.Capn.set(i, 10.0); }
						if (Math.random() < 0.5)
						{CUs.add(i); this.Capn.set(i, 10.0); }
						if (Math.random() < 0.5)
						{NGCs.add(i);this.Capn.set(i, 1000.0); }
	
					}
	
	
				}
		
				
		double RU_cost=1;
		double DU_cost=0.4;
		double CU_cost=0.2;
		double NGC_cost=0.5;
		
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
					}
					for (int j =0 ; j< NGCs.size();j++)
					{
						M_sun.set(new int[] { i, 3, NGCs.get(j)}, 1 );
					}
					K_su.set(i, 0, RU_cost);
					K_su.set(i, 1, DU_cost);
					K_su.set(i, 2, CU_cost);
					K_su.set(i, 3, NGC_cost);
					
					
					
					//set 3 links bandwidth
					b_suv.set(new int[] { i, 0, 1 }, fronthaul*slice_weight);
					//b_suv.set(new int[] { i, 1, 0 }, 20);//biderctionality
					b_suv.set(new int[] { i, 1, 2 }, midhaul*slice_weight);
					//b_suv.set(new int[] { i, 2, 1 }, 20);//biderctionality
					b_suv.set(new int[] { i, 2, 3 }, backhaul*slice_weight);
					
					//dedicated transport is only on the non common part
					D.add(0);
	
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
							K_su.set(i, 0, 0); //the cost is added only first time
							K_su.set(i, 1, DU_cost);
							K_su.set(i, 2, CU_cost);
							K_su.set(i, 3, NGC_cost);
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
						K_su.set(i, 0, RU_cost); //the cost is added only first time
						K_su.set(i, 1, DU_cost);
						K_su.set(i, 2, CU_cost);
						K_su.set(i, 3, NGC_cost);
					}
					
					//dedicated transport is only on the non common part
					D.add(0);
	
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
							K_su.set(i, 2, CU_cost);
							K_su.set(i, 3, NGC_cost);
							
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
						K_su.set(i, 2, CU_cost);
						K_su.set(i, 3, NGC_cost);
					}
					
					
				
					
					
					
					//dedicated transport is only on the non common part
					D.add(0);
	
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
							K_su.set(i, 3, NGC_cost);
							
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
						K_su.set(i, 3, NGC_cost);
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
							K_su.set(i, 3, NGC_cost);
							
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
						K_su.set(i, 3, NGC_cost);
					}
					
					
					//dedicated transport is only on the non common part
					D.add(0);
	
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
				bc_suv.set(new int[] {cseg.get(cs).getFirst(),0,1}, (max_b-cseg.get(cs).getThird())/max_b *fronthaul);
				
				Mc_sun.set(new int[] {cseg.get(cs).getFirst(),0,cseg.get(cs).getSecond() }, 1);
				
				for (int j =0 ; j< DUs.size();j++)
				{
					Mc_sun.set(new int[] { cseg.get(cs).getFirst(), 1, DUs.get(j)}, 1 );
				}
				
				Kc_su.set(cseg.get(cs).getFirst(), 0, RU_cost);
				Kc_su.set(cseg.get(cs).getFirst(), 1, DU_cost);
				//dedicated transport is only on the non common part
				cD.add(0);
	
				//reliability?
				cR.add(1);
	
			}
			
			break;
		case 3:
		case 4:
			for (int cs =0; cs<cseg.size();cs++)
			{
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
				
				
				Kc_su.set(cseg.get(cs).getFirst(), 0, RU_cost);
				Kc_su.set(cseg.get(cs).getFirst(), 1, DU_cost);
				Kc_su.set(cseg.get(cs).getFirst(), 2, CU_cost);
				//dedicated transport is only on the non common part
				cD.add(0);
	
				//reliability?
				cR.add(1);
	
			}
			
			break;
		case 5:
			for (int cs =0; cs<cseg.size();cs++)
			{
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
				
				Kc_su.set(cseg.get(cs).getFirst(), 0, RU_cost);
				Kc_su.set(cseg.get(cs).getFirst(), 1, DU_cost);
				Kc_su.set(cseg.get(cs).getFirst(), 2, CU_cost);
				Kc_su.set(cseg.get(cs).getFirst(), 3, NGC_cost);
				//dedicated transport is only on the non common part
				cD.add(0);
	
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


	public void createInput(NetPlan netPlan, int nslices, int system_bandwidth,int isolation,boolean grooming)
	{
		
		//Store indexes of RUs,DUs,CUs,NGCs
				//we will use for slice requests generation
				ArrayList<Integer> RUs = new ArrayList<Integer>(); 
				ArrayList<Integer> DUs = new ArrayList<Integer>(); 
				ArrayList<Integer> CUs = new ArrayList<Integer>(); 
				ArrayList<Integer> NGCs = new ArrayList<Integer>(); 
	
				//for (int i=0;i<N;i++)
				//this.Capn.add((double) Cn.getInt());
	
				for (int i = 0; i < N; i++) {
					this.Capn.add(0.0);
	
					try {
	
						if (netPlan.getNode(i).getAttribute("RU").toString().equals("1")) 
						{
							RUs.add(i);this.Capn.set(i, 1000000000.0); 
							antennas.add(i);
							activeNodes.add(i);
						}
	
						if (netPlan.getNode(i).getAttribute("DU").toString().equals("1")) 
						{DUs.add(i);this.Capn.set(i, 1000000000.0); }
	
						if (netPlan.getNode(i).getAttribute("CU").toString().equals("1")) 
						{CUs.add(i);this.Capn.set(i, 1000000000.0); }
	
						if (netPlan.getNode(i).getAttribute("NGC").toString().equals("1")) 
						{NGCs.add(i); this.Capn.set(i, 1000000000.0); }
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
						{NGCs.add(i);this.Capn.set(i, 1000000000.0); }
	
					}
	
	
				}
				
		Capn_native= new ArrayList<Double>(Capn);
		
				
		double RU_cost=57600.0;
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



}